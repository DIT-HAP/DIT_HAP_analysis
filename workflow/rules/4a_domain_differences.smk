# =============================================================================
# 4a_domain_differences.smk — Intra-gene DR segmentation (TGUH)
# =============================================================================
#
# `4a` is the <chapter><section> prefix convention (see the Includes block in
# Snakefile). It carries into results/4a_domain_differences/ and
# logs/4a_domain_differences/; the script and library directories keep the
# unprefixed name (workflow/scripts/domain_differences/, workflow/src/domain_differences/).
#
# Per dataset: position every in-gene insertion along its spliced CDS, then cut
# each gene's DR profile — ordered by CDS coordinate, not by genome coordinate —
# into pieces of constant mean with TGUH (tail-greedy unbalanced Haar). A boundary
# is called wherever two neighbouring pieces differ by more than the merge
# threshold
#
#     lambda = z * sigma_gene + sigma_add
#
# with sigma_gene the robust point-to-point noise of that gene's own profile.
# Both z and sigma_add are config (config.domain_differences). The scan is
# unbiased: it reads no domain annotation at all, and the hand-curated termini in
# resources/curated/ only ever serve as an outside check on what it found.
#
# This stage is a rewrite of the experiments in tmp/domain_scan/ — the README there
# documents every measurement behind the parameter choice (why by-site, why
# min_seg=1, why sigma is computed per gene) and tmp/domain_scan/explainer/ walks
# through the algorithm itself. That directory is deliberately frozen: it is the
# record of how the choices here were reached, not an input to this pipeline.
#
# Rules, in DAG order:
#   prepare_domain_insertions     -> _work/insertions.parquet (one row per
#                                    Forward/Reverse insertion, positioned on the CDS)
#   build_domain_hand_regions     -> domain_hand_regions.tsv (the curated spreadsheets
#                                    as residue intervals; an input to the plot and a
#                                    coverage log for the scan, never a filter)
#   compute_domain_segments       -> domain_segments.tsv (one row per segment),
#                                    domain_summary.tsv (one row per reported gene),
#                                    _work/segments.parquet (index ranges)
#   compute_domain_permutation_p  -> domain_permutation_p.tsv (OPTIONAL, see below)
#   plot_domain_segments          -> domain_segments.pdf + .review.png (each panel
#                                    titled with its gene, FYPOviability and gRNA DR)
#
# THE PERMUTATION RULE IS OFF BY DEFAULT. It is the stage's only expensive step
# and its two p columns are reference numbers — the analysis is filtered by the
# effect floor alone, so nothing downstream reads them and the figure never
# depends on this rule. config.domain_differences.permutation_p.enabled decides
# whether the rule is even defined, so with the default `false` a target naming
# domain_permutation_p.tsv fails with MissingRuleException instead of quietly
# computing something the config switched off.
#
# Data-format rule: per-stage intermediates are Parquet; the human-facing tables
# are TSV. The plot rule follows figures.py's house style, and because the gallery
# is one PDF page per group of genes it writes that PDF itself (PdfPages) rather
# than through save_dual(), whose PDF half would clobber it — the .review.png
# sibling is still declared as an output so `--delete-all-output` cleans it.
#
# DATA-PATH NOTE: the release tables come from upstream DIT_HAP_snakemake via
# DATASETS['datasets'][dataset]['release_dir']; pombase from
# resources/external/pombase/{DATASETS['reference']['pombase_version']}. The plot
# also reads 1c_annotate.smk's _ANNOT_REF — the annotation reference is where the
# panel headings' FYPOviability and gRNA_DR already meet, so this stage joins them
# there rather than re-deriving the join from the release table plus the curated
# gRNA spreadsheet.

_DOMAIN_CFG = config.get("domain_differences", {})
_DOMAIN = "results/4a_domain_differences/{dataset}"
_DOMAIN_WORK = f"{_DOMAIN}/_work"

# Read once here rather than as a `.get(key, default)` in every params block: the
# numbers ARE the analysis, and a threshold that differs between two rules of the
# same stage would not be a warning but a wrong result. Changing a key still
# rebuilds everything downstream — Snakemake's default rerun-triggers include
# `params`.
_DOMAIN_Z = _DOMAIN_CFG.get("z", 1.5)
_DOMAIN_SIGMA_ADD = _DOMAIN_CFG.get("sigma_add", 0.1)
_DOMAIN_MIN_SEG = _DOMAIN_CFG.get("min_seg", 1)
_DOMAIN_BY_SITE = _DOMAIN_CFG.get("by_site", True)
_DOMAIN_PERMUTATION = _DOMAIN_CFG.get("permutation_p", {})
_DOMAIN_PERMUTATION_ON = _DOMAIN_PERMUTATION.get("enabled", False)

# argparse's BooleanOptionalAction cannot parse the Python `True`/`False` a params
# value would render as, so the flag is spelled out; the two rules that take it
# must agree, or the permutation rule would re-derive a profile the segmentation
# rule never used and its index ranges would point at the wrong points.
_DOMAIN_SITE_FLAG = "--by-site" if _DOMAIN_BY_SITE else "--no-by-site"

# The hand annotation, and the protein lengths that turn "the last 20 residues"
# into a residue interval.
_DOMAIN_HAND_XLSX = "resources/curated/non_esential_domain_candidates.xlsx"
_DOMAIN_PEPTIDE_STATS = (
    f"resources/external/pombase/{DATASETS['reference']['pombase_version']}"
    f"/Protein_features/peptide_stats.tsv"
)


rule prepare_domain_insertions:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/insertion_level/fitting_results.tsv"
        ),
        annotations=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/insertion_level/annotations.tsv.gz"
        ),
        gene_level=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
    output:
        f"{_DOMAIN_WORK}/insertions.parquet",
    log:
        "logs/4a_domain_differences/prepare_insertions_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [domain_differences] Preparing in-gene insertions for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/domain_differences/prepare_domain_insertions.py \
            --fitting-results {input.fitting_results} \
            --annotations {input.annotations} \
            --gene-level {input.gene_level} \
            --output {output} &> {log}
        """


# Reads the spreadsheet with pandas' default Excel engine, which needs xlrd for
# .xls and openpyxl for .xlsx — hence the biopython env rather than either of the
# other two, which are the environments the stage's compute and plot rules use.
rule build_domain_hand_regions:
    input:
        xlsx=_DOMAIN_HAND_XLSX,
        peptide_stats=_DOMAIN_PEPTIDE_STATS,
        gene_level=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
    output:
        f"{_DOMAIN}/domain_hand_regions.tsv",
    log:
        "logs/4a_domain_differences/build_hand_regions_{dataset}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [domain_differences] Building hand-curated regions for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/domain_differences/build_domain_hand_regions.py \
            --xlsx {input.xlsx} \
            --peptide-stats {input.peptide_stats} \
            --gene-level {input.gene_level} \
            --output {output} &> {log}
        """


rule compute_domain_segments:
    input:
        insertions=f"{_DOMAIN_WORK}/insertions.parquet",
        hand_regions=f"{_DOMAIN}/domain_hand_regions.tsv",
        # the summary's coverage and context columns: protein lengths from PomBase,
        # and the 1c annotation reference the gallery titles also read
        peptide_stats=_DOMAIN_PEPTIDE_STATS,
        annotation=_ANNOT_REF,
    output:
        segments=f"{_DOMAIN}/domain_segments.tsv",
        summary=f"{_DOMAIN}/domain_summary.tsv",
        ranges=f"{_DOMAIN_WORK}/segments.parquet",
    params:
        min_seg=_DOMAIN_MIN_SEG,
        z=_DOMAIN_Z,
        sigma_add=_DOMAIN_SIGMA_ADD,
        site_flag=_DOMAIN_SITE_FLAG,
    log:
        "logs/4a_domain_differences/compute_segments_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [domain_differences] Segmenting DR profiles for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/domain_differences/compute_domain_segments.py \
            --insertions {input.insertions} \
            --hand-regions {input.hand_regions} \
            --peptide-stats {input.peptide_stats} \
            --annotation {input.annotation} \
            --min-seg {params.min_seg} \
            --z {params.z} \
            --sigma-add {params.sigma_add} \
            {params.site_flag} \
            --output-segments {output.segments} \
            --output-summary {output.summary} \
            --output-work {output.ranges} &> {log}
        """


# Defined only when the config asks for it — see the header. Run it with, e.g.:
#   snakemake --use-conda results/4a_domain_differences/HD_DIT_HAP/domain_permutation_p.tsv
if _DOMAIN_PERMUTATION_ON:
    rule compute_domain_permutation_p:
        input:
            insertions=f"{_DOMAIN_WORK}/insertions.parquet",
            ranges=f"{_DOMAIN_WORK}/segments.parquet",
        output:
            f"{_DOMAIN}/domain_permutation_p.tsv",
        params:
            n_permutations=_DOMAIN_PERMUTATION.get("n_permutations", 1000),
            min_seg=_DOMAIN_MIN_SEG,
            site_flag=_DOMAIN_SITE_FLAG,
        log:
            "logs/4a_domain_differences/compute_permutation_{dataset}.log",
        conda:
            "../envs/statistics_and_figure_plotting.yml"
        message:
            "*** [domain_differences] Permuting DR profiles for {wildcards.dataset}..."
        shell:
            """
            python workflow/scripts/domain_differences/compute_domain_permutation_p.py \
                --insertions {input.insertions} \
                --segments {input.ranges} \
                --n-permutations {params.n_permutations} \
                --min-seg {params.min_seg} \
                {params.site_flag} \
                --output {output} &> {log}
            """


# Depends on the three tables and the raw insertions, never on the permutation
# table: the p's are a separate artifact by design, so the figure stays buildable
# with the optional rule switched off.
rule plot_domain_segments:
    input:
        segments=f"{_DOMAIN}/domain_segments.tsv",
        summary=f"{_DOMAIN}/domain_summary.tsv",
        hand_regions=f"{_DOMAIN}/domain_hand_regions.tsv",
        insertions=f"{_DOMAIN_WORK}/insertions.parquet",
        # the project's annotation layer, for the panel headings' FYPOviability and
        # gRNA_DR — _ANNOT_REF is read from 1c_annotate.smk, which the Snakefile
        # includes ahead of this file
        annotation=_ANNOT_REF,
    output:
        pdf=f"{_DOMAIN}/domain_segments.pdf",
        review=f"{_DOMAIN}/domain_segments.review.png",
    log:
        "logs/4a_domain_differences/plot_segments_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [domain_differences] Drawing the segment gallery for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/domain_differences/plot_domain_segments.py \
            --segments {input.segments} \
            --summary {input.summary} \
            --hand-regions {input.hand_regions} \
            --insertions {input.insertions} \
            --annotation {input.annotation} \
            --output {output.pdf} &> {log}
        """
