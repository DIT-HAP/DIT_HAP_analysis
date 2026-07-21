# =============================================================================
# ml.smk — mljar AutoML analysis (target x mode fan-out)
# =============================================================================
#
# Per-dataset regression of a growth-fitness target (DR/DL) from gene features,
# in two mljar modes (Explain: fast hold-out; Perform: 5-fold CV, slow).
#
# Split: prepare_ml_data merges the feature matrix + the SELECTED finalize
# variant's final_clusters (config.clustering.selected_variant + per-dataset
# override) and applies the DR filter ONCE (shared spine); train_automl then reads
# that pickle for each target x mode instead of re-merging. Byte-faithful to the
# self-contained source notebook (its own train-only PowerTransform).

_MLWORK = "results/ml/models/{dataset}/{pombase_version}/_work"

wildcard_constraints:
    target="DR|DL",
    mode="Explain|Perform",


# --- Shared modeling-data spine (merge + DR filter, once per dataset) ---
rule prepare_ml_data:
    input:
        feature_matrix="results/features/{pombase_version}/pombe_coding_gene_protein_features.tsv",
        final_clusters=lambda wc: final_clusters_path(wc.dataset, selected_variant(wc.dataset)),
    output:
        modeling_data=f"{_MLWORK}/modeling_data.pkl",
    params:
        dr_filter=config.get("ml", {}).get("dr_filter", 0.3),
    log:
        "logs/ml/prepare_ml_data_{dataset}_{pombase_version}.log",
    conda:
        "../envs/machine_learning.yml"
    message:
        "*** [ml] Preparing shared modeling data for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/ml/prepare_ml_data.py \
            --feature-matrix {input.feature_matrix} \
            --final-clusters {input.final_clusters} \
            --output {output.modeling_data} \
            --dr-filter {params.dr_filter} &> {log}
        """


# --- One AutoML model (fanned out by target x mode) ---
rule train_automl:
    input:
        modeling_data=f"{_MLWORK}/modeling_data.pkl",
    output:
        metrics="results/ml/models/{dataset}/{pombase_version}/{target}_{mode}/metrics.tsv",
        importance="results/ml/models/{dataset}/{pombase_version}/{target}_{mode}/features_importance.csv",
        plot="results/ml/models/{dataset}/{pombase_version}/{target}_{mode}/prediction_and_residuals.pdf",
    params:
        output_dir="results/ml/models/{dataset}/{pombase_version}/{target}_{mode}",
        test_size=config.get("ml", {}).get("test_size", 0.2),
        random_state=config.get("ml", {}).get("random_state", 42),
        total_time_limit=config.get("ml", {}).get("total_time_limit", 14400),
    log:
        "logs/ml/train_automl_{dataset}_{pombase_version}_{target}_{mode}.log",
    conda:
        "../envs/machine_learning.yml"
    message:
        "*** Training AutoML for {wildcards.dataset} target={wildcards.target} mode={wildcards.mode}..."
    shell:
        """
        python workflow/scripts/ml/train_automl.py \
            --modeling-data {input.modeling_data} \
            --target {wildcards.target} \
            --mode {wildcards.mode} \
            --output-dir {params.output_dir} \
            --test-size {params.test_size} \
            --random-state {params.random_state} \
            --total-time-limit {params.total_time_limit} &> {log}
        """
