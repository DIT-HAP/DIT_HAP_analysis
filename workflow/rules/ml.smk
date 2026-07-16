# =============================================================================
# ml.smk — mljar AutoML analysis (target x mode fan-out)
# =============================================================================

# Per-dataset regression of a growth-fitness target (DR/DL) from gene features,
# in two mljar modes (Explain: fast hold-out; Perform: 5-fold CV, slow). Reads
# the raw feature matrix + curated final_clusters directly (byte-faithful to the
# self-contained source notebook — NOT the Task 6 transformed tables).
wildcard_constraints:
    target="DR|DL",
    mode="Explain|Perform",

rule train_automl:
    input:
        feature_matrix="results/features/{pombase_version}/pombe_coding_gene_protein_features.tsv",
        final_clusters="resources/curated/final_clusters.tsv",
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
            --feature-matrix {input.feature_matrix} \
            --final-clusters {input.final_clusters} \
            --target {wildcards.target} \
            --mode {wildcards.mode} \
            --output-dir {params.output_dir} \
            --test-size {params.test_size} \
            --random-state {params.random_state} \
            --total-time-limit {params.total_time_limit} &> {log}
        """
