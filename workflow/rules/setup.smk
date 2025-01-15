rule initial_stats:
    input:
        dataset=config["raw_dataset"],
    output:
        empty_metabolites="data/initial_stats/empty_metabolites.tsv",
        sample_stats="data/initial_stats/sample_stats.tsv",
        metabolite_stats="data/initial_stats/metabolite_stats.tsv",
    script:
        "../scripts/stats/raw_dataset_stats.py"


rule remove_outliers:
    input:
        dataset=config["raw_dataset"],
    output:
        missing_removed="data/initial_stats/missing_removed.tsv",
        dataset="data/imputation/input/clean_raw_dataset.xlsx",
    script:
        "../scripts/outliers/remove_outliers.py"
