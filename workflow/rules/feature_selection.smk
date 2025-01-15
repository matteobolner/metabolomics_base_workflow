rule boruta:
    input:
        dataset=rules.get_residuals.output.residuals,
    output:
        ranking="data/feature_selection/boruta/mice_seed_{mice_seed}/imputation_{imp_cycle}/boruta_seed_{boruta_seed}/ranking.tsv",
        long_df="data/feature_selection/boruta/mice_seed_{mice_seed}/imputation_{imp_cycle}/boruta_seed_{boruta_seed}/long_df.tsv",
    script:
        "../scripts/feature_selection/boruta/run_boruta.py"


rule merge_boruta_runs:
    input:
        long_counts=expand(
            "data/feature_selection/boruta/mice_seed_{mice_seed}/imputation_{imp_cycle}/boruta_seed_{boruta_seed}/long_df.tsv",
            mice_seed=mice_seeds,
            imp_cycle=imputation_cycles,
            boruta_seed=feature_selection_seeds,
        ),
    output:
        long_df="data/feature_selection/boruta/merged/long_data.tsv",
        summary="data/feature_selection/boruta/merged/summary.tsv",
    script:
        "../scripts/feature_selection/boruta/merge_boruta_runs.py"


rule prepare_CV_datasets:
    input:
        dataset=rules.get_residuals.output.residuals,
    output:
        directory=directory(
            "data/residuals/seed_{mice_seed}/imputation_{imp_cycle}_CV/"
        ),
        training_sets=expand(
            "data/residuals/seed_{{mice_seed}}/imputation_{{imp_cycle}}_CV/{fold}_train.xlsx",
            fold=[i for i in range(1, 11)],
        ),
    run:
        dataset = MetaboTK().io.from_excel(
            input.dataset,
            sample_id_column=config["sample_id_column"],
            metabolite_id_column=["metabolite_id_column"],
        )
        split_kfold = dataset.fs.stratified_kfold(
            n_splits=10, stratification_column="H_L", output_dir=output.directory
        )


rule boruta_CV:
    input:
        dataset="data/residuals/seed_{mice_seed}/imputation_{imp_cycle}_CV/{fold}_train.xlsx",
    output:
        ranking="data/feature_selection/boruta/mice_seed_{mice_seed}/imputation_{imp_cycle}_CV/fold_{fold}/boruta_seed_{boruta_seed}/ranking.tsv",
        long_df="data/feature_selection/boruta/mice_seed_{mice_seed}/imputation_{imp_cycle}_CV/fold_{fold}/boruta_seed_{boruta_seed}/long_df.tsv",
    script:
        "../scripts/feature_selection/boruta/run_boruta.py"


rule merge_boruta_runs_CV:
    input:
        long_counts=expand(
            "data/feature_selection/boruta/mice_seed_{mice_seed}/imputation_{imp_cycle}_CV/fold_{fold}/boruta_seed_{boruta_seed}/long_df.tsv",
            mice_seed=mice_seeds,
            imp_cycle=imputation_cycles,
            boruta_seed=feature_selection_seeds,
            fold=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        ),
    output:
        long_df="data/feature_selection/boruta/merged/long_data_CV.tsv",
        summary="data/feature_selection/boruta/merged/summary_CV.tsv",
    script:
        "../scripts/feature_selection/boruta/merge_boruta_runs.py"


rule summarize_boriuta:
    input:
        cv="data/feature_selection/boruta/merged/long_data_CV.tsv",
        long_df="data/feature_selection/boruta/merged/long_data.tsv",
        summary="data/feature_selection/boruta/merged/summary.tsv",
    output:
        long_df="data/feature_selection/boruta/merged_long_data.tsv",
        summary="tables/feature_selection/boruta/summary_with_cv.tsv",
    script:
        "../scripts/feature_selection/boruta/summarize.py"
