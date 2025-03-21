rule boruta:
    input:
        dataset=rules.get_residuals.output.residuals,
    output:
        long_df="data/feature_selection/boruta/mice_seed_{mice_seed}/imp_cycle_{imputation_cycle}/long_df.tsv",
        summary="data/feature_selection/boruta/mice_seed_{mice_seed}/imp_cycle_{imputation_cycle}/summary.tsv",
    script:
        "../scripts/feature_selection/boruta.py"


rule merge_boruta_across_imputations:
    input:
        datasets=expand(
            rules.boruta.output.summary,
            mice_seed=mice_seeds,
            imputation_cycle=imputation_cycles,
        ),
    output:
        summary="tables/feature_selection/boruta/summary.tsv",
    run:
        from functools import reduce

        dfs = [
            pd.read_table(i).set_index("metabolite").sort_index()
            for i in input.datasets
        ]
        summed = reduce(lambda a, b: a.add(b, fill_value=0), dfs)
        summed.to_csv(output.summary, index=True, sep="\t")


        """


rule annotate_dataset_with_boruta_results:
    input:
        summary=rules.summarize_boruta.output.summary,
        dataset=expand(
            rules.get_residuals.output.residuals,
            mice_seed=mice_seeds[0],
            imp_cycle=imputation_cycles[0],
        ),
    output:
        dataset="results/feature_selection/boruta/annotated_dataset.xlsx",
    run:
        dataset = setup_dataset(input.dataset[0])
        summary = pd.read_table(input.summary)
        summary.columns = ["boruta_" + i for i in summary.columns]
        summary = summary.rename(
            columns={"boruta_metabolite": config["metabolite_id_column"]}
        )
        dataset.chemical_annotation = dataset.chemical_annotation.merge(
            summary, left_index=True, right_on=config["metabolite_id_column"]
        ).set_index(config["metabolite_id_column"])
        dataset.io.save_excel(output.dataset)
"""
