rule initial_stats:
    input:
        dataset=config["raw_dataset"],
    output:
        empty_metabolites="data/initial_stats/empty_metabolites.tsv",
        sample_stats="data/initial_stats/sample_stats.tsv",
        metabolite_stats="data/initial_stats/metabolite_stats.tsv",
        pca_plot="figures/initial_stats/pca_grid.png",
    script:
        "../scripts/stats/raw_dataset_stats.py"


rule pca:
    input:
        dataset=rules.get_residuals.output.residuals,
    output:
        pca="tables/pca/selected_{selected}_{cv_selected}.tsv",
        pca_figure="figures/pca/selected_{selected}_{cv_selected}.png",
        pc1_pc2_figure="figures/pca/selected_{selected}_{cv_selected}_PC1_PC2.png",
        pc1_pc3_figure="figures/pca/selected_{selected}_{cv_selected}_PC1_PC3.png",
        explained_variance="tables/pca/selected_{selected}_{cv_selected}_explained_variance.tsv",
    params:
        hue_title=config["group_name"],
    script:
        "../scripts/pca/pca.py"


if config["paired_samples"]:
    rule summarize_results:
        input:
            dataset=expand(
                "data/normalization/seed_{mice_seed}/imputation_{imputation_cycle}.xlsx",
                mice_seed=mice_seeds[0],
                imputation_cycle=imputation_cycles[0],
            )[0],
            boruta=rules.merge_boruta_across_imputations.output.summary,
            wilcoxon=rules.wilcoxon.output.wilcoxon,
            auc=rules.get_ROC_AUC.output.auc,
        output:
            stats="results/metabolite_level_stats.tsv",
        script:
            "../scripts/merge_stats_paired.py"
else:
    rule summarize_results:
        input:
            dataset=expand(
                "data/normalization/seed_{mice_seed}/imputation_{imputation_cycle}.xlsx",
                mice_seed=mice_seeds[0],
                imputation_cycle=imputation_cycles[0],
            )[0],
            boruta=rules.merge_boruta_across_imputations.output.summary,
            mann_whitney=rules.mann_whitney.output.mann_whitney,
            auc=rules.get_ROC_AUC.output.auc,
        output:
            stats="results/metabolite_level_stats.tsv",
        script:
            "../scripts/merge_stats.py"
