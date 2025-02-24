rule get_ROC_AUC:
    input:
        dataset=expand(
            rules.get_residuals.output.residuals,
            mice_seed=mice_seeds[0],
            imputation_cycle=imputation_cycles[0],
        )[0],
    output:
        auc="tables/univariate_analyses/AUC.tsv",
    script:
        "../scripts/univariate_analyses/auc.py"


rule mann_whitney:
    input:
        dataset=expand(
            rules.get_residuals.output.residuals,
            mice_seed=mice_seeds[0],
            imputation_cycle=imputation_cycles[0],
        )[0],
    output:
        mann_whitney="tables/univariate_analyses/mann_whitney.tsv",
    script:
        "../scripts/univariate_analyses/mann_whitney.py"
