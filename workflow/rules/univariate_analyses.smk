rule get_ROC_AUC:
    input:
        dataset=rules.get_residuals.output.residuals,
    output:
        auc=expand(
            "tables/univariate_analyses/AUC/seed_{mice_seed}_imputation_{imputation_cycle}.tsv",
            mice_seed=mice_seeds[0],
            imputation_cycle=imputation_cycles[0],
        )[0],
    script:
        "../scripts/univariate_analyses/auc.py"


rule mann_whitney:
    input:
        dataset=rules.get_residuals.output.residuals,
    output:
        mann_whitney=expand(
            "tables/univariate_analyses/mann_whitney/seed_{mice_seed}_imputation_{imputation_cycle}.tsv",
            mice_seed=mice_seeds[0],
            imputation_cycle=imputation_cycles[0],
        )[0],
    script:
        "../scripts/univariate_analyses/mann_whitney.py"
