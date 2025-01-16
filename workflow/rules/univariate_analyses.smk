rule get_ROC_AUC:
    input:
        dataset=rules.get_residuals.output.residuals,
    output:
        auc="data/univariate_analyses/AUC/seed_{config['imputation_seeds'][0]}_imputation_1.tsv",
    script:
        "../scripts/univariate_analyses/auc.py"


rule mann_whitney:
    input:
        dataset="data/extremes/{trait}/seed_1000_imputation_1.xlsx",
    output:
        mann_whitney="data/univariate_analyses/{trait}/mann_whitney/seed_1000_imputation_1.tsv",
    script:
        "../scripts/univariate_analyses/mann_whitney.py"
