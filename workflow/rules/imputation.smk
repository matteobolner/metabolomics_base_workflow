module imputation_workflow_module:
    snakefile:
        github(
            "matteobolner/metabolomics_imputation_workflow",
            path="workflow/Snakefile",
            branch="main",
        )
    config:
        config


use rule * from imputation_workflow_module as imp_*


use rule setup_dataset from imputation_workflow_module as imp_setup_dataset with:
    input:
        dataset=config["raw_dataset"]


rule remove_metabolites_with_non_imputed_groups:
    input:
        #dataset="data/imputation/imputed/seed_{mice_seed}/imputation_{imputation_cycle}.xlsx",
        dataset=rules.imp_get_imputations.output.dataset

    output:
        dataset="data/imputation/imputed/seed_{mice_seed}/imputation_{imputation_cycle}_no_missing.xlsx",
    run:
        dataset = setup_dataset(input.dataset)
        dataset = dataset.ops.subset(
            what="metabolites",
            ids=dataset.stats.remove_missing(on="metabolites", threshold=0).columns,
        )
        dataset.io.save_excel(output.dataset)
