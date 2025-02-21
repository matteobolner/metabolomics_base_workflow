rule get_residuals:
    input:
        dataset=rules.normalize.output.normalized,
    output:
        residuals="data/residuals/seed_{mice_seed}/imputation_{imp_cycle}.xlsx",
    run:
        dataset = MetaboTK().io.from_excel(
            input.dataset,
            sample_id_column=config["sample_id_column"],
            metabolite_id_column=config["metabolite_id_column"],
        )
        dataset.data = dataset.lm.get_linear_model_residuals(
            formula=config["lm_formula"], models_path=None
        )
        dataset.io.save_excel(output.residuals)
