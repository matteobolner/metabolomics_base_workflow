rule normalize:
    input:
        dataset="data/imputation/imputed/seed_{mice_seed}/imputation_{imp_cycle}.xlsx",
    output:
        normalized="data/normalization/seed_{mice_seed}/imputation_{imp_cycle}.xlsx",
    run:
        dataset = MetaboTK.io.from_excel(
            file_path,
            sample_id_column=config["sample_id_column"],
            metabolite_id_column=config["metabolite_id_column"],
            sample_metadata_sheet=config["sample_metadata_sheed"],
            chemical_annotation_sheet=config["chemical_annotation_sheet"],
            data_sheet=config["data_sheet"],
        )

        normalization_method = config["normalization method"]
        if normalization_method == "None":
            dataset.io.save_excel(output.normalized)
        elif normalization_method == "log":
            dataset.data = dataset.data.apply(lambda x: np.log(x))
            dataset.io.save_excel(output.normalized)
        else:
            raise ValueError("Invalid or not yet implemented normalization method")
