rule normalize:
    input:
        dataset=rules.remove_metabolites_with_non_imputed_groups.output.dataset,
    output:
        normalized="data/normalization/seed_{mice_seed}/imputation_{imputation_cycle}.xlsx",
    run:
        import numpy as np

        dataset = setup_dataset(input.dataset)
        normalization_method = config["normalization_method"]
        if normalization_method == "None":
            dataset.io.save_excel(output.normalized)
        elif normalization_method == "log":
            dataset.data = dataset.data.apply(lambda x: np.log(x))
            dataset.io.save_excel(output.normalized)
        else:
            raise ValueError("Invalid or not yet implemented normalization method")
