rule prepare_dataset_for_imputation:
    input:
        dataset="data/imputation/input/clean_raw_dataset.xlsx",
    output:
        data_metadata="data/imputation/input/data_metadata.tsv",
        chemical_annotation="data/imputation/input/chemical_annotation.tsv",
        sample_metadata="data/imputation/input/samples.tsv",
    script:
        "../scripts/imputation/prepare_data_for_imputation.py"


rule download_and_setup_imputation_script:
    output:
        script="workflow/scripts/imputation/metabolomics_missing_data_imputation/impute.R" 
    shell:
        "curl https://raw.githubusercontent.com/matteobolner/metabolomics_missing_data_imputation/refs/heads/main/install_required_libraries.R -o setup_libraries.R && Rscript setup_libraries.R && rm setup_libraries.R && curl https://raw.githubusercontent.com/matteobolner/metabolomics_missing_data_imputation/refs/heads/main/impute.R -O {output.script}"

rule impute:
    input:
        data=rules.prepare_dataset_for_imputation.output.data_metadata,
        chemical_annotation=rules.prepare_dataset_for_imputation.output.chemical_annotation,
        script=rules.download_and_setup_imputation_script.output.script
    output:
        imputed="data/imputation/imputed/{mice_seed}.tsv",
    params:
        covariates=get_mice_covariates(),
    shell:
        "Rscript {input.script} -d {input.data} -c {input.chemical_annotation} -o {output.imputed} -s {wildcards.mice_seed} -m pmm -n 5 -r 0.25 -u 5 -a {params.covariates} --metabolite_id_column 'CHEMICAL ID' --super_pathway_column 'SUPER PATHWAY'"


rule split_imputed_datasets:
    input:
        imputed_data=rules.impute.output.imputed,
        samples=rules.prepare_dataset_for_imputation.output.sample_metadata,
        chemical_annotation=rules.prepare_dataset_for_imputation.output.chemical_annotation,
    output:
        imputations=expand(
            "data/imputation/imputed/seed_{{mice_seed}}/imputation_{imputation_cycle}.xlsx",
            imputation_cycle=[i for i in imputation_cycles],
        ),
    run:
        df = pd.read_table(input.imputed_data)
        df = df[df[".imp"] != 0]

        outputs = {
            i + 1: output.imputations[i] for i in range(0, len(config["mice_seeds"]))
        }

        for name, group in df.groupby(by=".imp"):
            group = group.drop(columns=[".id", ".imp"])
            dataset = MetaboTK().io.from_tables(
                sample_metadata=input.samples,
                chemical_annotation=input.chemical_annotation,
                data=group,
                sample_id_column=config["sample_id_column"],
                metabolite_id_column=config["metabolite_id_column"],
            )
            dataset.io.save_excel(outputs[name])
