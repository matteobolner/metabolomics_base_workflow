rule pca:
    input:
        dataset="results/{trait}/dataset_with_boruta_AUC_mw.xlsx",
    output:
        pca="tables/pca/selected_{selected}_{cv_selected}.tsv",
        pca_figure="figures/pca/selected_{selected}_{cv_selected}.png",
        pc1_pc2_figure="figures/pca/selected_{selected}_{cv_selected}_PC1_PC2.png",
        pc1_pc3_figure="figures/pca/selected_{selected}_{cv_selected}_PC1_PC3.png",
        explained_variance="tables/pca/selected_{selected}_{cv_selected}_explained_variance.tsv",
    params:
        group_column=config["group_column"],
    script:
        "../scripts/pca/pca.py"

