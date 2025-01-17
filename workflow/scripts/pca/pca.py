from metabotk import MetaboTK
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

group_column = snakemake.params.group_column

dataset = MetaboTK().io.from_excel(
    snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
    sample_metadata_sheet=snakemake.config["sample_metadata_sheet"],
    chemical_annotation_sheet=snakemake.config["chemical_annotation_sheet"],
    data_sheet=snakemake.config["data_sheet"],
)

selected = dataset.chemical_annotation[
    (
        dataset.chemical_annotation["boruta_Confirmed"]
        >= int(snakemake.wildcards.selected)
    )
    & (
        dataset.chemical_annotation["boruta_CV_confirmed"]
        >= int(snakemake.wildcards.cv_selected)
    )
]
dataset = dataset.ops.subset("metabolites", ids=selected.index)

new_title = snakemake.params.hue_title

if len(dataset.metabolites) > 3:
    pca, pca_object = dataset.dimred.get_pca(n_components=3, get_pca_object=True)
    pca.to_csv(snakemake.output.pca, sep="\t")
    explained_variance = pd.DataFrame(
        pca_object.explained_variance_ratio_ * 100,
        columns=["Explained variance (%)"],
        index=["PC1", "PC2", "PC3"],
    )
    explained_variance.index.name = "PC"
    explained_variance.to_csv(snakemake.output.explained_variance, sep="\t", index=True)
    pca_grid = dataset.viz.plot_pca_grid(
        hue=group_column, savepath=snakemake.output.pca_figure
    )
    plt.clf()
    pc1_pc2 = dataset.viz.plot_pca(
        hue=group_column,
        x="PC1",
        y="PC2",
    )
    pc1_pc2.legend().set_title(new_title)
    pc1_pc2.figure.savefig(snakemake.output.pc1_pc2_figure)
    pc1_pc2.figure.clf()
    pc1_pc3 = dataset.viz.plot_pca(
        hue=group_column,
        x="PC1",
        y="PC3",
    )
    pc1_pc3.legend().set_title(new_title)
    pc1_pc3.figure.savefig(snakemake.output.pc1_pc3_figure)
    pc1_pc3.figure.clf()
else:
    pd.DataFrame().to_csv(snakemake.output.pca)
    emptyfig = sns.scatterplot()
    emptyfig.figure.savefig(snakemake.output.pca_figure)
    emptyfig.figure.savefig(snakemake.output.pc1_pc2_figure)
    emptyfig.figure.savefig(snakemake.output.pc1_pc3_figure)
    pd.DataFrame().to_csv(snakemake.output.explained_variance)
