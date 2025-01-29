from snakemake.script import snakemake
from metabotk import MetaboTK

ds = MetaboTK().io.from_excel(
    file_path=snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
)
no_full_missing = ds.stats.remove_missing(on="metabolites", threshold=0.999999)
ds_all_missing = ds.ops.drop(what="metabolites", ids=list(no_full_missing.columns))
ds_all_missing.chemical_annotation.to_csv(
    snakemake.output.empty_metabolites, index=True, sep="\t"
)
ds = ds.ops.subset(what="metabolites", ids=list(no_full_missing.columns))
metabolite_stats = ds.stats.metabolite_stats()
sample_stats = ds.stats.sample_stats()

sample_stats.to_csv(snakemake.output.sample_stats, index=True, sep="\t")
metabolite_stats.to_csv(snakemake.output.metabolite_stats, index=True, sep="\t")


pca = ds.viz.plot_pca_grid(
    hue=snakemake.config["group_column"], savepath=snakemake.output.pca_plot
)
