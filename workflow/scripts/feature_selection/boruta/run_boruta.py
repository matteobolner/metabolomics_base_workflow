from metabotk import MetaboTK
import pandas as pd

dataset = MetaboTK().io.from_excel(
    snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
)

ranking = dataset.fs.boruta(
    y_column=snakemake.config["group_column"],
    kind="classifier",
    threads=1,
    random_state=int(snakemake.wildcards.borut_seed),
    max_depth=None,
    class_weight="balanced",
    n_estimators="auto",
    alpha=0.01,
    max_iterations=1000,
    output_dir=None,
)
ranking.to_csv(snakemake.output.ranking, index=False, sep="\t")

ranking["decision"] = "Rejected"
ranking["decision"] = ranking["decision"].where(ranking["rank"] > 2, "Tentative")
ranking["decision"] = ranking["decision"].where(ranking["rank"] > 1, "Confirmed")

for k, v in dict(snakemake.wildcards).items():
    ranking[k] = v

ranking = ranking.drop(columns=["rank"])
ranking.to_csv(snakemake.output.long_df, index=False, sep="\t")
