import pandas as pd
from metabotk import MetaboTK

dataset = MetaboTK().io.from_excel(
    snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
)

boruta = pd.read_table(snakemake.input.boruta).set_index("metabolite")
wilcoxon = pd.read_table(snakemake.input.wilcoxon).set_index(
    snakemake.config["metabolite_id_column"]
)
auc = pd.read_table(snakemake.input.auc).set_index(
    snakemake.config["metabolite_id_column"]
)

paired_sample_column = snakemake.config["paired_sample_column"]
# paired_sample_column = "CLIENT_SAMPLE_ID"
group_column = snakemake.config["group_column"]
# group_column = "TIME"


def stats_selected(data, group):
    stats = data.apply(lambda x: x.describe())
    stats = stats.transpose()
    stats.index.name = "Metabolite ID"
    stats = stats.drop(columns=["count", "25%", "75%"])
    stats = stats.rename(columns={"50%": "median"})
    stats.columns = ["mean", "std", "median", "min", "max"]
    col_tuples = [(group, i) for i in stats.columns]
    stats.columns = pd.MultiIndex.from_tuples(col_tuples)
    return stats


dataset = dataset.ops.sort(on="samples", by=paired_sample_column)
split_data = dataset.ops.split(by="samples", columns=group_column)

split_data = {
    k: v.ops.replace_sample_names_in_data(paired_sample_column)
    for k, v in split_data.items()
}
split_data_stats = {k: stats_selected(v.data, k) for k, v in split_data.items()}
data_values = list(split_data.values())
stat_values = list(split_data_stats.values())
stat_keys = list(split_data.keys())

stats = stat_values[0].merge(stat_values[1], left_index=True, right_index=True)
for k, v in split_data_stats.items():
    split_data_stats[k] = split_data_stats[k].reset_index(drop=True)

stats[("Comparison", "Mean delta %")] = (
    (
        split_data[stat_keys[snakemake.config["group_order"][0]]].data
        - split_data[stat_keys[snakemake.config["group_order"][1]]].data
    ).mean()
) * 100

stats[("Comparison", "Median delta %")] = (
    (
        split_data[stat_keys[snakemake.config["group_order"][0]]].data
        - split_data[stat_keys[snakemake.config["group_order"][1]]].data
    ).median()
) * 100
stats = stats.round(2).reset_index()

stats[("Comparison", "Confirmed")] = stats["Metabolite ID"].apply(
    lambda x: boruta["Confirmed"].to_dict()[x]
)
stats[("Comparison", "CV Confirmed")] = stats["Metabolite ID"].apply(
    lambda x: boruta["CV_confirmed"].to_dict()[x]
)

stats[("Comparison", "Wilcoxon test statistic")] = stats["Metabolite ID"].apply(
    lambda x: wilcoxon["wilcoxon_statistic"].to_dict()[x]
)

stats[("Comparison", "Wilcoxon test p-value")] = stats["Metabolite ID"].apply(
    lambda x: wilcoxon["wilcoxon_p_value"].to_dict()[x]
)
stats[("Comparison", "Wilcoxon test p-value (FDR corrected)")] = stats[
    "Metabolite ID"
].apply(lambda x: wilcoxon["wilcoxon_p_value_FDR_corrected"].to_dict()[x])

stats[("Comparison", "AUC")] = stats["Metabolite ID"].apply(
    lambda x: auc["AUC"].to_dict()[x]
)


metabolite_info = dataset.chemical_annotation
metabolite_info.index.name = "Metabolite ID"
metabolite_info.columns = pd.MultiIndex.from_tuples(
    [("Metabolite", i) for i in metabolite_info.columns]
)
metabolite_info = metabolite_info.fillna("Unknown")

stats = stats.set_index("Metabolite ID")
stats = metabolite_info.merge(stats, left_index=True, right_index=True)
stats = stats.reset_index()
stats.columns = pd.MultiIndex.from_tuples(
    [("Metabolite", "Metabolite ID")] + [i for i in stats.columns[1::]]
)


def code_from_super_pathway(input_string):
    input_string = input_string.replace("and ", "")
    input_string = input_string.split(" ")
    output = [i[0] for i in input_string]
    return "".join(output).upper()


stats[("Metabolite", "Super pathway code")] = stats[
    ("Metabolite", snakemake.config["super_pathway_column"])
].apply(code_from_super_pathway)

stats = stats.loc[
    stats[("Comparison", "Mean delta %")].abs().sort_values(ascending=False).index
]

stats = stats.sort_values(by=("Comparison", "Confirmed"), ascending=False)
stats.to_csv(snakemake.output.stats, index=False, sep="\t")
