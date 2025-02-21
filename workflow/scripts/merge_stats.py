import pandas as pd
from metabotk import MetaboTK

dataset = MetaboTK().io.from_excel(
    # "data/normalization/normalized.xlsx",
    snakemake.input.datasets[0],
    sample_id_column="SampleIdentification",
    metabolite_id_column="Metabolite",
)

# boruta = pd.read_table("tables/feature_selection/boruta/summary.tsv").set_index(
#    "metabolite"
# )
boruta = pd.read_table(snakemake.input.boruta).set_index("metabolite")
# dataset.import_excel(file_path=snakemake.input.dataset)

# raw_dataset.import_excel("data/extremes_raw/FCR/seed_1000_imputation_1.xlsx")
# raw_dataset.import_excel(snakemake.input.raw_dataset)
# dataset.data = raw_dataset.data

# group_column = "breed"
group_column = snakemake.config["group_column"]


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


split_data = dataset.ops.split(by="samples", columns=group_column)

split_data_stats = {k: stats_selected(v.data, k) for k, v in split_data.items()}
data_values = list(split_data.values())
stat_values = list(split_data_stats.values())
stat_keys = list(split_data.keys())

stats = stat_values[0].merge(stat_values[1], left_index=True, right_index=True)

for k, v in split_data_stats.items():
    split_data_stats[k] = split_data_stats[k].reset_index(drop=True)

stats[("Comparison", "Mean delta %")] = (
    (stat_values[0][(stat_keys[0], "mean")] - stat_values[1][(stat_keys[1], "mean")])
    / stat_values[0][(stat_keys[0], "mean")]
) * 100


stats[("Comparison", "Median delta %")] = (
    (
        stat_values[0][(stat_keys[0], "median")]
        - stat_values[1][(stat_keys[1], "median")]
    )
    / stat_values[0][(stat_keys[0], "median")]
) * 100

stats = stats.round(2).reset_index()

stats[("Comparison", "Confirmed")] = stats["Metabolite ID"].apply(
    lambda x: boruta["Confirmed"].to_dict()[x]
)
stats[("Comparison", "CV Confirmed")] = stats["Metabolite ID"].apply(
    lambda x: boruta["CV_confirmed"].to_dict()[x]
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


# stats[("Metabolite", "Super pathway")] = stats[("Metabolite", "Super pathway")].apply(
#    code_from_super_pathway
# )

stats = stats.loc[
    stats[("Comparison", "Mean delta %")].abs().sort_values(ascending=False).index
]

stats = stats.sort_values(by=("Comparison", "Confirmed"), ascending=False)

stats.to_csv(snakemake.output.stats, index=False, sep="\t")
