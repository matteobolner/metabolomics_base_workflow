import pandas as pd
from metabotk import MetaboTK

dataset = MetaboTK().io.from_excel(
    # "data/normalization/normalized.xlsx",
    snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
)


dataset.data = raw_dataset.data

# group_column='high_or_low'
group_column = snakemake.config["group_column"]

mw_tests = pd.read_table(snakemake.input.mann_whitney).set_index(
    snakemake.config["metabolite_id_column"]
)
mw_tests = mw_tests[["mann_whitney_u_statistic", "mann_whitney_p_value"]].copy()
mw_tests.columns = pd.MultiIndex.from_tuples(
    [
        ("Comparison", "Mann-Whitney U test statistic"),
        ("Comparison", "Mann-Whitney U test p-value"),
    ]
)
mw_tests.index.name = "Metabolite ID"


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


# dataset.sample_metadata[group_column]=dataset.sample_metadata[group_column]+"_FCR"

dataset.sample_metadata[group_column] = dataset.sample_metadata[group_column]

split_data = dataset.split_by_sample_column(group_column)
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

stats = stats.round(2)

stats = stats.merge(mw_tests, left_index=True, right_index=True)

stats = stats.reset_index()

stats[("Comparison", "AUC")] = stats["Metabolite ID"].apply(
    lambda x: dataset.chemical_annotation["AUC"].to_dict()[x]
)
stats[("Comparison", "Confirmed")] = stats["Metabolite ID"].apply(
    lambda x: dataset.chemical_annotation["Confirmed"].to_dict()[x]
)
stats[("Comparison", "CV Confirmed")] = stats["Metabolite ID"].apply(
    lambda x: dataset.chemical_annotation["CV_confirmed"].to_dict()[x]
)

metabolite_info = dataset.chemical_annotation[
    ["PLOT_NAME", "SUPER_PATHWAY", "SUB_PATHWAY"]
]
metabolite_info = metabolite_info.sort_values(
    by=["SUPER_PATHWAY", "SUB_PATHWAY", "CHEM_ID"]
)
metabolite_info.index.name = "Metabolite ID"
metabolite_info.columns = pd.MultiIndex.from_tuples(
    [("Metabolite", i) for i in ["Metabolite name", "Super pathway", "Sub pathway"]]
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

stats = stats.rename(
    columns={
        "Comparison": f"{snakemake.wildcards.trait}_Comparison",
    }
)


stats.to_csv(snakemake.output.stats, index=False, sep="\t")
