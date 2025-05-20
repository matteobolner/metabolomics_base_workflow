import pandas as pd
from metabotk import MetaboTK
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import mannwhitneyu

dataset = MetaboTK().io.from_excel(
    snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
)

split_data = dataset.ops.split(by="samples", columns=snakemake.config["group_column"])
data_values = list(split_data.values())

mw_tests = {}

for i in dataset.metabolites:
    mw_tests[i] = mannwhitneyu(data_values[0].data[i], data_values[1].data[i])

mw_tests = pd.DataFrame(mw_tests).transpose()
mw_tests.columns = ["mann_whitney_u_statistic", "mann_whitney_p_value"]
mw_tests["mann_whitney_p_value_FDR_corrected"] = fdrcorrection(
    mw_tests["mann_whitney_p_value"]
)[1]
mw_tests["mann_whitney_u_statistic"] = mw_tests["mann_whitney_u_statistic"].round(2)
mw_tests.index.name = snakemake.config["metabolite_id_column"]
mw_tests.to_csv(snakemake.output.mann_whitney, index=True, sep="\t")
