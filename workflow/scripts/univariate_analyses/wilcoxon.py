import pandas as pd
from metabotk import MetaboTK
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import wilcoxon

dataset = MetaboTK().io.from_excel(
    snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
)

dataset = dataset.ops.sort(on="samples", by=snakemake.config["paired_sample_column"])

split_data = dataset.ops.split(by="samples", columns=snakemake.config["group_column"])
data_values = list(split_data.values())

w_tests = {}

for i in dataset.metabolites:
    w_tests[i] = wilcoxon(x=data_values[0].data[i], y=data_values[1].data[i])


w_tests = pd.DataFrame(w_tests).transpose()
w_tests.columns = ["wilcoxon_statistic", "wilcoxon_p_value"]
w_tests["wilcoxon_p_value_FDR_corrected"] = fdrcorrection(w_tests["wilcoxon_p_value"])[
    1
]
w_tests.index.name = snakemake.config["metabolite_id_column"]
w_tests.to_csv(snakemake.output.wilcoxon, index=True, sep="\t")
