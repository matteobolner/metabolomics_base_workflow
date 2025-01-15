import pandas as pd

long_df = []
for long_counts in snakemake.input.long_counts:
    long_df.append(pd.read_table(long_counts))

long_df = pd.concat(long_df).reset_index(drop=True)

summary = long_df.groupby(by="metabolite")["decision"].value_counts()
summary = pd.DataFrame(summary)
summary = summary.reset_index()
summary = summary.pivot(index="metabolite", columns=["decision"])
summary.columns = summary.columns.droplevel()
summary = summary.fillna(0).astype(int)
summary.columns.name = None
summary = summary.reset_index()

summary.to_csv(snakemake.output.summary, index=False, sep="\t")
long_df.to_csv(snakemake.output.long_df, index=False, sep="\t")
