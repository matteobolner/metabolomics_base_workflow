import pandas as pd

cv = pd.read_table(snakemake.input.cv)
long_df = pd.read_table(snakemake.input.long_df)
summary = pd.read_table(snakemake.input.summary)

cv_df = pd.DataFrame()
cv_df["metabolite"] = cv["metabolite"].unique()

cv_long_df = []

for name, group in cv.groupby(by=["seed", "imputation_cycle", "seed_boruta"]):
    print(name)
    name = f"mice_{name[0]}_{name[1]}_boruta_{name[2]}_CV"
    counts = {}
    for metabolite, metgroup in group.groupby(by="metabolite"):
        counts[metabolite] = len(metgroup[metgroup["decision"] == "Confirmed"])
    cv_df[name] = cv_df["metabolite"].apply(lambda x: counts[x])

cv_counts = pd.DataFrame(
    cv.groupby(by=["metabolite", "seed", "imputation_cycle", "seed_boruta"])[
        "decision"
    ].value_counts()
).reset_index()
cv_counts = cv_counts[cv_counts["decision"] == "Confirmed"]
cv_counts = cv_counts.drop(columns="decision")
cv_counts = cv_counts.rename(columns={"count": "CV_confirmed"})

long_df = long_df.merge(
    cv_counts,
    on=["metabolite", "seed", "imputation_cycle", "seed_boruta"],
    how="left",
)
long_df["CV_confirmed"] = long_df["CV_confirmed"].fillna(0).astype(int)
long_df.to_csv(snakemake.output.long_df, index=False, sep="\t")

summary_counts_cv = pd.DataFrame(long_df.groupby(by="metabolite")["CV_confirmed"].sum())
summary_counts_cv = summary_counts_cv.reset_index()
summary = summary.merge(summary_counts_cv, on="metabolite")
summary.to_csv(snakemake.output.summary, index=False, sep="\t")
