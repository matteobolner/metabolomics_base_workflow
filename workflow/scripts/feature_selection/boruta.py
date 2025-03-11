from metabotk import MetaboTK
from tqdm import tqdm
from multiprocessing import Pool
import pandas as pd


# dataset = MetaboTK().io.from_excel(
#    "../../projects/anas_breed_comparisons/LW_LA_DU/datasets/Landrace_vs_Duroc/data/residuals/seed_1000/imputation_1.xlsx",
#    sample_id_column="PARENT_SAMPLE_NAME",
# )

ds = MetaboTK().io.from_excel(
    snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
    sample_metadata_sheet=snakemake.config["sample_metadata_sheet"],
    chemical_annotation_sheet=snakemake.config["chemical_annotation_sheet"],
    data_sheet=snakemake.config["data_sheet"],
)

ds.ops.subset(
    what="samples",
    ids=ds.sample_metadata[
        ds.sample_metadata[snakemake.config["group_column"]].isin(
            [snakemake.config["groups"]]
        )
    ].index,
)


def run_boruta(dataset, random_state):
    print(f"Running boruta with random state {random_state}", flush=True)
    ranking = dataset.fs.boruta(
        y_column=snakemake.config["group_column"],
        kind=snakemake.config["boruta_method"],
        threads=1,
        random_state=random_state,
        max_depth=None,
        class_weight="balanced",
        n_estimators="auto",
        alpha=0.01,
        max_iterations=1000,
        output_dir=None,
    )
    ranking["decision"] = "Rejected"
    ranking["decision"] = ranking["decision"].where(ranking["rank"] > 2, "Tentative")
    ranking["decision"] = ranking["decision"].where(ranking["rank"] > 1, "Confirmed")
    ranking["random_state"] = random_state
    print(f"Boruta run finished: \n {ranking['decision'].value_counts()}", flush=True)
    return ranking


def merge_boruta_results(results):
    long_df = pd.concat(results).reset_index(drop=True)
    summary = long_df.groupby(by="metabolite")["decision"].value_counts()
    summary = pd.DataFrame(summary)
    summary = summary.reset_index()
    summary = summary.pivot(index="metabolite", columns=["decision"])
    summary.columns = summary.columns.droplevel()
    summary = summary.fillna(0).astype(int)
    summary.columns.name = None
    summary = summary.reset_index()
    return long_df, summary


random_states = snakemake.config["feature_selection_seeds"]
# random_states = [10, 20, 30, 40, 50]
thread_number = snakemake.threads
# thread_number=1
with Pool(
    processes=thread_number,
) as p:
    results = list(
        p.starmap(
            run_boruta,
            [(ds, i) for i in random_states],
        )
    )
    p.close()
    p.join()


long_df, summary = merge_boruta_results(results)

if snakemake.config["CV_method"] == "stratified_kfold":
    cv_datasets = ds.fs.stratified_kfold(
        n_splits=snakemake.config["CV_folds_number"],
        stratification_column=snakemake.config["CV_stratification_column"],
    )["training_set"].values()
elif snakemake.config["CV_method"] == "LOO":
    cv_datasets = [ds.ops.drop(what="samples", ids=[i]) for i in ds.samples]
else:
    raise ValueError("Invalid CV method")

cv_combinations = []
for i in cv_datasets:
    for j in random_states:
        cv_combinations.append((i, j))

with Pool(
    processes=len(random_states),
) as p:
    cv_results = list(p.starmap(run_boruta, cv_combinations))
    p.close()
    p.join()

long_df_CV, summary_CV = merge_boruta_results(cv_results)

cv_df = pd.DataFrame()
cv_df["metabolite"] = long_df_CV["metabolite"].unique()

cv_long_df = []

for name, group in long_df_CV.groupby(by="random_state"):
    counts = {}
    for metabolite, metgroup in group.groupby(by="metabolite"):
        counts[metabolite] = len(metgroup[metgroup["decision"] == "Confirmed"])
    cv_df[name] = cv_df["metabolite"].apply(lambda x: counts[x])

cv_counts = pd.DataFrame(
    long_df_CV.groupby(by=["metabolite", "random_state"])["decision"].value_counts()
).reset_index()
cv_counts = cv_counts[cv_counts["decision"] == "Confirmed"]
cv_counts = cv_counts.drop(columns="decision")
cv_counts = cv_counts.rename(columns={"count": "CV_confirmed"})

long_df = long_df.merge(
    cv_counts,
    on=["metabolite", "random_state"],
    how="left",
)
long_df["CV_confirmed"] = long_df["CV_confirmed"].fillna(0).astype(int)
summary_counts_cv = pd.DataFrame(long_df.groupby(by="metabolite")["CV_confirmed"].sum())
summary_counts_cv = summary_counts_cv.reset_index()
summary = summary.merge(summary_counts_cv, on="metabolite")

long_df.to_csv(snakemake.output.long_df, index=False, sep="\t")
summary.to_csv(snakemake.output.summary, index=False, sep="\t")
