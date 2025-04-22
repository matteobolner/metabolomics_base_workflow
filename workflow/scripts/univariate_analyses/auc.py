from metabotk import MetaboTK
from sklearn import metrics

dataset = MetaboTK().io.from_excel(
    snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
)


def get_data_direction(data, metabolite, groupcol):
    """
    Determine the direction of difference between two groups based on their medians.

    Parameters:
    data (DataFrame): Input data
    metabolite (str): Name of the metabolite column
    groupcol (str): Name of the grouping column

    Returns:
    str: ">" if first group's median is lower than second group's, "<" if opposite
    """
    # Get unique groups and ensure there are exactly two
    groups = data[groupcol].unique()
    if len(groups) != 2:
        raise ValueError(f"Expected exactly 2 groups, got {len(groups)}: {groups}")

    # Sort groups to ensure consistent ordering
    group1, group2 = sorted(groups)

    # Calculate medians for each group
    medians = {}
    for name, group in data.groupby(by=groupcol):
        medians[name] = group[metabolite].median()

    if medians[group1] < medians[group2]:
        return ">"
    elif medians[group1] > medians[group2]:
        return "<"
    else:
        raise ValueError(f"MEDIANS: {medians}")


def get_roc_curve(data, metabolite, groupcol):
    """
    Calculate ROC curve for metabolite data between two groups.

    Parameters:
    data (DataFrame): Input data
    metabolite (str): Name of the metabolite column
    groupcol (str): Name of the grouping column

    Returns:
    tuple: (direction, false positive rate, true positive rate, thresholds)
    """
    # Get unique groups in sorted order
    groups = sorted(data[groupcol].unique())
    if len(groups) != 2:
        raise ValueError(f"Expected exactly 2 groups, got {len(groups)}: {groups}")

    group1, group2 = groups

    # Get direction using the modified function
    direction = get_data_direction(data, metabolite, groupcol)

    # Create binary classification based on direction
    if direction == ">":
        y = data[groupcol].replace({group1: 0, group2: 1})
    elif direction == "<":
        y = data[groupcol].replace({group1: 1, group2: 0})

    scores = data[metabolite]
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=1)
    return direction, fpr, tpr, thresholds


directions = {}
fprs = {}
tprs = {}
thresholds = {}
roc_aucs = {}
tempdata = dataset.ops.merge_sample_metadata_data()
for metabolite in dataset.metabolites:
    (
        directions[metabolite],
        fprs[metabolite],
        tprs[metabolite],
        thresholds[metabolite],
    ) = get_roc_curve(tempdata, metabolite, snakemake.config["group_column"])
    roc_aucs[metabolite] = metrics.auc(fprs[metabolite], tprs[metabolite])

dataset.chemical_annotation["median_direction"] = (
    dataset.chemical_annotation.reset_index()["CHEM_ID"]
    .apply(lambda x: directions[x])
    .tolist()
)

dataset.chemical_annotation["AUC"] = (
    dataset.chemical_annotation.reset_index()["CHEM_ID"]
    .apply(lambda x: roc_aucs[x])
    .tolist()
)

dataset.chemical_annotation[["median_direction", "AUC"]].to_csv(
    snakemake.output.auc, index=True, sep="\t"
)
