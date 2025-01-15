import pandas as pd
from metabotk import MetaboTK


configfile: "config/config.yaml"


def get_mice_covariates():
    return ",".join(config["mice_imputation_covariates"])


mice_seeds = config["mice_seeds"]
imputation_cycles = config["imputation_cycles"]
feature_selection_seeds = config["feature_selection_seeds"]


wildcard_constraints:
    mice_seed="[^_]+",
    fold="[^_]+",
    seed_boruta="[^_/]+",
    with_or_without="[^_/]+",
    date="[^_/]+",
