"""Capture the original NumPy AVX-512 GroupKFold partitions without fitting models."""

import argparse
import hashlib
import itertools
import json
import platform
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.model_selection import GroupKFold

from regional_code_paper.analysis.region_definition_sensitivity import (
    construct_definition,
    load_validated_sites,
)
from regional_code_paper.analysis.regional_code import build_tiles, load_inputs

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--archived-results", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
config_path = Path("config/config.yaml")
partitions = {}


def capture(groups, n_splits, name):
    unique, indices, counts = np.unique(groups, return_inverse=True, return_counts=True)
    cohort = list(zip(unique.tolist(), counts.tolist(), strict=True))
    key = hashlib.sha256(
        json.dumps([n_splits, cohort], separators=(",", ":")).encode()
    ).hexdigest()
    folds = np.full(len(groups), -1, dtype=int)
    splits = list(GroupKFold(n_splits).split(np.zeros(len(groups)), groups=groups))
    for fold, (_, test) in enumerate(splits):
        folds[test] = fold
    group_folds = []
    for index in range(len(unique)):
        values = np.unique(folds[indices == index])
        assert len(values) == 1
        group_folds.append(int(values[0]))
    item = {
        "cohort": [
            [protein, count, fold]
            for (protein, count), fold in zip(cohort, group_folds, strict=True)
        ],
        "n_splits": n_splits,
        "names": [name],
    }
    if key in partitions:
        assert partitions[key]["cohort"] == item["cohort"]
        partitions[key]["names"].append(name)
    else:
        partitions[key] = item
    print(name, len(groups), len(unique), key, flush=True)
    return splits, folds


with np.load("results/work/scanner/features.npz") as data:
    groups = data["groups"]
outer, folds = capture(groups, 5, "scanner.outer")
expected = pd.read_csv(args.archived_results / "analysis/scanner/outer_predictions.csv")
assert np.array_equal(folds, expected.fold), (
    "Reference CPU does not recover archived scanner folds"
)
for fold, (train, _) in enumerate(outer):
    capture(groups[train], 3, f"scanner.inner.{fold}")
with np.load("results/work/regional_code/prepared/tile_features.npz") as data:
    capture(data["groups"], 5, "regional_code")
config, _, site_map, _, _ = load_validated_sites(config_path)
settings = config.values["analysis"]
for gap, minimum in itertools.product(
    settings["region_sensitivity"]["core_gaps"], settings["region_sensitivity"]["minimum_sites"]
):
    regions, sites = construct_definition(
        site_map, core_gap=gap, minimum_sites=minimum, final_gap=settings["final_gap"]
    )
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        regions.to_csv(tmp / "regions.csv", index=False)
        sites.to_csv(tmp / "sites.csv", index=False)
        tiles = build_tiles(load_inputs(config_path, tmp / "regions.csv", tmp / "sites.csv"))
    capture(tiles.accession.to_numpy(str), 5, f"region_sensitivity.core{gap}.min{minimum}")
output = {
    "schema_version": 1,
    "source_commit": "c03b09bcb170767f8319bf08df18e1bedbed979e",
    "reference_runtime": {
        "numpy": np.__version__,
        "scikit_learn": sklearn.__version__,
        "platform": platform.platform(),
        "cpu": platform.processor(),
    },
    "partitions": partitions,
}
args.output.write_text(json.dumps(output, separators=(",", ":")) + "\n")
