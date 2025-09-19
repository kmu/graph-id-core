import argparse
import glob

# https://github.com/hackingmaterials/matminer_examples/blob/main/matminer_examples/machine_learning-nb/voronoi-ward-prb-2017.ipynb
import json
import os

import numpy as np
import pandas as pd
from matminer.featurizers.base import MultipleFeaturizer
from matminer.featurizers.composition import ElementProperty, IonProperty, Stoichiometry, ValenceOrbital
from matminer.featurizers.conversions import ConversionFeaturizer
from matminer.featurizers.structure import (
    ChemicalOrdering,
    MaximumPackingEfficiency,
    SiteStatsFingerprint,
    StructuralHeterogeneity,
    StructureComposition,
)
from pymatgen.core import Structure
import graph_id_cpp
# from graph_id.core.graph_id import GraphIDGenerator
# read json from a file
from monty.serialization import loadfn
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from tqdm import tqdm


def mean_absolute_error(y_true, y_pred):
    return np.mean(np.abs(np.array(y_true) - np.array(y_pred)))


random_seed = 1
cache_flag = False 

featurizer = MultipleFeaturizer(
  [
      SiteStatsFingerprint.from_preset("CoordinationNumber_ward-prb-2017"),
      StructuralHeterogeneity(),
      ChemicalOrdering(),
      MaximumPackingEfficiency(),
      SiteStatsFingerprint.from_preset("LocalPropertyDifference_ward-prb-2017"),
      StructureComposition(Stoichiometry()),
      StructureComposition(ElementProperty.from_preset("magpie")),
      StructureComposition(ValenceOrbital(props=["frac"])),
      StructureComposition(IonProperty(fast=True)),
  ],
)
feat_names = featurizer.feature_labels()

class CifToStructure(ConversionFeaturizer):  # type: ignore
  def __init__(self, target_col_id="cif", overwrite_data=False) -> None:
      self._target_col_id = target_col_id
      self._overwrite_data = overwrite_data
      super().__init__(target_col_id, overwrite_data)

  def citations(self):
      return []

  def featurize(self, string: str) -> list[Structure]:
      s = Structure.from_str(input_string=string, fmt="cif")
      return [s]

  def implementors(self) -> list[str]:
      return ["Koki Muraoka"]

pickle_path = "mp.2019.04.01.pickle"

def get_id_from_str(structure_str):
    structure = Structure.from_str(input_string=structure_str, fmt="cif")
    gid_gen = graph_id_cpp.GraphIDGenerator()
    gid = gid_gen.get_id(structure)

    return gid


if not os.path.exists(pickle_path):
  # Download from https://figshare.com/articles/dataset/Graphs_of_Materials_Project_20190401/8097992
  json_path = "mp.2019.04.01.json"
  data = loadfn(json_path)
  d1 = pd.DataFrame(data)
  print("Calculating Graph IDs...")
  d1["graph_id"] = d1["structure"].apply(get_id_from_str)

  c2s = CifToStructure()
  d2 = c2s.featurize_dataframe(d1, col_id="structure")

  d2["structure"] = d2["cif"]
  del d2["cif"]

  d3 = featurizer.featurize_dataframe(d2, "structure", ignore_errors=None)
#   d3["graph_id"] = gid_gen.get_many_ids(d3["structure"].values, parallel=True)

  print("Saving to pickle...")
  d3.to_pickle(pickle_path)
else:
  print("Reading from pickle...")
  d3 = pd.read_pickle(pickle_path)

if cache_flag:
    df = pd.read_csv("train_energy_df.csv")
    print(df)

else:
    if "descriptors" in d3:
        # if len(d3['descriptors']) == len(d3['feat_names']):
        df = pd.DataFrame([d3["descriptors"]], columns=d3["feat_names"])
        df["graph_id"] = d3["graph_id"]
        df["path"] = path
    else:
        df = d3.copy()
    df.to_csv("train_energy_df.csv", index=False)

    # if "mpid" not in df.columns:
    #     df["mpid"] = df.path.str.split("/").str[1]
    #     # print(df)

    # delete the column named 'Unnamed: 0' if exists in df
    if "Unnamed: 0" in df.columns:
        del df["Unnamed: 0"]

# 50% split
train_df, test_df = train_test_split(df, train_size=0.5, random_state=random_seed)
# train_df, test_df = train_test_split(d3, train_size=desired_train_rows, random_state=random_seed)

X_train = train_df[featurizer.feature_labels()].values
X_test = test_df[featurizer.feature_labels()].values

y_train = train_df["formation_energy_per_atom"].values
y_test = test_df["formation_energy_per_atom"].values

model = Pipeline(
    [
        ("imputer", SimpleImputer()),  # For the failed structures
        ("model", RandomForestRegressor(n_estimators=150, n_jobs=-1)),
    ],
)

model.fit(X_train, y_train)
y_test_pred = model.predict(X_test)
mae = mean_absolute_error(y_test, y_test_pred)
print(f"MAE of whole test data: {mae}")

obs_df = test_df[test_df.graph_id.isin(train_df.graph_id)]
novel_df = test_df[~test_df.graph_id.isin(train_df.graph_id)]

X_obs = obs_df[featurizer.feature_labels()].values
X_novel = novel_df[featurizer.feature_labels()].values

y_obs = obs_df["formation_energy_per_atom"].values
y_novel = novel_df["formation_energy_per_atom"].values

y_obs_pred = model.predict(X_obs)
obs_mae = mean_absolute_error(y_obs, y_obs_pred)

y_novel_pred = model.predict(X_novel)
novel_mae = mean_absolute_error(y_novel, y_novel_pred)

print(f"MAE of leaked test data: {obs_mae}")
print(f"MAE of unleaked test data: {novel_mae}")

obs_df = pd.DataFrame({"actual": y_obs, "pred": y_obs_pred})
obs_df["kind"] = "observed"

novel_df = pd.DataFrame({"actual": y_novel_pred, "pred": y_novel})
novel_df["kind"] = "novel"

leaked_unleaked_result_df = pd.concat([obs_df, novel_df])
leaked_unleaked_result_df.to_csv("leaked_unleaked_result_df.csv", index=False)

y_train_pred = model.predict(X_train)
train_mae = mean_absolute_error(y_train, y_train_pred)
print(f"MAE of train data: {train_mae}")

train_df = pd.DataFrame({"actual": y_train_pred, "pred": y_train})
train_df["kind"] = "train"

test_train_df = pd.concat([leaked_and_unleaked_result_df, train_df])
test_train_df.to_csv("test_train_df.csv", index=False)
