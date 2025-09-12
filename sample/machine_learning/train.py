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
import graph_id
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

if not os.path.exists(pickle_path):
  # Download from https://figshare.com/articles/dataset/Graphs_of_Materials_Project_20190401/8097992
  json_path = "mp.2019.04.01.json"
  data = loadfn(json_path)
  d1 = pd.DataFrame(data)

  c2s = CifToStructure()
  d2 = c2s.featurize_dataframe(d1, col_id="structure")

  d2["structure"] = d2["cif"]
  del d2["cif"]

  gid_gen = graph_id.GraphIDGenerator()
  print("Calculating Graph IDs...")

  d3 = featurizer.featurize_dataframe(d2, "structure", ignore_errors=None)
  d3["graph_id"] = gid_gen.get_many_ids(d3["structure"].values, parallel=True)

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
print(mae)

train_df.graph_id
test_df

obs_df = test_df[test_df.graph_id.isin(train_df.graph_id)]
novel_df = test_df[test_df.graph_id.isin(train_df.graph_id) is False]

X_obs = obs_df[featurizer.feature_labels()].values
X_novel = novel_df[featurizer.feature_labels()].values

y_obs = obs_df["formation_energy_per_atom"].values
y_novel = novel_df["formation_energy_per_atom"].values

y_obs_pred = model.predict(X_obs)
mae = mean_absolute_error(y_obs, y_obs_pred)
print("OBS")
print(mae)

df4 = pd.DataFrame({"actual": y_obs, "pred": y_obs_pred})
df4["kind"] = "observed"

y_novel_pred = model.predict(X_novel)
mae = mean_absolute_error(y_novel, y_novel_pred)
print("ACT")
print(mae)

df5 = pd.DataFrame({"actual": y_novel_pred, "pred": y_novel})
df5["kind"] = "novel"

df6 = pd.concat([df4, df5])
df6.to_csv("train3_df6.csv", index=False)

y_train_pred = model.predict(X_train)
mae = mean_absolute_error(y_train, y_train_pred)
print("ACT")
print(mae)

df7 = pd.DataFrame({"actual": y_train_pred, "pred": y_train})
df7["kind"] = "train"

df8 = pd.concat([df6, df7])
df8.to_csv("train3_df8.csv", index=False)
