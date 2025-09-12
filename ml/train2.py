# https://github.com/hackingmaterials/matminer_examples/blob/main/matminer_examples/machine_learning-nb/voronoi-ward-prb-2017.ipynb
import argparse
import os

import numpy as np
import pandas as pd
from chemsys.matminer.featurizers.conversion import CifToStructure
from graphid.core.graphid import GraphID
from matminer.featurizers.base import MultipleFeaturizer
from matminer.featurizers.composition import ElementProperty, IonProperty, Stoichiometry, ValenceOrbital
from matminer.featurizers.structure import (
    ChemicalOrdering,
    MaximumPackingEfficiency,
    SiteStatsFingerprint,
    StructuralHeterogeneity,
    StructureComposition,
)

# read json from a file
from monty.serialization import loadfn
from pymatgen.core import Structure
from senko.core.senkotask import SenkotaskBase
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline

from chemsys_matminer.featurizers.conversion import CifToStructure


# この記述子で、MEGNetのデータを用いる。


def mean_absolute_error(y_true, y_pred):
    return np.mean(np.abs(np.array(y_true) - np.array(y_pred)))


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", type=int, default=2019)
parser.add_argument("--debug", action="store_true")
parser.add_argument("-s", "--seed", type=int, default=0)

parser.add_argument("--senko", action="store_true")


args = parser.parse_args()


random_seed = args.seed

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

if args.dataset == 2018:
    json_path = "/home/mrok/data/mp.2018.6.1.json"

if args.dataset == 2019:
    # pickle_path = "mp2019.pickle"
    json_path = "mp.2019.04.01.json"

pickle_path = f"mp{args.dataset}-debug{args.debug}-senko{args.senko}.pickle"


class GraphIDSenkoTask(SenkotaskBase):
    def run_task(self, inputs):
        gid = GraphID()

        cif = inputs["cif"]

        structure = Structure.from_str(cif, fmt="cif")

        try:
            myid = gid.get_id(structure)
            descriptors = featurizer.featurize(structure)

            return {"graph_id": myid, "descriptors": descriptors, "feat_names": featurizer.feature_labels()}

        except Exception as e:
            return {"error": e}


if not os.path.exists(pickle_path):
    print(f"No {pickle_path}")
    print(f"Reading {json_path}")
    data = loadfn(json_path)
    print(f"The number of structures is {len(data)}")
    print("formation_energy_per_atom is {}".format(data[0]["formation_energy_per_atom"]))

    d1 = pd.DataFrame(data)
    if args.debug:
        d1 = d1.head(200)

    if args.senko:
        mpids = []
        senkoargs = []

        for entry in data:
            mpids.append(entry["material_id"])
            senkoargs.append({"cif": entry["structure"]})

        t = GraphIDSenkoTask()
        t.queue(names=mpids, args=senkoargs)
        del senkoargs
        t.launch()
        t.assimilate()

        d = loadfn(t.get_assimilated_path())
        d3 = pd.DataFrame(d)

    # dto = DictToObject(target_col_id='structure', overwrite_data=True)
    # d2 = dto.featurize_dataframe(d1, 'structure')
    else:
        c2s = CifToStructure()
        d2 = c2s.featurize_dataframe(d1, col_id="structure")

        d2["structure"] = d2["cif"]
        del d2["cif"]

        gidgen = GraphID()

        print("Calculating Graph IDs...")
        # print(d2['structure'])
        d3 = featurizer.featurize_dataframe(d2, "structure", ignore_errors=True)
        d3["graph_id"] = gidgen.get_many_ids(d3["structure"].values, parallel=True)

    print("Saving to pickle...")
    d3.to_pickle(pickle_path)
else:
    print("Reading from pickle...")
    d3 = pd.read_pickle(pickle_path)

# 50% split
train_df, test_df = train_test_split(d3, train_size=0.5, random_state=random_seed)
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
print(mae)


y_novel_pred = model.predict(X_novel)
mae = mean_absolute_error(y_novel, y_novel_pred)
print(mae)

# data = load_dataset("flla")
# print('Loaded {} entries'.format(len(data)))
# 3938

# dto = DictToObject(target_col_id='structure', overwrite_data=True)
# data = dto.featurize_dataframe(data, 'structure')
