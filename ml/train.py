import argparse
import os
import random

import numpy as np
import pandas as pd
from graphid.core.graphid import GraphID
from matminer.datasets.convenience_loaders import load_elastic_tensor
from matminer.featurizers.composition import ElementProperty, OxidationStates
from matminer.featurizers.conversions import CompositionToOxidComposition, StrToComposition
from matminer.featurizers.structure import DensityFeatures
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split


def mean_absolute_error(y_true, y_pred):
    return np.mean(np.abs(np.array(y_true) - np.array(y_pred)))


def train_and_test(train_df, test_df, dir_name):
    print(len(train_df))
    train_df.to_csv(f"{dir_name}/train.csv")
    test_df.to_csv(f"{dir_name}/test.csv")
    # y_train = train_df['formation_energy_per_atom']
    # y_test = test_df['formation_energy_per_atom']
    # # y = d3['formation_energy_per_atom']
    # X_train = train_df['structure']
    # X_test = test_df['structure']

    y_train = train_df["K_VRH"].values
    excluded = [
        "G_VRH",
        "K_VRH",
        "elastic_anisotropy",
        "formula",
        "material_id",
        "poisson_ratio",
        "structure",
        "composition",
        "composition_oxid",
    ]
    y_test = test_df["K_VRH"].values
    excluded = [
        "G_VRH",
        "K_VRH",
        "elastic_anisotropy",
        "formula",
        "material_id",
        "poisson_ratio",
        "structure",
        "composition",
        "composition_oxid",
        "graph_id",
    ]

    X_train = train_df.drop(excluded, axis=1)
    X_test = test_df.drop(excluded, axis=1)
    # print("There are {} possible descriptors:\n\n{}".format(X.shape[1], X.columns.values))

    rf = RandomForestRegressor(n_estimators=50, random_state=1)

    rf.fit(X_train, y_train)

    y_test_pred = rf.predict(X_test)

    mae = mean_absolute_error(y_test_pred, y_test)
    # print(mae)
    return mae
    # print("training R2 = " + str(round(rf.score(X, y), 3)))
    # print("training RMSE = %.3f" % np.sqrt(mean_squared_error(y_true=y, y_pred=rf.predict(X))))

    # # compute cross validation scores for random forest model
    # crossvalidation = KFold(n_splits=10, shuffle=True, random_state=1)
    # r2_scores = cross_val_score(rf, X, y, scoring="r2", cv=crossvalidation, n_jobs=-1)
    # scores = cross_val_score(rf, X, y, scoring="neg_mean_squared_error", cv=crossvalidation, n_jobs=-1)
    # rmse_scores = [np.sqrt(abs(s)) for s in scores]

    # print("Cross-validation results:")
    # print("Folds: %i, mean R2: %.3f" % (len(scores), np.mean(np.abs(r2_scores))))
    # print("Folds: %i, mean RMSE: %.3f" % (len(scores), np.mean(np.abs(rmse_scores))))


def split_ids_by_desired_rows(df, unique_ids, desired_train_rows):
    id_counts = df["graph_id"].value_counts().loc[unique_ids].sort_values(ascending=False)
    train_ids = []
    test_ids = []
    current_train_rows = 0

    for id, count in id_counts.items():
        if current_train_rows + count <= desired_train_rows:
            train_ids.append(id)
            current_train_rows += count
        else:
            test_ids.append(id)

        if current_train_rows == desired_train_rows:
            break

    # Add the remaining IDs to the test_ids list
    remaining_ids = set(id_counts.index) - set(train_ids)
    test_ids.extend(list(remaining_ids))

    return train_ids, test_ids


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seed", type=int, default=0)
parser.add_argument("-t", "--train", type=int, default=600)
parser.add_argument("-d", "--dedupe", action="store_true")

args = parser.parse_args()

print(args)


dirname = f"results/s_{args.seed}-r_{args.dedupe}-t_{args.train}"
random_seed = args.seed
os.makedirs(dirname, exist_ok=True)


if not os.path.exists("df.csv"):
    df = load_elastic_tensor()  # loads dataset in a pandas DataFrame object

    unwanted_columns = [
        "volume",
        "nsites",
        "compliance_tensor",
        "elastic_tensor",
        "elastic_tensor_original",
        "K_Voigt",
        "G_Voigt",
        "K_Reuss",
        "G_Reuss",
    ]
    df = df.drop(unwanted_columns, axis=1)

    df.head()

    print(df.describe())
    df = StrToComposition().featurize_dataframe(df, "formula")
    df.head()

    ep_feat = ElementProperty.from_preset(preset_name="magpie")
    df = ep_feat.featurize_dataframe(df, col_id="composition")  # input the "composition" column to the featurizer
    df.head()

    df = CompositionToOxidComposition().featurize_dataframe(df, "composition")

    os_feat = OxidationStates()
    df = os_feat.featurize_dataframe(df, "composition_oxid")
    df.head()

    df_feat = DensityFeatures()
    df = df_feat.featurize_dataframe(df, "structure")  # input the structure column to the featurizer

    gidgen = GraphID()
    df["graph_id"] = gidgen.get_many_ids(df["structure"].values, parallel=True)

    df.to_csv("df.csv", index=False)
    # df_feat.feature_labels()

else:
    df = pd.read_csv("df.csv")

desired_train_rows = args.train
# random case
# X_train_e, X_test_e, y_train_e, y_test_e = train_test_split(
# X, y, train_size=desired_train_rows, random_state=random_seed
# )

for random_seed in range(20):
    train_df, test_df = train_test_split(df, train_size=desired_train_rows, random_state=random_seed)

    mae_r = train_and_test(train_df, test_df, f"{dirname}")

    unique_ids = df["graph_id"].unique()

    print(len(df))
    print(len(unique_ids))
    random.seed(random_seed)
    random.shuffle(unique_ids)

    train_ids, test_ids = split_ids_by_desired_rows(df, unique_ids, desired_train_rows)

    train_df = df[df["graph_id"].isin(train_ids)]
    test_df = df[df["graph_id"].isin(test_ids)]

    mae_d = train_and_test(train_df, test_df, f"{dirname}")

    print(mae_r, mae_d)

# 人工的にデータを追加する
# Leakageの起こった状態を再現
