from __future__ import annotations

from typing import Any

import numpy as np
from pymatgen.analysis.local_env import CrystalNN, NearNeighbors
from pymatgen.core import IStructure, Structure
from sklearn.cluster import DBSCAN


def _get_original_site(structure, site):
    """Private convenience method for get_nn_info,
    gives original site index from ProvidedPeriodicSite.
    """

    if isinstance(structure, IStructure | Structure):
        site_fcoords = site.frac_coords
        strc_fcoords = structure.frac_coords
        tol = 1e-8  # threshold in Site.is_periodic_image
        # sort to reduce the iteration
        nearest_i = np.argsort(-(np.abs(strc_fcoords - site_fcoords) < tol).sum(axis=1))

        for i in nearest_i:
            if site.is_periodic_image(structure[i]):
                return i
    else:
        for i, s in enumerate(structure):
            if site == s:
                return i
    raise Exception("Site not found!")  # noqa: TRY002, TRY003, EM101


class DistanceClusteringNN(NearNeighbors):
    # 結合長のクラスタリングによって原子に近いクラスター順に番号を振る
    # その番号と元素記号を使ってGraph IDを計算する
    def __init__(self) -> None:
        """ """

    @property
    def structures_allowed(self) -> bool:
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    def get_nn_info(self, structure: Structure, n: int, rank_k: int, cutoff: float = 6.0) -> list[dict[str, Any]]:
        """
        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near
                neighbors.
            cutoff (float): distance cutoff parameter.
        Returns:
            siw (list[dict]): dicts with (Site, array, float) each one of which represents a
                neighbor site, its image location, and its weight.
        """

        site = structure[n]
        cutoff_cluster_list = self.get_cutoff_cluster(structure, n, cutoff)
        if len(cutoff_cluster_list) <= rank_k:
            return []

        neighs_dists = structure.get_neighbors(site, cutoff_cluster_list[rank_k])
        max_weight = round(cutoff_cluster_list[rank_k], 3)
        # is_periodic = isinstance(structure, Structure | IStructure) # Python 3.10 以降でのみサポート
        is_periodic = isinstance(structure, (IStructure, Structure))
        siw = []

        for nn in neighs_dists:
            weight = round(nn.nn_distance, 3)
            if (rank_k > 0 and weight <= max_weight and weight > round(cutoff_cluster_list[rank_k - 1], 3)) or (
                rank_k == 0 and weight <= max_weight
            ):
                siw.append(
                    {
                        "site": nn,
                        "image": self._get_image(structure, nn) if is_periodic else None,
                        "weight": weight,
                        "site_index": self._get_original_site(structure, nn),
                        "edge_properties": {"cluster_idx": rank_k + 1},
                    },
                )

        return siw

    def get_cutoff_cluster(self, structure: Structure, n: int, cutoff: float = 6.0) -> list:
        """
        DBSCANによって得られた距離のクラスターから結合判定の閾値を決定する
        """

        # # スーパーセルを作成し、6.0angまでの結合長を数え上げる
        # copy_structure = structure.copy()
        # supercell = copy_structure.make_supercell([3, 3, 3])
        # site_i = structure[n]

        # site_index = None
        # for idx, site in enumerate(supercell):
        #     # Siteのdistanceメソッドを使うとなぜか正しく距離が計算されない
        #     if float(np.linalg.norm(site_i.coords - site.coords)) < 0.01:
        #         site_index = idx
        #         break

        distance_list = []
        neighbors = structure.get_sites_in_sphere(structure[n].coords, cutoff)
        for neighbor in neighbors:
            dist = neighbor.nn_distance
            distance_list.append([dist, 0])

        dbscan = DBSCAN(eps=0.5, min_samples=2)
        dbscan.fit(distance_list)
        labels = dbscan.labels_

        max_dist_list = [0 for _ in range(max(labels) + 1)]
        for label_number in range(max(labels) + 1):
            max_dist = 0
            for label, distance in zip(labels, distance_list):
                if label == label_number:
                    max_dist = max(max_dist, distance[0])

            max_dist_list[label_number] = max_dist

        return sorted(max_dist_list)


class BondClusteringNN(CrystalNN):
    def __init__(
        self,
        weighted_cn=False,
        cation_anion=False,
        distance_cutoffs=(0.5, 1),
        x_diff_weight=3.0,
        porous_adjustment=True,
        search_cutoff=10,
        fingerprint_length=None,
    ) -> None:
        """
        Initialize CrystalNN with desired parameters. Default parameters assume
        "chemical bond" type behavior is desired. For geometric neighbor
        finding (e.g., structural framework), set (i) distance_cutoffs=None,
        (ii) x_diff_weight=0 and (optionally) (iii) porous_adjustment=False
        which will disregard the atomic identities and perform best for a purely
        geometric match.

        Args:
            weighted_cn (bool): if set to True, will return fractional weights
                for each potential near neighbor.
            cation_anion (bool): if set True, will restrict bonding targets to
                sites with opposite or zero charge. Requires an oxidation states
                on all sites in the structure.
            distance_cutoffs ([float, float]): - if not None, penalizes neighbor
                distances greater than sum of covalent radii plus
                distance_cutoffs[0]. Distances greater than covalent radii sum
                plus distance_cutoffs[1] are enforced to have zero weight.
            x_diff_weight (float): - if multiple types of neighbor elements are
                possible, this sets preferences for targets with higher
                electronegativity difference.
            porous_adjustment (bool): - if True, readjusts Voronoi weights to
                better describe layered / porous structures
            search_cutoff (float): cutoff in Angstroms for initial neighbor
                search; this will be adjusted if needed internally
            fingerprint_length (int): if a fixed_length CN "fingerprint" is
                desired from get_nn_data(), set this parameter
        """
        super().__init__(
            weighted_cn,
            cation_anion,
            distance_cutoffs,
            x_diff_weight,
            porous_adjustment,
            search_cutoff,
            fingerprint_length,
        )

    def get_nn_info(self, structure: Structure, n: int, cutoff: float = 10.0) -> list[dict[str, Any]]:
        """
        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near
                neighbors.
            cutoff (float): distance cutoff parameter.
        Returns:
            siw (list[dict]): dicts with (Site, array, float) each one of which represents a
                neighbor site, its image location, and its weight.
        """

        site = structure[n]
        nn_data = self.get_nn_data(structure, n)

        max_key = max(nn_data.cn_weights, key=lambda k: nn_data.cn_weights[k])
        nn = nn_data.cn_nninfo[max_key]

        weight_list = []
        for entry in nn_data.all_nninfo:
            weight = 0
            for cn in nn_data.cn_nninfo:
                for cn_entry in nn_data.cn_nninfo[cn]:
                    if entry["site"] == cn_entry["site"]:
                        weight += nn_data.cn_weights[cn]

            # entry["weight"] = weight
            weight_list.append([weight, 0])

        dbscan = DBSCAN(eps=0.5, min_samples=2)
        dbscan.fit(weight_list)
        labels = dbscan.labels_

        for entry, label in zip(nn_data.all_nninfo, labels):
            entry["weight"] = label + 1  # weightが0だとStructureGraphのedgeにweightを持たせてくれない

        return nn_data.all_nninfo
