from __future__ import annotations

from pathlib import Path

import orjson
from huggingface_hub import hf_hub_download
from huggingface_hub.utils import EntryNotFoundError

DB_PATH = Path(__file__).parent.parent / "raw/id_jsons"


class Finder:

    """Finder to find chemical structures from databases using Graph ID.

    Currently supports AFLOW, OQMD, and PCOD.
    If you want to find a structure in IZA and the Materials Project,
    you currently need to use the graph-id-db package.
    """

    def find(self, graph_id: str) -> list[dict[str, str]]:
        """Find chemical structures from databases using Graph ID.

        Args:
        ----
            graph_id(str): GraphID calculated using graph-id-core
        """
        ret_dict = []

        aflow_graph_ids = self.find_aflow_entries(graph_id)
        if aflow_graph_ids:
            ret_dict += aflow_graph_ids

        oqmd_graph_ids = self.find_oqmd_entries(graph_id)
        if oqmd_graph_ids:
            ret_dict += oqmd_graph_ids

        pcod_graph_ids = self.find_pcod_entries(graph_id)
        if pcod_graph_ids:
            ret_dict += pcod_graph_ids

        return ret_dict

    def find_aflow_entries(self, graph_id: str) -> list[dict[str, str]]:
        """Find only AFLOW entries."""
        ret_dict = []

        dir_name = graph_id[:2]
        file_name = graph_id[:4]

        try:
            local_path = hf_hub_download(
                repo_id="kamabata/aflow_graph_ids",
                filename=f"id_jsons/{dir_name}/{file_name}.json",
                repo_type="dataset",
            )

            with Path(local_path).open() as f:
                docs = orjson.loads(f.read())
                if docs.get(graph_id):
                    ret_dict = docs.get(graph_id)

        except EntryNotFoundError:
            return []

        return ret_dict

    def find_oqmd_entries(self, graph_id: str) -> list[dict[str, str]]:
        """Find only OQMD entries."""
        ret_dict = []

        dir_name = graph_id[:2]
        file_name = graph_id[:4]

        try:
            local_path = hf_hub_download(
                repo_id="kamabata/oqmd_graph_ids",
                filename=f"id_jsons/{dir_name}/{file_name}.json",
                repo_type="dataset",
            )

            with Path(local_path).open() as f:
                docs = orjson.loads(f.read())
                if docs.get(graph_id):
                    ret_dict = docs.get(graph_id)

        except EntryNotFoundError:
            return []

        return ret_dict

    def find_pcod_entries(self, graph_id: str) -> list[dict[str, str]]:
        """Find only PCOD entries."""
        ret_dict = []

        dir_name = graph_id[:2]
        file_name = graph_id[:4]

        try:
            local_path = hf_hub_download(
                repo_id="kamabata/pcod_graph_ids",
                filename=f"id_jsons/{dir_name}/{file_name}.json",
                repo_type="dataset",
            )

            with Path(local_path).open() as f:
                docs = orjson.loads(f.read())
                if docs.get(graph_id):
                    ret_dict = docs.get(graph_id)

        except EntryNotFoundError:
            return []

        return ret_dict
