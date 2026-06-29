from __future__ import annotations

from pathlib import Path

import orjson
from huggingface_hub import hf_hub_download
from huggingface_hub.utils import EntryNotFoundError

DB_PATH = Path(__file__).parent.parent / "raw/id_jsons"


class Finder:
    def find(self, graph_id: str, is_fast: bool = False) -> list[dict[str, str]]:
        """Args:
        ----
            graph_id(str): GraphID calculated using graph-id-core
            is_fast(bool): If True, find only on-memory entries (from IZA and Materials Project)
        """
        if is_fast:
            return self.find_fast(graph_id)

        else:
            ret_dict = []
            fast_graph_ids = self.find_fast(graph_id)
            if fast_graph_ids:
                ret_dict += fast_graph_ids

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

    def find_fast(self, graph_id: str) -> list[dict[str, str]]:
        """Find only on-memory entries."""
        ret_dict = []

        dir_name = graph_id[:2]
        file_name = graph_id[:4]

        db_path = DB_PATH / dir_name / f"{file_name}.json"
        if db_path.exists():
            with open(db_path) as f:
                docs = orjson.loads(f.read())
                if docs.get(graph_id):
                    ret_dict = docs.get(graph_id)

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

            with open(local_path) as f:
                docs = orjson.loads(f.read())
                if docs.get(graph_id):
                    ret_dict = docs.get(graph_id)

            return ret_dict
        except EntryNotFoundError:
            return []

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

            with open(local_path) as f:
                docs = orjson.loads(f.read())
                if docs.get(graph_id):
                    ret_dict = docs.get(graph_id)

            return ret_dict
        except EntryNotFoundError:
            return []

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

            with open(local_path) as f:
                docs = orjson.loads(f.read())
                if docs.get(graph_id):
                    ret_dict = docs.get(graph_id)

            return ret_dict
        except EntryNotFoundError:
            return []
