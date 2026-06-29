from __future__ import annotations

import logging
import re
from pathlib import Path

import orjson
from huggingface_hub import hf_hub_download
from huggingface_hub.utils import EntryNotFoundError

DB_PATH = Path(__file__).parent.parent / "raw/id_jsons"

logger = logging.getLogger(__name__)

# A canonical Graph ID is a blake2b hex digest (lowercase hex, >= 4 chars so the
# graph_id[:2]/graph_id[:4] path components are well-formed). Restricting to hex
# also prevents path separators or traversal sequences from reaching the lookup
# path builders.
GRAPH_ID_PATTERN = re.compile(r"[0-9a-f]{4,}")


class Finder:

    """Finder to find chemical structures from databases using Graph ID.

    Currently supports AFLOW, OQMD, and PCOD.
    If you want to find a structure in IZA and the Materials Project,
    you currently need to use the graph-id-db package.
    """

    def __init__(self, backends: list[str] | None = None):
        """Initialize the Finder.

        Args:
        ----
            backends(list[str]): List of backends to use.
        """
        if backends is None:
            backends = ["aflow", "oqmd", "pcod"]
        self.backends = backends
        self.backends_map = {backend: getattr(self, f"find_{backend}_entries") for backend in backends}

    def find(self, graph_id: str) -> list[dict[str, str]]:
        """Find chemical structures from databases using Graph ID.

        Args:
        ----
            graph_id(str): GraphID calculated using graph-id-core
        """
        ret_dict: list[dict[str, str]] = []

        if not GRAPH_ID_PATTERN.fullmatch(graph_id):
            return ret_dict

        for backend in self.backends:
            try:
                graph_ids = self.backends_map[backend](graph_id)
            except Exception:  # noqa: BLE001
                # A transient failure in one source must not discard results
                # already collected from the other sources.
                logger.warning("Lookup via %s failed", backend, exc_info=True)
                continue
            if graph_ids:
                ret_dict += graph_ids

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
