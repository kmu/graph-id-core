import argparse
import csv
from datetime import datetime

from graph_id import GraphIDGenerator
from pymatgen.core import Structure
from tabulate import tabulate

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Graph ID: graph-based ID for materails")

    parser.add_argument(
        dest="filenames",
        metavar="filename",
        nargs="+",
        help="List of structure files.",
        default=[],
    )

    parser.add_argument("-p", "--parallel", help="parallel execution", action="store_true")

    gid = GraphIDGenerator()
    args = parser.parse_args()

    table = []

    for fname in args.filenames:
        s = Structure.from_file(fname)
        s.merge_sites(mode="delete")

        my_id = gid.get_id(s)

        table.append([my_id, fname])

    t_headers = ["GraphIDGenerator", "Filename"]

    print(tabulate(table, headers=t_headers))

    now = datetime.now()
    now_str = now.strftime("%Y%m%d-%H%M%S")

    with open(f"graph_id_{now_str}.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(table)
