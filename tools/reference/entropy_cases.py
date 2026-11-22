"""Score fixed partitions with graph-tool 3.9.

Usage: python entropy_cases.py edges.csv blocks.csv directed
edges.csv: 0-based "source,target" rows; blocks.csv: one 0-based label per line.
Prints dl for (deg_corr, degree_dl_kind) combinations as CSV.
"""
import sys
import numpy as np
import graph_tool.all as gt

edges = np.loadtxt(sys.argv[1], delimiter=",", dtype=np.int64, ndmin=2)
blocks = np.loadtxt(sys.argv[2], dtype=np.int64, ndmin=1)
directed = sys.argv[3].lower() in ("1", "true")
g = gt.Graph(directed=directed)
g.add_vertex(len(blocks))
g.add_edge_list(edges)
b = g.new_vp("int", vals=blocks)
print("degree_corrected,degree_prior,description_length")
for dc in (False, True):
    state = gt.BlockState(g, b=b, deg_corr=dc)
    kinds = ("distributed", "uniform") if dc else ("distributed",)
    for kind in kinds:
        S = state.entropy(degree_dl_kind=kind)
        print(f"{dc},{kind},{S:.12f}")
