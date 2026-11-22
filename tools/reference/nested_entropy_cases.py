"""Score a fixed nested hierarchy with graph-tool 3.9.

Usage: python nested_entropy_cases.py edges.csv hierarchy.json directed
edges.csv: 0-based "source,target" rows. hierarchy.json: a list of levels with
1-based contiguous labels; level 1 labels vertices, level l labels the groups
of level l - 1. Prints the description length for each model as CSV.
"""
import json
import sys
import numpy as np
import graph_tool.all as gt

edges = np.loadtxt(sys.argv[1], delimiter=",", dtype=np.int64, ndmin=2)
bs = [np.asarray(b, dtype=np.int64) - 1 for b in json.load(open(sys.argv[2]))]
directed = sys.argv[3].lower() in ("1", "true")
g = gt.Graph(directed=directed)
g.add_vertex(len(bs[0]))
g.add_edge_list(edges)
print("degree_corrected,degree_prior,description_length")
for dc in (False, True):
    kinds = ("distributed", "uniform") if dc else ("distributed",)
    for kind in kinds:
        state = gt.NestedBlockState(g, bs=bs, base_state_args=dict(deg_corr=dc))
        S = state.entropy(degree_dl_kind=kind)
        print(f"{dc},{kind},{S:.12f}")
