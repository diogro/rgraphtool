"""Fit the nested SBM to a Netzschleuder network with graph-tool 3.9.

Usage: pixi run --manifest-path tools/reference/pixi.toml \
         python tools/reference/nested_fits.py nseeds [outdir [network]]
Fits a Netzschleuder network (default "football").
Prints CSV: seed, degree_corrected, description_length, level sizes. With
outdir, also writes each fitted hierarchy as JSON (1-based contiguous labels,
element l labels the groups of level l - 1) for fixed-hierarchy score checks.
"""
import json
import os
import sys
import numpy as np
import graph_tool.all as gt

nseeds = int(sys.argv[1]) if len(sys.argv) > 1 else 10
outdir = sys.argv[2] if len(sys.argv) > 2 else None
name = sys.argv[3] if len(sys.argv) > 3 else "football"
g = gt.collection.ns[name]


def contiguous(bs):
    out, lower = [], None
    for l, b in enumerate(bs):
        b = np.asarray(b)
        items = b if l == 0 else b[lower]
        _, first = np.unique(items, return_index=True)
        order = items[np.sort(first)]
        index = {r: i + 1 for i, r in enumerate(order)}
        out.append([int(index[r]) for r in items])
        lower = order
    return out


print("seed,degree_corrected,description_length,B")
for seed in range(nseeds):
    for dc in (False, True):
        gt.seed_rng(seed)
        state = gt.minimize_nested_blockmodel_dl(g, state_args=dict(deg_corr=dc))
        bs = contiguous(state.get_bs())
        B = "-".join(str(max(b)) for b in bs)
        print(f"{seed},{dc},{state.entropy():.10f},{B}")
        if outdir:
            with open(os.path.join(outdir, f"nested_{seed}_{dc}.json"), "w") as f:
                json.dump({"bs": bs, "entropy": state.entropy()}, f)
