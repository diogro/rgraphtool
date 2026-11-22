"""Fit the football network with graph-tool 3.9 for several seeds.

Usage: pixi run --manifest-path tools/reference/pixi.toml \
         python tools/reference/football_fits.py [nseeds]
Prints CSV: seed, degree_corrected, B, description_length.
"""
import sys
import graph_tool.all as gt

nseeds = int(sys.argv[1]) if len(sys.argv) > 1 else 10
g = gt.collection.ns["football"]
print("seed,degree_corrected,B,description_length")
for seed in range(nseeds):
    for dc in (False, True):
        gt.seed_rng(seed)
        state = gt.minimize_blockmodel_dl(g, state_args=dict(deg_corr=dc))
        print(f"{seed},{dc},{len(set(state.get_blocks().a))},{state.entropy():.10f}")
