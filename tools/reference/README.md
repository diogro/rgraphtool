# graph-tool reference values

Scripts here produce the reference numbers used in the package tests and in
`ROADMAP.md`. They run against the official graph-tool 3.9 conda-forge build in
a project-local [pixi](https://pixi.sh) environment. Python is only needed to
regenerate references, never to use the package.

```sh
cd tools/reference
pixi install                     # graph-tool-base 3.9, pinned by pixi.lock
pixi run python entropy_cases.py edges.csv blocks.csv false
pixi run python football_fits.py 20 > football_upstream.csv
pixi run python nested_entropy_cases.py edges.csv hierarchy.json false
pixi run python nested_fits.py 20 hierarchies/ > nested_upstream.csv
pixi run python nested_fits.py 10 hierarchies/ celegansneural > celegans_nested_upstream.csv
```

- `entropy_cases.py` scores a fixed partition (0-based edge list and labels)
  with `BlockState.entropy()` for the non-degree-corrected model and for the
  degree-corrected model with both degree priors.
- `football_fits.py` runs `minimize_blockmodel_dl()` on the cookbook football
  network for seeds `0..n-1`, for both models. `football_upstream.csv` holds the
  output for 20 seeds.
- `nested_entropy_cases.py` scores a fixed nested hierarchy (a JSON list of
  1-based level label vectors) with `NestedBlockState.entropy()`.
- `nested_fits.py` runs `minimize_nested_blockmodel_dl()` on a Netzschleuder
  network for seeds `0..n-1` and both models. It can also write each fitted
  hierarchy as JSON, so that `nested_block_state()` can be checked against
  graph-tool's score for the same hierarchy. `nested_upstream.csv` (football,
  20 seeds) and `celegans_nested_upstream.csv` (C. elegans, 10 seeds) hold its
  output.

The drawing module is not installed (no cairo), which graph-tool reports as a
warning on import; inference does not need it.
