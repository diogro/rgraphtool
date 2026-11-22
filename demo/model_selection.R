# Model selection between degree-corrected and non-degree-corrected nested SBMs
# (graph-tool cookbook, "Model selection"). Data are downloaded from
# Netzschleuder. The political blogs fits take about 15 s each.
library(graphtool)

polblogs <- netzschleuder_network("polblogs")
set.seed(1)
print(compare_blockmodels(polblogs$graph, nested = TRUE))

football <- football_network()
set.seed(1)
comparison <- compare_blockmodels(football$graph, nested = TRUE, nruns = 3)
print(comparison)
comparison$comparison
