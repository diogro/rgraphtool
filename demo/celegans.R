# Hierarchical partition of the C. elegans neural network (graph-tool cookbook,
# "Hierarchical partitions"). The data are downloaded from Netzschleuder.
library(graphtool)

celegans <- netzschleuder_network("celegansneural")
g <- celegans$graph
set.seed(1)
state <- minimize_nested_blockmodel_dl(g)
print(state)                       # equivalent of state.print_summary()
nested_summary(state)              # per-level sizes and description lengths

# Group memberships of one vertex at every level (vertex 47 here is vertex 46
# in graph-tool's zero-based numbering)
sapply(seq_len(nrow(nested_summary(state))) - 1L,
       function(l) nested_project(state, l)[47])

# The level-0 partition as a flat block state
level0 <- nested_level_state(state, 0)
block_sizes(level0)

plot_hierarchy(state, main = "C. elegans: nested SBM")
