library(graphtool)

football <- football_network()
g <- football$graph
set.seed(43)
state <- minimize_blockmodel_dl(g)
print(block_description_length(state))
print(block_sizes(state))

pos <- football$positions
edges <- graph_edges(g)
groups <- match(block_membership(state), sort(unique(block_membership(state))))
colors <- grDevices::hcl.colors(max(groups), "Dynamic")
graphics::plot(pos, type = "n", xlab = "", ylab = "", axes = FALSE, asp = 1)
graphics::segments(pos[edges[, 1], 1], pos[edges[, 1], 2],
                   pos[edges[, 2], 1], pos[edges[, 2], 2], col = "grey80")
graphics::points(pos, pch = 19, col = colors[groups], cex = 0.7)
graphics::title("Football network: inferred blocks")

counts <- block_edge_counts(state)
graphics::image(seq_len(nrow(counts)), seq_len(ncol(counts)), counts,
                xlab = "Block", ylab = "Block", main = "Block edge counts")
