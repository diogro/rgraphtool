# Cookbook demo: Geographic network plotting on projected world and regional maps
# R-native counterpart of graph-tool's Basemap/Cartopy map-integration example.
# Runs 100% offline using built-in coastlines, cartographic projections, and
# spherical great-circle geodesic arcs (no Python or external GIS dependencies).

library(rgraphtool)

# 18 global aviation/research hub cities with (longitude, latitude)
cities <- data.frame(
  code = c("NYC", "LAX", "ORD", "MIA", "GRU", "EZE",
           "LON", "PAR", "FRA", "MAD", "IST", "CAI",
           "JNB", "DXB", "BOM", "SIN", "TYO", "SYD"),
  lon = c(-74.0, -118.2, -87.6, -80.2, -46.6, -58.4,
          -0.1,   2.35,   8.68, -3.70, 28.98, 31.24,
          28.05,  55.27,  72.88, 103.82, 139.69, 151.21),
  lat = c( 40.7,  34.0,  41.9,  25.8, -23.5, -34.6,
           51.5,  48.9,  50.1,  40.4,  41.0,  30.0,
          -26.2,  25.2,  19.1,   1.35,  35.68, -33.87),
  stringsAsFactors = FALSE
)

# Intercontinental and regional routes (including trans-Pacific routes)
routes <- rbind(
  c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(1, 7), c(1, 8), c(1, 9), c(1, 10),
  c(2, 3), c(2, 17), c(2, 18), c(3, 7), c(3, 9), c(3, 17),
  c(4, 5), c(4, 6), c(4, 10), c(5, 6), c(5, 7), c(5, 10), c(5, 13),
  c(7, 8), c(7, 9), c(7, 10), c(7, 11), c(7, 13), c(7, 14), c(7, 16), c(7, 17),
  c(8, 9), c(8, 10), c(8, 12), c(8, 14), c(9, 11), c(9, 14), c(9, 15), c(9, 16),
  c(11, 12), c(11, 14), c(12, 13), c(12, 14), c(13, 14), c(13, 18),
  c(14, 15), c(14, 16), c(14, 17), c(14, 18),
  c(15, 16), c(16, 17), c(16, 18), c(17, 18)
)

g <- graph_create(nrow(cities), directed = FALSE)
graph_add_edges(g, routes)

deg <- graph_degree(g)
vcol <- map_colors(deg, palette = "YlOrRd", rev = TRUE)
vsiz <- prop_to_size(deg, mi = 1.0, ma = 2.4, power = 0.6)

# 1. Global Robinson projection with geodesic great-circle arcs (vector PDF + raster PNG)
f_rob_pdf <- file.path(tempdir(), "rgraphtool_map_robinson.pdf")
f_rob_png <- file.path(tempdir(), "rgraphtool_map_robinson.png")

plot_graph_map(
  g, lon = cities$lon, lat = cities$lat,
  projection = "robinson",
  geodesic = TRUE,
  vertex_color = vcol,
  vertex_size = vsiz,
  vorder = deg,
  vertex_label = cities$code,
  vertex_label_cex = 0.65,
  edge_color = grDevices::adjustcolor("#1D4ED8", 0.45),
  edge_lwd = 1.1,
  main = "Global Hub Network (Robinson Projection, Great-Circle Routes)",
  output = f_rob_pdf,
  width = 10, height = 5.4
)

plot_graph_map(
  g, lon = cities$lon, lat = cities$lat,
  projection = "robinson",
  geodesic = TRUE,
  vertex_color = vcol,
  vertex_size = vsiz,
  vorder = deg,
  vertex_label = cities$code,
  vertex_label_cex = 0.65,
  edge_color = grDevices::adjustcolor("#1D4ED8", 0.45),
  edge_lwd = 1.1,
  main = "Global Hub Network (Robinson Projection, Great-Circle Routes)",
  output = f_rob_png,
  width = 10, height = 5.4, dpi = 150
)

# 2. Orthographic 3D globe view centered on the Atlantic
f_ortho_png <- file.path(tempdir(), "rgraphtool_map_orthographic.png")
plot_graph_map(
  g, lon = cities$lon, lat = cities$lat,
  projection = "orthographic",
  lon0 = -20, lat0 = 25,
  geodesic = TRUE,
  vertex_color = vcol,
  vertex_size = vsiz,
  vertex_label = cities$code,
  vertex_label_cex = 0.7,
  edge_color = grDevices::adjustcolor("#0F766E", 0.55),
  main = "Atlantic Hemisphere (Orthographic Globe)",
  output = f_ortho_png,
  width = 6, height = 6, dpi = 150
)

# 3. Regional Mercator projection zoomed on North Atlantic & Europe
f_reg_pdf <- file.path(tempdir(), "rgraphtool_map_regional.pdf")
plot_graph_map(
  g, lon = cities$lon, lat = cities$lat,
  projection = "mercator",
  xlim = c(-100, 45), ylim = c(18, 62),
  geodesic = TRUE,
  vertex_color = vcol,
  vertex_size = vsiz,
  vertex_label = cities$code,
  vertex_label_cex = 0.75,
  main = "Trans-Atlantic & European Corridor (Mercator)",
  output = f_reg_pdf,
  width = 8.5, height = 5.2
)

cat(sprintf("Wrote Robinson PDF (%d bytes): %s\n", file.info(f_rob_pdf)$size, f_rob_pdf))
cat(sprintf("Wrote Robinson PNG (%d bytes): %s\n", file.info(f_rob_png)$size, f_rob_png))
cat(sprintf("Wrote Orthographic PNG (%d bytes): %s\n", file.info(f_ortho_png)$size, f_ortho_png))
cat(sprintf("Wrote Regional Mercator PDF (%d bytes): %s\n", file.info(f_reg_pdf)$size, f_reg_pdf))
