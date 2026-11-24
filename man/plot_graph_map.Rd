\name{plot_graph_map}
\alias{plot_graph_map}
\alias{project_coords}
\alias{great_circle_arc}
\alias{world_coastlines}
\title{Geographic map projections, great-circle arcs, and spatial network plotting}
\usage{
project_coords(lon, lat,
  projection = c("robinson", "mercator", "equirectangular", "mollweide",
  "orthographic"), lon0 = 0, lat0 = 0)

great_circle_arc(lon1, lat1, lon2, lat2, n = 24L)

world_coastlines()

plot_graph_map(graph, lon, lat,
  projection = c("robinson", "mercator", "equirectangular", "mollweide",
  "orthographic"), lon0 = 0, lat0 = 0, xlim = NULL, ylim = NULL,
  geodesic = TRUE, n_arc = 24L, coastlines = TRUE, polygons = NULL,
  graticule = TRUE, land_col = "#ECEAE4", border_col = "#94A3B8",
  ocean_col = "#E0F2FE",
  graticule_col = grDevices::adjustcolor("#64748B", 0.30),
  vertex_color = "#D62828", vertex_border = "#1F2937", vertex_size = 1.2,
  vertex_label = NULL, vertex_label_cex = 0.7, vertex_label_pos = 3L,
  edge_color = grDevices::adjustcolor("#1D4ED8", 0.55), edge_lwd = 1.1,
  eorder = NULL, vorder = NULL, arrows = NULL, arrow_length = 0.06,
  main = NULL, output = NULL, width = 9, height = 5, dpi = 150, ...)
}
\arguments{
  \item{lon, lat}{Numeric vectors of longitudes (degrees East) and latitudes
  (degrees North).}
  \item{projection}{Cartographic map projection: \code{"robinson"},
  \code{"mercator"}, \code{"equirectangular"}, \code{"mollweide"}, or
  \code{"orthographic"}.}
  \item{lon0, lat0}{Central meridian and central parallel (degrees) for the
  projection.}
  \item{lon1, lat1, lon2, lat2}{Start and end geographic coordinates (degrees)
  for \code{great_circle_arc()}.}
  \item{n}{Number of interpolation points along the spherical great-circle arc.}
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{xlim, ylim}{Optional geographic bounding box (\code{c(lon_min, lon_max)}
  and \code{c(lat_min, lat_max)}) for regional maps.}
  \item{geodesic}{Logical; if \code{TRUE}, draw edges along spherical
  great-circle routes split cleanly across the antimeridian.}
  \item{n_arc}{Number of points per geodesic edge arc.}
  \item{coastlines}{Logical; whether to draw coastline land polygons.}
  \item{polygons}{Optional custom data frame with columns \code{lon},
  \code{lat}, and optional \code{group} (defaults to built-in offline
  \code{world_coastlines()}).}
  \item{graticule}{Logical; whether to draw latitude/longitude grid lines.}
  \item{land_col, border_col, ocean_col, graticule_col}{Basemap colours.}
  \item{vertex_color, vertex_border, vertex_size, vertex_label,
  vertex_label_cex, vertex_label_pos}{Vertex aesthetics and label placement.}
  \item{edge_color, edge_lwd, eorder, vorder, arrows, arrow_length}{Edge
  aesthetics, z-ordering, and arrowheads.}
  \item{main}{Optional title.}
  \item{output}{Optional destination file path exported via
  \code{\link{save_graph_plot}} (PDF, SVG, EPS, PNG, JPEG, TIFF).}
  \item{width, height, dpi}{Output dimensions (inches) and raster DPI when
  \code{output} is supplied.}
  \item{...}{Additional arguments forwarded to \code{\link[graphics]{points}}.}
}
\description{
R-native equivalent of graph-tool's Basemap/Cartopy geographic plotting
integration. Runs 100\% offline with zero external GIS or Python dependencies
using pure-R cartographic projections, spherical great-circle interpolation,
and built-in low-resolution world coastlines (with support for custom
user-supplied polygon data frames).
}
\value{
\code{project_coords()} returns a data frame with columns \code{x}, \code{y},
and \code{visible}. \code{great_circle_arc()} returns a data frame with columns
\code{lon}, \code{lat}, and \code{segment}. \code{world_coastlines()} returns a
data frame with columns \code{lon}, \code{lat}, \code{group}, and
\code{region}. \code{plot_graph_map()} invisibly returns a list containing
\code{projected}, \code{projection}, \code{lon0}, and \code{lat0}.
}
\examples{
g <- graph_create(4L, directed = FALSE)
graph_add_edges(g, rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L), c(4L, 1L)))
lon <- c(-74.0, -0.1, 139.7, -118.2)
lat <- c(40.7, 51.5, 35.7, 34.0)
plot_graph_map(g, lon = lon, lat = lat, projection = "robinson",
               vertex_label = c("NYC", "LON", "TYO", "LAX"),
               output = tempfile(fileext = ".pdf"))
}
