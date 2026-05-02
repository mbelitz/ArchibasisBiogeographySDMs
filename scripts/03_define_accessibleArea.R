# =============================================================================
# 03_define_accessibleArea.R — Step 2: Define the Accessible Area (M Region)
# =============================================================================
#
# WHAT THIS STEP DOES (conceptually)
# ------------------------------------
# MaxEnt models the relationship between occurrence locations and environment.
# A foundational choice is: which background environment should the model
# contrast against occurrences? Using the entire world as background means
# most background cells are climatically dissimilar to the occurrence region
# simply because of geography, not ecology — making the model trivially easy
# and ecologically uninformative.
#
# The solution is to restrict background sampling to the "accessible area"
# (sometimes called the M region in the BAM framework): the geographic area
# that the species *could* have colonised given dispersal ability and time.
# This script defines M as an alpha hull drawn around occurrence points,
# buffered outward by a distance derived from how spread-out the occurrences
# are (plus a user-specified minimum buffer).
#
# Alpha hulls are better than simple bounding boxes or convex hulls because
# they can represent non-convex ranges (e.g., a species that follows a river
# system or archipelago) without enclosing large unsuitable gaps.
#
# WHAT THE USER CAN CHANGE
# -------------------------
#   minBuff      — Minimum buffer distance in metres around the alpha hull.
#                  The actual buffer is max(80th-percentile nearest-neighbour
#                  distance, minBuff). Increase minBuff if you want a wider
#                  accessible area even when records are tightly clustered
#                  (e.g., 150000 for 150 km). Decrease if occurrences already
#                  span a large region and you do not want to pull in distant
#                  background environments. Default: 75000 (75 km).
#
#   fraction / partCount / initialAlpha (getDynamicAlphaHull)
#                — These control alpha-hull geometry. fraction = 1 means
#                  all points must be enclosed. Lower values (e.g., 0.95)
#                  allow a small fraction of outlier points to fall outside
#                  the hull (useful when a few outlier records would otherwise
#                  create a huge hull). partCount = 1 enforces a single
#                  contiguous polygon; increase if the species has genuinely
#                  disjunct island populations that should each get their own
#                  M region. initialAlpha controls granularity of the hull
#                  — smaller values produce tighter, more jagged hulls.
#
#   clipToCoast  — "terrestrial" removes marine portions of the hull, which
#                  is appropriate for dragonflies. Change to "aquatic" or
#                  omit for marine species.
# =============================================================================

library(rangeBuilder)
library(sf)
library(terra)
library(tidyverse)
library(rnaturalearth)

#' Define accessible area for SDMs using an alpha hull approach.
#'
#' The buffer radius is the larger of:
#'   (a) the 80th percentile of pairwise nearest-neighbour distances among
#'       occurrence records (captures the typical dispersal gap in the data),
#'   (b) minBuff (a hard lower bound so the M region never collapses to
#'       just the occurrence points for very dense datasets).
#'
#' All distance calculations use the Equal-Area Cylindrical (CEA) projection
#' to ensure buffer distances are in metres, not degrees.
#'
#' @param species_df     Data frame with columns LONG, LAT, GEN_DB, SP_DB.
#' @param minBuff        Minimum buffer in metres (default 75000 = 75 km).
#' @param imagedir       Directory for the optional preview PNG.
#' @param shapefiledir   Directory for the optional shapefile.
#' @param saveImage      Logical. Save a preview map? Default FALSE.
#' @param saveShapefile  Logical. Save the M polygon as a shapefile? Default FALSE.
#' @return  sf polygon (WGS84 / EPSG:4326) representing the accessible area.

define_accessibleArea <- function(species_df, minBuff = 75000,
                                   imagedir     = "NotSaving",
                                   shapefiledir = "NotSaving",
                                   saveImage    = FALSE,
                                   saveShapefile = FALSE) {

  temp <- as.data.frame(species_df)

  # ── Distance-based buffer calculation ───────────────────────────────────────
  # Project to CEA (cylindrical equal-area) so st_distance() returns metres,
  # not degrees. Degrees are not a meaningful unit for buffering because 1°
  # of longitude shrinks toward the poles.
  temp_sf   <- st_as_sf(temp, coords = c("LONG", "LAT"), crs = 4326)
  tempTrans <- st_transform(temp_sf, crs = "+proj=cea +lat_ts=0 +lon_0=0")

  # Compute the 80th percentile of nearest-neighbour distances.
  # The 80th percentile (rather than mean or max) is robust to a few very
  # isolated records that would otherwise inflate the buffer.
  dist_mat  <- as.matrix(st_distance(tempTrans))
  diag(dist_mat) <- NA
  nn_dists  <- apply(dist_mat, 2, min, na.rm = TRUE)
  buffDist  <- as.numeric(quantile(nn_dists, probs = 0.80, na.rm = TRUE))
  buffDist  <- max(buffDist, minBuff)   # enforce the user-specified floor

  # ── Alpha hull ───────────────────────────────────────────────────────────────
  # getDynamicAlphaHull iteratively tightens the alpha until a single polygon
  # (partCount = 1) encloses fraction = 1 of all points, then clips to
  # terrestrial land. The result is a SpatialPolygons object; we convert to sf.
  shape <- getDynamicAlphaHull(
    x            = temp[, c("LONG", "LAT")],
    fraction     = 1,          # all points must be inside the hull
    partCount    = 1,          # enforce one contiguous polygon
    initialAlpha = 20,         # starting alpha value (degrees)
    clipToCoast  = "terrestrial",
    coordHeaders = c("LONG", "LAT")
  )

  # ── Buffer + dissolve ────────────────────────────────────────────────────────
  # Buffer the alpha hull by buffDist + 25000 m (a fixed 25 km padding added
  # on top of the data-driven distance, giving a bit more background context),
  # dissolve any overlapping sub-polygons into one, then reproject to WGS84.
  shape_sf   <- st_as_sf(shape[[1]])
  shapeTrans <- st_transform(shape_sf, crs = "+proj=cea +lat_ts=0 +lon_0=0")
  shape2     <- st_buffer(shapeTrans, dist = buffDist + 25000) |>
                  st_union() |>
                  st_as_sf()
  shape2_geo <- st_transform(shape2, crs = 4326)

  # ── Visualisation ────────────────────────────────────────────────────────────
  world  <- ne_countries(scale = "medium", returnclass = "sf")
  plotz  <- ggplot() +
    geom_sf(data = world,      fill = NA) +
    geom_sf(data = shape2_geo, fill = "orange", alpha = 0.5) +
    geom_sf(data = temp_sf,    colour = "black", size = 0.8) +
    coord_sf(
      xlim = c(min(temp$LONG) - 10, max(temp$LONG) + 10),
      ylim = c(min(temp$LAT)  - 10, max(temp$LAT)  + 10)
    ) +
    theme_classic()

  if (saveImage) {
    ggsave(
      filename = paste0(imagedir,
                        paste0(species_df$GEN_DB[1], "_", species_df$SP_DB[1]),
                        ".png"),
      plot = plotz
    )
  } else {
    message("Alpha hull plot not saved. Set saveImage = TRUE to save.")
  }

  if (saveShapefile) {
    st_write(shape2_geo,
             dsn          = paste0(shapefiledir,
                                   paste0(species_df$GEN_DB[1], "_",
                                          species_df$SP_DB[1]),
                                   "_AA.shp"),
             delete_layer = TRUE)
  } else {
    message("Accessible area shapefile not saved. Set saveShapefile = TRUE to save.")
  }

  return(shape2_geo)
}
