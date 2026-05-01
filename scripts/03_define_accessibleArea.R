library(rangeBuilder)
library(sf)
library(terra)
library(tidyverse)
library(rnaturalearth)

#' Define accessible area for SDMs using an alpha hull approach.
#' Buffer is the larger of (80th percentile nearest-neighbor distance) or minBuff.
#' Returns an sf polygon in WGS84 (EPSG:4326).

define_accessibleArea <- function(species_df, minBuff = 75000,
                                   imagedir    = "NotSaving",
                                   shapefiledir = "NotSaving",
                                   saveImage    = FALSE,
                                   saveShapefile = FALSE) {

  temp <- as.data.frame(species_df)

  # --- Distance-based buffer calculation (CEA projection, metres) ---
  temp_sf  <- st_as_sf(temp, coords = c("LONG", "LAT"), crs = 4326)
  tempTrans <- st_transform(temp_sf, crs = "+proj=cea +lat_ts=0 +lon_0=0")

  # 80th percentile of nearest-neighbour distances
  dist_mat  <- as.matrix(st_distance(tempTrans))
  diag(dist_mat) <- NA
  nn_dists  <- apply(dist_mat, 2, min, na.rm = TRUE)
  buffDist  <- as.numeric(quantile(nn_dists, probs = 0.80, na.rm = TRUE))
  buffDist  <- max(buffDist, minBuff)

  # --- Alpha hull (rangeBuilder returns SpatialPolygons; convert to sf) ---
  shape <- getDynamicAlphaHull(
    x            = temp[, c("LONG", "LAT")],
    fraction     = 1,
    partCount    = 1,
    initialAlpha = 20,
    clipToCoast  = "terrestrial",
   # proj         = "+proj=longlat +datum=WGS84",
    coordHeaders = c("LONG", "LAT")
  )

  # Transform to CEA, buffer, dissolve, return to WGS84
  shape_sf   <- st_as_sf(shape[[1]])
  shapeTrans <- st_transform(shape_sf, crs = "+proj=cea +lat_ts=0 +lon_0=0")
  shape2     <- st_buffer(shapeTrans, dist = buffDist + 25000) |> st_union() |> st_as_sf()
  shape2_geo <- st_transform(shape2, crs = 4326)

  # --- Visualisation ---
  world  <- ne_countries(scale = "medium", returnclass = "sf")
  plotz  <- ggplot() +
    geom_sf(data = world,     fill = NA) +
    geom_sf(data = shape2_geo, fill = "orange", alpha = 0.5) +
    geom_sf(data = temp_sf,    colour = "black", size = 0.8) +
    coord_sf(
      xlim = c(min(temp$LONG) - 10, max(temp$LONG) + 10),
      ylim = c(min(temp$LAT)  - 10, max(temp$LAT)  + 10)
    ) +
    theme_classic()

  if (saveImage) {
    ggsave(filename = paste0(imagedir, paste0(species_df$GEN_DB[1], "_", species_df$SP_DB[1]), ".png"), plot = plotz)
  } else {
    message("Alpha hull plot not saved. Set saveImage = TRUE to save.")
  }

  if (saveShapefile) {
    st_write(shape2_geo,
             dsn          = paste0(shapefiledir, paste0(species_df$GEN_DB[1], "_", species_df$SP_DB[1]), "_AA.shp"),
             delete_layer = TRUE)
  } else {
    message("Accessible area shapefile not saved. Set saveShapefile = TRUE to save.")
  }

  return(shape2_geo)
}
