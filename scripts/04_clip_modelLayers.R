library(terra)
library(sf)

#' Crop and mask a directory of raster layers to an accessible area polygon.
#'
#' @param layerDir   Path to directory containing .tif environmental rasters.
#' @param accessibleArea  sf polygon (WGS84) returned by define_accessibleArea().
#' @return  terra SpatRaster masked to the accessible area.

clip_variableLayers <- function(layerDir, accessibleArea) {

  filelist <- list.files(layerDir, full.names = TRUE, pattern = "\\.tif$")
  rstack   <- terra::rast(filelist)

  # Reproject accessible area to match raster CRS, then convert to SpatVector
  aa_proj <- st_transform(accessibleArea, crs = terra::crs(rstack))
  aa_vect <- terra::vect(aa_proj)

  r_crop <- terra::crop(rstack, aa_vect)
  r_mask <- terra::mask(r_crop, aa_vect)

  return(r_mask)
}
