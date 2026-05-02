# =============================================================================
# 04_clip_modelLayers.R — Step 3: Clip Environmental Layers to the M Region
# =============================================================================
#
# WHAT THIS STEP DOES (conceptually)
# ------------------------------------
# MaxEnt needs environmental predictor rasters as both:
#   1. The "foreground" — values extracted at occurrence locations.
#   2. The "background" — values sampled across the study area to characterise
#      what environments are available to the species.
#
# If you pass in global rasters, MaxEnt draws background from the whole world.
# That violates the accessible-area principle (Step 2) and makes the model far
# too easy to discriminate — it will detect trivial climate gradients (e.g.,
# tropical vs Arctic) rather than the subtle environmental gradients that
# differentiate occupied from unoccupied habitat within the species' range.
#
# This step crops and masks every environmental raster to the accessible area
# polygon, so all subsequent model fitting and variable selection operates
# only within M.
#
# WHAT THE USER CAN CHANGE
# -------------------------
#   layerDir     — Path to your raster directory. All .tif files found there
#                  will be loaded and clipped. Make sure the directory contains
#                  only the variables you want to consider (or use
#                  clip_projectionLayers with varNames to subset later).
#
#   accessibleArea — The sf polygon from define_accessibleArea(). If you want
#                  to experiment with a different M definition (e.g., a country
#                  boundary or convex hull), you can substitute any sf polygon
#                  here and the rest of the pipeline is unaffected.
#
# NOTE ON RESOLUTION
# ------------------
# After clipping, the pipeline aggregates rasters by a factor of 5
# (terra::aggregate in sdm_pipeline.R). This coarsens the resolution to
# reduce computation time and smooth local noise. Increase the aggregation
# factor to speed up model fitting on large areas; decrease (or remove) it
# if you have sparse occurrences and need fine-grained predictions.
# =============================================================================

library(terra)
library(sf)

#' Crop and mask a directory of raster layers to an accessible area polygon.
#'
#' Crops to the bounding box first (fast), then masks to the exact polygon
#' boundary (precise). Any raster cell whose centre falls outside the polygon
#' is set to NA.
#'
#' @param layerDir       Path to directory containing .tif environmental rasters.
#' @param accessibleArea sf polygon (WGS84) returned by define_accessibleArea().
#' @return  terra SpatRaster with all input layers masked to the accessible area.

clip_variableLayers <- function(layerDir, accessibleArea) {

  filelist <- list.files(layerDir, full.names = TRUE, pattern = "\\.tif$")
  rstack   <- terra::rast(filelist)

  # Reproject the accessible-area polygon to match the raster's CRS before
  # cropping. This is necessary if your climate layers are in a projected CRS
  # rather than WGS84. terra::vect converts sf → SpatVector (terra's format).
  aa_proj <- st_transform(accessibleArea, crs = terra::crs(rstack))
  aa_vect <- terra::vect(aa_proj)

  r_crop <- terra::crop(rstack, aa_vect)   # trim to bounding box
  r_mask <- terra::mask(r_crop, aa_vect)   # set cells outside polygon to NA

  return(r_mask)
}
