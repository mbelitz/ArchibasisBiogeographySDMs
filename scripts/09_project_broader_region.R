library(terra)
library(sf)
library(dplyr)
library(dismo)   # for predict.MaxEnt
library(raster)  # bridge: dismo::predict requires RasterStack input
library(rnaturalearth)
library(ggplot2)
library(stringr)

# Countries included in the projection extent (Australia + Southeast Asia,
# spanning both sides of Wallace's Line)
.PROJ_COUNTRIES <- c(
  "Australia", "Indonesia", "Malaysia", "Philippines",
  "Papua New Guinea", "Singapore", "Brunei", "Vietnam",
  "Thailand", "Myanmar", "Cambodia", "Laos", "Timor-Leste",
  "Solomon Islands", "Vanuatu", "New Caledonia"
)

#' Build an sf polygon covering the AU+SEA projection extent.
#'
#' @return  sf MULTIPOLYGON in WGS84 representing the union of target countries.

define_projectionExtent <- function() {
  world <- ne_countries(scale = "medium", returnclass = "sf")
  proj_area <- world |>
    dplyr::filter(name %in% .PROJ_COUNTRIES |
                  admin %in% .PROJ_COUNTRIES |
                  sovereignt %in% .PROJ_COUNTRIES) |>
    sf::st_union() |>
    sf::st_as_sf()
  return(proj_area)
}

#' Clip global WorldClim (or equivalent) layers to the projection extent.
#'
#' @param layerDir        Directory of global .tif rasters (one per variable).
#' @param projectionExtent  sf polygon returned by define_projectionExtent().
#' @param varNames        Character vector of layer names to keep (subset of
#'                        all rasters in layerDir). NULL keeps all layers.
#' @return  terra SpatRaster masked to projectionExtent.

clip_projectionLayers <- function(layerDir, projectionExtent, varNames = NULL) {
  filelist    <- list.files(layerDir, full.names = TRUE, pattern = "\\.tif$")
  global_vars <- terra::rast(filelist)

  if (!is.null(varNames)) {
    global_vars <- global_vars[[varNames]]
  }

  proj_vect <- terra::vect(st_transform(projectionExtent, crs = terra::crs(global_vars)))
  r_crop    <- terra::crop(global_vars, proj_vect)
  r_mask    <- terra::mask(r_crop, proj_vect)

  return(r_mask)
}

#' Project the best MaxEnt model to the AU+SEA region.
#'
#' MaxEnt's native clamping is applied by default during dismo::predict():
#' predictor values outside the training range are clamped to [min, max]
#' before prediction, preventing unconstrained extrapolation.
#'
#' A MESS (Multivariate Environmental Similarity Surface) raster is also
#' saved alongside the prediction. Negative MESS values flag cells whose
#' environment is outside the multivariate range of the training occurrences;
#' use this as a companion layer when interpreting the projection.
#'
#' @param ENMeval_output  ENMevaluation object.
#' @param training_vars   terra SpatRaster used to fit the model (accessible area).
#' @param occ_df          Data frame with columns LONG / LAT
#'                        (thinned occurrences used in modelling).
#' @param projection_vars terra SpatRaster for the full projection extent
#'                        (must contain all layers present in training_vars).
#' @param tss_threshold   Numeric cloglog threshold value returned by
#'                        save_SDM_results() — applied to produce the PA map.
#' @param resultDir       Output directory.
#' @param spp             Species binomial.
#' @param mapDir          Optional directory for saving a PNG map. NULL skips.
#' @return  Invisible named list: projection (SpatRaster), mess (SpatRaster),
#'          proj_pa (SpatRaster).

project_toRegion <- function(ENMeval_output, training_vars, occ_df,
                              projection_vars, tss_threshold,
                              resultDir, spp, mapDir = NULL) {

  spp_bw <- stringr::str_replace(spp, " ", "_")

  # --- Identify best model ---
  bestmod     <- ENMeval_output@results |> dplyr::filter(delta.AICc == 0)
  bestmod     <- bestmod[1, ]

  # Fall back to highest AUC if best AICc model has low AUC
  if (!"auc.train" %in% names(bestmod) || is.na(bestmod$auc.train) ||
      bestmod$auc.train < 0.7) {
    bestmod <- ENMeval_output@results |>
      dplyr::filter(auc.train == max(auc.train, na.rm = TRUE))
    bestmod <- bestmod[1, ]
  }

  maxent_args <- as.character(bestmod$tune.args)

  # Use numeric-index lookup for the same robustness reason as in 08_save_SDMoutputs_TSS.R
  model_names <- as.character(names(ENMeval_output@models))
  model_idx   <- which(model_names == maxent_args)
  if (length(model_idx) == 0) {
    stop(sprintf("Model '%s' not found in @models.\nAvailable: %s",
                 maxent_args, paste(model_names, collapse = ", ")))
  }
  best_model <- ENMeval_output@models[[model_idx[1]]]

  # --- Subset projection layers to the variables used in training ---
  # Drive the subset from the model's own stored variable names (colnames of
  # @presence) rather than names(training_vars): ENMeval passes data to
  # maxent.jar via temp ASC files whose names can differ slightly from the
  # SpatRaster layer names, and dismo::predict.MaxEnt checks @presence colnames.
  model_var_names <- colnames(best_model@presence)
  proj_subset     <- projection_vars[[model_var_names]]

  missing_vars <- setdiff(model_var_names, names(proj_subset))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "Projection raster missing variables required by model: %s\nAvailable layers: %s",
      paste(missing_vars, collapse = ", "),
      paste(names(projection_vars), collapse = ", ")
    ))
  }

  # --- Continuous cloglog projection (MaxEnt native clamping enabled by default) ---
  # dismo::predict.MaxEnt runs maxent.jar in project mode; clamping is on by
  # default (doclamp=true in the saved model settings).
  # The raster::stack() call is a bridge required by dismo — it does not imply
  # continued use of the raster package for any other operations.
  proj_raster_raw <- dismo::predict(best_model, raster::stack(proj_subset))
  proj_raster     <- terra::rast(proj_raster_raw)
  names(proj_raster) <- "cloglog"

  terra::writeRaster(proj_raster,
                     filename  = file.path(resultDir, paste0(spp_bw, "_projection_cloglog.tif")),
                     overwrite = TRUE)

  # --- Presence-absence map using threshold from training area ---
  pa_mat   <- matrix(c(0, tss_threshold, 0, tss_threshold, 1, 1), ncol = 3, byrow = TRUE)
  proj_pa  <- terra::classify(proj_raster, pa_mat)
  names(proj_pa) <- "presence"

  terra::writeRaster(proj_pa,
                     filename  = file.path(resultDir, paste0(spp_bw, "_projection_PA.tif")),
                     overwrite = TRUE)

  # --- MESS surface ---
  # Reference values: environmental conditions at training occurrence localities.
  # Negative MESS = novel environment; positive = within training range.
  occ_coords <- as.matrix(occ_df[, c("LONG", "LAT")])
  occ_env    <- terra::extract(training_vars, occ_coords, ID = FALSE) |> na.omit()

  mess_raster <- terra::mess(x = proj_subset, v = occ_env)
  names(mess_raster) <- "MESS"

  terra::writeRaster(mess_raster,
                     filename  = file.path(resultDir, paste0(spp_bw, "_MESS.tif")),
                     overwrite = TRUE)

  pct_novel <- round(100 * mean(terra::values(mess_raster) < 0, na.rm = TRUE), 1)
  message(sprintf("[%s] Projection complete. Novel-environment cells (MESS < 0): %.1f%%",
                  spp, pct_novel))

  # --- Optional map ---
  if (!is.null(mapDir)) {
    world <- ne_countries(scale = "medium", returnclass = "sf")

    proj_df <- as.data.frame(proj_raster, xy = TRUE) |> na.omit()
    pa_df   <- as.data.frame(proj_pa,    xy = TRUE) |> na.omit() |>
                 dplyr::mutate(presence = as.character(presence))
    mess_df <- as.data.frame(mess_raster, xy = TRUE) |> na.omit()

    p_cont <- ggplot() +
      geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.2) +
      geom_tile(data = proj_df, aes(x = x, y = y, fill = cloglog)) +
      scale_fill_viridis_c(name = "Suitability") +
      coord_sf(xlim = range(proj_df$x), ylim = range(proj_df$y)) +
      ggtitle(paste(spp, "— projected suitability (clamped)")) +
      theme_classic()

    p_pa <- ggplot() +
      geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.2) +
      geom_tile(data = pa_df, aes(x = x, y = y, fill = presence)) +
      scale_fill_viridis_d(name = "Presence") +
      coord_sf(xlim = range(pa_df$x), ylim = range(pa_df$y)) +
      ggtitle(paste(spp, "— projected presence/absence")) +
      theme_classic()

    p_mess <- ggplot() +
      geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.2) +
      geom_tile(data = mess_df, aes(x = x, y = y, fill = MESS)) +
      scale_fill_gradient2(name = "MESS", low = "red", mid = "white", high = "blue",
                           midpoint = 0) +
      coord_sf(xlim = range(mess_df$x), ylim = range(mess_df$y)) +
      ggtitle(paste(spp, "— MESS (red = novel environment)")) +
      theme_classic()

    combined <- egg::ggarrange(p_cont, p_pa, p_mess, nrow = 2, draw = FALSE)
    ggsave(plot     = combined,
           filename = file.path(mapDir, paste0(spp_bw, "_projectionMap.png")),
           width    = 12, height = 8)
  }

  return(invisible(list(projection = proj_raster, mess = mess_raster, proj_pa = proj_pa)))
}
