# =============================================================================
# 09_project_broader_region.R — Step 7: Project Model to AU + Southeast Asia
# =============================================================================
#
# WHAT THIS STEP DOES (conceptually)
# ------------------------------------
# The MaxEnt model was trained on the accessible area (M region). That gives
# us a habitat-suitability surface for the species' *known* range. The next
# question is: "Where else in Australia and Southeast Asia might this species
# be able to survive if it could get there (or if it is already there but
# undetected)?"
#
# This step re-fits the best model parameters using maxnet (a pure-R MaxEnt
# implementation) and projects the model onto a broader environmental canvas
# covering AU + SEA. Three outputs are produced:
#
#   1. Continuous suitability raster (cloglog scale, 0–1).
#   2. Binary presence/absence raster using the TSS threshold from step 6.
#   3. MESS surface (optional) — identifies where the projection environment
#      falls *outside* the range seen during training (novel environments).
#
# WHY REFIT WITH maxnet?
# ----------------------
# ENMeval uses maxent.jar (Java), which writes temporary files. After the
# script finishes, those temp files are deleted. To project later (or on a
# different machine), we'd need to re-run ENMeval. Instead, we re-fit using
# maxnet (Elith et al. 2017), which is a pure-R reimplementation of MaxEnt
# that produces identical predictions. We simply parse the best rm and fc
# from the ENMeval results table and pass them to maxnet.
#
# CLAMPING
# --------
# When projecting to a new region, some cells may have predictor values
# outside the range seen during training (e.g., drier than anywhere in M).
# Without clamping, MaxEnt extrapolates using response curves fitted on
# training data, which can produce artefactual predictions (high or low)
# far outside the training envelope. With clamp = TRUE, predictor values
# are constrained to [training min, training max] before prediction —
# the model "plateaus" at its edge rather than extrapolating.
#   TUNE: set clamp = FALSE if you specifically want to explore extrapolation
#         behaviour, but interpret results cautiously.
#
# MESS (Multivariate Environmental Similarity Surface)
# -----------------------------------------------------
# MESS quantifies, for each projection cell, how similar the environmental
# conditions are to the range of conditions at occurrence locations.
#   • Positive MESS: conditions are within the training envelope → reliable
#     predictions.
#   • Negative MESS: at least one variable is outside the training range →
#     novel environment → model is extrapolating, predictions are uncertain.
# MESS is computationally intensive (scales with raster size × number of
# predictors). Use calc_mess = FALSE during exploratory runs to skip it.
#
# WHAT THE USER CAN CHANGE
# -------------------------
#   .PROJ_COUNTRIES  — The list of countries defining the projection extent.
#                      Add or remove countries to expand/shrink the region.
#                      Uses rnaturalearth `name`, `admin`, or `sovereignt`.
#
#   bg sample size   — terra::spatSample(size = 10000) draws 10 000 background
#                      points for refitting the maxnet model. Increase to 50000
#                      for more stable background characterisation (slower);
#                      decrease for quick exploratory runs.
#
#   clamp            — TRUE (default) restricts predictor values to training
#                      range; FALSE allows extrapolation (see above).
#
#   calc_mess        — TRUE (default) computes and saves the MESS surface.
#                      Set to FALSE to skip MESS and save substantial time on
#                      large projection extents.
#
#   mapDir           — If not NULL, saves PNG maps to this directory. Set to
#                      NULL to suppress plot output (e.g. when the pipeline
#                      builds its own composite figure).
# =============================================================================

library(terra)
library(sf)
library(dplyr)
library(maxnet)
library(rnaturalearth)
library(ggplot2)
library(stringr)

# Countries included in the projection extent (Australia + Southeast Asia,
# spanning both sides of Wallace's Line). Edit this vector to change the
# projection region — names must match rnaturalearth's name/admin/sovereignt.
.PROJ_COUNTRIES <- c(
  "Australia", "Indonesia", "Malaysia", "Philippines",
  "Papua New Guinea", "Singapore", "Brunei", "Vietnam",
  "Thailand", "Myanmar", "Cambodia", "Laos", "Timor-Leste",
  "Solomon Islands", "Vanuatu", "New Caledonia"
)

#' Build an sf polygon covering the AU+SEA projection extent.
#'
#' Unions the boundaries of all countries in .PROJ_COUNTRIES into a single
#' multipolygon. This polygon is used to clip the projection rasters and
#' set map extents.
#'
#' @return  sf MULTIPOLYGON in WGS84 representing the union of target countries.

define_projectionExtent <- function() {
  world     <- ne_countries(scale = "medium", returnclass = "sf")
  proj_area <- world |>
    dplyr::filter(name      %in% .PROJ_COUNTRIES |
                  admin     %in% .PROJ_COUNTRIES |
                  sovereignt %in% .PROJ_COUNTRIES) |>
    sf::st_union() |>
    sf::st_as_sf()
  return(proj_area)
}

#' Clip global WorldClim (or equivalent) layers to the projection extent.
#'
#' Called in sdm_pipeline.R to prepare projection-extent rasters that contain
#' only the variables retained after VIF selection (varNames).
#'
#' @param layerDir          Directory of global .tif rasters (one per variable).
#' @param projectionExtent  sf polygon returned by define_projectionExtent().
#' @param varNames          Character vector of layer names to keep.
#'                          NULL keeps all layers (use when no selection needed).
#' @return  terra SpatRaster masked to projectionExtent.

clip_projectionLayers <- function(layerDir, projectionExtent, varNames = NULL) {
  filelist    <- list.files(layerDir, full.names = TRUE, pattern = "\\.tif$")
  global_vars <- terra::rast(filelist)

  if (!is.null(varNames)) {
    global_vars <- global_vars[[varNames]]
  }

  proj_vect <- terra::vect(st_transform(projectionExtent,
                                         crs = terra::crs(global_vars)))
  r_crop    <- terra::crop(global_vars, proj_vect)
  r_mask    <- terra::mask(r_crop, proj_vect)

  return(r_mask)
}

#' Project the best MaxEnt model to the AU+SEA region.
#'
#' Refits the best ENMeval parameters using maxnet, projects onto the full
#' AU+SEA raster with clamping, and optionally computes a MESS surface to
#' flag novel environments.
#'
#' @param ENMeval_output  ENMevaluation object from step 5.
#' @param training_vars   terra SpatRaster of variables used to fit the model
#'                        (accessible area extent, after VIF selection).
#' @param occ_df          Data frame with columns LONG / LAT (occurrence records).
#' @param projection_vars terra SpatRaster for the full projection extent
#'                        (same variable names as training_vars).
#' @param tss_threshold   Numeric cloglog threshold from save_SDM_results().
#'                        Cells with suitability >= threshold → predicted present.
#' @param resultDir       Directory for output rasters and CSVs.
#' @param spp             Species binomial string (spaces allowed).
#' @param mapDir          Directory for optional PNG maps. NULL suppresses output.
#' @param calc_mess       Logical. Compute and save the MESS surface? Default TRUE.
#'                        Set to FALSE to skip MESS and save computation time.
#' @return  Invisible named list with elements:
#'            projection — continuous cloglog SpatRaster (projection extent)
#'            mess       — MESS SpatRaster, or NULL if calc_mess = FALSE
#'            proj_pa    — binary presence/absence SpatRaster (projection extent)

project_toRegion <- function(ENMeval_output, training_vars, occ_df,
                              projection_vars, tss_threshold,
                              resultDir, spp, mapDir = NULL,
                              calc_mess = TRUE) {

  spp_bw <- stringr::str_replace(spp, " ", "_")

  # ── Identify best model parameters from ENMeval results ─────────────────────
  # Primary criterion: lowest AICc (delta.AICc == 0).
  # Fallback: highest training AUC (used if AICc-best model has AUC < 0.7,
  # which can happen when cross-validation blocks are very unbalanced).
  bestmod <- ENMeval_output@results |> dplyr::filter(delta.AICc == 0)
  bestmod <- bestmod[1, ]

  if (is.na(bestmod$auc.train) || bestmod$auc.train < 0.7) {
    bestmod <- ENMeval_output@results |>
      dplyr::filter(auc.train == max(auc.train, na.rm = TRUE))
    bestmod <- bestmod[1, ]
  }

  # Parse the tune.args string (format: "rm.0.5_fc.LQ") to extract RM and FC.
  tune_str  <- as.character(bestmod$tune.args)
  rm_val    <- as.numeric(sub("rm\\.([^_]+)_fc\\..*", "\\1", tune_str))
  fc_str    <- sub(".*_fc\\.(.+)", "\\1", tune_str)
  fc_maxnet <- tolower(fc_str)   # maxnet uses lowercase feature class strings

  message(sprintf("[%s] Refitting best model: rm = %s, fc = %s", spp, rm_val, fc_str))

  # ── Refit best model with maxnet ─────────────────────────────────────────────
  # Extract environmental values at occurrence and background locations.
  # Background is sampled from the training (accessible-area) rasters, not
  # the projection extent — the model must "learn" from within M only.
  # TUNE: increase bg sample size (currently 10000) for more stable fits on
  #       large accessible areas; the cost is proportionally longer fitting.
  occ_coords <- as.matrix(occ_df[, c("LONG", "LAT")])
  occ_env    <- terra::extract(training_vars, occ_coords) |> na.omit()

  bg_env <- terra::spatSample(training_vars, size = 10000,
                               method = "random", na.rm = TRUE,
                               values = TRUE, xy = FALSE) |>
    as.data.frame()

  p    <- c(rep(1L, nrow(occ_env)), rep(0L, nrow(bg_env)))
  data <- rbind(occ_env, bg_env)

  f          <- maxnet::maxnet.formula(p, data, classes = fc_maxnet)
  maxnet_mod <- maxnet::maxnet(p, data, f, regmult = rm_val)

  # ── Subset projection layers to training variable names ──────────────────────
  # Ensure the projection raster contains exactly the variables the model was
  # trained on — no more, no less. Stops with an informative error if a
  # variable is missing from the projection rasters.
  proj_subset  <- projection_vars[[names(training_vars)]]
  missing_vars <- setdiff(names(training_vars), names(proj_subset))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "Projection raster missing variables: %s\nAvailable: %s",
      paste(missing_vars, collapse = ", "),
      paste(names(projection_vars), collapse = ", ")
    ))
  }

  # ── Continuous cloglog projection with clamping ──────────────────────────────
  # terra::predict() passes raster chunks as data frames to predict.maxnet().
  # clamp = TRUE restricts each predictor to [training min, training max]
  # before prediction, preventing extrapolation artefacts in novel environments.
  # type = "cloglog" scales output to a 0–1 probability-like suitability index.
  proj_raster <- terra::predict(proj_subset, maxnet_mod,
                                 fun   = predict,
                                 clamp = TRUE,
                                 type  = "cloglog",
                                 na.rm = TRUE)
  names(proj_raster) <- "cloglog"

  terra::writeRaster(proj_raster,
                     filename  = file.path(resultDir,
                                           paste0(spp_bw, "_projection_cloglog.tif")),
                     overwrite = TRUE)

  # ── Presence-absence map using TSS threshold from training area ───────────────
  # Apply the same cloglog threshold derived in step 6 (save_SDM_results) to
  # classify each projected cell as predicted-present (1) or predicted-absent (0).
  pa_mat  <- matrix(c(0, tss_threshold, 0, tss_threshold, 1, 1), ncol = 3, byrow = TRUE)
  proj_pa <- terra::classify(proj_raster, pa_mat)
  names(proj_pa) <- "presence"

  terra::writeRaster(proj_pa,
                     filename  = file.path(resultDir,
                                           paste0(spp_bw, "_projection_PA.tif")),
                     overwrite = TRUE)

  # ── MESS surface (optional — can be skipped to save time) ────────────────────
  # dismo::mess() computes, for each cell in the projection extent, the minimum
  # similarity score across all predictors relative to the distribution of values
  # at occurrence points (v). Negative = novel; positive = interpolation.
  # Set calc_mess = FALSE to skip; mess_raster will be NULL in the return value.
  if (calc_mess) {
    occ_env_mess <- terra::extract(training_vars, occ_coords) |>
      na.omit()
    mess_raster  <- dismo::mess(x = raster::stack(proj_subset), v = occ_env_mess)
    names(mess_raster) <- "MESS"

    terra::writeRaster(mess_raster,
                       filename  = file.path(resultDir,
                                             paste0(spp_bw, "_MESS.tif")),
                       overwrite = TRUE)

    pct_novel <- round(100 * mean(terra::values(mess_raster) < 0, na.rm = TRUE), 1)
    message(sprintf(
      "[%s] Projection complete. Novel-environment cells (MESS < 0): %.1f%%",
      spp, pct_novel
    ))
  } else {
    mess_raster <- NULL
    message(sprintf(
      "[%s] Projection complete. MESS calculation skipped (calc_mess = FALSE).",
      spp
    ))
  }

  # ── Optional maps ─────────────────────────────────────────────────────────────
  # When mapDir is provided, saves a PNG with continuous, PA, and (if available)
  # MESS panels. In sdm_pipeline.R the four-panel composite figure is built
  # separately; pass mapDir = NULL there to avoid duplicate PNG output.
  if (!is.null(mapDir)) {
    world   <- ne_countries(scale = "medium", returnclass = "sf")
    proj_df <- as.data.frame(proj_raster, xy = TRUE) |> na.omit()
    pa_df   <- as.data.frame(proj_pa,    xy = TRUE) |> na.omit() |>
                 dplyr::mutate(presence = as.character(presence))

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

    if (!is.null(mess_raster)) {
      mess_df <- as.data.frame(mess_raster, xy = TRUE) |> na.omit()
      p_mess  <- ggplot() +
        geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.2) +
        geom_tile(data = mess_df, aes(x = x, y = y, fill = MESS)) +
        scale_fill_gradient2(name = "MESS", low = "red", mid = "white",
                             high = "blue", midpoint = 0) +
        coord_sf(xlim = range(mess_df$x), ylim = range(mess_df$y)) +
        ggtitle(paste(spp, "— MESS (red = novel environment)")) +
        theme_classic()
      combined <- egg::ggarrange(p_cont, p_pa, p_mess, nrow = 2, draw = FALSE)
    } else {
      combined <- egg::ggarrange(p_cont, p_pa, nrow = 1, draw = FALSE)
    }

    ggsave(plot     = combined,
           filename = file.path(mapDir, paste0(spp_bw, "_projectionMap.png")),
           width    = 12, height = 8)
  }

  return(invisible(list(projection = proj_raster,
                        mess       = mess_raster,
                        proj_pa    = proj_pa)))
}
