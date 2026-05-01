library(terra)
library(sf)
library(dplyr)
library(maxnet)
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
#' @param layerDir          Directory of global .tif rasters (one per variable).
#' @param projectionExtent  sf polygon returned by define_projectionExtent().
#' @param varNames          Character vector of layer names to keep. NULL keeps all.
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
#' The best ENMeval parameter combination (fc + rm) is re-fitted using maxnet
#' (pure R) so projection works regardless of whether the original maxent.jar
#' temp files still exist. Clamping is applied via maxnet's clamp = TRUE
#' argument, restricting predictor values to the training range before
#' prediction. A MESS surface is saved alongside the continuous prediction.
#'
#' @param ENMeval_output  ENMevaluation object.
#' @param training_vars   terra SpatRaster used to fit the model (accessible area).
#' @param occ_df          Data frame with columns LONG / LAT.
#' @param projection_vars terra SpatRaster for the full projection extent.
#' @param tss_threshold   Numeric cloglog threshold from save_SDM_results().
#' @param resultDir       Output directory.
#' @param spp             Species binomial.
#' @param mapDir          Optional directory for PNG maps. NULL skips.
#' @return  Invisible named list: projection, mess, proj_pa (all SpatRasters).

project_toRegion <- function(ENMeval_output, training_vars, occ_df,
                              projection_vars, tss_threshold,
                              resultDir, spp, mapDir = NULL) {

  spp_bw <- stringr::str_replace(spp, " ", "_")

  # --- Identify best model parameters from ENMeval results ---
  bestmod <- ENMeval_output@results |> dplyr::filter(delta.AICc == 0)
  bestmod <- bestmod[1, ]

  if (is.na(bestmod$auc.train) || bestmod$auc.train < 0.7) {
    bestmod <- ENMeval_output@results |>
      dplyr::filter(auc.train == max(auc.train, na.rm = TRUE))
    bestmod <- bestmod[1, ]
  }

  # Parse tune.args string (format: "rm.0.5_fc.LQ")
  tune_str <- as.character(bestmod$tune.args)
  rm_val   <- as.numeric(sub("rm\\.([^_]+)_fc\\..*", "\\1", tune_str))
  fc_str   <- sub(".*_fc\\.(.+)",                    "\\1", tune_str)
  fc_maxnet <- tolower(fc_str)   # maxnet uses lowercase: "l","q","h","p","t"

  message(sprintf("[%s] Refitting best model: rm = %s, fc = %s", spp, rm_val, fc_str))

  # --- Refit best model with maxnet (no temp files needed) ---
  # Background sampled from the accessible area (training extent + mask)
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

  # --- Subset projection layers to training variables ---
  proj_subset <- projection_vars[[names(training_vars)]]

  missing_vars <- setdiff(names(training_vars), names(proj_subset))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "Projection raster missing variables: %s\nAvailable: %s",
      paste(missing_vars, collapse = ", "),
      paste(names(projection_vars), collapse = ", ")
    ))
  }

  # --- Continuous cloglog projection with clamping ---
  # terra::predict passes each chunk as a data frame to predict.maxnet;
  # clamp = TRUE restricts values to [training min, training max] per variable.
  proj_raster <- terra::predict(proj_subset, maxnet_mod,
                                 fun    = predict,
                                 clamp  = TRUE,
                                 type   = "cloglog",
                                 na.rm  = TRUE)
  names(proj_raster) <- "cloglog"

  terra::writeRaster(proj_raster,
                     filename  = file.path(resultDir, paste0(spp_bw, "_projection_cloglog.tif")),
                     overwrite = TRUE)

  # --- Presence-absence map using threshold from training area ---
  pa_mat  <- matrix(c(0, tss_threshold, 0, tss_threshold, 1, 1), ncol = 3, byrow = TRUE)
  proj_pa <- terra::classify(proj_raster, pa_mat)
  names(proj_pa) <- "presence"

  terra::writeRaster(proj_pa,
                     filename  = file.path(resultDir, paste0(spp_bw, "_projection_PA.tif")),
                     overwrite = TRUE)

  # --- MESS surface ---
  # Negative MESS = novel environment outside the training occurrence range.
  occ_env_mess <- terra::extract(training_vars, occ_coords) |> 
  #  dplyr::select(-ID) |> 
    na.omit()
  mess_raster  <- dismo::mess(x = raster::stack(proj_subset), v = occ_env_mess)
  names(mess_raster) <- "MESS"

  terra::writeRaster(mess_raster,
                     filename  = file.path(resultDir, paste0(spp_bw, "_MESS.tif")),
                     overwrite = TRUE)

  pct_novel <- round(100 * mean(terra::values(mess_raster) < 0, na.rm = TRUE), 1)
  message(sprintf("[%s] Projection complete. Novel-environment cells (MESS < 0): %.1f%%",
                  spp, pct_novel))

  # --- Optional maps ---
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
