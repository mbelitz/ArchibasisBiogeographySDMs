library(terra)
library(sf)
library(dplyr)
library(rnaturalearth)
library(ggplot2)
library(stringr)
library(dismo)   # maxent.jar via ENMeval; also used in 06 / 09
library(raster)  # bridge for dismo (dismo::maxent requires RasterStack)
library(ENMeval)

source("scripts/03_define_accessibleArea.R")
source("scripts/04_clip_modelLayers.R")
source("scripts/05_select_modelVariables.R")
source("scripts/08_save_SDMoutputs_TSS.R")
source("scripts/09_project_broader_region.R")

# World boundaries used in all maps
world <- ne_countries(scale = "medium", returnclass = "sf")

# Pre-build the projection extent once (avoids repeated ne_countries() calls)
proj_extent <- define_projectionExtent()

# ── Read occurrence data ──────────────────────────────────────────────────────

occs <- data.table::fread('data/Combined_data_AUS-PHI-IND.csv') %>% 
  filter(GEN_DB == "Archibasis")

# ── Helper: spatial thinning ──────────────────────────────────────────────────
# Returns a thinned data frame with one record per raster cell.
thin_occurrences <- function(occ_df, ref_layer, agg_factor) {
  ref_agg  <- terra::aggregate(ref_layer, agg_factor)
  cell_ids <- terra::cellFromXY(ref_agg,
                                 as.matrix(occ_df[, c("decimalLongitude",
                                                       "decimalLatitude")]))
  occ_df |>
    dplyr::mutate(cell_id = cell_ids) |>
    dplyr::group_by(cell_id) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(-cell_id)
}

# ── Main pipeline function ────────────────────────────────────────────────────

# test species is Archibasis crucigera

  cleanedOccs <- occs |>
    dplyr::filter(GEN_DB == "Archibasis", SP_DB == "crucigera") |>
    na.omit()

  # 1. Accessible area (M region) ------------------------------------------------
  aa_shp <- define_accessibleArea(species_df  = cleanedOccs,
                                   minBuff      = 75000,
                                   saveImage    = FALSE,
                                   saveShapefile = FALSE)

  # 2. Environmental variables clipped to M region --------------------------------
  mod_vars <- clip_variableLayers(layerDir      = "ClimateOnly/",
                                   accessibleArea = aa_shp)
  mod_vars <- terra::aggregate(mod_vars, 5)

  # 3. Spatial thinning (adaptive to accessible-area size) -----------------------
  area_sqkm <- as.numeric(sf::st_area(aa_shp)) / 1e6

  spp_df <- if (area_sqkm < 100000) {
    cleanedOccs
  } else if (area_sqkm < 250000) {
    thin_occurrences(cleanedOccs, mod_vars[[2]], agg_factor = 5)
  } else if (area_sqkm < 1000000) {
    thin_occurrences(cleanedOccs, mod_vars[[2]], agg_factor = 10)
  } else if (area_sqkm < 2500000) {
    thin_occurrences(cleanedOccs, mod_vars[[2]], agg_factor = 20)
  } else {
    thin_occurrences(cleanedOccs, mod_vars[[2]], agg_factor = 40)
  }

  occ_matrix <- as.matrix(spp_df[, c("decimalLongitude", "decimalLatitude")])

  # 4. Variable selection (VIF < 5) ----------------------------------------------
  max_model <- dismo::maxent(x        = raster::stack(mod_vars),
                              p        = occ_matrix,
                              progress = "text")

  predictors <- select_sdmVariables(pred_vars  = mod_vars,
                                    maxent_mod  = max_model,
                                    maxVIF      = 5)

  # 5. ENMevaluate — spatial block CV, tuned RM & FC ----------------------------
  eval1 <- ENMeval::ENMevaluate(
    occs       = occ_matrix,
    envs       = predictors,
    tune.args  = list(rm = c(0.5, 1, 2, 3, 4),
                      fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT")),
    partitions  = "block",
    algorithm   = "maxent.jar",
    parallel    = TRUE,
    numCores    = 5
  )

  # 6. Set up output directories -------------------------------------------------
  bw  <- stringr::str_replace(binomial, " ", "_")
  fp  <- "Fixed_Sdms/Nearctic_sdmOutputs_R2/"
  fp2 <- "Fixed_Sdms/Nearctic_sdmMaps_R2/"
  dir.create(file.path(fp, bw), showWarnings = FALSE, recursive = TRUE)

  # 7. Save training-area SDM outputs; get TSS threshold for projection ----------
  tss_out <- save_SDM_results(ENMeval_output = eval1,
                               AUCmin         = 0.7,
                               resultDir      = file.path(fp, bw),
                               spp            = binomial,
                               occ_df         = spp_df)

  save(eval1, file = file.path(fp, bw, paste0(bw, "_ENMeval.RData")))

  # 8. Project to AU+SEA with MaxEnt clamping + MESS ----------------------------
  proj_vars <- clip_projectionLayers(
    layerDir        = "ClimateOnly/",
    projectionExtent = proj_extent,
    varNames         = names(predictors)
  )

  project_toRegion(
    ENMeval_output  = eval1,
    training_vars   = predictors,
    occ_df          = spp_df,
    projection_vars  = proj_vars,
    tss_threshold    = tss_out$threshold,
    resultDir        = file.path(fp, bw),
    spp              = binomial,
    mapDir           = fp2
  )

  # 9. Training-area maps --------------------------------------------------------
  r  <- terra::rast(file.path(fp, bw, paste0(bw, "_SDM.tif")))
  r2 <- terra::rast(file.path(fp, bw, paste0(bw, "_SDM_PA.tif")))

  xlim <- c(min(cleanedOccs$decimalLongitude) - 5,  max(cleanedOccs$decimalLongitude) + 5)
  ylim <- c(min(cleanedOccs$decimalLatitude)  - 5,  max(cleanedOccs$decimalLatitude)  + 5)

  p <- ggplot() +
    geom_sf(data = world, fill = NA) +
    geom_sf(data = aa_shp, fill = "orange", alpha = 0.5) +
    geom_point(data = cleanedOccs,
               aes(x = decimalLongitude, y = decimalLatitude),
               shape = 1, size = 0.75) +
    coord_sf(xlim = xlim, ylim = ylim) +
    ggtitle(paste(binomial, "n =", nrow(cleanedOccs))) +
    theme_classic()

  p2 <- ggplot() +
    geom_sf(data = world, fill = NA) +
    geom_tile(data = as.data.frame(r, xy = TRUE) |> na.omit() |> dplyr::rename(ClogLog = 3),
              aes(x = x, y = y, fill = ClogLog)) +
    scale_fill_viridis_c() +
    coord_sf(xlim = xlim, ylim = ylim) +
    theme_classic()

  p3 <- ggplot() +
    geom_sf(data = world, fill = NA) +
    geom_tile(data = as.data.frame(r2, xy = TRUE) |> na.omit() |>
                dplyr::rename(ClogLog = 3) |>
                dplyr::mutate(ClogLog = as.character(ClogLog)),
              aes(x = x, y = y, fill = ClogLog)) +
    scale_fill_viridis_d() +
    labs(fill = "Presence") +
    coord_sf(xlim = xlim, ylim = ylim) +
    theme_classic()

  e <- egg::ggarrange(p2, p3, p, nrow = 2, heights = c(1, 0.65), draw = FALSE)
  ggsave(plot     = e,
         filename = file.path(fp2, paste0(bw, "_sdmMap.png")),
         width    = 8, height = 6)

  bp <- ggplot() +
    geom_sf(data = world, fill = NA) +
    geom_tile(data = as.data.frame(r2, xy = TRUE) |> na.omit() |>
                dplyr::rename(ClogLog = 3) |>
                dplyr::mutate(ClogLog = as.character(ClogLog)),
              aes(x = x, y = y, fill = ClogLog), alpha = 0.8) +
    geom_point(data = cleanedOccs,
               aes(x = decimalLongitude, y = decimalLatitude),
               size = 0.5, alpha = 0.75, shape = 1) +
    scale_fill_viridis_d() +
    labs(fill = "Presence") +
    coord_sf(xlim = c(xlim[1] - 5, xlim[2] + 5),
             ylim = c(ylim[1] - 5, ylim[2] + 5)) +
    ggtitle(paste(binomial, "n =", nrow(cleanedOccs))) +
    theme_classic()

  ggsave(plot     = bp,
         filename = file.path(fp2, paste0(bw, "_sdmMapOccPoints.png")))

  # 10. Cleanup ------------------------------------------------------------------
  rm(p, r, p2, r2, p3, e, bp)
  gc()
  files <- list.files(tempdir(), full.names = TRUE, all.files = TRUE, recursive = TRUE)
  file.remove(files)
  gc()
}
