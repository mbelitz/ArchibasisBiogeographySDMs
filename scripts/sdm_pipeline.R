# =============================================================================
# sdm_pipeline.R — Main Orchestration Script
# =============================================================================
#
# OVERVIEW OF THE MODELLING WORKFLOW
# ------------------------------------
# This script runs a complete MaxEnt species distribution model (SDM) for a
# single species, from raw occurrence records to a projected suitability surface.
# The numbered steps below correspond to the helper scripts in /scripts/.
#
# Step 0 — Load and filter occurrences from the combined database.
# Step 1 — Define the accessible area (M region): the geographic area the
#           species could plausibly have reached, used to sample background.
# Step 2 — Clip environmental rasters to M, coarsen resolution, and reproject
#           to equal-area CRS so background sampling and thinning are unbiased.
# Step 3 — Spatially thin occurrences to reduce autocorrelation.
# Step 4 — Fit an initial MaxEnt model for variable importance, then select
#           non-collinear predictors via VIF iteration.
# Step 5 — Tune MaxEnt via ENMevaluate: test all combinations of regularisation
#           multipliers and feature classes using spatial block cross-validation.
# Step 6 — Select the best model (AICc), threshold it with a modified TSS,
#           and save training-area rasters + variable importance.
# Step 7 — Project the best model to the AU+SEA extent using maxnet,
#           with clamping and optionally a MESS surface.
# Step 8 — Produce a four-panel summary figure (training + projection).
#
# KEY USER-TUNABLE SETTINGS (search for TUNE: throughout this file)
# -----------------------------------------------------------------
#   genus / species      — Change to run a different species.
#   proj4_aea            — Equal-area CRS for model fitting and projection.
#   minBuff              — Minimum accessible-area buffer in metres (Step 1).
#   aggregation factor   — How aggressively rasters are coarsened (Step 2).
#   thinning thresholds  — Thinning intensity by accessible-area size (Step 3).
#   maxVIF               — Collinearity threshold for variable selection (Step 4).
#   rm / fc              — Regularisation multipliers and feature classes (Step 5).
#   AUCmin               — Minimum acceptable AUC for model selection (Step 6).
#   calc_mess            — Whether to compute the MESS surface (Step 7; slow).
#   numCores             — Parallel cores for ENMevaluate (Step 5).
# =============================================================================

library(terra)
library(sf)
library(dplyr)
library(rnaturalearth)
library(ggplot2)
library(stringr)
library(cowplot)         # plot_grid() for the four-panel summary figure
library(dismo)           # maxent.jar via ENMeval; also used in 06 / 09
library(raster)          # bridge for dismo (dismo::maxent requires RasterStack)
library(ENMeval)

source("scripts/03_define_accessibleArea.R")
source("scripts/04_clip_modelLayers.R")
source("scripts/05_select_modelVariables.R")
source("scripts/08_save_SDMoutputs_TSS.R")
source("scripts/09_project_broader_region.R")

# ── Equal-area projection CRS ─────────────────────────────────────────────────
# Custom Albers Equal Area centred on the AU+SEA study region.
#   lat_1 =  7°N  }  standard parallels — areas between these are most accurate;
#   lat_2 = 36°S  }  chosen to bracket the full AU+SEA latitudinal range
#   lat_0 = 15°S     latitude of the projection origin (centre of region)
#   lon_0 = 132°E    central meridian bisecting the region longitudinally
#
# Using equal area ensures that:
#   • background cells each represent the same land area (unbiased sampling),
#   • spatial thinning removes one record per equal-area cell (not degree cell),
#   • output maps do not exaggerate high-latitude areas.
#
# TUNE: to switch to Australian Albers (EPSG:3577, published standard):
#   proj4_aea <- "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132
#                 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
proj4_aea <- paste("+proj=aea +lat_1=7 +lat_2=-36 +lat_0=-15 +lon_0=132",
                   "+x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

# ── Species identity ──────────────────────────────────────────────────────────
# TUNE: change these two lines to model a different species.
genus    <- "Archibasis"
species  <- "crucigera"
binomial <- paste(genus, species)
bw       <- paste(genus, species, sep = "_")   # filename-safe version

# ── Output directories ────────────────────────────────────────────────────────
fp  <- "SDM_outputs/"
fp2 <- "SDM_maps/"
dir.create(file.path(fp, bw), showWarnings = FALSE, recursive = TRUE)
dir.create(fp2,              showWarnings = FALSE, recursive = TRUE)

# Build world boundaries and projection extent once (reused across steps).
world       <- ne_countries(scale = "medium", returnclass = "sf")
proj_extent <- define_projectionExtent()

# ── Helper: spatial thinning ──────────────────────────────────────────────────
# Returns a thinned data frame with at most one occurrence per raster cell.
# agg_factor controls the cell size used for thinning: larger = more aggressive.
# If ref_layer is in a projected (non-longlat) CRS, occurrence LONG/LAT
# coordinates are automatically transformed before cell assignment.
thin_occurrences <- function(occ_df, ref_layer, agg_factor) {
  ref_agg <- terra::aggregate(ref_layer, agg_factor)
  coords  <- as.matrix(occ_df[, c("LONG", "LAT")])

  # Project coordinates to match the raster CRS when it is not longlat.
  # terra::is.lonlat() returns FALSE for any projected (metre-based) CRS.
  if (!terra::is.lonlat(ref_layer)) {
    coords <- terra::project(coords,
                              from = "EPSG:4326",
                              to   = terra::crs(ref_layer, proj = TRUE))
  }

  cell_ids <- terra::cellFromXY(ref_agg, coords)
  occ_df |>
    dplyr::mutate(cell_id = cell_ids) |>
    dplyr::group_by(cell_id) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(-cell_id)
}

# =============================================================================
# Step 0: Load and filter occurrence records
# =============================================================================

occs <- data.table::fread("data/Combined_data_AUS-PHI-IND.csv")

cleanedOccs <- occs |>
  dplyr::filter(GEN_DB == genus, SP_DB == species) |>
  dplyr::filter(!is.na(LONG), !is.na(LAT))

message(sprintf("Records for %s: %d", binomial, nrow(cleanedOccs)))

# =============================================================================
# Step 1: Define the accessible area (M region)
# =============================================================================
# TUNE: minBuff — minimum buffer radius in metres around the alpha hull.

aa_shp <- define_accessibleArea(species_df    = cleanedOccs,
                                 minBuff       = 75000,
                                 saveImage     = FALSE,
                                 saveShapefile = FALSE)

# =============================================================================
# Step 2: Clip, coarsen, and reproject environmental layers
# =============================================================================
# Rasters are clipped to M in their native CRS, coarsened by aggregation,
# then reprojected to the equal-area CRS defined above.
#
# Reprojection happens AFTER aggregation so that terra::project() operates on
# fewer cells (faster) and the aggregation smoothing precedes reprojection.
# method = "bilinear" is appropriate for continuous climate variables; use
# method = "near" for categorical layers.
#
# TUNE: the aggregation factor (currently 5).

mod_vars <- clip_variableLayers(layerDir       = "ClimateOnly/",
                                 accessibleArea = aa_shp)
mod_vars <- terra::aggregate(mod_vars, 5)
mod_vars <- terra::project(mod_vars, proj4_aea, method = "bilinear")

# =============================================================================
# Step 3: Spatial thinning (adaptive to accessible-area size)
# =============================================================================
# TUNE: area thresholds (km²) and agg_factor values.

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

message(sprintf("Records after thinning: %d", nrow(spp_df)))

# Project occurrence coordinates from WGS84 to the equal-area CRS.
# All terra::extract() calls on projected rasters require coordinates in the
# same CRS as the raster — these are passed to dismo, ENMevaluate, and maxnet.
occ_matrix <- terra::project(
  as.matrix(spp_df[, c("LONG", "LAT")]),
  from = "EPSG:4326",
  to   = proj4_aea
)

# =============================================================================
# Step 4: Variable selection (VIF < 5)
# =============================================================================
# TUNE: maxVIF in select_sdmVariables().

max_model <- dismo::maxent(x        = raster::stack(mod_vars),
                            p        = occ_matrix,
                            progress = "text")

predictors <- select_sdmVariables(pred_vars  = mod_vars,
                                   maxent_mod = max_model,
                                   maxVIF     = 5)

message(sprintf("Variables retained after VIF selection: %s",
                paste(names(predictors), collapse = ", ")))

# =============================================================================
# Step 5: Model tuning via ENMevaluate (spatial block cross-validation)
# =============================================================================
# TUNE: rm, fc, numCores.

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

# =============================================================================
# Step 6: Save training-area SDM outputs and derive the TSS threshold
# =============================================================================
# TUNE: AUCmin.

tss_out <- save_SDM_results(ENMeval_output = eval1,
                             AUCmin         = 0.7,
                             resultDir      = file.path(fp, bw),
                             spp            = binomial,
                             occ_df         = spp_df)

save(eval1, file = file.path(fp, bw, paste0(bw, "_ENMeval.RData")))

# =============================================================================
# Step 7: Clip, reproject, and project to AU+SEA
# =============================================================================
# proj_vars is clipped in native CRS then reprojected to proj4_aea — the same
# CRS as the training rasters — so the projection step uses a consistent grid.
#
# mapDir = NULL suppresses the per-step PNG; Step 8 builds the figure instead.
# TUNE: calc_mess = FALSE to skip MESS and save computation time.

proj_vars <- clip_projectionLayers(
  layerDir         = "ClimateOnly/",
  projectionExtent = proj_extent,
  varNames         = names(predictors)
)
proj_vars <- terra::project(proj_vars, proj4_aea, method = "bilinear")

proj_out <- project_toRegion(
  ENMeval_output  = eval1,
  training_vars   = predictors,
  occ_df          = spp_df,
  projection_vars = proj_vars,
  tss_threshold   = tss_out$threshold,
  resultDir       = file.path(fp, bw),
  spp             = binomial,
  mapDir          = NULL,
  calc_mess       = TRUE
)

# =============================================================================
# Step 8: Four-panel summary figure
# =============================================================================
# All spatial objects are projected to proj4_aea for a consistent equal-area
# display. Occurrence points and the world basemap are transformed here;
# raster data frames already carry AEA coordinates from terra::rast().
#
# Panel layout:
#   A (top-left)  — Training-area suitability, zoomed to M region.
#                   Orange outline = accessible area; circles = occurrences.
#   B (top-right) — Projected suitability across AU+SEA.
#   C (bot-left)  — Projected binary presence/absence (TSS threshold).
#   D (bot-right) — MESS surface; blank if calc_mess = FALSE above.
#
# TUNE: pad_m (padding around the training-area panel in metres);
#       ggsave width / height for figure dimensions.

r_train  <- terra::rast(file.path(fp, bw, paste0(bw, "_SDM.tif")))
train_df <- as.data.frame(r_train, xy = TRUE) |> na.omit() |>
              dplyr::rename(ClogLog = 3)

# Project world basemap and accessible-area polygon to equal-area CRS.
world_aea  <- sf::st_transform(world, crs = proj4_aea)
aa_shp_aea <- sf::st_transform(aa_shp, crs = proj4_aea)

# Project occurrence points for the training-area panel.
occ_aea_mat <- terra::project(
  as.matrix(cleanedOccs[, c("LONG", "LAT")]),
  from = "EPSG:4326",
  to   = proj4_aea
)
occ_aea_df <- data.frame(X = occ_aea_mat[, 1], Y = occ_aea_mat[, 2])

# Bounding box of the accessible area in AEA metres; add padding on each side.
aa_bbox <- sf::st_bbox(aa_shp_aea)
pad_m   <- 500000   # 500 km padding — reduce for tightly distributed species
xlim_aa <- c(aa_bbox["xmin"] - pad_m, aa_bbox["xmax"] + pad_m)
ylim_aa <- c(aa_bbox["ymin"] - pad_m, aa_bbox["ymax"] + pad_m)

# Panel A: training cloglog + accessible-area outline + occurrences
p_tl <- ggplot() +
  geom_sf(data = world_aea, fill = "grey90", colour = "grey60", linewidth = 0.2) +
  geom_tile(data = train_df, aes(x = x, y = y, fill = ClogLog)) +
  scale_fill_viridis_c(name = "Suitability") +
  geom_sf(data = aa_shp_aea, fill = NA, colour = "orange", linewidth = 0.6) +
  geom_point(data = occ_aea_df, aes(x = X, y = Y),
             shape = 1, size = 0.75, colour = "black") +
  coord_sf(crs = proj4_aea, xlim = xlim_aa, ylim = ylim_aa) +
  ggtitle(paste(binomial, "— training area")) +
  theme_classic()

# Panel B: projected suitability (full AU+SEA extent)
proj_df <- as.data.frame(proj_out$projection, xy = TRUE) |> na.omit()

p_tr <- ggplot() +
  geom_sf(data = world_aea, fill = "grey90", colour = "grey60", linewidth = 0.2) +
  geom_tile(data = proj_df, aes(x = x, y = y, fill = cloglog)) +
  scale_fill_viridis_c(name = "Suitability") +
  coord_sf(crs = proj4_aea,
           xlim = range(proj_df$x), ylim = range(proj_df$y)) +
  ggtitle("Projected suitability (AU+SEA)") +
  theme_classic()

# Panel C: projected binary presence/absence
pa_df <- as.data.frame(proj_out$proj_pa, xy = TRUE) |> na.omit() |>
           dplyr::mutate(presence = as.character(presence))

p_bl <- ggplot() +
  geom_sf(data = world_aea, fill = "grey90", colour = "grey60", linewidth = 0.2) +
  geom_tile(data = pa_df, aes(x = x, y = y, fill = presence)) +
  scale_fill_viridis_d(name = "Presence") +
  coord_sf(crs = proj4_aea,
           xlim = range(pa_df$x), ylim = range(pa_df$y)) +
  ggtitle("Projected presence/absence (AU+SEA)") +
  theme_classic()

# Panel D: MESS surface — NULL panel if calc_mess = FALSE
if (!is.null(proj_out$mess)) {
  mess_df <- as.data.frame(proj_out$mess, xy = TRUE) |> na.omit()
  p_br <- ggplot() +
    geom_sf(data = world_aea, fill = "grey90", colour = "grey60", linewidth = 0.2) +
    geom_tile(data = mess_df, aes(x = x, y = y, fill = MESS)) +
    scale_fill_gradient2(name = "MESS", low = "red", mid = "white",
                         high = "blue", midpoint = 0) +
    coord_sf(crs = proj4_aea,
             xlim = range(mess_df$x), ylim = range(mess_df$y)) +
    ggtitle("MESS (red = novel environment)") +
    theme_classic()
} else {
  p_br <- NULL
}

panel_fig <- cowplot::plot_grid(p_tl, p_tr, p_bl, p_br,
                                 ncol   = 2, nrow = 2,
                                 labels = c("A", "B", "C", "D"))

ggsave(plot     = panel_fig,
       filename = file.path(fp2, paste0(bw, "_sdmMap_4panel.png")),
       width    = 14, height = 10)

# =============================================================================
# Step 9: Cleanup temporary files
# =============================================================================

rm(p_tl, p_tr, p_bl, p_br, panel_fig, proj_df, pa_df, train_df,
   occ_aea_df, occ_aea_mat, world_aea, aa_shp_aea)
gc()
files <- list.files(tempdir(), full.names = TRUE, all.files = TRUE, recursive = TRUE)
file.remove(files)
gc()
