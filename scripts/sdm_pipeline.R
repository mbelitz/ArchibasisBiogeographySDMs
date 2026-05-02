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
# Step 2 — Clip environmental rasters to M and coarsen resolution.
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
#   minBuff              — Minimum accessible-area buffer in metres (Step 1).
#   agg_factor thresholds — How aggressively rasters are coarsened by area (Step 2).
#   thin_occurrences     — Spatial thinning intensity (Step 3).
#   maxVIF               — Collinearity threshold for variable selection (Step 4).
#   rm / fc              — Regularisation multipliers and feature classes to tune (Step 5).
#   AUCmin               — Minimum acceptable AUC for best-model selection (Step 6).
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

# ── Species identity ──────────────────────────────────────────────────────────
# TUNE: change these two lines to model a different species.
# The data file (below) must contain matching GEN_DB / SP_DB entries.
genus   <- "Archibasis"
species <- "crucigera"
binomial <- paste(genus, species)
bw       <- paste(genus, species, sep = "_")   # filename-safe version

# ── Output directories ────────────────────────────────────────────────────────
fp  <- "SDM_outputs/"     # rasters, CSVs, and RData objects
fp2 <- "SDM_maps/"        # PNG figures
dir.create(file.path(fp, bw), showWarnings = FALSE, recursive = TRUE)
dir.create(fp2,              showWarnings = FALSE, recursive = TRUE)

# Build world boundaries and projection extent once (reused across steps).
world       <- ne_countries(scale = "medium", returnclass = "sf")
proj_extent <- define_projectionExtent()

# ── Helper: spatial thinning ──────────────────────────────────────────────────
# Returns a thinned data frame with at most one occurrence per raster cell.
# agg_factor controls the cell size used for thinning: larger values = coarser
# grid = more aggressive thinning. This reduces spatial autocorrelation among
# training points, which can inflate cross-validation AUC.
thin_occurrences <- function(occ_df, ref_layer, agg_factor) {
  ref_agg  <- terra::aggregate(ref_layer, agg_factor)
  cell_ids <- terra::cellFromXY(ref_agg, as.matrix(occ_df[, c("LONG", "LAT")]))
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
# The CSV contains records for all species in the study; filter to the target.
# LONG and LAT are required — rows with missing coordinates are dropped.

occs <- data.table::fread("data/Combined_data_AUS-PHI-IND.csv")

cleanedOccs <- occs |>
  dplyr::filter(GEN_DB == genus, SP_DB == species) |>
  dplyr::filter(!is.na(LONG), !is.na(LAT))

message(sprintf("Records for %s: %d", binomial, nrow(cleanedOccs)))

# =============================================================================
# Step 1: Define the accessible area (M region)
# =============================================================================
# See 03_define_accessibleArea.R for full documentation.
#
# TUNE: minBuff — minimum buffer in metres around the alpha hull.
#   Increase (e.g. 150000) for a wider M region that samples more
#   background climate; decrease for range-restricted species where
#   a tight M is more defensible ecologically.

aa_shp <- define_accessibleArea(species_df    = cleanedOccs,
                                 minBuff       = 75000,
                                 saveImage     = FALSE,
                                 saveShapefile = FALSE)

# =============================================================================
# Step 2: Clip and coarsen environmental layers
# =============================================================================
# Rasters are clipped to M (Step 1) and then aggregated to reduce computation
# time. terra::aggregate(, fact = 5) increases cell size by a factor of 5
# (e.g. 1 km → 5 km; 2.5 arcmin → ~12.5 arcmin).
#
# TUNE: the aggregation factor (currently 5).
#   Reduce to 2–3 if you have many records and need finer resolution.
#   Increase to 10+ for very large accessible areas or slow machines.

mod_vars <- clip_variableLayers(layerDir       = "ClimateOnly/",
                                 accessibleArea = aa_shp)
mod_vars <- terra::aggregate(mod_vars, 5)

# =============================================================================
# Step 3: Spatial thinning (adaptive to accessible-area size)
# =============================================================================
# Occurrences clustered in small areas give MaxEnt an inflated signal for
# that habitat type. Thinning to one record per grid cell reduces this bias.
# The thinning intensity scales with the accessible area (in km²) to keep
# the model tractable and avoid over-thinning small ranges.
#
# TUNE: the area thresholds (e.g. 100000, 250000, …) and the agg_factor
#   values (5, 10, 20, 40) to match the resolution of your study extent.
#   For a species with very few records (< 30), consider skipping thinning
#   (use cleanedOccs directly) to avoid losing too many points.

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

occ_matrix <- as.matrix(spp_df[, c("LONG", "LAT")])

# =============================================================================
# Step 4: Variable selection (VIF < 5)
# =============================================================================
# First an initial MaxEnt model is fitted on all variables to extract
# permutation importances. These importances are used as a tiebreaker in the
# VIF loop: when two variables are similarly collinear, the less-important one
# is removed. See 05_select_modelVariables.R for full documentation.
#
# TUNE: maxVIF in select_sdmVariables().
#   5 is the standard threshold. Lower (e.g. 3) produces a leaner set;
#   higher (e.g. 10) retains more variables but risks collinearity artefacts.

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
# Tests all combinations of regularisation multiplier (rm) and feature class
# (fc). Models are ranked by AICc. See 06_generate_ENMevals.R for full
# documentation of what rm and fc control ecologically/statistically.
#
# TUNE: rm (regularisation multiplier values) and fc (feature class strings).
#   Add finer rm steps (e.g. 1.5) for more resolution around the optimum.
#   Remove "LQHP" and "LQHPT" for small datasets (< 50 records) to avoid
#   overfitting from product and threshold features.
#
# TUNE: numCores — set to (number of physical CPU cores − 1) on your machine.
#   On Windows, parallel = TRUE requires the doParallel package.

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
# Selects the best model, thresholds it via modified TSS, and saves:
#   *_SDM.tif          — continuous cloglog raster (training extent)
#   *_SDM_PA.tif       — binary presence/absence raster (training extent)
#   *_bestModel.csv    — model selection summary
#   *_variableImportance.csv
#
# Returns: tss_out$threshold (used in Step 7 to threshold the projection).
# See 08_save_SDMoutputs_TSS.R for full documentation of AUCmin and TSS.
#
# TUNE: AUCmin — minimum acceptable AUC for the best AICc model.
#   Default 0.7. Raise to 0.8 for stricter quality control; lower to 0.6
#   for data-sparse species where cross-validated AUC is inherently noisy.

tss_out <- save_SDM_results(ENMeval_output = eval1,
                             AUCmin         = 0.7,
                             resultDir      = file.path(fp, bw),
                             spp            = binomial,
                             occ_df         = spp_df)

save(eval1, file = file.path(fp, bw, paste0(bw, "_ENMeval.RData")))

# =============================================================================
# Step 7: Project to AU+SEA with clamping and optional MESS
# =============================================================================
# Clips global climate rasters to the AU+SEA projection extent, refits the
# best model with maxnet (no Java temp files needed), projects with clamping,
# and optionally computes a MESS surface.
#
# mapDir = NULL suppresses per-step PNG output; the four-panel figure (Step 8)
# is the primary visual summary and is built from the returned rasters.
#
# TUNE: calc_mess — set to FALSE to skip MESS computation and save time on
#   large projection extents. The MESS panel in the Step 8 figure will be
#   omitted automatically if calc_mess = FALSE.

proj_vars <- clip_projectionLayers(
  layerDir         = "ClimateOnly/",
  projectionExtent = proj_extent,
  varNames         = names(predictors)
)

proj_out <- project_toRegion(
  ENMeval_output  = eval1,
  training_vars   = predictors,
  occ_df          = spp_df,
  projection_vars = proj_vars,
  tss_threshold   = tss_out$threshold,
  resultDir       = file.path(fp, bw),
  spp             = binomial,
  mapDir          = NULL,      # suppress per-step PNG; Step 8 builds the figure
  calc_mess       = TRUE
)

# =============================================================================
# Step 8: Four-panel summary figure
# =============================================================================
# Combines the key model outputs into one publication-ready figure:
#
#   Top-left  (A): Training-area suitability (cloglog), zoomed to M region.
#                  Accessible-area boundary shown as an orange outline.
#                  Occurrence points overlaid.
#   Top-right (B): Projected suitability across AU+SEA (cloglog).
#   Bot-left  (C): Projected binary presence/absence map (AU+SEA),
#                  thresholded by the TSS value from Step 6.
#   Bot-right (D): MESS surface (AU+SEA) — red cells are novel environments
#                  where the model is extrapolating; interpret with caution.
#                  Omitted (blank panel) if calc_mess = FALSE above.
#
# TUNE: ggsave width / height below to change figure dimensions.

r_train   <- terra::rast(file.path(fp, bw, paste0(bw, "_SDM.tif")))
train_df  <- as.data.frame(r_train, xy = TRUE) |> na.omit() |>
               dplyr::rename(ClogLog = 3)

xlim_aa   <- c(min(cleanedOccs$LONG) - 5, max(cleanedOccs$LONG) + 5)
ylim_aa   <- c(min(cleanedOccs$LAT)  - 5, max(cleanedOccs$LAT)  + 5)

# Panel A: training cloglog + accessible area outline + occurrences
p_tl <- ggplot() +
  geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.2) +
  geom_tile(data = train_df, aes(x = x, y = y, fill = ClogLog)) +
  scale_fill_viridis_c(name = "Suitability") +
  geom_sf(data = aa_shp, fill = NA, colour = "orange", linewidth = 0.6) +
  geom_point(data = cleanedOccs, aes(x = LONG, y = LAT),
             shape = 1, size = 0.75, colour = "black") +
  coord_sf(xlim = xlim_aa, ylim = ylim_aa) +
  ggtitle(paste(binomial, "— training area")) +
  theme_classic()

# Panel B: projected suitability (full AU+SEA extent)
proj_df <- as.data.frame(proj_out$projection, xy = TRUE) |> na.omit()

p_tr <- ggplot() +
  geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.2) +
  geom_tile(data = proj_df, aes(x = x, y = y, fill = cloglog)) +
  scale_fill_viridis_c(name = "Suitability") +
  coord_sf(xlim = range(proj_df$x), ylim = range(proj_df$y)) +
  ggtitle("Projected suitability (AU+SEA)") +
  theme_classic()

# Panel C: projected binary presence/absence (full AU+SEA extent)
pa_df <- as.data.frame(proj_out$proj_pa, xy = TRUE) |> na.omit() |>
           dplyr::mutate(presence = as.character(presence))

p_bl <- ggplot() +
  geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.2) +
  geom_tile(data = pa_df, aes(x = x, y = y, fill = presence)) +
  scale_fill_viridis_d(name = "Presence") +
  coord_sf(xlim = range(pa_df$x), ylim = range(pa_df$y)) +
  ggtitle("Projected presence/absence (AU+SEA)") +
  theme_classic()

# Panel D: MESS surface — NULL panel (blank) if calc_mess = FALSE
if (!is.null(proj_out$mess)) {
  mess_df <- as.data.frame(proj_out$mess, xy = TRUE) |> na.omit()
  p_br <- ggplot() +
    geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.2) +
    geom_tile(data = mess_df, aes(x = x, y = y, fill = MESS)) +
    scale_fill_gradient2(name = "MESS", low = "red", mid = "white",
                         high = "blue", midpoint = 0) +
    coord_sf(xlim = range(mess_df$x), ylim = range(mess_df$y)) +
    ggtitle("MESS (red = novel environment)") +
    theme_classic()
} else {
  p_br <- NULL   # blank bottom-right cell when MESS was skipped
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

rm(p_tl, p_tr, p_bl, p_br, panel_fig, proj_df, pa_df, train_df)
gc()
files <- list.files(tempdir(), full.names = TRUE, all.files = TRUE, recursive = TRUE)
file.remove(files)
gc()
