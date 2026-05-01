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

# ── Settings ──────────────────────────────────────────────────────────────────

genus   <- "Archibasis"
species <- "crucigera"
binomial <- paste(genus, species)
bw       <- paste(genus, species, sep = "_")   # filename-safe version

fp  <- "SDM_outputs/"
fp2 <- "SDM_maps/"
dir.create(file.path(fp, bw), showWarnings = FALSE, recursive = TRUE)
dir.create(fp2, showWarnings = FALSE, recursive = TRUE)

# World boundaries and projection extent built once
world       <- ne_countries(scale = "medium", returnclass = "sf")
proj_extent <- define_projectionExtent()

# ── Helper: spatial thinning ──────────────────────────────────────────────────
# Returns a thinned data frame with one record per raster cell.
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

# ── 0. Load & filter occurrences ─────────────────────────────────────────────

occs <- data.table::fread("data/Combined_data_AUS-PHI-IND.csv")

cleanedOccs <- occs |>
  dplyr::filter(GEN_DB == genus, SP_DB == species) |>
  dplyr::filter(!is.na(LONG), !is.na(LAT))

message(sprintf("Records for %s: %d", binomial, nrow(cleanedOccs)))

# ── 1. Accessible area (M region) ────────────────────────────────────────────

aa_shp <- define_accessibleArea(species_df    = cleanedOccs,
                                 minBuff       = 75000,
                                 saveImage     = FALSE,
                                 saveShapefile = FALSE)

# ── 2. Environmental variables clipped to M region ───────────────────────────

mod_vars <- clip_variableLayers(layerDir       = "ClimateOnly/",
                                 accessibleArea = aa_shp)
mod_vars <- terra::aggregate(mod_vars, 5)

# ── 3. Spatial thinning (adaptive to accessible-area size) ───────────────────

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

# ── 4. Variable selection (VIF < 5) ──────────────────────────────────────────

max_model <- dismo::maxent(x        = raster::stack(mod_vars),
                            p        = occ_matrix,
                            progress = "text")

predictors <- select_sdmVariables(pred_vars  = mod_vars,
                                   maxent_mod = max_model,
                                   maxVIF     = 5)

message(sprintf("Variables retained after VIF selection: %s",
                paste(names(predictors), collapse = ", ")))

# ── 5. ENMevaluate — spatial block CV, tuned RM & FC ─────────────────────────

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

# ── 6. Save training-area SDM outputs; get TSS threshold for projection ───────

tss_out <- save_SDM_results(ENMeval_output = eval1,
                             AUCmin         = 0.7,
                             resultDir      = file.path(fp, bw),
                             spp            = binomial,
                             occ_df         = spp_df)

save(eval1, file = file.path(fp, bw, paste0(bw, "_ENMeval.RData")))

# ── 7. Project to AU+SEA with MaxEnt clamping + MESS ─────────────────────────

proj_vars <- clip_projectionLayers(
  layerDir         = "ClimateOnly/",
  projectionExtent = proj_extent,
  varNames         = names(predictors)
)

project_toRegion(
  ENMeval_output  = eval1,
  training_vars   = predictors,
  occ_df          = spp_df,
  projection_vars = proj_vars,
  tss_threshold   = tss_out$threshold,
  resultDir       = file.path(fp, bw),
  spp             = binomial,
  mapDir          = fp2
)

# ── 8. Training-area maps ─────────────────────────────────────────────────────

r  <- terra::rast(file.path(fp, bw, paste0(bw, "_SDM.tif")))
r2 <- terra::rast(file.path(fp, bw, paste0(bw, "_SDM_PA.tif")))

xlim <- c(min(cleanedOccs$LONG) - 5,  max(cleanedOccs$LONG) + 5)
ylim <- c(min(cleanedOccs$LAT)  - 5,  max(cleanedOccs$LAT)  + 5)

p <- ggplot() +
  geom_sf(data = world, fill = NA) +
  geom_sf(data = aa_shp, fill = "orange", alpha = 0.5) +
  geom_point(data = cleanedOccs,
             aes(x = LONG, y = LAT),
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
             aes(x = LONG, y = LAT),
             size = 0.5, alpha = 0.75, shape = 1) +
  scale_fill_viridis_d() +
  labs(fill = "Presence") +
  coord_sf(xlim = c(xlim[1] - 5, xlim[2] + 5),
           ylim = c(ylim[1] - 5, ylim[2] + 5)) +
  ggtitle(paste(binomial, "n =", nrow(cleanedOccs))) +
  theme_classic()

ggsave(plot     = bp,
       filename = file.path(fp2, paste0(bw, "_sdmMapOccPoints.png")))

# ── 9. Cleanup ────────────────────────────────────────────────────────────────

rm(p, r, p2, r2, p3, e, bp)
gc()
files <- list.files(tempdir(), full.names = TRUE, all.files = TRUE, recursive = TRUE)
file.remove(files)
gc()
