# =============================================================================
# 02_clean_occurrences.R — Step 1: Occurrence Record Cleaning
# =============================================================================
#
# WHAT THIS STEP DOES (conceptually)
# ------------------------------------
# Before building a species distribution model, occurrence records must be
# filtered to remove records that are erroneous, imprecise, or duplicated.
# Garbage-in → garbage-out applies especially strongly to MaxEnt, which uses
# the spatial pattern of occurrences to learn what environments a species
# occupies. A single record placed at a museum headquarters or in the ocean
# can noticeably distort that signal.
#
# This script uses the CoordinateCleaner package, which applies a battery of
# spatial heuristics to flag and remove suspect records.
#
# WHAT THE USER CAN CHANGE
# -------------------------
# See individual parameter notes inside coord_clean() below. The most commonly
# tuned parameters are:
#   buffer (cc_cen / cc_gbif) — distance in metres around known problematic
#       coordinates (institution centroids, GBIF HQ). Larger = more aggressive.
#   tdi (cc_outl) — distance-based outlier threshold in km. Decrease to remove
#       records far from the species' core range.
#   mltpl (cc_outl, second pass) — IQR multiplier for the quantile-based
#       outlier filter. Decrease to be more aggressive about removing
#       geographically peripheral records.
# =============================================================================

library(dplyr)
library(ggplot2)
library(CoordinateCleaner)
library(rnaturalearth)
library(sf)
library(stringr)


#' Clean occurrence coordinates for a single species using CoordinateCleaner.
#'
#' Applies a sequential pipeline of spatial quality filters and returns the
#' cleaned data frame. Optionally saves the cleaned records as a CSV and/or
#' a before/after comparison map.
#'
#' @param species      Species name string matching the `scientific_name` column.
#' @param df           Full occurrence data frame containing all species.
#' @param imagedir     Directory to write the comparison map PNG (if saveimage=TRUE).
#' @param occdir       Directory to write the cleaned CSV (if saveocc=TRUE).
#' @param saveimage    Logical. Save a before/after map? Default FALSE.
#' @param saveocc      Logical. Save cleaned records to CSV? Default FALSE.
#' @return  Data frame of cleaned occurrence records.

coord_clean <- function(species, df,
                        imagedir = "NotNeeded",
                        occdir   = "NotNeeded",
                        saveimage = FALSE,
                        saveocc   = FALSE) {

  cs <- filter(df, scientific_name == species)

  # ── Step 1: Remove records with missing or non-numeric coordinates ──────────
  # NA coordinates cannot be used; non-numeric values suggest data entry errors.
  cs1 <- cs %>%
    filter(!is.na(decimalLatitude)) %>%
    filter(!is.na(decimalLongitude)) %>%
    mutate(decimalLatitude  = as.numeric(decimalLatitude),
           decimalLongitude = as.numeric(decimalLongitude)) %>%
    cc_val(lon = "decimalLongitude", lat = "decimalLatitude", value = "clean")

  cs2 <- cs1 %>%

    # ── Step 2: Remove records where lat == lon (e.g. 10°N, 10°E) ─────────────
    # A common data-entry mistake; these coordinates are almost certainly wrong.
    cc_equ(lon = "decimalLongitude", lat = "decimalLatitude",
           value = "clean") %>%

    # ── Step 3: Remove records too close to country/province centroids ─────────
    # Occurrences pinned to administrative centroids (often used as a
    # placeholder when exact coordinates are unknown) create artificial
    # hotspots. buffer = 500 m removes records within 500 m of a centroid.
    # TUNE: increase buffer (e.g. 1000–5000 m) if you suspect many centroid-
    #       pinned records; decrease if occurrences legitimately cluster near
    #       capital cities.
    cc_cen(lon = "decimalLongitude", lat = "decimalLatitude",
           species = "validName", buffer = 500, value = "clean") %>%

    # ── Step 4: Remove exact duplicates ────────────────────────────────────────
    # Duplicate lat/lon pairs for the same species add false weight to certain
    # locations. ENMeval's spatial thinning (in the pipeline) handles spatial
    # autocorrelation further, but removing true duplicates here is cleaner.
    cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude",
            value = "clean", species = "validName") %>%

    # ── Step 5: Remove records near the GBIF headquarters ──────────────────────
    # GBIF's office in Copenhagen (55.7°N, 12.6°E) is a common erroneous
    # default location in aggregated databases. buffer = 500 m is a safe zone.
    cc_gbif(lon = "decimalLongitude", lat = "decimalLatitude",
            value = "clean", species = "validName", buffer = 500) %>%

    # ── Step 6: Remove records at or very near (0°, 0°) ───────────────────────
    # The intersection of the prime meridian and the equator is the most common
    # placeholder coordinate in global databases. buffer = 0.05° ≈ 5 km.
    cc_zero(lon = "decimalLongitude", lat = "decimalLatitude",
            buffer = 0.05, value = "clean") %>%

    # ── Step 7: Distance-based outlier removal (first pass) ────────────────────
    # cc_outl identifies records whose nearest-neighbour distance exceeds a
    # threshold (tdi), suggesting they may be misidentifications or georeference
    # errors rather than true range extensions.
    # tdi = 1000 km — records > 1000 km from their nearest neighbour are flagged.
    # TUNE: lower tdi (e.g. 500) to be more aggressive; raise (e.g. 2000) if the
    #       species has known disjunct populations separated by large distances.
    cc_outl(lon = "decimalLongitude", lat = "decimalLatitude",
            value = "clean", species = "validName", method = "distance",
            tdi = 1000)

  # ── Step 8: Remove marine/coastal records (with scale fallback) ─────────────
  # Freshwater insects like dragonflies should not have marine occurrences.
  # cc_sea uses a land polygon; records falling in the ocean are removed.
  # The tryCatch tries the default (high-resolution) coastline first and falls
  # back to scale = 50 if that errors (e.g. package installation issues).
  cs2 <-
    tryCatch(cc_sea(x = cs2, lon = "decimalLongitude", lat = "decimalLatitude",
                    value = "clean"),
             error = function(e)
               cc_sea(x = cs2, lon = "decimalLongitude", lat = "decimalLatitude",
                      value = "clean", scale = 50))

  # ── Step 9: Adaptive quantile-based outlier removal (second pass) ───────────
  # After the distance filter above, this pass uses a quantile-based approach
  # to catch clusters of outliers that the distance method can miss.
  #
  # The IQR of the coordinate distribution informs how strict to be:
  #   • Narrow IQR (< 0.5°) → species highly clustered → use a very lenient
  #     multiplier (mltpl = 50) to avoid removing valid records near the core.
  #   • Moderate IQR (0.5–1°) → mltpl = 25 (moderately strict).
  #   • Wider IQR (≥ 1°) → mltpl = 10 (standard strictness).
  #
  # TUNE: lower mltpl values remove more records; raise them to be more
  #       conservative. A species known to have a tight range may warrant
  #       mltpl = 5–10 even with a narrow IQR if outliers are clearly erroneous.
  if (IQR(cs2$decimalLatitude) < 0.5 | IQR(cs2$decimalLongitude) < 0.5) {
    cs3 <- cc_outl(x = cs2, lon = "decimalLongitude", lat = "decimalLatitude",
                   value = "clean", species = "validName", method = "quantile",
                   mltpl = 50)
  } else if (IQR(cs2$decimalLatitude) >= 0.5 & IQR(cs2$decimalLatitude) < 1 |
             IQR(cs2$decimalLongitude) >= 0.5 & IQR(cs2$decimalLongitude) < 1) {
    cs3 <- cc_outl(x = cs2, lon = "decimalLongitude", lat = "decimalLatitude",
                   value = "clean", species = "validName", method = "quantile",
                   mltpl = 25)
  } else {
    cs3 <- cc_outl(x = cs2, lon = "decimalLongitude", lat = "decimalLatitude",
                   value = "clean", species = "validName", method = "quantile",
                   mltpl = 10)
  }

  # ── Optional outputs ─────────────────────────────────────────────────────────

  if (saveocc == TRUE) {
    if (nrow(cs3) > 10) {
      write.csv(x = cs3, file = paste(occdir, "/", species, ".csv", sep = ""),
                row.names = FALSE)
    } else {
      print(paste(species, "does not have enough data points, no csv produced"))
    }
  } else {
    print("Decision made to not save occurrence points to directory")
  }

  world <- ne_countries(scale = "medium", returnclass = "sf")

  a <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_point(cs1, mapping = aes(x = decimalLongitude, y = decimalLatitude),
               color = "purple") +
    coord_sf(xlim = c(min(cs1$decimalLongitude) - 10, max(cs1$decimalLongitude) + 10),
             ylim = c(min(cs1$decimalLatitude)  - 10, max(cs1$decimalLatitude)  + 10)) +
    ggtitle("Uncleaned") +
    theme_bw()

  b <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_point(cs3, mapping = aes(x = decimalLongitude, y = decimalLatitude),
               color = "purple") +
    coord_sf(xlim = c(min(cs1$decimalLongitude) - 10, max(cs1$decimalLongitude) + 10),
             ylim = c(min(cs1$decimalLatitude)  - 10, max(cs1$decimalLatitude)  + 10)) +
    ggtitle("Cleaned") +
    theme_bw()

  cp <- cowplot::plot_grid(a, b)

  if (saveimage == TRUE) {
    if (nrow(cs3) > 10) {
      ggsave(filename = paste(imagedir, "/", species, "_occs", ".png", sep = ""),
             plot = cp)
    } else {
      print(paste(species, "does not have enough data points, no image produced"))
    }
  } else {
    print("Decision made not to save image of occurrence points to image directory")
  }

  return(cs3)
}
