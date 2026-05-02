# =============================================================================
# 08_save_SDMoutputs_TSS.R — Step 6: Save Best Model & Derive TSS Threshold
# =============================================================================
#
# WHAT THIS STEP DOES (conceptually)
# ------------------------------------
# ENMevaluate produces dozens of candidate models (one per rm × fc combination).
# This script:
#   1. Selects the best model using AICc (information-criterion approach).
#   2. Validates that the best AICc model meets a minimum AUC requirement.
#      If not, it falls back to the model with the highest AUC.
#   3. Saves the continuous (cloglog) prediction raster for the training area.
#   4. Converts the continuous prediction to a binary presence/absence (PA) map
#      by finding the cloglog threshold that maximises a modified TSS.
#   5. Saves variable importance to CSV.
#   6. Returns the TSS threshold so that project_toRegion() can apply the
#      same threshold to the projected surface.
#
# MODEL SELECTION: ΔAIC c == 0
# ----------------------------
# AICc (Akaike Information Criterion, corrected for small samples) penalises
# model complexity. The model with the lowest AICc (delta.AICc == 0) has the
# best trade-off between goodness-of-fit and parsimony. Among multiple models
# tied at delta.AICc == 0, the first row is taken.
#
# AUC MINIMUM (AUCmin)
# --------------------
# AUC (Area Under the ROC Curve) measures the model's ability to rank occupied
# cells above background cells. AUC = 0.5 is random; AUC = 1.0 is perfect.
#   • AUCmin = 0.7: the default "acceptable" threshold. Models below this are
#     considered poorly discriminating and the fallback (highest-AUC model) is
#     used instead.
#   • Both training AUC (on the full dataset) and validation AUC (cross-
#     validated, more conservative) must meet the threshold.
#   TUNE: raise to 0.8 for stricter quality control; lower to 0.6 for rare
#         species with few records where cross-validation AUC is noisy.
#
# TSS THRESHOLD
# -------------
# The True Skill Statistic (TSS = sensitivity + specificity − 1) measures
# classification accuracy accounting for chance agreement. Here we use a
# *modified* TSS that down-weights specificity (Sp × 0.33) to favour
# sensitivity (Se), reflecting the preference in conservation planning to
# minimise false negatives (missing suitable habitat) over false positives
# (predicting unsuitable habitat as suitable).
#
# The threshold is determined by testing several quantile levels applied to
# occurrence predictions:
#   quant_levels = c(0, 0.01, 0.025, 0.05, 0.1)
#   — 0 means "threshold = minimum predicted value at occurrences" (all
#     occurrence cells predicted present; very liberal).
#   — 0.1 means "threshold = 10th percentile of occurrence predictions"
#     (allows 10% of occurrences to fall in 'absent' cells; more conservative).
#
#   TUNE: Add 0.2 or 0.25 to the vector to test more conservative thresholds.
#         Remove 0 if you never want every occurrence cell classified as present.
#
# WHAT THE USER CAN CHANGE
# -------------------------
#   AUCmin       — Minimum AUC for a model to be accepted (see above).
#   quant_levels — The quantile levels tested when optimising the TSS threshold.
# =============================================================================

library(terra)
library(dplyr)

# ── Internal helper: extract a prediction raster by model name ────────────────
# Uses numeric-position fallback because name-based [[]] lookup on a SpatRaster
# can silently return an empty raster when layer names contain dots or differ
# slightly in numeric formatting (e.g. "rm.1" vs "rm.1.0").
.extract_prediction <- function(enm_obj, tune_args_str) {
  preds      <- enm_obj@predictions
  pred_names <- as.character(names(preds))

  r_raw <- tryCatch(preds[[tune_args_str]], error = function(e) NULL)

  if (is.null(r_raw) ||
      (inherits(r_raw, "SpatRaster") && terra::ncell(r_raw) == 0)) {
    idx <- which(pred_names == tune_args_str)
    if (length(idx) == 0) {
      stop(sprintf(
        "Prediction for '%s' not found.\nAvailable models: %s",
        tune_args_str, paste(pred_names, collapse = ", ")
      ))
    }
    r_raw <- preds[[idx[1]]]
  }

  if (!inherits(r_raw, "SpatRaster")) r_raw <- terra::rast(r_raw)
  r_raw
}

#' Save best-model SDM outputs and return the optimal TSS threshold.
#'
#' @param ENMeval_output  ENMevaluation object from ENMeval::ENMevaluate().
#' @param AUCmin          Minimum acceptable training and validation AUC.
#' @param resultDir       Directory for output files.
#' @param spp             Species binomial (spaces OK; converted internally).
#' @param occ_df          Data frame with columns LONG / LAT.
#' @return  Named list: threshold (cloglog value used for PA map) and
#'          quant_level (quantile level that produced that threshold).

save_SDM_results <- function(ENMeval_output, AUCmin, resultDir, spp, occ_df) {

  spp_bw  <- stringr::str_replace(spp, " ", "_")
  bestmod <- ENMeval_output@results |> dplyr::filter(delta.AICc == 0)
  bestmod <- bestmod[1, ]

  # ── Inner helper: TSS at a given quantile level ────────────────────────────
  # p = quantile of occurrence predictions used as the cut-off threshold.
  # Se (sensitivity) = proportion of occurrences correctly classified as present.
  # Sp (specificity) = proportion of background correctly classified as absent.
  # Modified TSS = (Se + 0.33 * Sp) − 1; down-weighting Sp favours sensitivity.
  calculate_tss <- function(p, r_best, spp_pts, back_pts) {
    thresh <- quantile(terra::extract(r_best, spp_pts)[[2]], probs = p, na.rm = TRUE)
    pa_mat <- matrix(c(0, thresh, 0, thresh, 1, 1), ncol = 3, byrow = TRUE)
    r_pa   <- terra::classify(r_best, pa_mat)

    pres_vals <- terra::extract(r_pa, spp_pts) |>
      dplyr::select(-ID) |> dplyr::rename(pa = 1) |> na.omit()
    abs_vals  <- terra::extract(r_pa, back_pts) |>
      dplyr::select(-ID) |> dplyr::rename(pa = 1) |> na.omit()

    Se  <- sum(pres_vals$pa)    / nrow(pres_vals)
    Sp  <- sum(1 - abs_vals$pa) / nrow(abs_vals)
    (Se + 0.33 * Sp) - 1
  }

  # ── Inner helper: build PA raster and save all outputs ────────────────────
  build_pa_outputs <- function(r_best, maxent_args) {

    # Project occurrence coordinates to the raster's CRS if it is not longlat.
    # This handles the case where the rasters have been reprojected to an
    # equal-area CRS (e.g. AU+SEA Albers) before model fitting.
    spp_coords <- as.matrix(occ_df[, c("LONG", "LAT")])
    if (!terra::is.lonlat(r_best)) {
      spp_coords <- terra::project(spp_coords,
                                    from = "EPSG:4326",
                                    to   = terra::crs(r_best, proj = TRUE))
    }
    spp_pts <- as.data.frame(spp_coords)
    colnames(spp_pts) <- c("x", "y")

    # Background points sampled from the accessible-area raster (same extent as
    # the model was trained on). Size = min(nOcc × 1000, 100000) then subsampled
    # to nOcc to keep the comparison balanced.
    back_pts <- terra::spatSample(r_best,
                                  size   = min(nrow(occ_df) * 1000, 100000),
                                  method = "random", na.rm = TRUE, xy = TRUE) |>
      dplyr::slice_sample(n = nrow(occ_df)) |>
      dplyr::select(x, y)

    # Test each quantile level; pick the one maximising modified TSS.
    quant_levels <- c(0, 0.01, 0.025, 0.05, 0.1)
    tss_vals <- sapply(quant_levels, calculate_tss,
                       r_best = r_best, spp_pts = spp_pts, back_pts = back_pts)

    tss_df  <- data.frame(tss = tss_vals, quant = quant_levels)
    best_q  <- dplyr::filter(tss_df, tss == max(tss))$quant
    if (length(best_q) > 1) best_q <- stats::median(best_q)

    occ_pred   <- terra::extract(r_best, spp_pts)[[2]]
    lpt_thresh <- quantile(occ_pred, probs = best_q, na.rm = TRUE)

    # Classify the continuous cloglog raster into 0/1 using the TSS threshold.
    pa_mat <- matrix(c(0, lpt_thresh, 0, lpt_thresh, 1, 1), ncol = 3, byrow = TRUE)
    r_pa   <- terra::classify(r_best, pa_mat)
    terra::writeRaster(r_pa,
                       filename  = file.path(resultDir, paste0(spp_bw, "_SDM_PA.tif")),
                       overwrite = TRUE)

    varimp <- ENMeval_output@variable.importance[[maxent_args]]
    write.csv(varimp,
              file      = file.path(resultDir, paste0(spp_bw, "_variableImportance.csv")),
              row.names = FALSE)

    list(threshold = as.numeric(lpt_thresh), quant_level = best_q)
  }

  # ── Path A: best AICc model meets AUC threshold ────────────────────────────
  if (bestmod$auc.train >= AUCmin && bestmod$auc.val.avg >= AUCmin) {

    maxent_args <- as.character(bestmod$tune.args)
    write.csv(bestmod,
              file      = file.path(resultDir, paste0(spp_bw, "_bestModel.csv")),
              row.names = FALSE)

    r_best <- .extract_prediction(ENMeval_output, maxent_args)
    terra::writeRaster(r_best,
                       filename  = file.path(resultDir, paste0(spp_bw, "_SDM.tif")),
                       overwrite = TRUE)

    return(invisible(build_pa_outputs(r_best, maxent_args)))

  } else {
  # ── Path B: best AICc model fails AUC check → fall back to highest AUC ─────
  # This can happen if the best-complexity model happens to be in the AICc
  # sweet spot but a spatial block happens to be poorly sampled. The fallback
  # ensures we always save a raster.
    bestmod <- ENMeval_output@results |>
      dplyr::filter(auc.train   == max(auc.train)   |
                    auc.val.avg == max(auc.val.avg))
    bestmod <- bestmod[1, ]

    maxent_args <- as.character(bestmod$tune.args)
    write.csv(bestmod,
              file      = file.path(resultDir, paste0(spp_bw, "_bestModel.csv")),
              row.names = FALSE)

    r_best <- .extract_prediction(ENMeval_output, maxent_args)
    terra::writeRaster(r_best,
                       filename  = file.path(resultDir, paste0(spp_bw, "_SDM.tif")),
                       overwrite = TRUE)

    return(invisible(build_pa_outputs(r_best, maxent_args)))
  }
}
