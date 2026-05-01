library(terra)
library(dplyr)

#' Save best-model SDM outputs and return the optimal TSS threshold.
#'
#' Selects the best ENMeval model (delta.AICc == 0, subject to AUC >= AUCmin),
#' saves a continuous cloglog prediction raster and a presence-absence raster
#' thresholded by the quantile level that maximises a modified TSS, and writes
#' variable importance to CSV.
#'
#' @param ENMeval_output  ENMevaluation object from ENMeval::ENMevaluate().
#' @param AUCmin          Minimum acceptable AUC (training and validation).
#' @param resultDir       Directory path for output files.
#' @param spp             Species binomial (spaces OK; converted internally).
#' @param occ_df          Data frame with columns LONG / LAT.
#' @return  Named list: threshold (cloglog value used for PA map) and
#'          quant_level (quantile level that produced that threshold).

# Extract a prediction layer from an ENMevaluation object by tune.args string.
# Uses numeric-position fallback because name-based [[]] lookup on a SpatRaster
# can silently return an empty raster when layer names contain dots or differ
# slightly in numeric formatting (e.g. "rm.1" vs "rm.1.0").
.extract_prediction <- function(enm_obj, tune_args_str) {
  preds <- enm_obj@predictions
  pred_names <- as.character(names(preds))

  # Try name lookup first
  r_raw <- tryCatch(preds[[tune_args_str]], error = function(e) NULL)

  # Fall back to numeric index if name lookup returned empty or failed
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

save_SDM_results <- function(ENMeval_output, AUCmin, resultDir, spp, occ_df) {

  spp_bw <- stringr::str_replace(spp, " ", "_")

  bestmod <- ENMeval_output@results |> dplyr::filter(delta.AICc == 0)
  bestmod <- bestmod[1, ]

  # Helper: calculate modified TSS at a given quantile level p
  # TSS = (Se + 0.33*Sp) - 1  where Se = sensitivity, Sp = specificity
  calculate_tss <- function(p, r_best, spp_pts, back_pts) {
    thresh  <- quantile(terra::extract(r_best, spp_pts, ID = FALSE)[[1]], probs = p, na.rm = TRUE)
    pa_mat  <- matrix(c(0, thresh, 0, thresh, 1, 1), ncol = 3, byrow = TRUE)
    r_pa    <- terra::classify(r_best, pa_mat)

    pres_vals <- terra::extract(r_pa, spp_pts, ID = FALSE) |> dplyr::rename(pa = 1) |> na.omit()
    abs_vals  <- terra::extract(r_pa, back_pts, ID = FALSE) |> dplyr::rename(pa = 1) |> na.omit()

    Se  <- sum(pres_vals$pa)    / nrow(pres_vals)
    Sp  <- sum(1 - abs_vals$pa) / nrow(abs_vals)
    (Se + 0.33 * Sp) - 1
  }

  # Shared PA-threshold logic used regardless of which model is selected
  build_pa_outputs <- function(r_best, maxent_args) {

    spp_pts <- occ_df |>
      dplyr::rename(x = LONG, y = LAT) |>
      dplyr::select(x, y)

    # Sample background proportional to occurrence count (max 100 k)
    back_pts <- terra::spatSample(r_best, size = min(nrow(occ_df) * 1000, 100000),
                                   method = "random", na.rm = TRUE, xy = TRUE) |>
      dplyr::slice_sample(n = nrow(occ_df)) |>
      dplyr::select(x, y)

    quant_levels <- c(0, 0.01, 0.025, 0.05, 0.1)
    tss_vals <- sapply(quant_levels, calculate_tss,
                       r_best = r_best, spp_pts = spp_pts, back_pts = back_pts)

    tss_df <- data.frame(tss = tss_vals, quant = quant_levels)
    best_q <- dplyr::filter(tss_df, tss == max(tss))$quant
    if (length(best_q) > 1) best_q <- stats::median(best_q)

    occ_pred   <- terra::extract(r_best, spp_pts, ID = FALSE)[[1]]
    lpt_thresh <- quantile(occ_pred, probs = best_q, na.rm = TRUE)

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

  # ---- Path A: best AICc model meets AUC threshold ----
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
  # ---- Path B: fall back to highest-AUC model ----

    bestmod <- ENMeval_output@results |>
      dplyr::filter(auc.train  == max(auc.train) |
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
