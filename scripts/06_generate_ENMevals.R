# =============================================================================
# 06_generate_ENMevals.R — Step 5: Model Tuning via ENMevaluate
# =============================================================================
#
# WHAT THIS STEP DOES (conceptually)
# ------------------------------------
# MaxEnt has two main "complexity knobs" that interact to control how closely
# the model fits the training data. Running ENMevaluate tests all combinations
# of these knobs and uses cross-validation to pick the settings that generalise
# best — neither underfitting (missing real patterns) nor overfitting (learning
# noise in the training occurrences).
#
# REGULARISATION MULTIPLIER (RM / reg_mult)
# -----------------------------------------
# MaxEnt is a maximum-entropy model with a regularisation penalty that prevents
# feature weights from becoming arbitrarily large. The multiplier scales that
# penalty:
#   • Low RM (e.g. 0.5): weak penalty → tighter fit to occurrences → more
#     complex, potentially overfitted model with many fine-grained range
#     boundaries. Can be appropriate for well-sampled, widespread species.
#   • High RM (e.g. 4): strong penalty → smoother predictions → more
#     conservative, potentially under-fitted model. Useful for data-sparse
#     species or where sampling is geographically biased.
#   Default grid: c(0.5, 1, 2, 3, 4). Remove extreme values to speed up tuning.
#
# FEATURE CLASSES (FC / fcs)
# --------------------------
# Feature classes are mathematical transformations applied to predictor
# variables before fitting:
#   L  = Linear           (straight-line response; simplest)
#   Q  = Quadratic        (unimodal response; captures optima)
#   H  = Hinge            (piecewise linear; flexible threshold responses)
#   P  = Product          (interaction terms between two variables)
#   T  = Threshold        (step functions; sharp transitions)
#
# More complex FC combinations (e.g. LQHPT) allow the model to fit more
# intricate response curves but also increase overfitting risk. Start with
# simpler combinations (L, LQ) for small datasets (< 50 records); richer
# combinations are appropriate with > 100 records.
#   Default: c("L", "LQ", "H", "LQH", "LQHP", "LQHPT")
#
# CROSS-VALIDATION / PARTITION METHOD
# ------------------------------------
# "block" partitioning divides the study area into four spatial blocks
# (quadrants) and uses each block as a hold-out in turn. This is more
# informative than random k-fold for evaluating transferability because
# each fold tests the model on a geographically distinct region. However,
# it requires ≥ 4 occurrences in each block — for very few records, switch
# to partitionMethod = "randomkfold" with k = 4 or k = 2.
#
# WHAT THE USER CAN CHANGE
# -------------------------
#   reg_mult         — RM values to test (see above). Add finer resolution
#                      (e.g. 1.5) or a wider range if default results cluster
#                      at the boundary.
#   fcs              — FC combinations to test. Remove H, P, T for small
#                      datasets to reduce overfitting.
#   partitionMethod  — "block" (default), "randomkfold", "jackknife", etc.
#   parallel / numCores — ENMevaluate can parallelise across models. Set
#                      numCores to the number of physical cores on your
#                      machine minus 1. On Windows, parallel = TRUE requires
#                      the doParallel backend.
# =============================================================================

library(terra)
library(dismo)   # for initial maxent() call (variable importance)
library(raster)  # bridge for dismo (dismo::maxent requires RasterStack)
library(ENMeval)

#' Fit a tuned MaxEnt model using ENMevaluate (ENMeval >= 2.0).
#'
#' Two-stage approach:
#'   1. An initial full-variable dismo::maxent() model is fitted solely to
#'      obtain permutation importances for VIF-based variable selection.
#'   2. ENMevaluate then exhaustively tests all (rm × fc) combinations using
#'      spatial block cross-validation and ranks them by AICc.
#'
#' NOTE: This function is provided as a standalone helper. In sdm_pipeline.R
#' the variable-selection step and ENMevaluate call are inlined for clarity,
#' but the logic is identical.
#'
#' @param occ_df          Data frame with columns decimalLongitude / decimalLatitude.
#' @param model_vars      terra SpatRaster of environmental predictors (M region).
#' @param partitionMethod ENMeval partition method string. Default "block".
#' @param reg_mult        Numeric vector of regularisation multiplier values.
#' @param fcs             Character vector of feature class combinations.
#' @return  ENMevaluation object containing predictions, results table,
#'          variable importances, and the optimal model.

generate_ENMeval <- function(occ_df, model_vars,
                              partitionMethod,
                              reg_mult,
                              fcs) {

  occ_matrix <- as.matrix(occ_df[, c("decimalLongitude", "decimalLatitude")])

  # Stage 1: initial full-variable model for permutation importance.
  # dismo::maxent() requires a RasterStack (not a terra SpatRaster), so
  # raster::stack() is used as a bridge. This model is NOT used for
  # prediction — only its @results table (variable importances) is read.
  max_model <- dismo::maxent(
    x        = raster::stack(model_vars),
    p        = occ_matrix,
    progress = "text"
  )

  # Stage 2: VIF-based variable selection (see 05_select_modelVariables.R).
  predictors <- select_sdmVariables(pred_vars = model_vars,
                                    maxent_mod = max_model,
                                    maxVIF     = 5)

  # Stage 3: exhaustive tuning across all (rm × fc) combinations.
  eval1 <- ENMeval::ENMevaluate(
    occs       = occ_matrix,
    envs       = predictors,
    tune.args  = list(rm = reg_mult, fc = fcs),
    partitions = partitionMethod,
    algorithm  = "maxent.jar",
    parallel   = TRUE,
    numCores   = 5
  )

  return(eval1)
}
