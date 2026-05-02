# =============================================================================
# 05_select_modelVariables.R — Step 4: Variable Selection via VIF
# =============================================================================
#
# WHAT THIS STEP DOES (conceptually)
# ------------------------------------
# MaxEnt can technically fit a model with any number of predictors, but
# including highly correlated variables causes two problems:
#
#   1. Collinearity inflates model complexity and makes coefficients unstable —
#      swap one correlated variable for another and you get a very different
#      model, even though they carry nearly the same ecological information.
#
#   2. When projecting to a new region or future climate, correlated variables
#      may diverge. A model that relied on both temperature and a temperature-
#      derived index will be confused if those variables no longer track each
#      other in the projection space.
#
# The Variance Inflation Factor (VIF) quantifies how much each variable's
# variance is explained by the other predictors. A VIF of 5 means ~80% of
# that variable's variance is explained by the rest — a common threshold for
# "too collinear." We iteratively remove the highest-VIF variable (using
# permutation importance to break ties) until all VIFs are below the threshold.
#
# WHAT THE USER CAN CHANGE
# -------------------------
#   maxVIF       — The VIF threshold above which a variable is considered
#                  redundant. Lower values (e.g. 3) produce a leaner variable
#                  set and more conservative models. Higher values (e.g. 10)
#                  retain more variables but risk collinearity issues.
#                  Default: 5 (standard in SDM literature).
#
# ABOUT THE TIE-BREAKING LOGIC
# -----------------------------
# Among the top two highest-VIF candidates, the one with lower permutation
# importance (i.e., less useful to the MaxEnt model) is removed. This
# prioritises retaining ecologically informative variables even when two
# variables are similarly collinear. The permutation importance comes from
# an initial "full" MaxEnt model fitted before ENMevaluate (see sdm_pipeline.R
# step 4 for that call).
# =============================================================================

library(terra)
library(dplyr)
library(stringr)
library(usdm)   # VIF calculation; requires usdm >= 2.1-4 for terra SpatRaster

#' Remove the single most-collinear variable from a SpatRaster.
#'
#' Among variables with VIF > th, the two with the highest VIF are identified.
#' The one with lower permutation importance is dropped. This preserves
#' ecologically informative variables over redundant ones.
#'
#' @param pred_vars   terra SpatRaster of candidate predictor layers.
#' @param maxent_mod  dismo MaxEnt model object supplying permutation importances.
#' @param th          VIF threshold; only variables above this are candidates.
#' @return  terra SpatRaster with the least-informative collinear variable removed.

remove_singleVariables <- function(pred_vars, maxent_mod, th) {

  pVIF <- usdm::vif(pred_vars)

  # Extract permutation importance from the dismo MaxEnt results table.
  # Permutation importance reflects how much model AUC drops when that
  # variable's values are randomly permuted — higher = more important.
  m_results <- as.data.frame(as.table(maxent_mod@results)) |>
    dplyr::rename(variables = 1, rem = 2, permutation.importance = 3) |>
    dplyr::select(variables, permutation.importance)

  vIMP <- m_results |>
    dplyr::filter(stringr::str_detect(variables, "\\.permutation\\.importance")) |>
    dplyr::mutate(Variables = stringr::word(variables, sep = stringr::fixed(".")))

  jdf <- dplyr::left_join(pVIF, vIMP, by = "Variables")

  # Among the two highest-VIF variables, remove the one with lower permutation
  # importance. Using the top two (rather than just the single highest) avoids
  # oscillating between two equally collinear variables across iterations.
  lowVar <- jdf |>
    dplyr::filter(VIF > th) |>
    dplyr::filter(VIF == sort(VIF, decreasing = TRUE)[1] |
                  VIF == sort(VIF, decreasing = TRUE)[2]) |>
    dplyr::filter(permutation.importance == min(permutation.importance))

  keep <- names(pred_vars)[!names(pred_vars) %in% as.character(lowVar$Variables)]
  return(pred_vars[[keep]])
}

#' Iteratively remove collinear variables until all VIFs are below maxVIF.
#'
#' Calls remove_singleVariables() in a loop, removing one variable per
#' iteration, until max(VIF) < maxVIF.
#'
#' @param pred_vars   terra SpatRaster of predictor layers.
#' @param maxent_mod  dismo MaxEnt model object (for permutation importance).
#' @param maxVIF      Maximum acceptable VIF. Default 5.
#' @return  terra SpatRaster containing only low-collinearity variables.

select_sdmVariables <- function(pred_vars, maxent_mod, maxVIF = 5) {
  vv <- pred_vars
  while (max(usdm::vif(vv)$VIF) >= maxVIF) {
    vv <- remove_singleVariables(pred_vars = vv, maxent_mod = maxent_mod, th = maxVIF)
  }
  return(vv)
}
