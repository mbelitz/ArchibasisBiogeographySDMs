library(terra)
library(dismo)   # for initial maxent() call (variable importance)
library(raster)  # bridge for dismo (dismo::maxent requires RasterStack)
library(ENMeval)

#' Fit a tuned MaxEnt model using ENMevaluate (ENMeval >= 2.0).
#'
#' An initial full-variable MaxEnt model is built first solely to obtain
#' permutation importances for VIF-based variable selection. ENMevaluate
#' then tests all combinations of regularisation multipliers and feature
#' classes using spatial block cross-validation.
#'
#' @param occ_df          Data frame with columns decimalLongitude / decimalLatitude.
#' @param model_vars      terra SpatRaster of environmental predictors (accessible area).
#' @param partitionMethod ENMeval partition method (e.g. "block", "randomkfold").
#' @param reg_mult        Numeric vector of regularisation multiplier values.
#' @param fcs             Character vector of feature class combinations.
#' @return  ENMevaluation object.

generate_ENMeval <- function(occ_df, model_vars,
                              partitionMethod,
                              reg_mult,
                              fcs) {

  occ_matrix <- as.matrix(occ_df[, c("decimalLongitude", "decimalLatitude")])

  # Initial model for permutation importance used in VIF selection.
  # dismo::maxent() requires a RasterStack — bridge via raster::stack().
  max_model <- dismo::maxent(
    x        = raster::stack(model_vars),
    p        = occ_matrix,
    progress = "text"
  )

  predictors <- select_sdmVariables(pred_vars = model_vars,
                                    maxent_mod = max_model,
                                    maxVIF     = 5)

  eval1 <- ENMeval::ENMevaluate(
    occs      = occ_matrix,
    envs      = predictors,
    tune.args = list(rm = reg_mult, fc = fcs),
    partitions = partitionMethod,
    algorithm  = "maxent.jar",
    parallel   = TRUE,
    numCores   = 5
  )

  return(eval1)
}
