library(terra)
library(dplyr)
library(stringr)
library(usdm)   # requires usdm >= 2.1-4 for terra SpatRaster support

#' Drop the variable with the highest VIF (or lowest permutation importance
#' when two variables tie on VIF) from a terra SpatRaster.
#'
#' @param pred_vars   terra SpatRaster of predictor layers.
#' @param maxent_mod  dismo MaxEnt model object (used for permutation importance).
#' @param th          VIF threshold above which variables are candidates for removal.
#' @return  terra SpatRaster with one layer removed.

remove_singleVariables <- function(pred_vars, maxent_mod, th) {

  pVIF <- usdm::vif(pred_vars)

  m_results <- as.data.frame(as.table(maxent_mod@results)) |>
    dplyr::rename(variables = 1, rem = 2, permutation.importance = 3) |>
    dplyr::select(variables, permutation.importance)

  vIMP <- m_results |>
    dplyr::filter(stringr::str_detect(variables, "\\.permutation\\.importance")) |>
    dplyr::mutate(Variables = stringr::word(variables, sep = stringr::fixed(".")))

  jdf <- dplyr::left_join(pVIF, vIMP, by = "Variables")

  # Among the two highest-VIF variables, remove the one with lower permutation importance
  lowVar <- jdf |>
    dplyr::filter(VIF > th) |>
    dplyr::filter(VIF == sort(VIF, decreasing = TRUE)[1] |
                  VIF == sort(VIF, decreasing = TRUE)[2]) |>
    dplyr::filter(permutation.importance == min(permutation.importance))

  keep <- names(pred_vars)[!names(pred_vars) %in% as.character(lowVar$Variables)]
  return(pred_vars[[keep]])
}

#' Iteratively remove variables until all VIFs are below maxVIF.
#'
#' @param pred_vars   terra SpatRaster of predictor layers.
#' @param maxent_mod  dismo MaxEnt model object.
#' @param maxVIF      Maximum acceptable VIF (default 5).
#' @return  terra SpatRaster containing only low-collinearity variables.

select_sdmVariables <- function(pred_vars, maxent_mod, maxVIF = 5) {
  vv <- pred_vars
  while (max(usdm::vif(vv)$VIF) >= maxVIF) {
    vv <- remove_singleVariables(pred_vars = vv, maxent_mod = maxent_mod, th = maxVIF)
  }
  return(vv)
}
