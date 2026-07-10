# =============================================================================
# Darrell Kaufman's GLS/softplus two-acid age calculator (reference method)
# =============================================================================
#
# A faithful R implementation of the spreadsheet age calculator described in
# the in-prep Np AAR calibration study: generalized least squares (GLS)
# combination of independent TDK (Asp) and SPK (Glu) age estimates, using a
# softplus-smoothed weight floor to keep both weights positive even when the
# raw GLS weight would go negative (see `vignette("two-acid-age-inversion")`
# for the scientific background). This is not a Bayesian model -- it is
# included so the two Bayesian approaches built elsewhere in this package
# (`R/age_calibration.R` and `R/age_calibration_forward.R`) can be compared
# directly against the original point-estimate-plus-prediction-interval
# method, using its own published formulas and calibration constants.

#' Calibration constants for the GLS age calculator
#'
#' Table 2 / calibration constants from the in-prep Np AAR calibration study
#' (n = 128 Neogloboquadrina pachyderma samples), reproduced from the
#' author's spreadsheet calculator's `Constants` sheet.
#'
#' @format A named list with elements `a_Asp`, `e_Asp`, `MSE_Asp`, `n_Asp`,
#'   `xbar_Asp`, `Sxx_Asp`, `a_Glu`, `e_Glu`, `MSE_Glu`, `n_Glu`, `xbar_Glu`,
#'   `Sxx_Glu`, `rho`, `SE_rho`, `t50`, `t66`, `t90`.
#'
#' @export
gls_calibration_constants <- list(
  a_Asp    = 2872.4972617303,
  e_Asp    = 2.7585198455,
  MSE_Asp  = 0.1647421056,
  n_Asp    = 128,
  xbar_Asp = -1.3427301286,
  Sxx_Asp  = 36.1405601873,
  a_Glu    = 9581.1294167325,
  e_Glu    = 2.2005812166,
  MSE_Glu  = 0.2207267401,
  n_Glu    = 128,
  xbar_Glu = -2.2305749821,
  Sxx_Glu  = 55.333375227,
  rho      = 0.824940751,
  SE_rho   = 0.0306425681,
  t50      = 0.676441857,
  t66      = 0.9577956685,
  t90      = 1.657036982
)


#' GLS-combined two-acid age with softplus weight smoothing
#'
#' Reference implementation of Darrell Kaufman's spreadsheet age calculator:
#' independent TDK (Asp) and SPK (Glu) age predictions are combined by
#' generalized least squares, with the (possibly negative) raw GLS weights
#' passed through a softplus floor before combination. See
#' \code{\link{gls_calibration_constants}} for the underlying calibration
#' constants.
#'
#' @param DL_Asp,DL_Glu Observed D/L ratios (sample means).
#' @param SD_Asp,SD_Glu Standard deviation across subsamples.
#' @param N_Asp,N_Glu Number of subsamples used for the mean.
#' @param constants List of calibration constants, see
#'   \code{\link{gls_calibration_constants}}. Default uses the package's
#'   built-in Np constants.
#'
#' @return A one-row data frame with columns `t_Asp`, `t_Glu`, `V_Asp`,
#'   `V_Glu`, `C`, `w_Asp`, `w_Glu`, `D`, `V_int`, `V_excess`, `V_total`,
#'   `sigma_comb`, `t_comb`, `lo50`, `hi50`, `lo66`, `hi66`, `lo90`, `hi90`.
#'
#' @examples
#' # Reproduces "Example 1" from the spreadsheet calculator
#' predict_age_gls(DL_Asp = 0.350, SD_Asp = 0.01, N_Asp = 5,
#'                 DL_Glu = 0.153, SD_Glu = 0.015, N_Glu = 5)
#'
#' @export
predict_age_gls <- function(DL_Asp, SD_Asp, N_Asp, DL_Glu, SD_Glu, N_Glu,
                            constants = gls_calibration_constants) {
  cst <- constants
  SE_Asp <- SD_Asp / sqrt(N_Asp)
  SE_Glu <- SD_Glu / sqrt(N_Glu)

  t_Asp <- cst$a_Asp * atanh(DL_Asp)^cst$e_Asp
  t_Glu <- cst$a_Glu * DL_Glu^cst$e_Glu

  # Predictor-derivative terms propagate D/L measurement uncertainty
  # through the log(f(D/L)) transform via the chain rule.
  slope_Asp <- cst$e_Asp / (atanh(DL_Asp) * (1 - DL_Asp^2))
  slope_Glu <- cst$e_Glu / DL_Glu

  V_Asp <- cst$MSE_Asp * (1 + 1 / cst$n_Asp +
                           (log(atanh(DL_Asp)) - cst$xbar_Asp)^2 / cst$Sxx_Asp) +
    (slope_Asp * SE_Asp)^2
  V_Glu <- cst$MSE_Glu * (1 + 1 / cst$n_Glu +
                           (log(DL_Glu) - cst$xbar_Glu)^2 / cst$Sxx_Glu) +
    (slope_Glu * SE_Glu)^2

  C <- cst$rho * sqrt(V_Asp * V_Glu)
  s <- sqrt(V_Asp * V_Glu) * cst$SE_rho

  .softplus_floor <- function(raw, s) {
    pmax(raw, 0) + s * log(1 + exp(-abs(raw / s)))
  }
  w_Asp <- .softplus_floor(V_Glu - C, s)
  w_Glu <- .softplus_floor(V_Asp - C, s)
  D <- w_Asp + w_Glu

  V_int    <- (V_Asp * V_Glu - C^2) / D
  V_excess <- w_Asp * w_Glu * (log(t_Asp) - log(t_Glu))^2 / D^2
  V_total  <- V_int + V_excess
  sigma_comb <- sqrt(V_total)

  t_comb <- t_Asp^(w_Asp / D) * t_Glu^(w_Glu / D)

  data.frame(
    t_Asp = t_Asp, t_Glu = t_Glu, V_Asp = V_Asp, V_Glu = V_Glu, C = C,
    w_Asp = w_Asp, w_Glu = w_Glu, D = D, V_int = V_int, V_excess = V_excess,
    V_total = V_total, sigma_comb = sigma_comb, t_comb = t_comb,
    lo50 = t_comb * exp(-cst$t50 * sigma_comb),
    hi50 = t_comb * exp(+cst$t50 * sigma_comb),
    lo66 = t_comb * exp(-cst$t66 * sigma_comb),
    hi66 = t_comb * exp(+cst$t66 * sigma_comb),
    lo90 = t_comb * exp(-cst$t90 * sigma_comb),
    hi90 = t_comb * exp(+cst$t90 * sigma_comb)
  )
}
