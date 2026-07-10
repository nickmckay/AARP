# =============================================================================
# Stage 2: forward-model two-acid age inversion for a single sample
# =============================================================================
#
# Companion to R/age_inversion.R, but forward-simulating D/L for a candidate
# age (via racemize_one_depth(), using the shared kinetic model fit in
# R/age_calibration_forward.R) instead of converting D/L to a point-age
# estimate through an invertible empirical curve.
#
# For a sample of unknown age, there is no real age-depth model (age is the
# unknown). Since only a single depth/age pair is ever needed, a trivial
# two-point age model -- a straight line from the core top (0 cm, 0 ka)
# through (depth_cm, candidate age) -- is built fresh for each MCMC
# proposal via rate_to_age_model(). The absolute value of depth_cm is
# arbitrary and cancels out: it only fixes the (depth, age) line whose sole
# free quantity is the candidate age, so samples with no recorded depth
# (e.g. Darrell's spreadsheet examples) can use any placeholder value.

#' Log-likelihood for forward-model two-acid age inversion
#'
#' Forward-simulates D/L for both acids at a candidate age -- via a trivial
#' single-point age-depth model and a simulated placeholder temperature
#' history (see \code{\link{simulate_temp_history}}) -- and evaluates a
#' correlated bivariate-normal likelihood on the D/L residuals, using the
#' shared kinetic model fit in Stage 1
#' (\code{\link{log_posterior_age_calibration_fwd}}).
#'
#' @param log_age Candidate value of \eqn{\log(\text{age in ka})}.
#' @param calib_params_fwd List with elements `Ae_Asp`, `Ea_Asp`, `x_Asp`,
#'   `sigma_Asp`, `Ae_Glu`, `Ea_Glu`, `x_Glu`, `sigma_Glu`, `rho` (see
#'   \code{\link{summarize_calibration_fwd}}).
#' @param DL_Asp_obs,DL_Glu_obs Observed D/L ratios for the sample.
#' @param temp_C Representative bottom-water temperature for the sample.
#' @param depth_cm Placeholder depth (cm) used only to define the trivial
#'   age-depth model; the result does not depend on its value. Default `100`.
#' @param dT_yr Integration timestep (yr). Default `2000`.
#'
#' @return Scalar log-likelihood.
#'
#' @export
log_lik_two_acid_age_fwd <- function(log_age, calib_params_fwd,
                                     DL_Asp_obs, DL_Glu_obs, temp_C,
                                     depth_cm = 100, dT_yr = 2000) {
  age_ka <- exp(log_age)
  am     <- rate_to_age_model(rate_cm_per_ka = depth_cm / age_ka,
                              max_depth_cm   = depth_cm)
  tm     <- simulate_temp_history(temp_C, age_ka)

  params_Asp <- list(Ae = calib_params_fwd$Ae_Asp, Ea = calib_params_fwd$Ea_Asp,
                     R = 0.001987, x = calib_params_fwd$x_Asp)
  params_Glu <- list(Ae = calib_params_fwd$Ae_Glu, Ea = calib_params_fwd$Ea_Glu,
                     R = 0.001987, x = calib_params_fwd$x_Glu)

  DL_Asp_pred <- racemize_one_depth(depth_cm, am, tm, dT_yr = dT_yr,
                                    AAR_params = params_Asp)$DL
  DL_Glu_pred <- racemize_one_depth(depth_cm, am, tm, dT_yr = dT_yr,
                                    AAR_params = params_Glu)$DL

  if (!is.finite(DL_Asp_pred) || !is.finite(DL_Glu_pred)) return(-Inf)

  resid_Asp <- DL_Asp_obs - DL_Asp_pred
  resid_Glu <- DL_Glu_obs - DL_Glu_pred

  .dbvnorm_log(resid_Asp, resid_Glu, calib_params_fwd$sigma_Asp,
              calib_params_fwd$sigma_Glu, calib_params_fwd$rho)
}


#' Log-posterior for forward-model two-acid age inversion (Stage 2)
#'
#' @param params Named numeric vector with element `log_age`.
#' @param calib_params_fwd List of calibration constants, see
#'   \code{\link{log_lik_two_acid_age_fwd}}.
#' @param DL_Asp_obs,DL_Glu_obs Observed D/L ratios for the sample.
#' @param temp_C Representative bottom-water temperature for the sample.
#' @param depth_cm Placeholder depth (cm). Default `100`.
#' @param dT_yr Integration timestep (yr). Default `2000`.
#' @param range_ka Numeric length-2 vector passed to
#'   \code{\link{log_prior_age}}.
#'
#' @return Scalar log-posterior.
#'
#' @export
log_posterior_two_acid_age_fwd <- function(params, calib_params_fwd,
                                           DL_Asp_obs, DL_Glu_obs, temp_C,
                                           depth_cm = 100, dT_yr = 2000,
                                           range_ka = c(0.5, 3000)) {
  log_age <- params[["log_age"]]
  lp <- log_prior_age(log_age, range_ka = range_ka)
  if (!is.finite(lp)) return(-Inf)

  lp + log_lik_two_acid_age_fwd(log_age, calib_params_fwd, DL_Asp_obs,
                                DL_Glu_obs, temp_C, depth_cm = depth_cm,
                                dT_yr = dT_yr)
}


#' Run the forward-model two-acid age inversion for one sample
#'
#' Convenience wrapper around \code{\link{run_mcmc}} for
#' \code{\link{log_posterior_two_acid_age_fwd}}. Produces the same
#' single-parameter (`log_age`) chain shape as
#' \code{\link{invert_two_acid_age}}, so \code{\link{summarize_age_posterior}}
#' and \code{\link{plot_age_posterior}} apply unchanged.
#'
#' @param DL_Asp_obs,DL_Glu_obs Observed D/L ratios for the sample.
#' @param calib_params_fwd List of calibration constants, see
#'   \code{\link{log_lik_two_acid_age_fwd}}.
#' @param temp_C Representative bottom-water temperature for the sample.
#' @param depth_cm Placeholder depth (cm). Default `100`.
#' @param n_iter Integer. Total MCMC iterations. Default `20000`.
#' @param proposal_sd Proposal SD for `log_age`. Default `0.3`.
#' @param dT_yr Integration timestep (yr). Default `2000`.
#' @param range_ka Numeric length-2 vector passed to
#'   \code{\link{log_prior_age}}.
#'
#' @return List returned by \code{\link{run_mcmc}}.
#'
#' @export
invert_two_acid_age_fwd <- function(DL_Asp_obs, DL_Glu_obs, calib_params_fwd,
                                    temp_C, depth_cm = 100, n_iter = 20000,
                                    proposal_sd = 0.3, dT_yr = 2000,
                                    range_ka = c(0.5, 3000)) {
  # Initialise from a crude age guess: whichever acid's forward curve at a
  # broad grid of ages comes closest to the observed D/L
  age_grid <- exp(seq(log(range_ka[1]), log(range_ka[2]), length.out = 200))
  grid_ll  <- vapply(age_grid, function(a) {
    log_lik_two_acid_age_fwd(log(a), calib_params_fwd, DL_Asp_obs, DL_Glu_obs,
                             temp_C, depth_cm = depth_cm, dT_yr = dT_yr)
  }, numeric(1))
  init <- c(log_age = log(age_grid[which.max(grid_ll)]))

  run_mcmc(log_posterior_two_acid_age_fwd, init = init, n_iter = n_iter,
           proposal_sd = c(log_age = proposal_sd),
           calib_params_fwd = calib_params_fwd,
           DL_Asp_obs = DL_Asp_obs, DL_Glu_obs = DL_Glu_obs, temp_C = temp_C,
           depth_cm = depth_cm, dT_yr = dT_yr, range_ka = range_ka)
}


# =============================================================================
# Single-acid age inversion (for comparing against the two-acid joint model)
# =============================================================================

#' Log-likelihood for single-acid forward-model age inversion
#'
#' Same forward model as \code{\link{log_lik_two_acid_age_fwd}}, but using
#' only one acid's observation and its own (univariate) residual SD --
#' useful for quantifying how much a second acid narrows the age estimate.
#'
#' @param log_age Candidate value of \eqn{\log(\text{age in ka})}.
#' @param calib_params_fwd List of calibration constants, see
#'   \code{\link{log_lik_two_acid_age_fwd}}.
#' @param DL_obs Observed D/L ratio for the chosen acid.
#' @param acid Character, `"Asp"` or `"Glu"`.
#' @param temp_C,depth_cm,dT_yr As in \code{\link{log_lik_two_acid_age_fwd}}.
#'
#' @return Scalar log-likelihood.
#'
#' @importFrom stats dnorm
#' @export
log_lik_one_acid_age_fwd <- function(log_age, calib_params_fwd, DL_obs,
                                     acid = c("Asp", "Glu"), temp_C,
                                     depth_cm = 100, dT_yr = 2000) {
  acid   <- match.arg(acid)
  age_ka <- exp(log_age)
  am     <- rate_to_age_model(rate_cm_per_ka = depth_cm / age_ka,
                              max_depth_cm   = depth_cm)
  tm     <- simulate_temp_history(temp_C, age_ka)

  params <- list(Ae = calib_params_fwd[[paste0("Ae_", acid)]],
                 Ea = calib_params_fwd[[paste0("Ea_", acid)]],
                 R  = 0.001987,
                 x  = calib_params_fwd[[paste0("x_", acid)]])
  sigma  <- calib_params_fwd[[paste0("sigma_", acid)]]

  DL_pred <- racemize_one_depth(depth_cm, am, tm, dT_yr = dT_yr,
                                AAR_params = params)$DL
  if (!is.finite(DL_pred)) return(-Inf)

  stats::dnorm(DL_obs, DL_pred, sigma, log = TRUE)
}


#' Log-posterior for single-acid forward-model age inversion
#'
#' @param params Named numeric vector with element `log_age`.
#' @inheritParams log_lik_one_acid_age_fwd
#' @param range_ka Numeric length-2 vector passed to
#'   \code{\link{log_prior_age}}.
#'
#' @return Scalar log-posterior.
#'
#' @export
log_posterior_one_acid_age_fwd <- function(params, calib_params_fwd, DL_obs,
                                           acid = c("Asp", "Glu"), temp_C,
                                           depth_cm = 100, dT_yr = 2000,
                                           range_ka = c(0.5, 3000)) {
  acid    <- match.arg(acid)
  log_age <- params[["log_age"]]
  lp <- log_prior_age(log_age, range_ka = range_ka)
  if (!is.finite(lp)) return(-Inf)

  lp + log_lik_one_acid_age_fwd(log_age, calib_params_fwd, DL_obs, acid,
                                temp_C, depth_cm = depth_cm, dT_yr = dT_yr)
}


#' Run the single-acid forward-model age inversion for one sample
#'
#' @inheritParams log_lik_one_acid_age_fwd
#' @param n_iter Integer. Total MCMC iterations. Default `20000`.
#' @param proposal_sd Proposal SD for `log_age`. Default `0.3`.
#' @param range_ka Numeric length-2 vector passed to
#'   \code{\link{log_prior_age}}.
#'
#' @return List returned by \code{\link{run_mcmc}}.
#'
#' @export
invert_one_acid_age_fwd <- function(DL_obs, acid = c("Asp", "Glu"),
                                    calib_params_fwd, temp_C, depth_cm = 100,
                                    n_iter = 20000, proposal_sd = 0.3,
                                    dT_yr = 2000, range_ka = c(0.5, 3000)) {
  acid <- match.arg(acid)
  age_grid <- exp(seq(log(range_ka[1]), log(range_ka[2]), length.out = 200))
  grid_ll  <- vapply(age_grid, function(a) {
    log_lik_one_acid_age_fwd(log(a), calib_params_fwd, DL_obs, acid, temp_C,
                             depth_cm = depth_cm, dT_yr = dT_yr)
  }, numeric(1))
  init <- c(log_age = log(age_grid[which.max(grid_ll)]))

  run_mcmc(log_posterior_one_acid_age_fwd, init = init, n_iter = n_iter,
           proposal_sd = c(log_age = proposal_sd),
           calib_params_fwd = calib_params_fwd, DL_obs = DL_obs, acid = acid,
           temp_C = temp_C, depth_cm = depth_cm, dT_yr = dT_yr,
           range_ka = range_ka)
}


# =============================================================================
# Age inversion with temperature-history uncertainty folded in
# =============================================================================

#' Log-likelihood for two-acid age inversion with an uncertain temperature
#'
#' As \code{\link{log_lik_two_acid_age_fwd}}, but the temperature history is
#' anchored on `temp_C_center + temp_offset` rather than a fixed value --
#' `temp_offset` is a free parameter representing genuine uncertainty in the
#' sample's effective bottom-water temperature (we know it's cold and in a
#' narrow range, but not precisely where in that range). Sampling
#' `temp_offset` jointly with `log_age` marginalizes over plausible
#' temperature histories instead of conditioning on one guess.
#'
#' @param log_age Candidate value of \eqn{\log(\text{age in ka})}.
#' @param temp_offset Candidate temperature offset from `temp_C_center`
#'   (deg C).
#' @param calib_params_fwd List of calibration constants, see
#'   \code{\link{log_lik_two_acid_age_fwd}}.
#' @param DL_Asp_obs,DL_Glu_obs Observed D/L ratios for the sample.
#' @param temp_C_center Central/recorded bottom-water temperature (deg C).
#' @param depth_cm,dT_yr As in \code{\link{log_lik_two_acid_age_fwd}}.
#'
#' @return Scalar log-likelihood.
#'
#' @export
log_lik_two_acid_age_tempunc <- function(log_age, temp_offset, calib_params_fwd,
                                         DL_Asp_obs, DL_Glu_obs, temp_C_center,
                                         depth_cm = 100, dT_yr = 2000) {
  log_lik_two_acid_age_fwd(log_age, calib_params_fwd, DL_Asp_obs, DL_Glu_obs,
                           temp_C = temp_C_center + temp_offset,
                           depth_cm = depth_cm, dT_yr = dT_yr)
}


#' Log-posterior for two-acid age inversion with temperature uncertainty
#'
#' @param params Named numeric vector with elements `log_age`, `temp_offset`.
#' @param calib_params_fwd List of calibration constants, see
#'   \code{\link{log_lik_two_acid_age_fwd}}.
#' @param DL_Asp_obs,DL_Glu_obs Observed D/L ratios for the sample.
#' @param temp_C_center Central/recorded bottom-water temperature (deg C).
#' @param depth_cm,dT_yr As in \code{\link{log_lik_two_acid_age_fwd}}.
#' @param range_ka Numeric length-2 vector passed to
#'   \code{\link{log_prior_age}}.
#' @param temp_sd Prior SD on `temp_offset` (deg C). Default `1`, reflecting
#'   a plausible cold/narrow-but-not-precisely-known range.
#'
#' @return Scalar log-posterior.
#'
#' @importFrom stats dnorm
#' @export
log_posterior_two_acid_age_tempunc <- function(params, calib_params_fwd,
                                               DL_Asp_obs, DL_Glu_obs,
                                               temp_C_center, depth_cm = 100,
                                               dT_yr = 2000,
                                               range_ka = c(0.5, 3000),
                                               temp_sd = 1) {
  log_age     <- params[["log_age"]]
  temp_offset <- params[["temp_offset"]]

  lp <- log_prior_age(log_age, range_ka = range_ka) +
    stats::dnorm(temp_offset, 0, temp_sd, log = TRUE)
  if (!is.finite(lp)) return(-Inf)

  lp + log_lik_two_acid_age_tempunc(log_age, temp_offset, calib_params_fwd,
                                    DL_Asp_obs, DL_Glu_obs, temp_C_center,
                                    depth_cm = depth_cm, dT_yr = dT_yr)
}


#' Run the two-acid age inversion with temperature uncertainty folded in
#'
#' @param DL_Asp_obs,DL_Glu_obs Observed D/L ratios for the sample.
#' @param calib_params_fwd List of calibration constants, see
#'   \code{\link{log_lik_two_acid_age_fwd}}.
#' @param temp_C_center Central/recorded bottom-water temperature (deg C).
#' @param depth_cm Placeholder depth (cm). Default `100`.
#' @param n_iter Integer. Total MCMC iterations. Default `20000`.
#' @param proposal_sd Named numeric vector of proposal SDs for `log_age` and
#'   `temp_offset`. Default `c(log_age = 0.3, temp_offset = 0.3)`.
#' @param dT_yr Integration timestep (yr). Default `2000`.
#' @param range_ka Numeric length-2 vector passed to
#'   \code{\link{log_prior_age}}.
#' @param temp_sd Prior SD on `temp_offset` (deg C). Default `1`.
#'
#' @return List returned by \code{\link{run_mcmc}}.
#'
#' @export
invert_two_acid_age_tempunc <- function(DL_Asp_obs, DL_Glu_obs, calib_params_fwd,
                                        temp_C_center, depth_cm = 100,
                                        n_iter = 20000,
                                        proposal_sd = c(log_age = 0.3,
                                                        temp_offset = 0.3),
                                        dT_yr = 2000, range_ka = c(0.5, 3000),
                                        temp_sd = 1) {
  init <- c(log_age = log(sqrt(range_ka[1] * range_ka[2])), temp_offset = 0)
  ll0 <- log_lik_two_acid_age_fwd(init[["log_age"]], calib_params_fwd,
                                  DL_Asp_obs, DL_Glu_obs, temp_C_center,
                                  depth_cm = depth_cm, dT_yr = dT_yr)
  if (!is.finite(ll0)) {
    age_grid <- exp(seq(log(range_ka[1]), log(range_ka[2]), length.out = 200))
    grid_ll  <- vapply(age_grid, function(a) {
      log_lik_two_acid_age_fwd(log(a), calib_params_fwd, DL_Asp_obs,
                               DL_Glu_obs, temp_C_center, depth_cm = depth_cm,
                               dT_yr = dT_yr)
    }, numeric(1))
    init[["log_age"]] <- log(age_grid[which.max(grid_ll)])
  }

  run_mcmc(log_posterior_two_acid_age_tempunc, init = init, n_iter = n_iter,
           proposal_sd = proposal_sd, calib_params_fwd = calib_params_fwd,
           DL_Asp_obs = DL_Asp_obs, DL_Glu_obs = DL_Glu_obs,
           temp_C_center = temp_C_center, depth_cm = depth_cm, dT_yr = dT_yr,
           range_ka = range_ka, temp_sd = temp_sd)
}