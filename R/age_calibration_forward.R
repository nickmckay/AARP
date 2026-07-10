# =============================================================================
# Forward-model two-acid age calibration: shared kinetics for Asp and Glu
# =============================================================================
#
# Unlike R/age_calibration.R (which fits separate empirical curve forms --
# TDK for Asp, SPK for Glu -- by regressing log(age) on log(f(D/L))), this
# file fits ONE shared Arrhenius + Bada power-law kinetic model (the same
# form already used for isoleucine elsewhere in the package) to each acid
# independently, forward-simulated through the package's existing
# age-depth-model / depth-series machinery (`make_age_model()`,
# `rate_to_age_model()`, `racemize_one_depth()`) rather than a single-shot
# constant-temperature shortcut. Asp and Glu get their own `Ae`, `Ea`, `x`,
# but are evaluated with the identical forward-model code path.

# =============================================================================
# Per-core age-depth models
# =============================================================================

#' Build per-sample age-depth models from grouped core data
#'
#' Groups calibration rows by `(study, core)` and constructs an age-depth
#' model for each row. Cores with two or more distinct depths get a real
#' interpolated \code{\link{make_age_model}} (duplicate depths are collapsed
#' by averaging their ages first); cores with only one usable depth fall
#' back to \code{\link{rate_to_age_model}}, a straight line from the core
#' top (0 cm, 0 ka) through that single point.
#'
#' @param calibration_data Data frame with columns `study`, `core`,
#'   `depth_mbsf`, `age_ka` (one row per calibration sample).
#'
#' @return A list, one element per row of `calibration_data` (same order),
#'   each an age-depth model as returned by \code{\link{make_age_model}} /
#'   \code{\link{rate_to_age_model}}.
#'
#' @export
build_core_age_models <- function(calibration_data) {
  core_key <- paste(calibration_data$study, calibration_data$core)
  age_models <- vector("list", nrow(calibration_data))

  for (key in unique(core_key)) {
    idx <- which(core_key == key)
    depth_cm <- calibration_data$depth_mbsf[idx] * 100
    age_ka   <- calibration_data$age_ka[idx]

    # Collapse duplicate depths by averaging their ages
    agg      <- stats::aggregate(age_ka, by = list(depth_cm = depth_cm), FUN = mean)
    depth_u  <- agg$depth_cm
    age_u    <- agg$x
    ord      <- order(depth_u)
    depth_u  <- depth_u[ord]
    age_u    <- age_u[ord]

    if (length(depth_u) >= 2) {
      if (depth_u[1] > 0) {
        depth_u <- c(0, depth_u)
        age_u   <- c(0, age_u)
      }
      am <- make_age_model(depth_cm = depth_u, age_ka = age_u)
    } else {
      am <- rate_to_age_model(rate_cm_per_ka = depth_u[1] / age_u[1],
                              max_depth_cm   = depth_u[1])
    }

    age_models[idx] <- list(am)
  }

  age_models
}


# =============================================================================
# Simulated placeholder temperature history
# =============================================================================

#' Simulate a placeholder glacial-interglacial temperature history
#'
#' The calibration dataset records only a single "representative" bottom
#' water temperature per sample -- not an actual paleotemperature
#' reconstruction. This generates a mild synthetic oscillation around that
#' recorded value, giving the depth-series time integration a temperature
#' *history* to work with rather than a constant. The amplitude and period
#' are placeholders (order-of-magnitude plausible for glacial-interglacial
#' deep polar ocean variability), not a real reconstruction, and are meant
#' to be replaced with actual per-core paleotemperature histories when
#' available.
#'
#' @param temp_C Recorded (present-day / representative) bottom water
#'   temperature (deg C).
#' @param max_age_ka Oldest age (ka BP) the history needs to cover.
#' @param amplitude Oscillation amplitude (deg C). Default `1`.
#' @param period_ka Oscillation period (ka). Default `100` (approximate
#'   Pleistocene glacial-interglacial pacing).
#' @param resolution_ka Spacing of control points (ka). Default `5`.
#'
#' @return A data frame with columns `age_ka`, `temp_C`, suitable for
#'   \code{\link{get_temp_at_time}} / \code{\link{racemize_one_depth}}.
#'
#' @export
simulate_temp_history <- function(temp_C, max_age_ka, amplitude = 1,
                                  period_ka = 100, resolution_ka = 5) {
  age_ka <- seq(0, max(max_age_ka, resolution_ka), by = resolution_ka)
  data.frame(
    age_ka = age_ka,
    temp_C = temp_C + amplitude * sin(2 * pi * age_ka / period_ka)
  )
}


# =============================================================================
# Precompute (age model + temperature history) once per calibration sample
# =============================================================================

#' Precompute per-sample forward-model inputs for calibration
#'
#' Builds the (depth, age-model, temperature-history) triple for each
#' calibration sample once, so the expensive parts of the forward model
#' (age-depth interpolation, temperature-history construction) are not
#' repeated on every MCMC iteration -- only the kinetic parameters change
#' between iterations.
#'
#' @param calibration_data Data frame with columns `study`, `core`,
#'   `depth_mbsf`, `age_ka`, `temp_C`, `DL_Asp`, `DL_Glu`. Rows with missing
#'   `temp_C` are dropped.
#' @param amplitude,period_ka,resolution_ka Passed to
#'   \code{\link{simulate_temp_history}}.
#'
#' @return A list with elements `depth_cm`, `age_model` (list), `temp_model`
#'   (list), `DL_Asp`, `DL_Glu`, one entry/row per retained sample.
#'
#' @export
precompute_forward_inputs <- function(calibration_data, amplitude = 1,
                                      period_ka = 100, resolution_ka = 5) {
  cal <- calibration_data[!is.na(calibration_data$temp_C), ]
  age_models <- build_core_age_models(cal)
  temp_models <- lapply(seq_len(nrow(cal)), function(i) {
    simulate_temp_history(cal$temp_C[i], cal$age_ka[i],
                          amplitude = amplitude, period_ka = period_ka,
                          resolution_ka = resolution_ka)
  })

  list(
    depth_cm   = cal$depth_mbsf * 100,
    age_model  = age_models,
    temp_model = temp_models,
    DL_Asp     = cal$DL_Asp,
    DL_Glu     = cal$DL_Glu
  )
}


# =============================================================================
# Stage 1: forward-model calibration fit (shared kinetics, Asp + Glu)
# =============================================================================

#' Forward-model D/L for a set of precomputed calibration samples
#'
#' Runs \code{\link{racemize_one_depth}} once per sample for a given set of
#' kinetic parameters, reusing the precomputed age-depth models and
#' temperature histories from \code{\link{precompute_forward_inputs}}.
#'
#' @param precomputed List returned by \code{\link{precompute_forward_inputs}}.
#' @param AAR_params Named list with elements `Ae`, `Ea`, `R`, `x`.
#' @param dT_yr Integration timestep (yr). Default `1500`; coarse but
#'   accurate here because the abiotic accumulation is linear in time within
#'   each step (see package vignettes for a convergence check).
#'
#' @return Numeric vector of predicted D/L, one per precomputed sample.
#'
#' @export
predict_forward_DL <- function(precomputed, AAR_params, dT_yr = 1500) {
  n   <- length(precomputed$depth_cm)
  out <- numeric(n)
  for (i in seq_len(n)) {
    out[i] <- racemize_one_depth(precomputed$depth_cm[i],
                                 precomputed$age_model[[i]],
                                 precomputed$temp_model[[i]],
                                 dT_yr      = dT_yr,
                                 AAR_params = AAR_params)$DL
  }
  out
}


#' Log-likelihood for the forward-model two-acid calibration
#'
#' Forward-simulates D/L for both acids at every calibration sample using a
#' single shared kinetic model form (Arrhenius + Bada power-law, via
#' \code{\link{racemize_one_depth}}), each acid with its own `Ae`, `Ea`,
#' `x`, and evaluates a correlated bivariate-normal likelihood on the D/L
#' residuals (paralleling \code{\link{log_lik_age_calibration}}, but in D/L
#' space rather than log-age space, since prediction now happens forward
#' from age to D/L rather than via an invertible empirical curve).
#'
#' @param log_Ae_Asp,Ea_Asp,x_Asp,log_sigma_Asp Aspartic acid kinetic
#'   parameters and log residual SD (D/L scale).
#' @param log_Ae_Glu,Ea_Glu,x_Glu,log_sigma_Glu Glutamic acid kinetic
#'   parameters and log residual SD.
#' @param z_rho Fisher-z transform of the residual correlation.
#' @param precomputed List returned by \code{\link{precompute_forward_inputs}}.
#' @param dT_yr Integration timestep (yr). Default `1500`.
#'
#' @return Scalar log-likelihood.
#'
#' @export
log_lik_age_calibration_fwd <- function(log_Ae_Asp, Ea_Asp, x_Asp, log_sigma_Asp,
                                        log_Ae_Glu, Ea_Glu, x_Glu, log_sigma_Glu,
                                        z_rho, precomputed, dT_yr = 1500) {
  if (x_Asp <= 0 || x_Glu <= 0 || Ea_Asp <= 0 || Ea_Glu <= 0) return(-Inf)

  sigma_Asp <- exp(log_sigma_Asp)
  sigma_Glu <- exp(log_sigma_Glu)
  rho       <- tanh(z_rho)

  params_Asp <- list(Ae = exp(log_Ae_Asp), Ea = Ea_Asp, R = 0.001987, x = x_Asp)
  params_Glu <- list(Ae = exp(log_Ae_Glu), Ea = Ea_Glu, R = 0.001987, x = x_Glu)

  DL_Asp_pred <- predict_forward_DL(precomputed, params_Asp, dT_yr = dT_yr)
  DL_Glu_pred <- predict_forward_DL(precomputed, params_Glu, dT_yr = dT_yr)

  if (any(!is.finite(c(DL_Asp_pred, DL_Glu_pred)))) return(-Inf)

  resid_Asp <- precomputed$DL_Asp - DL_Asp_pred
  resid_Glu <- precomputed$DL_Glu - DL_Glu_pred

  sum(.dbvnorm_log(resid_Asp, resid_Glu, sigma_Asp, sigma_Glu, rho))
}


#' Log-prior for the forward-model two-acid calibration
#'
#' `log_Ae_*` and `x_*` are weakly informative, centred loosely on the
#' package's existing (isoleucine) defaults (\code{\link{default_AAR_params}}).
#' `Ea_Asp`/`Ea_Glu` use an **informed** prior centred in the ~25-35 kcal/mol
#' range typical of published racemization activation energies for
#' calcite-hosted amino acids -- deliberately more informative than the
#' other parameters, because the calibration dataset's narrow natural
#' temperature range (~3.5 deg C) cannot identify Ea on its own (see
#' `vignette("two-acid-age-forward-model")`).
#'
#' @inheritParams log_lik_age_calibration_fwd
#'
#' @return Scalar log-prior.
#'
#' @export
log_prior_age_calibration_fwd <- function(log_Ae_Asp, Ea_Asp, x_Asp, log_sigma_Asp,
                                          log_Ae_Glu, Ea_Glu, x_Glu, log_sigma_Glu,
                                          z_rho) {
  if (x_Asp <= 0 || x_Glu <= 0 || Ea_Asp <= 0 || Ea_Glu <= 0) return(-Inf)
  stats::dnorm(log_Ae_Asp,    42,        5,   log = TRUE) +
  stats::dnorm(Ea_Asp,        30,        3,   log = TRUE) +
  stats::dnorm(x_Asp,         3,         1,   log = TRUE) +
  stats::dnorm(log_sigma_Asp, log(0.03), 1.5, log = TRUE) +
  stats::dnorm(log_Ae_Glu,    42,        5,   log = TRUE) +
  stats::dnorm(Ea_Glu,        30,        3,   log = TRUE) +
  stats::dnorm(x_Glu,         3,         1,   log = TRUE) +
  stats::dnorm(log_sigma_Glu, log(0.03), 1.5, log = TRUE) +
  stats::dnorm(z_rho,         0,         1,   log = TRUE)
}


#' Log-posterior for the forward-model two-acid calibration (Stage 1)
#'
#' @param params Named numeric vector: `log_Ae_Asp`, `Ea_Asp`, `x_Asp`,
#'   `log_sigma_Asp`, `log_Ae_Glu`, `Ea_Glu`, `x_Glu`, `log_sigma_Glu`,
#'   `z_rho`.
#' @param precomputed List returned by \code{\link{precompute_forward_inputs}}.
#' @param dT_yr Integration timestep (yr). Default `1500`.
#'
#' @return Scalar log-posterior.
#'
#' @export
log_posterior_age_calibration_fwd <- function(params, precomputed, dT_yr = 1500) {
  lp <- log_prior_age_calibration_fwd(
    log_Ae_Asp    = params[["log_Ae_Asp"]],
    Ea_Asp        = params[["Ea_Asp"]],
    x_Asp         = params[["x_Asp"]],
    log_sigma_Asp = params[["log_sigma_Asp"]],
    log_Ae_Glu    = params[["log_Ae_Glu"]],
    Ea_Glu        = params[["Ea_Glu"]],
    x_Glu         = params[["x_Glu"]],
    log_sigma_Glu = params[["log_sigma_Glu"]],
    z_rho         = params[["z_rho"]]
  )
  if (!is.finite(lp)) return(-Inf)

  lp + log_lik_age_calibration_fwd(
    log_Ae_Asp    = params[["log_Ae_Asp"]],
    Ea_Asp        = params[["Ea_Asp"]],
    x_Asp         = params[["x_Asp"]],
    log_sigma_Asp = params[["log_sigma_Asp"]],
    log_Ae_Glu    = params[["log_Ae_Glu"]],
    Ea_Glu        = params[["Ea_Glu"]],
    x_Glu         = params[["x_Glu"]],
    log_sigma_Glu = params[["log_sigma_Glu"]],
    z_rho         = params[["z_rho"]],
    precomputed   = precomputed,
    dT_yr         = dT_yr
  )
}


#' Extract posterior-mean forward-model calibration parameters
#'
#' @param mcmc_out List returned by \code{\link{run_mcmc}} using
#'   \code{\link{log_posterior_age_calibration_fwd}}.
#' @param burnin Integer. Burn-in iterations to discard. Default `1000`.
#'
#' @return A named list with elements `Ae_Asp`, `Ea_Asp`, `x_Asp`,
#'   `sigma_Asp`, `Ae_Glu`, `Ea_Glu`, `x_Glu`, `sigma_Glu`, `rho`.
#'
#' @export
summarize_calibration_fwd <- function(mcmc_out, burnin = 1000) {
  post <- utils::tail(mcmc_out$samples, nrow(mcmc_out$samples) - burnin)
  m    <- colMeans(post)
  list(
    Ae_Asp    = exp(m[["log_Ae_Asp"]]),
    Ea_Asp    = m[["Ea_Asp"]],
    x_Asp     = m[["x_Asp"]],
    sigma_Asp = exp(m[["log_sigma_Asp"]]),
    Ae_Glu    = exp(m[["log_Ae_Glu"]]),
    Ea_Glu    = m[["Ea_Glu"]],
    x_Glu     = m[["x_Glu"]],
    sigma_Glu = exp(m[["log_sigma_Glu"]]),
    rho       = tanh(m[["z_rho"]])
  )
}


#' Plot forward-model predicted vs. observed D/L for the calibration set
#'
#' @param mcmc_out List returned by \code{\link{run_mcmc}} using
#'   \code{\link{log_posterior_age_calibration_fwd}}.
#' @param precomputed List returned by \code{\link{precompute_forward_inputs}}.
#' @param burnin Integer. Burn-in iterations to discard. Default `1000`.
#' @param dT_yr Integration timestep (yr). Default `1500`.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline facet_wrap labs
#'   theme_bw
#' @importFrom rlang .data
#' @export
plot_age_calibration_fit_fwd <- function(mcmc_out, precomputed, burnin = 1000,
                                         dT_yr = 1500) {
  cp <- summarize_calibration_fwd(mcmc_out, burnin = burnin)

  params_Asp <- list(Ae = cp$Ae_Asp, Ea = cp$Ea_Asp, R = 0.001987, x = cp$x_Asp)
  params_Glu <- list(Ae = cp$Ae_Glu, Ea = cp$Ea_Glu, R = 0.001987, x = cp$x_Glu)

  df <- rbind(
    data.frame(acid = "Aspartic acid", observed = precomputed$DL_Asp,
               predicted = predict_forward_DL(precomputed, params_Asp, dT_yr = dT_yr)),
    data.frame(acid = "Glutamic acid", observed = precomputed$DL_Glu,
               predicted = predict_forward_DL(precomputed, params_Glu, dT_yr = dT_yr))
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$observed, y = .data$predicted)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey50",
                         linetype = "dashed") +
    ggplot2::geom_point(size = 1.5, alpha = 0.6, colour = "steelblue") +
    ggplot2::facet_wrap(~acid) +
    ggplot2::labs(x = "Observed D/L", y = "Forward-model predicted D/L",
                  title = "Forward-model calibration fit") +
    ggplot2::theme_bw()
}


#' Plot the Ae-Ea posterior compensation ridge
#'
#' Amino acid racemization rate constants show a classic Arrhenius
#' "compensation effect": with a narrow natural temperature range, `Ae` and
#' `Ea` become strongly (positively) correlated in the posterior, because
#' many `(Ae, Ea)` combinations give nearly the same rate constant over the
#' small span of temperatures actually observed. This plot makes that
#' identifiability limitation visible directly.
#'
#' @param mcmc_out List returned by \code{\link{run_mcmc}} using
#'   \code{\link{log_posterior_age_calibration_fwd}}.
#' @param burnin Integer. Burn-in iterations to discard. Default `1000`.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap labs theme_bw
#' @importFrom rlang .data
#' @export
plot_ae_ea_posterior <- function(mcmc_out, burnin = 1000) {
  post <- as.data.frame(utils::tail(mcmc_out$samples,
                                    nrow(mcmc_out$samples) - burnin))
  df <- rbind(
    data.frame(acid = "Aspartic acid", log_Ae = post$log_Ae_Asp, Ea = post$Ea_Asp),
    data.frame(acid = "Glutamic acid", log_Ae = post$log_Ae_Glu, Ea = post$Ea_Glu)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$Ea, y = .data$log_Ae)) +
    ggplot2::geom_point(size = 0.6, alpha = 0.3, colour = "steelblue") +
    ggplot2::facet_wrap(~acid, scales = "free") +
    ggplot2::labs(x = "Ea (kcal/mol)", y = "log(Ae)",
                  title = "Ae-Ea posterior compensation ridge",
                  subtitle = "Narrow natural temperature range weakly identifies Ea") +
    ggplot2::theme_bw()
}


#' Plot prior vs. posterior densities for all Stage 1 forward-model parameters
#'
#' Overlays each parameter's prior density (as specified in
#' \code{\link{log_prior_age_calibration_fwd}}) with a posterior density
#' estimated from the post-burn-in MCMC samples. A posterior much narrower
#' than its prior means the data are informative for that parameter; a
#' posterior that closely tracks the prior (as for `Ea_Asp`/`Ea_Glu`, see
#' \code{\link{plot_ae_ea_posterior}}) means the data alone cannot move it
#' far from the assumed prior belief.
#'
#' @param mcmc_out List returned by \code{\link{run_mcmc}} using
#'   \code{\link{log_posterior_age_calibration_fwd}}.
#' @param burnin Integer. Burn-in iterations to discard. Default `1000`.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_density geom_area geom_line facet_wrap
#'   labs theme_bw
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @export
plot_prior_posterior_fwd <- function(mcmc_out, burnin = 1000) {
  prior_pars <- list(
    log_Ae_Asp    = c(mean = 42,        sd = 5),
    Ea_Asp        = c(mean = 30,        sd = 3),
    x_Asp         = c(mean = 3,         sd = 1),
    log_sigma_Asp = c(mean = log(0.03), sd = 1.5),
    log_Ae_Glu    = c(mean = 42,        sd = 5),
    Ea_Glu        = c(mean = 30,        sd = 3),
    x_Glu         = c(mean = 3,         sd = 1),
    log_sigma_Glu = c(mean = log(0.03), sd = 1.5),
    z_rho         = c(mean = 0,         sd = 1)
  )

  post <- as.data.frame(utils::tail(mcmc_out$samples,
                                    nrow(mcmc_out$samples) - burnin))
  post_long <- tidyr::pivot_longer(post, dplyr::everything(),
                                   names_to = "parameter", values_to = "value")
  post_long$parameter <- factor(post_long$parameter, levels = names(prior_pars))

  prior_df <- do.call(rbind, lapply(names(prior_pars), function(nm) {
    p <- prior_pars[[nm]]
    v <- seq(p["mean"] - 4 * p["sd"], p["mean"] + 4 * p["sd"], length.out = 400)
    data.frame(parameter = nm, value = v,
               density   = stats::dnorm(v, p["mean"], p["sd"]))
  }))
  prior_df$parameter <- factor(prior_df$parameter, levels = names(prior_pars))

  ggplot2::ggplot() +
    ggplot2::geom_area(data = prior_df,
                       ggplot2::aes(x = .data$value, y = .data$density),
                       fill = "darkorange", alpha = 0.3) +
    ggplot2::geom_line(data = prior_df,
                       ggplot2::aes(x = .data$value, y = .data$density),
                       colour = "darkorange", linewidth = 0.8) +
    ggplot2::geom_density(data = post_long,
                          ggplot2::aes(x = .data$value),
                          fill = "steelblue", alpha = 0.4) +
    ggplot2::facet_wrap(~parameter, scales = "free", ncol = 3) +
    ggplot2::labs(x = "Value", y = "Density",
                  title = "Stage 1 forward model: prior (orange) vs. posterior (blue)") +
    ggplot2::theme_bw()
}
