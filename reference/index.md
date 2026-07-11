# Package index

## Forward model

Core Arrhenius/power-law racemization model and thermal integration
utilities.

- [`arrhenius()`](https://nickmckay.github.io/AARP/reference/arrhenius.md)
  : Arrhenius rate constant for amino acid racemization
- [`dRacPL()`](https://nickmckay.github.io/AARP/reference/dRacPL.md) :
  Incremental racemization in power-law space
- [`racemize()`](https://nickmckay.github.io/AARP/reference/racemize.md)
  : Simulate D/L accumulation over a full temperature history
- [`racemizeSample()`](https://nickmckay.github.io/AARP/reference/racemizeSample.md)
  : Simulate racemization for a single downcore sample (with burial
  damping)
- [`racemize_one_depth()`](https://nickmckay.github.io/AARP/reference/racemize_one_depth.md)
  : Simulate AAR racemization for a single sediment sample at a given
  depth
- [`racemize_depth_series()`](https://nickmckay.github.io/AARP/reference/racemize_depth_series.md)
  : Simulate AAR racemization for a depth series of sediment samples
- [`predict_heating_DL()`](https://nickmckay.github.io/AARP/reference/predict_heating_DL.md)
  : Predict D/L for lab heating experiments
- [`get_temp_at_time()`](https://nickmckay.github.io/AARP/reference/get_temp_at_time.md)
  : Return lake bottom water temperature at one or more times
- [`diffusion_attenuation()`](https://nickmckay.github.io/AARP/reference/diffusion_attenuation.md)
  : Thermal diffusion attenuation factor for a periodic temperature
  signal
- [`dampSedTemps()`](https://nickmckay.github.io/AARP/reference/dampSedTemps.md)
  : Apply sediment burial temperature damping

## Age–depth and TOC models

Construct and query age–depth relationships and TOC profiles for
sediment cores.

- [`make_age_model()`](https://nickmckay.github.io/AARP/reference/make_age_model.md)
  : Build an age-depth model from tie-point depths and ages
- [`rate_to_age_model()`](https://nickmckay.github.io/AARP/reference/rate_to_age_model.md)
  : Build an age-depth model from a constant sedimentation rate
- [`make_toc_model()`](https://nickmckay.github.io/AARP/reference/make_toc_model.md)
  : Build a TOC-depth model for bacterial resetting

## Simulation

Generate synthetic heating-experiment and downcore D/L datasets for
testing and vignette demonstrations.

- [`sim_heating_data()`](https://nickmckay.github.io/AARP/reference/sim_heating_data.md)
  : Simulate lab heating experiment D/L data
- [`sim_downcore_data()`](https://nickmckay.github.io/AARP/reference/sim_downcore_data.md)
  : Simulate downcore D/L observations from a known-temperature lake
- [`default_AAR_params`](https://nickmckay.github.io/AARP/reference/default_AAR_params.md)
  : Default AAR kinetic parameters

## Bayesian inference

Log-likelihood, log-prior, and log-posterior functions for the two-stage
MCMC framework, plus the MCMC sampler itself.

- [`log_lik_heating()`](https://nickmckay.github.io/AARP/reference/log_lik_heating.md)
  : Log-likelihood for lab heating experiment data
- [`log_lik_downcore()`](https://nickmckay.github.io/AARP/reference/log_lik_downcore.md)
  : Log-likelihood for downcore D/L data at known temperature
- [`log_lik_temperature()`](https://nickmckay.github.io/AARP/reference/log_lik_temperature.md)
  : Log-likelihood for downcore D/L data at unknown temperature
- [`log_prior_kinetics()`](https://nickmckay.github.io/AARP/reference/log_prior_kinetics.md)
  : Log-prior for kinetic parameters
- [`log_prior_temperature()`](https://nickmckay.github.io/AARP/reference/log_prior_temperature.md)
  : Log-prior for temperature history control points
- [`log_posterior_kinetics()`](https://nickmckay.github.io/AARP/reference/log_posterior_kinetics.md)
  : Log-posterior for kinetic parameters (Stage 1)
- [`log_posterior_temperature()`](https://nickmckay.github.io/AARP/reference/log_posterior_temperature.md)
  : Log-posterior for temperature history (Stage 2)
- [`run_mcmc()`](https://nickmckay.github.io/AARP/reference/run_mcmc.md)
  : Random-walk Metropolis-Hastings MCMC sampler

## Visualisation & diagnostics

Trace plots, posterior density plots, temperature reconstructions, and
posterior predictive checks.

- [`plot_mcmc_chains()`](https://nickmckay.github.io/AARP/reference/plot_mcmc_chains.md)
  : Trace plots for MCMC output
- [`plot_posterior_kinetics()`](https://nickmckay.github.io/AARP/reference/plot_posterior_kinetics.md)
  : Posterior density plots for kinetic parameters
- [`plot_temperature_reconstruction()`](https://nickmckay.github.io/AARP/reference/plot_temperature_reconstruction.md)
  : Temperature history reconstruction plot
- [`posterior_predictive_DL()`](https://nickmckay.github.io/AARP/reference/posterior_predictive_DL.md)
  : Posterior predictive D/L depth profiles
- [`plotInputs()`](https://nickmckay.github.io/AARP/reference/plotInputs.md)
  : Plot temperature inputs for a damped sediment record

## Two-acid AAR age dating

Bayesian calibration and age inversion combining aspartic (TDK) and
glutamic (SPK) acid D/L ratios, as an alternative to GLS-combined
point-estimate age calculators.

- [`f_tdk()`](https://nickmckay.github.io/AARP/reference/f_tdk.md) :
  Predictor transform for the time-dependent kinetics (TDK) model
- [`f_spk()`](https://nickmckay.github.io/AARP/reference/f_spk.md) :
  Predictor transform for the simple power-law kinetics (SPK) model
- [`predict_age_asp()`](https://nickmckay.github.io/AARP/reference/predict_age_asp.md)
  : Predict age from D/L using the TDK (Asp) calibration curve
- [`predict_age_glu()`](https://nickmckay.github.io/AARP/reference/predict_age_glu.md)
  : Predict age from D/L using the SPK (Glu) calibration curve
- [`predict_DL_asp()`](https://nickmckay.github.io/AARP/reference/predict_DL_asp.md)
  : Predict D/L from age using the TDK (Asp) calibration curve
- [`predict_DL_glu()`](https://nickmckay.github.io/AARP/reference/predict_DL_glu.md)
  : Predict D/L from age using the SPK (Glu) calibration curve
- [`log_lik_age_calibration()`](https://nickmckay.github.io/AARP/reference/log_lik_age_calibration.md)
  : Log-likelihood for the joint two-acid age calibration
- [`log_prior_age_calibration()`](https://nickmckay.github.io/AARP/reference/log_prior_age_calibration.md)
  : Log-prior for the joint two-acid age calibration
- [`log_posterior_age_calibration()`](https://nickmckay.github.io/AARP/reference/log_posterior_age_calibration.md)
  : Log-posterior for the joint two-acid age calibration (Stage 1)
- [`summarize_calibration()`](https://nickmckay.github.io/AARP/reference/summarize_calibration.md)
  : Extract posterior-mean calibration parameters from a Stage 1 MCMC
  run
- [`plot_age_calibration_fit()`](https://nickmckay.github.io/AARP/reference/plot_age_calibration_fit.md)
  : Plot fitted Asp and Glu calibration curves against calibration data
- [`log_lik_two_acid_age()`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age.md)
  : Log-likelihood for two-acid age inversion
- [`log_prior_age()`](https://nickmckay.github.io/AARP/reference/log_prior_age.md)
  : Log-prior for sample age
- [`log_posterior_two_acid_age()`](https://nickmckay.github.io/AARP/reference/log_posterior_two_acid_age.md)
  : Log-posterior for two-acid age inversion (Stage 2)
- [`invert_two_acid_age()`](https://nickmckay.github.io/AARP/reference/invert_two_acid_age.md)
  : Run the two-acid Bayesian age inversion for one sample
- [`summarize_age_posterior()`](https://nickmckay.github.io/AARP/reference/summarize_age_posterior.md)
  : Posterior summary for a two-acid age inversion
- [`plot_age_posterior()`](https://nickmckay.github.io/AARP/reference/plot_age_posterior.md)
  : Plot the posterior age distribution for a two-acid age inversion

## Two-acid AAR age dating — shared forward model

Alternative to the empirical TDK/SPK curve fit above: fits one shared
Arrhenius + Bada power-law kinetic model per acid, forward-simulated
through real age-depth models and a placeholder temperature history, via
the same depth-series machinery used for lake-sediment paleothermometry.

- [`build_core_age_models()`](https://nickmckay.github.io/AARP/reference/build_core_age_models.md)
  : Build per-sample age-depth models from grouped core data
- [`simulate_temp_history()`](https://nickmckay.github.io/AARP/reference/simulate_temp_history.md)
  : Simulate a placeholder glacial-interglacial temperature history
- [`precompute_forward_inputs()`](https://nickmckay.github.io/AARP/reference/precompute_forward_inputs.md)
  : Precompute per-sample forward-model inputs for calibration
- [`predict_forward_DL()`](https://nickmckay.github.io/AARP/reference/predict_forward_DL.md)
  : Forward-model D/L for a set of precomputed calibration samples
- [`log_lik_age_calibration_fwd()`](https://nickmckay.github.io/AARP/reference/log_lik_age_calibration_fwd.md)
  : Log-likelihood for the forward-model two-acid calibration
- [`log_prior_age_calibration_fwd()`](https://nickmckay.github.io/AARP/reference/log_prior_age_calibration_fwd.md)
  : Log-prior for the forward-model two-acid calibration
- [`log_posterior_age_calibration_fwd()`](https://nickmckay.github.io/AARP/reference/log_posterior_age_calibration_fwd.md)
  : Log-posterior for the forward-model two-acid calibration (Stage 1)
- [`summarize_calibration_fwd()`](https://nickmckay.github.io/AARP/reference/summarize_calibration_fwd.md)
  : Extract posterior-mean forward-model calibration parameters
- [`plot_age_calibration_fit_fwd()`](https://nickmckay.github.io/AARP/reference/plot_age_calibration_fit_fwd.md)
  : Plot forward-model predicted vs. observed D/L for the calibration
  set
- [`plot_ae_ea_posterior()`](https://nickmckay.github.io/AARP/reference/plot_ae_ea_posterior.md)
  : Plot the Ae-Ea posterior compensation ridge
- [`plot_prior_posterior_fwd()`](https://nickmckay.github.io/AARP/reference/plot_prior_posterior_fwd.md)
  : Plot prior vs. posterior densities for all Stage 1 forward-model
  parameters
- [`log_lik_two_acid_age_fwd()`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md)
  : Log-likelihood for forward-model two-acid age inversion
- [`log_posterior_two_acid_age_fwd()`](https://nickmckay.github.io/AARP/reference/log_posterior_two_acid_age_fwd.md)
  : Log-posterior for forward-model two-acid age inversion (Stage 2)
- [`invert_two_acid_age_fwd()`](https://nickmckay.github.io/AARP/reference/invert_two_acid_age_fwd.md)
  : Run the forward-model two-acid age inversion for one sample
- [`log_lik_one_acid_age_fwd()`](https://nickmckay.github.io/AARP/reference/log_lik_one_acid_age_fwd.md)
  : Log-likelihood for single-acid forward-model age inversion
- [`log_posterior_one_acid_age_fwd()`](https://nickmckay.github.io/AARP/reference/log_posterior_one_acid_age_fwd.md)
  : Log-posterior for single-acid forward-model age inversion
- [`invert_one_acid_age_fwd()`](https://nickmckay.github.io/AARP/reference/invert_one_acid_age_fwd.md)
  : Run the single-acid forward-model age inversion for one sample
- [`log_lik_two_acid_age_tempunc()`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_tempunc.md)
  : Log-likelihood for two-acid age inversion with an uncertain
  temperature
- [`log_posterior_two_acid_age_tempunc()`](https://nickmckay.github.io/AARP/reference/log_posterior_two_acid_age_tempunc.md)
  : Log-posterior for two-acid age inversion with temperature
  uncertainty
- [`invert_two_acid_age_tempunc()`](https://nickmckay.github.io/AARP/reference/invert_two_acid_age_tempunc.md)
  : Run the two-acid age inversion with temperature uncertainty folded
  in

## Two-acid AAR age dating — reference GLS calculator

Faithful R port of the original spreadsheet age calculator (GLS
combination of independent Asp/Glu age estimates with softplus weight
smoothing), included so the two Bayesian approaches above can be
compared directly against the point-estimate method they build on.

- [`gls_calibration_constants`](https://nickmckay.github.io/AARP/reference/gls_calibration_constants.md)
  : Calibration constants for the GLS age calculator
- [`predict_age_gls()`](https://nickmckay.github.io/AARP/reference/predict_age_gls.md)
  : GLS-combined two-acid age with softplus weight smoothing
