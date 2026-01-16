# This script contains afunction that runs all other functions for one replication of the study
# this function should contain all the arguments that are defined in the previous sections
# ------------------------------------------------------------------------------------------------

run_one_rep_study <- function(
    rep_id,                                                 # replication index (set by outer loop)
    N,                                                      # sample size
    T,                                                      # number of waves
    k,                                                      # number of confounders
    R2_1,                                                   # confounder R^2 at wave 1
    target_sd,                                              # SD of time-variation in B
    scenarios,                                              # character vector: c("constant","stepwise",...)
    B_scenarios = NULL,                                     # optional: user-defined list of B trajectories
    A,                                                      # 2×2 autoregressive + cross-lag matrix
    Psi,                                                    # k×k confounder covariance
    rho_extra,                                              # extra covariance added to X,Y each wave
    models_to_run,                                          # e.g. c("clpm","riclpm","dpm","lbca","adj")
    base_seed = 1234                                        # base seed
){
  
  # set seed for this replication, use rep_id so it varies
  set.seed(base_seed + rep_id)

  # sample baseline confounder effects matrix B1
  B1 <- sample_B_linear(
    k    = k,
    R2_1 = R2_1
  )

  # prepare output list for each scenario
  out_list <- vector("list", length(scenarios))

  # loop over the scenarios
  for (j in seq_along(scenarios)) {
    
    scen <- scenarios[j]

    # choose trajectory generator
    # if user provided trajectories, use those
    if (!is.null(B_scenarios)) {

      # check that the scenario exists in the list
      if (!scen %in% names(B_scenarios)) {
        stop("Scenario not found in B_scenarios: ", scen)
      }

      # take the pre-defined trajectory
      B_list <- B_scenarios[[scen]]

    } else {

      # otherwise, fall back to generating trajectories from the sampled B1
      if (scen == "constant") {
        B_list <- generate_B_constant(B1, T)
      } else if (scen == "stepwise") {
        B_list <- generate_B_stepwise(B1, T, target_sd = target_sd)
      } else {
        stop("Unknown scenario: ", scen)
      }
    }

    # extract mean betas
    beta_X_vec <- sapply(B_list, function(Bt) mean(Bt[1, ]))
    beta_Y_vec <- sapply(B_list, function(Bt) mean(Bt[2, ]))
    beta_vec   <- beta_X_vec

    # try to simulate panel data
    df <- tryCatch(
      simulate_panel_data(
        N         = N,
        T         = T,
        A         = A,
        B_list    = B_list,
        Psi       = Psi,
        rho_extra = rho_extra
      ),
      error = function(e) NULL
    )

    if (is.null(df)) {
      
      # simulation failed return NA
      out_list[[j]] <- data.frame(
        run      = rep(rep_id, T),
        occasion = 1:T,
        scenario = scen,
        beta     = beta_vec,
        beta_X   = beta_X_vec,
        beta_Y   = beta_Y_vec,

        estXY_CLPM      = NA,
        estXY_RI_CLPM   = NA,
        estXY_DPM       = NA,
        estXY_CLPM_Adj  = NA,
        estXY_CLPM_LBCA = NA,

        estYX_CLPM      = NA,
        estYX_RI_CLPM   = NA,
        estYX_DPM       = NA,
        estYX_CLPM_Adj  = NA,
        estYX_CLPM_LBCA = NA,

        estA_CLPM      = NA,
        estA_RI_CLPM   = NA,
        estA_DPM       = NA,
        estA_CLPM_Adj  = NA,
        estA_CLPM_LBCA = NA,

        # AR Y
        estAY_CLPM      = NA,
        estAY_RI_CLPM   = NA,
        estAY_DPM       = NA,
        estAY_CLPM_Adj  = NA,
        estAY_CLPM_LBCA = NA,

        estRho_CLPM = NA,

        fail_CLPM      = TRUE,
        fail_RI_CLPM   = TRUE,
        fail_DPM       = TRUE,
        fail_CLPM_Adj  = TRUE,
        fail_CLPM_LBCA = TRUE,

        err_CLPM      = "sim failed",
        err_RI_CLPM   = "sim failed",
        err_DPM       = "sim failed",
        err_CLPM_Adj  = "sim failed",
        err_CLPM_LBCA = "sim failed",

        is_na_run = 1L
      )

      next
    }

    # build model strings
    model_clpm         <- build_clpm(T)
    model_riclpm       <- build_riclpm(T)
    model_dpm          <- build_dpm(T)
    model_clpm_with_Cs <- build_clpm_with_Cs(T, k)

    # fit the models safely
    res_clpm <- if ("clpm"   %in% models_to_run) safe_fit_clpm(model_clpm, df) else list(fit=NULL, err=NA)
    res_ric  <- if ("riclpm" %in% models_to_run) safe_fit_riclpm(model_riclpm, df) else list(fit=NULL, err=NA)
    res_dpm0 <- if ("dpm"    %in% models_to_run) safe_fit_dpm(model_dpm, df) else list(fit=NULL, err=NA)
    res_adj  <- if ("adj"    %in% models_to_run) safe_fit_clpm_C(model_clpm_with_Cs, df) else list(fit=NULL, err=NA)
    res_lbca <- if ("lbca"   %in% models_to_run) safe_fit_clpm_resid(model_clpm, df) else list(fit=NULL, err=NA)

    fit_clpm_raw <- res_clpm$fit
    fit_ric      <- res_ric$fit
    fit_dpm0     <- res_dpm0$fit
    fit_adj      <- res_adj$fit
    fit_lbca     <- res_lbca$fit

    # extract lagged parameters
    lag_raw  <- extract_lagged_parameters(fit_clpm_raw, T, "clpm")
    lag_ric  <- extract_lagged_parameters(fit_ric,       T, "riclpm")
    lag_dpm0 <- extract_lagged_parameters(fit_dpm0,      T, "dpm")
    lag_adj  <- extract_lagged_parameters(fit_adj,       T, "clpm")
    lag_lbca <- extract_lagged_parameters(fit_lbca,      T, "clpm")

    # residual correlations
    rho_clpm <- extract_rho_vec(fit_clpm_raw, T, "clpm")

    # assemble output row
    out_list[[j]] <- data.frame(

      run      = rep(rep_id, T),
      occasion = 1:T,
      scenario = scen,

      beta     = beta_vec,
      beta_X   = beta_X_vec,
      beta_Y   = beta_Y_vec,

      # cross-lag XY
      estXY_CLPM      = c(NA, lag_raw$xy),
      estXY_RI_CLPM   = c(NA, lag_ric$xy),
      estXY_DPM       = c(NA, lag_dpm0$xy),
      estXY_CLPM_Adj  = c(NA, lag_adj$xy),
      estXY_CLPM_LBCA = c(NA, lag_lbca$xy),

      # cross-lag YX
      estYX_CLPM      = c(NA, lag_raw$yx),
      estYX_RI_CLPM   = c(NA, lag_ric$yx),
      estYX_DPM       = c(NA, lag_dpm0$yx),
      estYX_CLPM_Adj  = c(NA, lag_adj$yx),
      estYX_CLPM_LBCA = c(NA, lag_lbca$yx),

      # autoregressive X
      estA_CLPM      = c(NA, lag_raw$ar_x),
      estA_RI_CLPM   = c(NA, lag_ric$ar_x),
      estA_DPM       = c(NA, lag_dpm0$ar_x),
      estA_CLPM_Adj  = c(NA, lag_adj$ar_x),
      estA_CLPM_LBCA = c(NA, lag_lbca$ar_x),

      # autoregressive Y 
      estAY_CLPM      = c(NA, lag_raw$ar_y),
      estAY_RI_CLPM   = c(NA, lag_ric$ar_y),
      estAY_DPM       = c(NA, lag_dpm0$ar_y),
      estAY_CLPM_Adj  = c(NA, lag_adj$ar_y),
      estAY_CLPM_LBCA = c(NA, lag_lbca$ar_y),

      # residual correlation
      estRho_CLPM = rho_clpm,

      # failure indicators
      fail_CLPM      = is.null(fit_clpm_raw),
      fail_RI_CLPM   = is.null(fit_ric),
      fail_DPM       = is.null(fit_dpm0),
      fail_CLPM_Adj  = is.null(fit_adj),
      fail_CLPM_LBCA = is.null(fit_lbca),

      # error messages
      err_CLPM      = rep(res_clpm$err,   T),
      err_RI_CLPM   = rep(res_ric$err,    T),
      err_DPM       = rep(res_dpm0$err,   T),
      err_CLPM_Adj  = rep(res_adj$err,   T),
      err_CLPM_LBCA = rep(res_lbca$err,  T),

      # NA run marker
      is_na_run = as.integer(all(is.na(c(
        lag_raw$xy, lag_ric$xy, lag_dpm0$xy, lag_adj$xy, lag_lbca$xy
      ))))
    )
  }

  dplyr::bind_rows(out_list)
}
