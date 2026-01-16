# this script runs the wrapper function that executes one replication of the simulation study
# adding a progress bar and parallelization
# ------------------------------------------------------------------------------------------------

run_simulation_study1 <- function(
    reps,                                                                    # number of replications
    N,                                                                       # sample size
    T,                                                                       # number of waves
    k,                                                                       # number of confounders
    R2_1,                                                                    # confounder R^2 at wave 1
    target_sd,                                                               # SD of time-varying B
    scenarios,                                                               # e.g., c("constant","stepwise")
    B_scenarios = NULL,                                                      # optional: user-defined list of B trajectories
    A,                                                                       # 2×2 AR + cross-lag matrix
    Psi,                                                                     # k×k confounder covariance
    rho_extra,                                                               # extra covariance added to observations
    models_to_run,                                                           # c("clpm","riclpm","dpm","adj","lbca")
    cores = NULL,                                                            # default is detectCores()/2
    base_seed = 1234                                                         # master seed for reproducible reps
) {

  # if the number of cores is not specified, detect and use half of available cores
  if (is.null(cores)) {

    # detect and use half of available cores
    cores <- max(1, floor(parallel::detectCores() / 2))
  }

  # if cores is 1, run sequentially without parallelization
  if (cores == 1L) {

    # run sequentially
    results_list <- lapply(
      X = 1:reps,
      FUN = function(rep_id) {
        run_one_rep_study(
          rep_id        = rep_id,
          N             = N,
          T             = T,
          k             = k,
          R2_1          = R2_1,
          target_sd     = target_sd,
          scenarios     = scenarios,
          B_scenarios   = B_scenarios,
          A             = A,
          Psi           = Psi,
          rho_extra     = rho_extra,
          models_to_run = models_to_run,
          base_seed     = base_seed
        )
      }
    )
    return(dplyr::bind_rows(results_list))
  }

  # make the cluster
  cl <- parallel::makeCluster(cores)

  # set RNG streams on the cluster for reproducibility
  parallel::clusterSetRNGStream(cl, iseed = base_seed)

  # load required packages on each worker
  parallel::clusterEvalQ(cl, {
    library(lavaan)
    library(mvtnorm)
    NULL
  })

  # export all necessary functions and variables to the cluster
  parallel::clusterExport(
    cl,
    c(
      "sample_B_linear",
      "generate_B_constant",
      "generate_B_stepwise",
      "simulate_panel_data",
      "build_clpm",
      "build_riclpm",
      "build_dpm",
      "build_clpm_with_Cs",
      "safe_fit_clpm",
      "safe_fit_riclpm",
      "safe_fit_dpm",
      "safe_fit_clpm_C",
      "safe_fit_clpm_resid",
      "extract_lagged_parameters",
      "extract_rho_vec",
      "residualise_panel_linearC",
      "run_one_rep_study",
      "N","T","k","R2_1","target_sd","scenarios","B_scenarios",
      "A","Psi","rho_extra","models_to_run","base_seed"
    ),
    envir = environment()
  )

  # run the simulation with a progress bar
  results_list <- pbapply::pblapply(
    X = 1:reps,
    cl = cl,
    FUN = function(rep_id) {
      run_one_rep_study(
        rep_id        = rep_id,
        N             = N,
        T             = T,
        k             = k,
        R2_1          = R2_1,
        target_sd     = target_sd,
        scenarios     = scenarios,
        B_scenarios   = B_scenarios,
        A             = A,
        Psi           = Psi,
        rho_extra     = rho_extra,
        models_to_run = models_to_run,
        base_seed     = base_seed
      )
    }
  )

  # stop the cluster
  parallel::stopCluster(cl)

  # return the results
  dplyr::bind_rows(results_list)
}
