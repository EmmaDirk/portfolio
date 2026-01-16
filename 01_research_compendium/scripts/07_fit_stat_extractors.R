# since all fits are slightly different, we need to figure out where in the model object
# our parameters of interest: autoregressive, cross-lagged, and residual correlation (rho) are located
# this file contains functions to extract those parameters from the fitted model objects
# ------------------------------------------------------------------------------------------------

# lagged parameters extractor
extract_lagged_parameters <- function(
    fit,                                                      # lavaan model object
    T,                                                        # number of time points
    model_type = c("clpm", "riclpm", "dpm")                   # model type
){

  # match model type
  model_type <- match.arg(model_type)

  # if the model fit failed, return NAs
  if (is.null(fit)) {
    return(list(
      ar_x = rep(NA, T-1),    # autoregressive X_t ← X_{t-1}
      ar_y = rep(NA, T-1),    # autoregressive Y_t ← Y_{t-1}
      xy   = rep(NA, T-1),    # cross-lag Y_t ← X_{t-1}   (X → Y)
      yx   = rep(NA, T-1)     # cross-lag X_t ← Y_{t-1}   (Y → X)
    ))
  }

  # try to extract parameter table
  pe <- tryCatch(lavaan::parameterEstimates(fit), error=function(e) NULL)

  # if extraction failed, return NAs
  if (is.null(pe)) {
    return(list(
      ar_x = rep(NA, T-1),
      ar_y = rep(NA, T-1),
      xy   = rep(NA, T-1),
      yx   = rep(NA, T-1)
    ))
  }

  # RI-CLPM uses latent within-person variables wx, wy
  if (model_type == "riclpm") {
    xvar <- "wx"
    yvar <- "wy"
  } else {
    xvar <- "x"
    yvar <- "y"
  }

  # helper to grab a single parameter
  grab <- function(lhs, rhs) {
    ix <- which(pe$lhs == lhs & pe$rhs == rhs)
    if (length(ix) == 0) return(NA_real_)
    pe$est[ix[1]]
  }

  # containers
  ar_x <- numeric(T-1)
  ar_y <- numeric(T-1)
  xy   <- numeric(T-1)
  yx   <- numeric(T-1)

  # extract all lagged parameters
  for (t in 2:T) {

    # autoregressive
    ar_x[t-1] <- grab(paste0(xvar, t), paste0(xvar, t-1))
    ar_y[t-1] <- grab(paste0(yvar, t), paste0(yvar, t-1))

    # cross-lag (X → Y)
    xy[t-1]   <- grab(paste0(yvar, t), paste0(xvar, t-1))

    # cross-lag (Y → X)
    yx[t-1]   <- grab(paste0(xvar, t), paste0(yvar, t-1))
  }

  list(
    ar_x = ar_x,
    ar_y = ar_y,
    xy   = xy,
    yx   = yx
  )
}

# residual correclations extractor
extract_rho_vec <- function(
    fit,                                                           # lavaan model object
    T,                                                             # number of time points
    model_type = c("clpm","riclpm","dpm")                          # model type
){

  # match model type
  model_type <- match.arg(model_type)

  # if the model fit failed, return NAs
  if (is.null(fit)) return(rep(NA_real_, T))

  # try to extract parameter estimates
  pe <- tryCatch(lavaan::parameterEstimates(fit), error=function(e) NULL)

  # if extraction failed, return NAs
  if (is.null(pe)) return(rep(NA_real_, T))

  # determine variable names based on model type, default is x, otherwise is wx
  if (model_type == "riclpm") {
    xvar <- "wx"
    yvar <- "wy"
  } else {
    xvar <- "x"
    yvar <- "y"
  }

  # prepare container
  rho <- numeric(T)

  # extract rho
  for (t in 1:T) {

    # the left hand side of the correlation equation
    lhs_xy <- paste0(xvar, t)

    # the right hand side of the correlation equation
    lhs_yx <- paste0(yvar, t)

    # find the covariance estimate between x_t and y_t
    ix <- which(pe$lhs == lhs_xy & pe$rhs == lhs_yx)

    # if not found, try the other direction
    if (length(ix) == 0) {

      # find the covariance estimate between y_t and x_t
      ix <- which(pe$lhs == lhs_yx & pe$rhs == lhs_xy)
    }

    # if not found, return NA
    if (length(ix) == 0) {
      rho[t] <- NA_real_
      next
    }

    # covariance estimate
    cov_xy <- pe$est[ix[1]]

    # find the variance estimates for x_t and y_t
    vx_idx <- which(pe$lhs == lhs_xy & pe$rhs == lhs_xy)
    vy_idx <- which(pe$lhs == lhs_yx & pe$rhs == lhs_yx)

    # if not found, return NA
    if (length(vx_idx) == 0 || length(vy_idx) == 0) {
      rho[t] <- NA_real_
      next
    }

    # variance estimates
    vx <- pe$est[vx_idx[1]]
    vy <- pe$est[vy_idx[1]]

    # compute rho: cov_xy / sqrt(vx * vy)
    if (is.na(vx) || is.na(vy) || vx <= 0 || vy <= 0) {
      rho[t] <- NA_real_
    } else {
      rho[t] <- cov_xy / sqrt(vx * vy)
    }
  }

  # return the residual correlations vector
  rho
}
