# this script contains a function to sample beta coefficients for confounders
# such that the total R^2 of confounders at time t = 1 is equal to a specified value
# ---------------------------------------------------------------------

sample_B_linear <- function(
    k,                                                             # number of confounders
    R2_1,                                                          # total confounder R^2 at t = 1 
    min_abs   = 0.01,                                              # minimum absolute value for each beta
    max_abs   = 0.60,                                              # maximum absolute value for each beta
    max_tries = 100000                                             # maximum sampling attempts
) {

  # split R2_1 equally across X and Y
  target_X <- R2_1 
  target_Y <- R2_1 

  # start the loop: for i in the max number of tries
  for (i in seq_len(max_tries)) {

    # sample k random numbers from a normal distribution
    u_x <- rnorm(k)

    # normalize so the sum of squares of this vector is 1
    u_x <- u_x / sqrt(sum(u_x^2))

    # scale to the target variance 
    b_x <- sqrt(target_X) * u_x

    # repeat for Y
    u_y <- rnorm(k)
    u_y <- u_y / sqrt(sum(u_y^2))
    b_y <- sqrt(target_Y) * u_y

    # check if all absolute values are within the specified bounds
    if (all(abs(b_x) >= min_abs,
            abs(b_x) <= max_abs,
            abs(b_y) >= min_abs,
            abs(b_y) <= max_abs)) {

      # if so, return the B matrix
      B1 <- rbind(b_x, b_y)

      # set row and column names
      rownames(B1) <- c("X", "Y")
      colnames(B1) <- paste0("c", 1:k)
      return(B1)
    }
  }

  # if the loop doesn't return a valid B matrix, throw an error
  stop("Failed to sample a valid B matrix within max_tries.")
}
