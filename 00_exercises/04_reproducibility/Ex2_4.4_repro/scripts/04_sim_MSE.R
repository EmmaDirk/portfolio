# this script provides the simulation study for measurement error
# ---------------------------------------------------------------------

# simulation of measurement error:
# store the reference estimate
ref <- lm(bp ~ HbA1C + bmi + age + as.factor(sex), data=dc)$coef[2]

# number of simulations
n.sim <- 1e3

# define percentages of measurement error for exposure and confounder (i.e. measurement error levels)
perc.me.exp <- seq(0,.5,.1)
perc.me.conf<- seq(0,.5,.1)

# for all combinations of measurement error levels across the two variables
scenarios <- expand.grid(perc.me.exp,perc.me.conf)

# calculate variances of exposure and confounder for the original data
var.exp <- var(dc$HbA1C)
var.conf <- var(dc$bmi)

# sample size
n <- dim(dc)[1]

# matrix for runs x scenarios (1000 x 36)
beta.hat <- matrix(ncol=dim(scenarios)[1], nrow=n.sim)

# run simulations
for (k in 1:n.sim){

  # print progress: turned off now to reduce output clutter
  # print(k)

  # set seed for reproducibility, but different across runs
  set.seed(k)

  # for each scenario
  for (i in 1:dim(scenarios)[1]){

    # calculate measurement error variances for exposure and confounder
    var.me.exp <- var.exp*scenarios[i,1]/(1-scenarios[i,1])
    var.me.conf <- var.conf*scenarios[i,2]/(1-scenarios[i,2])

    # create measured variables by adding random normal error to true variables
    dc$HbA1C.me <- dc$HbA1C + rnorm(dim(dc)[1], 0, sqrt(var.me.exp) )
    dc$bmi.me <- dc$bmi + rnorm(dim(dc)[1], 0, sqrt(var.me.conf) )

    # estimate beta and store in the matrix
    beta.hat[k,i] <- lm(bp ~ HbA1C.me + age + bmi.me + as.factor(sex), data=dc)$coef[2]
  }}