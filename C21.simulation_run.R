#### C21.simulation_body.R is to be run by source("C21.simulation_run.R")

### 0. Get the setting number.
arg <- commandArgs(trailingOnly = TRUE)
if (length(arg) < 4) {
  warning("commandArgs was not provided. Being set as c(1,1,1,1).")
  arg = c(1, 1, 1, 1) # by default
}
names(arg)[1:4] = c("beta", "propensity", "size", "crit")
arg1 <- as.numeric(arg[1]) # 1..4
arg2 <- as.numeric(arg[2]) # 1..2
arg3 <- as.numeric(arg[3]) # 1..2
arg4 <- as.numeric(arg[4]) # 1..3
arg.date <- if (is.na(arg[5]) | arg[5] == "") Sys.Date() else as.character(arg[5])


### 1. setup

# default settings
default <- list(n.eval = 10000, n.sim = 200, tau = 10, n.stages = 3)

# arg4 crit
crit <- list(crit1 = list(criterion = "mean", crit.value = NULL, value = "truncated mean E[T]"),
             crit2 = list(criterion = "surv.mean", crit.value = 5, value = "S(5)"))

# arg3 size
size <- list(small.sample.size = list(n = 300),
             large.sample.size = list(n = 1000))

# arg2 propensity
propensity <-   # (int), surv.prev, covariate (1~5)
  list(obs = list(beta.propensity = function(p) c(0, 1, -rep(0.5, p))),  # no unmeasured confounder
       rct  = list(beta.propensity = function(p) c(0, 0, rep(0, p))))    # RCT

# arg1 beta
betas <- list(
  beta1 =
    list(beta.hazard0 = c(0, -1, c(1,1,1, 1, 1)),              # (int), surv.prev, covariate (1~5)
         beta.hazard1 = c(0, -3, c(2,2,1,-1,-1)),              # (int), surv.prev, covariate (1~5)
         beta.censor0 = c(-3, 0.1, 0.2 * c(1,1,1, 1, 1)),      # (int), surv.prev, covariate (1~5)
         beta.censor1 = c(-3, 0.2, 0.2 * c(1,1,1,-1,-1))),     # (int), surv.prev, covariate (1~5)
  beta2 =  # moderate censoring rate
    list(beta.hazard0 = c(0, -1, c(1,1,1, 1, 1)),              # (int), surv.prev, covariate (1~5)
         beta.hazard1 = c(0, -3, c(2,2,1,-1,-1)),              # (int), surv.prev, covariate (1~5)
         beta.censor0 = c(-2, 0.1, 0.2 * c(1,1,1, 1, 1)),      # (int), surv.prev, covariate (1~5)
         beta.censor1 = c(-2, 0.2, 0.2 * c(1,1,1,-1,-1))),     # (int), surv.prev, covariate (1~5)
  beta3 =  # common.beta1 + smaller p
    list(beta.hazard0 = c(0, -1, 1, 1),                        # (int), surv.prev, covariate (1~2)
         beta.hazard1 = c(0, -2, 2, 0),                        # (int), surv.prev, covariate (1~2)
         beta.censor0 = c(-3, 0.1, 0.2, 0.2),                  # (int), surv.prev, covariate (1~2)
         beta.censor1 = c(-3, 0.2, 0.2, -0.2)),                 # (int), surv.prev, covariate (1~2)
  beta4 =  # common.beta1 + larger p
    list(beta.hazard0 = c(0, -1, rep(1,8),  rep(0,2)),    # (int), surv.prev, covariate (1~10)
         beta.hazard1 = c(0, -2, rep(2,4), rep(-1,4), rep(-2,2)),    # (int), surv.prev, covariate (1~10)
         beta.censor0 = c(-3, 0.1, 0.2 * c(rep(1,4), rep(1,3), rep(0,3))),   # (int), surv.prev, covariate (1~10)
         beta.censor1 = c(-3, 0.2, 0.2 * c(rep(1,4), rep(-1,3), rep(0,3))))   # (int), surv.prev, covariate (1~10)
)

p.list <- lapply(betas, function(x) length(x$beta.hazard0) - 2)

    # setting1 n.boot 50 / n 300
    setting = c(arg = list(arg), default, betas[[arg1]], p = p.list[[arg1]],
                propensity[[arg2]], size[[arg3]], crit[[arg4]])
    
### 2. put a selected setting into the global environment    
    cat("setting (beta, propensity, size, crit) ", arg, "\n")
    list2env(setting, envir = globalenv())
    beta.propensity <- beta.propensity(p)
    filename = paste0("output/simResult_", arg.date, "_beta", arg1, "_prop", arg2, 
                      "_n", arg3, "_crit", arg4, ".rds")
    
    # complete the scores
    predPropensityFn1 <- function(surv.previous, rho, omega, covariate) {
      # rho, omega are there for compatibility only.
      plogis(cbind(1, log(surv.previous + 1), covariate) %*% beta.propensity)
    }
    predHazardFn1 <- function(surv.previous, rho, omega, action, covariate) {
      # rho, omega are there for compatibility only.
      ifelse(action == 1,
             cbind(1, log(surv.previous + 1), covariate) %*% beta.hazard1,
             cbind(1, log(surv.previous + 1), covariate) %*% beta.hazard0)
    }
    predCensorFn1 <- function(surv.previous, rho, omega, action, covariate) {
      # rho, omega are there for compatibility only.
      ifelse(action == 1,
             cbind(1, log(surv.previous + 1), covariate) %*% beta.censor1,
             cbind(1, log(surv.previous + 1), covariate) %*% beta.censor0)
    }
    noCensorFn1 <- function(surv.previous, rho, omega, action, covariate) {
      matrix(-10, nrow = length(surv.previous), ncol = 1)
    }
    predCensorFn <- predHazardFn <- predPropensityFn <- noCensorFn <- vector("list", n.stages)
    for (i in 1:n.stages) {
      predPropensityFn[[i]] <- predPropensityFn1
      predHazardFn[[i]]     <- predHazardFn1
      predCensorFn[[i]]     <- predCensorFn1
      noCensorFn[[i]]       <- noCensorFn1
    }

### 3. Run the simulation
    skip.csk <- skip.gk <- skip.dw <- skip.zom <- FALSE
    cv.nodesize = FALSE
    
    source("C21.simulation_body.R")


### More description
  # n       # training sample size
  # n.eval  # evaluation sample size
  # n.sim   # number of simulation replicates
  # p = 5   # number of covariates (should agree with beta coefficients' dimension)
  # tau = 10 # total study length
  # n.stages = 3 # number of stages

