library(MASS); library(ggplot2); library(dplyr); library(Rcpp)

#' @return event.time observed time
#' @return gamma 1(T <= U); T = failure time, U = treatment time, T_tilde = min(T, U)
#' @return delta 1(T_tilde <= C); C = censoring time
#' @return rho tumor size
#' @return omega wellness
dynamics <- function(time.max = 3, tick = 0.01, rho.plus, omega.plus, 
                     pred.hazard = 0, pred.censor = 0, admin.censor = Inf,
                     at.risk = 1, corr = -0.5, multiplier = 2, full_data = FALSE, terminal.stage = FALSE,
                     omega.threshold = 0.1) {
  # at.risk  : a vector of who's eligible for this stage. (patients who already died or censored will not be at risk.)
  # full_data: present all hidden event times (failure, treatment, censoring) and trajectories
  #            if (full_data == 2) present all hidden event times without trajectories.
  # admin.censor: remaining study length. For stage 1, = tau, but for later stages, = tau - time survived.
  # corr: correlatino of OU process (rho, omega)
  # if rho is greater than the upper bound, the patient is in serious state and cannot make it to the next stage.
  
  ## outputs :
  # event.time = X = min(T, D), gamma = 1(D <= T), delta = 1(not censored)
  # rho.x / omega.x = tumor size / wellness measure at the event time.
  if (!at.risk) {
    output = c(event.time = NA, gamma = NA, delta = NA,  rho.x = NA, omega.x = NA)
    if (full_data == 1) {
      return(list(statistics = output, 
                  times = c(failure.time = NA, treatment.time = NA,
                            censor.time = NA),
                  rho = NA, omega = NA, surv = NA))
    } else if (full_data == 2) {
      return(c(output, failure.time = NA, treatment.time = NA,
               censor.time = NA))
    } else {
      return(output)
    }
  }
  
  timeframe = seq(0, time.max, by = tick)
  n.time = length(timeframe)
  
  ## error generation
  # brownian_motion <- mvrnorm(n.time, rep(0,2), matrix(c(1,corr, corr, 1), 2) / n.time)
  brownian_motion <- mvrnorm(n.time, c(0.1, 0.05) * tick, matrix(c(1, corr, corr, 1), 2) * tick^2)
  brownian_motion <- apply(brownian_motion, 2, cum_ou, reverting_const = 2 * tick) * multiplier
  
  # ggplot(data.frame(brownian_motion, col = 1:n.time), aes(X1, X2, col=col)) + geom_path()
  
  ## trajectory of rho and omega
  rho = pmax(0, brownian_motion[, 1] + 2 * rho.plus * timeframe + rho.plus)
  omega = pmax(0.1, omega.plus + (1 - omega.plus) * (1 - exp(- timeframe/2)) + brownian_motion[, 2])
# tmp.bm <<- brownian_motion
  ## event times
  # surv = exp(-cumsum(timeframe * tick *4 /3 * exp(-brownian_motion[, 2]/3 + pred.hazard) * rho.plus))
  hazard = (rho / omega * exp(pred.hazard) + (omega <= omega.threshold)) * tick / 5  
  surv = exp(-cumsum(hazard))
# print(surv)
# tmp.haz <<- hazard
  failure.time = which(surv <= runif(1))[1] * tick # NA if administratively censored
    if (is.na(failure.time)) failure.time <- Inf
  # critical <- which(omega <= omega.threshold)[1]
  # if (!is.na(critical)) if (critical * tick < failure.time) failure.time = critical * tick
  if (terminal.stage) {
    trt.time = Inf
  } else {
    trt.time = which(rho >= 1)[1] * tick           # NA if administratively censored
    if (is.na(trt.time)) trt.time <- Inf
  }
  
  ## summary statistics 
  X = min(failure.time, trt.time)
  gamma = failure.time <= trt.time
  x.index = which(timeframe == X)
  # censoring
  censor.time = min(admin.censor, rexp(1, exp(pred.censor)))
# tmp.cens <<- pred.censor
# tmp.cens.time <<- censor.time
  XX = min(X, censor.time)
  delta = (X <= censor.time)
  
  if (!delta) { #if censored (when delta = 0)
    output = c(event.time = XX, gamma = NA, delta = delta, rho.x = NA, omega.x = NA)
    if (full_data == 1) {
      return(list(statistics = output, 
                  times = c(failure.time = failure.time, treatment.time = trt.time,
                            censor.time = censor.time),
                  rho = rho, omega = omega, surv = surv))
    } else if (full_data == 2) {
      return(c(output, failure.time = failure.time, treatment.time = trt.time,
                censor.time = censor.time))
    } else {
      return(output)
    }
  }
  
  
  # cat(failure.time, " ", trt.time, " ",  gamma,"\n")  
  if (!length(x.index)) {
    rho.x = NA
    omega.x = NA
  } else {
    rho.x = rho[x.index]
    omega.x = omega[x.index]
  }
  output = c(event.time = X, gamma = gamma, delta = delta, rho.x = rho.x, omega.x = omega.x)
  
  if (full_data == 1) {
    return(list(statistics = output, 
                times = c(failure.time = failure.time, treatment.time = trt.time,
                          censor.time = censor.time),
                rho = rho, omega = omega, surv = surv))
  } else if (full_data == 2) {
    return(c(output, failure.time = failure.time, treatment.time = trt.time,
             censor.time = censor.time))
  } else {
    return(output)
  }
}
dynamics.vec <- Vectorize(dynamics, vectorize.args = c("rho.plus", "omega.plus", "pred.hazard", 
                                                       "pred.censor", "at.risk", "admin.censor"))

## OU process with parameter theta -> replaced to an external Rcpp function
  # dx_t = -theta x_t dt + sigma dW_t
  # cum_ou <- function(vec, reverting_const) { #reverting_const = theta dt
  #   output = vec
  #   for (k in  2:length(vec)) {
  #     output[k] = output[k] + output[k-1] * (1 - reverting_const)
  #     # x_k = r_k + x_{k-1} (1-theta dt)
  #   }
  #   output
  # }
  # brownian_motion <- apply(brownian_motion, 2, cumsum)

  #cum_ou in Rcpp
  cppFunction(
    "NumericVector cum_ou(NumericVector vec, double reverting_const) {
    const int n = vec.size();
    NumericVector y(n);
    y[0] = vec[0];
    for (int i=1; i < n; ++i) {
    y[i] = y[i] + y[i-1] * (1 - reverting_const);
    }
    return y;
    }")
