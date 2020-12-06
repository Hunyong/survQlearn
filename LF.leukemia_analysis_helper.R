
weights_leuk <- function(data, lvls1) {
  # stage 1 was randomized
  propensity1 = sapply(lvls1, function(x) mean(data$A.1 == x))
  propensity1 = sapply(data$A.1, function(s) propensity1[lvls1 == s])
  prop2.model =
    glm(A.2 ~ age + cyto + prevTime.2 + prevA.2 + response.2, 
        data = data %>% filter(!is.na(A.2)), family = binomial)
  propensity2 = predict(prop2.model, type = "response")
  propensity2 = ifelse(is.na(data$A.2), 1, ifelse(data$A.2 == 1, propensity2, 1 - propensity2))
  propensity  = propensity1 * propensity2
  
  # DWSurv can only hanlde binary Tx.
  propensity1.DW = sapply(0:1, function(x) mean(data$binaryA1 == x))
  propensity1.DW = sapply(as.numeric(data$binaryA1), function(s) propensity1.DW[0:1 == s])
  propensity.DW = propensity1.DW * propensity2
  
  # IPCW
  Sc.hat1 <- coxph(Surv(T.1, 1 - d.1) ~ age, data = data) # Only 1 censored, so regress on only one var.
  Sc.hat1 <- exp( - predict(Sc.hat1, type = "expected"))
  Sc.hat2 <- coxph(Surv(T.2, 1 - d.2) ~ age, data = data) # Only 1 censored, so regress on only one var.
  Sc.hat2 <- exp( - predict(Sc.hat2, type = "expected"))
  Sc.hat <- Sc.hat1
  Sc.hat[!is.na(data$d.2)] = Sc.hat[!is.na(data$d.2)] * Sc.hat2 # only those who made to the second stage
  weight.censor = data$delta/Sc.hat
  
  return(list(propensity = propensity, propensity.DW = propensity.DW, weight.censor = weight.censor))
}
evaluate <- function(testY, weight, criterion, tau) {
  if (criterion[1] == "mean") {
    mean(pmin(tau, testY) * weight)/ mean(weight)
  } else {
    mean(as.numeric(testY >= as.numeric(criterion[2])) * weight)/ mean(weight)
  }
}
getValue <- function(test, actual, estimated, propensity, weight.censor, criterion, tau) {
  weight = apply(actual == estimated, 1, all, na.rm = TRUE) / propensity
  weight = weight * weight.censor
  evaluate (test[, "T.0"], weight = weight, criterion = criterion, tau = tau)
}
ET <- function(surv, time, tau) {
  tau.index = max(which(time <= tau))
  time.last = time[tau.index]
  surv = surv[1:tau.index]
  time = time[1:tau.index]
  surv = c(1, surv)
  time.diff = c(time, tau) - c(0, time)
  trunc.mean = sum(surv * time.diff)
  trunc.mean
}
St <- function(surv, time, t.eval, tau) {
  # tau is a placeholder for compatibility
  t.eval = as.numeric(t.eval)
  t.index.inf = max(which(time <= t.eval))
  t.index.sup = min(which(time >= t.eval))
  t.inf = time[t.index.inf]
  t.sup = time[t.index.sup]
  s.inf = surv[t.index.inf]
  s.sup = surv[t.index.sup]
  if (t.sup == t.inf) {
    S.t = s.inf
  } else {
    S.t = s.inf -  (s.inf - s.sup) * (t.eval - t.inf)/(t.sup - t.inf)
  }
  S.t
}