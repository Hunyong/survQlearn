
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

weights_aric <- function(data, weight.formula.list, weight.formula.bin.list) {
  require(dplyr) # %>%
  require(randomForest)
  require(randomForestSRC)
  
  n = dim(data)[1]
  Tx.nm = sapply(weight.formula.list, function(x) x[2] %>% 
                   as.character %>% gsub("(factor\\()(.*)(\\))", "\\2", .) )
  prop1.model =
    randomForest(weight.formula.list[[1]], data = data, )
  propensity1 = predict(prop1.model, data, type = "prob")                   # n x n.levels matrix
  if (is.null(dim(propensity1))) {
    # in case it is binary
    propensity1 = data.frame(1 - propensity1, propensity1)
    names(propensity1) = 0:1
  }
  missing1 = is.na(propensity1[, 1]) & !is.na(data$d.1)             # Bookkeeping covariate missing subjects (neither censored nor failed)
  propensity1 = propensity1[cbind(rownames(data), data[, Tx.nm[1]] %>% as.character)]  # A vector of predicted pi(A|X) with NA's
  propensity1 = ifelse(is.na(propensity1), 1, propensity1)             # Replace NA's with 1's.
  
  avail2 = !is.na(data$d.2)
  prop2.model =
    randomForest(weight.formula.list[[2]], 
                 data = data[avail2, ] %>% mutate(CENTER = CENTER %>% droplevels))
  Tx2.level = nlevels(data[, Tx.nm[2]])
  if (Tx2.level <= 2) {Tx2.level = 2}
  propensity2 = matrix(NA, n, Tx2.level, dimnames = list(rownames(data), levels(data[, Tx.nm[2]]) %>% as.character))
  if (Tx2.level > 2) {
    propensity2[avail2, ] = predict(prop2.model, data[avail2, ], type = "prob")                   # n x n.levels matrix
  } else {
    colnames(propensity2) = 0:1
    propensity2[avail2, "1"] = predict(prop2.model, data[avail2, ], type = "prob")[, "1"]                   # n x n.levels matrix
    propensity2[avail2, "0"] = 1 - propensity2[avail2, "1"]
  }
  propensity2 = propensity2[cbind(rownames(data), data[, Tx.nm[2]] %>% as.character)]  # A vector of predicted pi(A|X) with NA's
  propensity2 = ifelse(is.na(propensity2), 1, propensity2)             # Replace NA's with 1's.
  
  propensity  = propensity1 * propensity2
  
  prop1.DW =
    randomForest (weight.formula.bin.list[[1]], data = data)
  Tx.var = 
    weight.formula.bin.list[[1]][2]  %>% 
    as.character %>% gsub("(factor\\()(.*)(\\))", "\\2", .) 
  propensity1.DW = predict(prop1.DW, data, type = "prob")[, "1"]
  propensity1.DW = ifelse(is.na(data[, Tx.var]), 1, 
                          ifelse(data[, Tx.var] == 1, propensity1.DW, 1 - propensity1.DW))
  prop2.DW =
    randomForest (weight.formula.bin.list[[2]], 
                  data = data[avail2, ] %>% mutate(CENTER = CENTER %>% droplevels))
  Tx.var2 = 
    weight.formula.bin.list[[2]][2] %>% 
    as.character %>% gsub("(factor\\()(.*)(\\))", "\\2", .) 
  propensity2.DW = matrix(NA, n, 1)
  propensity2.DW[avail2] = predict(prop2.DW, data[avail2, ], type = "prob")[, "1"]
  propensity2.DW = ifelse(is.na(data[, Tx.var2]), 1, 
                          ifelse(data[, Tx.var2] == 1, propensity2.DW, 1 - propensity2.DW))
  
  propensity.DW  = propensity1.DW * propensity2.DW
  
  # IPCW
  Sc.hat1 <- rfsrc(Surv(V.1, 1 - d.1) ~ 
                     PREVHF01 + PRVCHD05 + hf + CENTER + gender + race + age + 
                     cig_yrs + AA.1 + AC.1 + AS.1 + 
                     bmi.1 + wth.1 + drink.1 + hypert.1 + glucose.1 + smoke.1 + hdl.1, 
                   data = data)
  Sc.hat1 <- sapply(1:dim(data)[1], 
                    function(i) St2(Sc.hat1$survival[i, ], Sc.hat1$time.interest, 
                                    t.eval = data$V.1[i]))
  
  Sc.hat2.res <- rfsrc(Surv(V.2, 1 - d.2) ~ 
                         PREVHF01 + PRVCHD05 + hf + CENTER + gender + race + age + 
                         cig_yrs + AA.1 + AC.1 + AS.1 + V.1 + 
                         AA.2 + AC.2 + AS.2 + V.2 + 
                         bmi.2 + wth.2 + drink.2 + hypert.2 + glucose.2 + smoke.2 + hdl.2, 
                       data = data[avail2, ] %>% mutate(CENTER = CENTER %>% droplevels))
  Sc.hat2 <- rep(NA, n)
  Sc.hat2[avail2] <- 
    sapply(1:sum(avail2), 
           function(i) St2(Sc.hat2.res$survival[i, ], Sc.hat2.res$time.interest,
                           t.eval = data$V.2[avail2][i]))
  
  Sc.hat <- Sc.hat1
  Sc.hat[!is.na(data$d.2)] = Sc.hat[!is.na(data$d.2)] * Sc.hat2[!is.na(data$d.2)] # only those who made to the second stage
  # clipping
  Sc.hat <- Sc.hat %>% pmax(0.05) %>% pmin(0.95)
  weight.censor = data$delta/Sc.hat
  
  # p1 <<- propensity1
  # pd1 <<- propensity1.DW
  # c1 <<- Sc.hat1
  
  return(list(propensity = propensity, propensity.DW = propensity.DW, weight.censor = weight.censor))
}

all2 = function(x) {
  # if everything is NA, return NA. Otherwise, TRUE only if all is TRUE.
  # The naive all() returns NA even if there are only TRUEs except NAs.
  na.index = is.na(x)
  if (all(na.index)) return(NA)
  all(x[!na.index])
}
evaluate <- function(testY, weight, criterion, tau) {
  if (criterion[1] == "mean") {
    mean(pmin(tau, testY) * weight)/ mean(weight)
  } else {
    mean(as.numeric(testY >= as.numeric(criterion[2])) * weight)/ mean(weight)
  }
}
getValue <- function(test, actual, estimated, propensity, weight.censor, criterion, tau) {
  # weight = apply(actual == estimated, 1, all, na.rm = TRUE) / propensity
  weight = 
    apply(actual == estimated, 1, all2) %>%   
    # When there is at least one NA,
    # 1. all() does not return TRUE                     all(c(NA, NA, T)) = NA; all(c(T, T, T)) = TRUE
    # 2. all() returns NA if there is no FALSE          all(c(NA, NA, F, T)) = FALSE;  all(c(NA, NA)) = NA
    # 3. all() returns FALSE if there is at least one FALSE
    # When the second stage is not available in test set, the NA-match cases should still be counted. NA => 1.
    {ifelse(is.na(.), 1, as.numeric(.))} %>% 
    "/" (propensity)
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

St2 <- function(surv, time, t.eval, tau = Inf, exponential.tail = TRUE) {
  # tau is a placeholder for compatibility
  time = c(0, time)
  surv =  c(1, surv)
  max.t = max(time)
  min.s = min(surv)
  if (t.eval > max.t) {
    if (exponential.tail) {
      S.t = 1 - (1-min.s)^(t.eval / max.t)
    } else {
      S.t = min.s
    }
  } else {
    S.t = approxfun(time, surv)(t.eval)
  }
  return(S.t)
}

