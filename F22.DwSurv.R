#### JASA 2019 paper: Estimating Optimal Dynamic Treatment Regimes With Survival Outcomes
#### Gabrielle Simoneau, Erica E. M. Moodie, Jagtar S. Nijjar, Robert W. Platt & the Scottish Early Rheumatoid Arthritis Inception Cohort Investigators
source("F00.generic.R")
library(DTRreg)

## transform the array data into a wide matrix.
dwTrans <- function(dat.array, n.stages = dim(dat.array)[3], p = attr(dat.array, "p")) {
  # overall censoring status
  delta <- dat.array[, "delta", ]
  delta[is.na(delta)] <- Inf
  delta <- apply(delta, 1, min)
  
  # outcomes
  Y <- dat.array[, "event.time", ]
  Y[is.na(Y)] <- 0
  colnames(Y) <- paste0("Y.", 1:n.stages)
  
  # treatment
  A <- dat.array[, "action", ]
  colnames(A) <- paste0("A.", 1:n.stages)
  # effect coding into ordinary coding
  lvls <- A %>% as.vector() %>% unique %>% sort
  if (all(lvls == c(-1, 1))) A[!is.na(A) & A==-1] <- 0
  # A[is.na(A)] <- 0
  
  # log(B + 1)
  # lB <- log(dat.array[, "surv.previous", ] + 1)
  lB <- dat.array[, "lB", ]
  colnames(lB) <- paste0("lB.", 1:n.stages)
  lB <- lB[, -1] # remove B.1 to avoid singulairty problem
    
  Z <- dat.array[, paste0("Z", 1:p), ] 
    # putting NA for non at.risk subjects # This acts as a flag later on. (Some subjects still have Z values even when they are not available.)
    na.index = dat.array[,"at.risk",] == 0
    Z[, "Z1", ][na.index] <- NA
  
  Z <- matrix(Z, nrow = dim(dat.array)[1])
  colnames(Z) <- paste0("Z", 1:p, ".", rep(1:n.stages, each = p))
  
  subj.id <- dat.array[, "subj.id", 1]
  
  data.frame(subj.id, Y, delta, A, lB, Z)
}

dwFormula <- function(dv = "Y", iv = c("Z1", "Z2"), global.dv = FALSE, n.stages) {
  if (dv == "") global.dv = TRUE
  lapply(1:n.stages, function(s) {
    LHS = paste0(dv, ifelse(global.dv, "", paste0(".", s)), " ~ ")
    RHS = paste0(iv, ".", s) %>% paste(collapse = " + ")
    RHS <- gsub("lB\\.1 \\+", "", RHS)  # remove B.1 to avoid singularity
    as.formula(paste0(LHS, RHS))
    })
}

# Goldberg Kosorok with linear regression Pi = {H_i \beta} at stage k
dwSurv <- function(data, n.stages = dim(data)[3], tau = attr(data, "tau"), p = attr(data, "p"), 
                   formula = ...) {
  # transform original data into the auxiliary form
  dat <- dwTrans(data)
  iv <- c("lB", paste0("Z", 1:p))
  
  form.time <- dwFormula(dv = "", iv = "Y", n.stages = n.stages)
  form.blip <- dwFormula(dv = "", iv = iv, n.stages = n.stages)
  form.tf   <- dwFormula(dv = "", iv = iv, n.stages = n.stages)
  form.treat<- dwFormula(dv = "A", iv = iv, n.stages = n.stages)
  form.cens <- dwFormula(dv = "delta", iv = iv, n.stages = n.stages, global.dv = TRUE)

# DTRreg
tmp.dw <<- list(time = form.time, blip.mod = form.blip, tf.mod = form.tf,
                treat.mod = form.treat, cens.mod = form.cens, data = dat)
# dw.est <-
#   DWSurv2(time = form.time, blip.mod = form.blip, tf.mod = form.tf,
#          treat.mod = form.treat, cens.mod = form.cens, data = dat)
  
  dw.est <-
    DWSurv(time = form.time, blip.mod = form.blip, tf.mod = form.tf,
           treat.mod = form.treat, cens.mod = form.cens, data = dat)
  
  attr(dw.est, "class") <- "dwSurv"
  ITR <- lapply(1:n.stages, function(stage)
    decision.rule.dw(dw.est, stage = stage, return.optimal.Q = TRUE))

  list(dw.est = dw.est, action.est = ITR)
  
}





decision.rule.dw <- function(dw.est, newdata = NULL, stage = 1,
                             return.optimal.Q = FALSE, return.val.only = FALSE) {
  # dw.est is estimate for all stages.
  # newdata is the transformed data (dwTrans())
  
  # if (!is.null(newdata) & is.null(newdata$subj.id)) stop("subj.id is needed for new data.")
  if (is.null(newdata)) newdata = dw.est$data
  
  subj.id <- if (!is.null(newdata$subj.id)) {
    newdata$subj.id
  } else {
    rownames(newdata)
  }
  
  # subsetting the newdata by 
  full.tf   <- complete.cases(get_all_vars(dw.est$tf.mod[[stage]], newdata))
  full.blip <- complete.cases(get_all_vars(dw.est$blip.mod[[stage]], newdata))
  full <- full.tf & full.blip
tmp <<- list(newdata = newdata, full = full, newdata2 = newdata[full, ])
  newdata <- newdata[full, ]
  n.sample <- dim(newdata)[1]
  
  beta <- dw.est$beta[[stage]] # tf model coefficient
  psi  <- dw.est$psi[[stage]]  # blip model coefficient
  
  tf   <- model.matrix(dw.est$tf.mod[[stage]], newdata) %*% beta
  blip <- model.matrix(dw.est$blip.mod[[stage]], newdata) %*% psi

  pred0 <- pred1 <- rule <- rep(0, n.sample)
  pred0[full] <- tf
  pred1[full] <- tf + blip
  pred0[is.na(pred0)] <- 0
  pred1[is.na(pred1)] <- 0
  blip[is.na(blip)] <- 0
  
  attr(pred1, "criterion") <- attr(pred0, "criterion") <- list(criterion = "mean", crit.value = NULL)
  
  if (return.val.only) {
    return(list(val1 = pred1, val0 = pred0))
  }
  
  rule[full] = sign(blip)
  rule.fix = ifelse(rule == 0, rbinom(n.sample, 1, 0.5) * 2 -1, rule)
  rule = data.frame(subj.id = subj.id, rule = rule, rule.fix = rule.fix)
  
  if (return.optimal.Q) {
    pred.star <- pmax(pred1, pred0)
    names(pred.star) <- subj.id
    output <- list(rule = rule, 
                   Q.optim = pred.star)
    return(output)
  }
  return(rule)
}
# decision.rule.dw(dw.tmp)
# decision.rule.dw(dw.tmp, stage=2)
