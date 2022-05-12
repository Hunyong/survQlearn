library(caret); library(purrr); library(dplyr); library(survival)
library(DTRreg); library(randomForestSRC); library(dtrSurv)
source("LF.helper.R")  # All the library and source files
source("F00.generic.R")
source("F21.GK1.R")
source("F21.GK2.R")
source("F22.DwSurv.R")


### 0. parameters
value.criterion = "mean"

if (value.criterion[1] != "mean") value.criterion[2] = 180

## time points
tau = 450
timepoints = seq(0, sqrt(tau), length.out = 100)^2

## tuning parameters
Ntree = 100
nodesize = 10
mindeath = round(sqrt(c(nodesize)), 0)
rule = if (value.criterion[1] == "mean") "mean" else "logrank"
ert = TRUE; rs = 0.2 # randomSplit = 0.2

cat("value.criterion =",value.criterion[1],
    "value.s =", value.criterion[2], "ntree =", Ntree, 
    "tau =", tau, "nodesize =", nodesize, "mindeath =", mindeath, "\n")


## models
form1 <- Surv(T, delta) ~ A + age + cyto + prevTime + prevA + response



### 1. data preprocessing
leukemia <- sas7bdat::read.sas7bdat("../data_Leukemia/leukemia_xh.sas7bdat")
dat.leuk <-
  leukemia %>% 
  transmute(T.0 = T0,
            T.1 = ifelse(resp == "FAIL", T0, T1),
            T.1 = ifelse(is.na(T.1), T.0, T.1),
            T.2 = ifelse(is.na(T22), 0, T22),
            delta = status1, 
            d.1 = ifelse(T.2 == 0 & delta == 0, 0, 1),
            d.2 = ifelse(d.1 == 0 | T.2 == 0, NA, delta),
            A.1 = frontline_treatment,
            A.2 = ifelse(HDAC == -1, NA, HDAC),
            binaryA1 = ifelse(A.1 %in% c("FAI+ATRA", "FAI+G+ATRA"), 1, 0), # For dwSurv only
            A.1 = as.factor(A.1),
            prevTime.2 = T.1,
            prevA.2 = A.1,
            response.2 = as.factor(resp),
            age = age,
            cyto = as.factor(cyto3))
lvls1 = levels(dat.leuk$A.1)
if (is.null(lvls1)) lvls1 = 0:1



## cross-validation
K = 100
set.seed(100, kind="Mersenne-Twister")
cv.insample  = createDataPartition(dat.leuk$A.1, p = .8, list = FALSE, times = K)
cv.outsample = map_dfc(1:K, ~which(!(1:210 %in% cv.insample[, .]))) %>% as.matrix

# skeleton
values <- matrix(NA, K, 11, 
                 dimnames = list(1:K, c("CSK", "GKRF", "GKLM", "DW", "ZOM", "observed", "cens1", "cens2",
                                        "ns.CSK", "ns.GK", "ns.GKLM")))

# common parameters
stg = 1:2  # stages

# common file name
nm = paste0(value.criterion[1], "_", 
            round(as.numeric(value.criterion[2]), 1), "_", 
            rule, "Split_", "tau_", tau)
fnm = paste0("output_leuk/", Sys.Date(), "/") # folder name
if (!dir.exists("output_leuk")) dir.create("output_leuk")
if (!dir.exists(fnm)) dir.create(fnm)
rds = paste0("_", Sys.Date(), ".rds")

# weights
weight = weights_leuk (data = dat.leuk, lvls1 = lvls1)


for (cv in 1:K) {
  cat(cv, "th cv.\n")
  set.seed(cv)
  in.cv = cv.insample[, cv]   # insample index
  out.cv = cv.outsample[, cv] # outsample index
  train = dat.leuk[in.cv, ]
  test  = dat.leuk[out.cv,]
  values[cv, "cens1"] = mean(train$d.1==0); values[cv, "cens2"] = mean(train$d.2==0, na.rm = T)
  values[cv,]
  
  
  ############################################################################################################
  #### A. rule estimation ####################################################################################
  ############################################################################################################
  
  ### A1. The proposed method
  
  args.CSK <- list(data = train, 
                  txName = paste("A", 1:2, sep = "."),
                  models = list(Surv(T.1, d.1) ~ A.1 + age + cyto,
                                Surv(T.2, d.2) ~ A.2 + age + cyto + prevTime.2 + prevA.2 + response.2),
                  usePrevTime = FALSE, tau = tau, timePoints = timepoints,
                  criticalValue = value.criterion[1], evalTime = as.numeric(value.criterion[2]), 
                  splitRule = ifelse(value.criterion[1] == "mean", "mean", "logrank"),
                  ERT = ert, uniformSplit = ert, replace = !ert,
                  randomSplit = rs, nTree = Ntree, mTry = c(3, 6),
                  pooled = FALSE, stratifiedSplit = FALSE)
  
  # actual fitting
  values[cv, "ns.CSK"] = nodesize
  set.seed(cv)
  CSK.leuk.i <- 
    try(do.call(dtrSurv, c(args.CSK, list(nodeSize = nodesize, minEvent = mindeath))))
  err.CSK = class(CSK.leuk.i)[1] == "try-error"
  
  
  ### A2. Goldberg-Kosorok RF
  args.GK = list(common.formula = form1,
                 common.Tx.label = "A", stage.label = 1:2, tau = tau,
                 data = train, stage.sep = ".", regress.prev.time = F)
  
  # actual fitting
  values[cv, "ns.GK"] = nodesize
  set.seed(cv)
  GKRF.leuk.i <-
    try(do.call(gk.separate, c(args.GK, list(nodesize = nodesize, method = "rf"))))
  err.GKRF = class(GKRF.leuk.i)[1] == "try-error"
  
  ### A3. Goldberg-Kosorok linear
  set.seed(cv)
  GKLM.leuk.i <-
    try(gk.separate(common.formula = form1,
                    common.Tx.label = "A", stage.label = 1:2, tau = tau,
                    data = train, stage.sep = ".", method = "lm", regress.prev.time = F))
  err.GKLM = class(GKLM.leuk.i)[1] == "try-error"
  
  
  set.seed(cv)
  DW.leuk.i <-
    try(DWSurv(time = list(~ T.1, ~ T.2), 
               blip.mod = list(~ age + cyto, ~ age + cyto),
               treat.mod = list(binaryA1 ~ 1, A.2 ~ age + cyto), # stage 1 is RCT
               tf.mod = list(~ age + cyto, ~ age + cyto),
               cens.mod = list(delta ~ age + cyto, delta ~ age + cyto),
               data = train))
  err.DW = class(DW.leuk.i)[1] == "try-error"
  
  ### A5. zero-order model
  zom.leuk.i <- 
    try(do.call(dtrSurv, c(args.CSK, list(nodeSize = 1e+4, minEvent = 1e+4))))
  err.zom = class(zom.leuk.i)[1] == "try-error"
  
  if (!err.zom) {
    zom.pred = 
      data.frame(
        zom.leuk.i@stageResults[[1]]@optimal@optimalTx[1],
        zom.leuk.i@stageResults[[2]]@optimal@optimalTx[1])
  }
  
  ############################################################################################################
  #### B. Rules applied to a test set ########################################################################
  ############################################################################################################
  
  ### B0. skeletons for the predicted A's
  opt.rule.pred = matrix(NA, dim(cv.outsample)[1], 2, dimnames = list(NULL, stg))
  opt.rule.CSK <- opt.rule.GKRF <- opt.rule.GKLM <- opt.rule.DW <-  opt.rule.zom <- opt.rule.pred
  
  
  for (q in seq_along(stg)) {
    elig = !test[, paste0("A.", q)] %>% is.na
    
    ### B1. the proposed method
    
    if (!err.CSK) {
      opt.CSK =
        predict(CSK.leuk.i, 
                # newdata = test[elig, ],
                newdata = test[elig, ] %>% mutate(A.2 = ifelse(is.na(A.2), -1, A.2), # temporary code for alpha version.
                                                  d.2 = ifelse(is.na(d.2), 0, d.2)), 
                stage = q)
      if (q==1) {
        opt.rule.CSK[elig, q] = factor(opt.CSK$optimal@optimalTx, levels = lvls1) %>% as.numeric()
      } else {
        opt.rule.CSK[elig, q] = opt.CSK$optimal@optimalTx
      }
    }
    ### B2. Goldberg-Kosorok RF
    if (!err.GKRF) {
      opt.GKRF = predict.opt.sep(GKRF.leuk.i$survRF[[q]],
                                 newdata = test[elig, ] %>% dplyr::select(-T.2, -d.2) %>% mutate(A.2 = ifelse(is.na(A.2), -1, A.2)),
                                 Tx.label = paste0("A.", q))
      opt.rule.GKRF[elig, q] = opt.GKRF$optimal.Tx
    }
    ### B3. Goldberg-Kosorok linear
    if (!err.GKLM) {
      opt.GKLM = predict.opt.sep(GKLM.leuk.i$survRF[[q]],
                                 newdata = test[elig, ] %>% dplyr::select(-T.2, -d.2) %>% mutate(A.2 = ifelse(is.na(A.2), -1, A.2)),
                                 Tx.label = paste0("A.", q))
      opt.rule.GKLM[elig, q] = opt.GKLM$optimal.Tx
    }
    
    ### B4. Simoneau et al.
    # DW (effect coded -> dummy)
    if (!err.DW) {
      opt.DW = decision.rule.dw(DW.leuk.i, newdata = test[elig, ] %>% dplyr::select(-T.2, -d.2) %>% mutate(A.2 = ifelse(is.na(A.2), -1, A.2)), stage = q)$rule.fix
      opt.rule.DW[elig, q] = opt.DW*.5 + .5
    }
    
    ### B5. zero-order model
    if (!err.zom) {
      if (q==1) {
        opt.rule.zom[elig, q] = factor(zom.pred[1, q], levels = lvls1) %>% as.numeric()
      } else {
        opt.rule.zom[elig, q] = zom.pred[1, q]
      }
    }
  }
  
  
  ############################################################################################################
  #### C. Value estimation ###################################################################################
  ############################################################################################################
  
  
  test.tmp = test %>% transmute(A.1 = as.numeric(A.1), A.2 = ifelse(is.na(A.2), -1, A.2))
  test.tmp.DW = test %>% transmute(A.1 = as.numeric(as.numeric(binaryA1)), A.2 = ifelse(is.na(A.2), -1, A.2))
  if (!err.CSK)  opt.tmp.CSK  = opt.rule.CSK %>% as.data.frame %>% transmute(A.1 = `1`, A.2 = ifelse(is.na(`2`), -1, `2`))
  if (!err.GKRF) opt.tmp.GKRF = opt.rule.GKRF %>% as.data.frame %>% transmute(A.1 = `1`, A.2 = ifelse(is.na(`2`), -1, `2`))
  if (!err.GKLM) opt.tmp.GKLM = opt.rule.GKLM %>% as.data.frame %>% transmute(A.1 = `1`, A.2 = ifelse(is.na(`2`), -1, `2`))
  if (!err.DW)   opt.tmp.DW   = opt.rule.DW %>% as.data.frame %>% transmute(A.1 = `1`, A.2 = ifelse(is.na(`2`), -1, `2`))
  if (!err.zom)  opt.tmp.zom  = opt.rule.zom %>% as.data.frame %>% transmute(A.1 = `1`, A.2 = ifelse(is.na(`2`), -1, `2`))
  
  arg.val = list(test = test, actual = test.tmp, propensity = weight$propensity[out.cv],
                 weight.censor = weight$weight.censor[out.cv], criterion = value.criterion,
                 tau = tau)
  if (!err.CSK)  values[cv, "CSK"] = do.call(getValue, c(arg.val, list(estimated = opt.tmp.CSK)))
  if (!err.GKRF) values[cv, "GKRF"] = do.call(getValue, c(arg.val, list(estimated = opt.tmp.GKRF)))
  if (!err.GKLM) values[cv, "GKLM"] = do.call(getValue, c(arg.val, list(estimated = opt.tmp.GKLM)))
                 values[cv, "observed"] = do.call(getValue, c(arg.val, list(estimated = test.tmp)))
  if (!err.zom)  values[cv, "ZOM"] = do.call(getValue, c(arg.val, list(estimated = opt.tmp.zom)))
  
  arg.val.DW = list(test = test, actual = test.tmp.DW, 
                    propensity = weight$propensity.DW[out.cv],
                    weight.censor = weight$weight.censor[out.cv], 
                    criterion = value.criterion,
                    tau = tau, estimated = opt.tmp.DW)
  if (!err.DW) values[cv, "DW"] = do.call(getValue, arg.val.DW)
  
  print(values[cv, ])
  print(apply(values, 2, mean, na.rm = TRUE))
  
  ############################################################################################################
  #### C. saving the results #################################################################################
  ############################################################################################################
  
  
  attr(values, "spec") = data.frame(criterion = value.criterion[1], criterion.s = as.numeric(value.criterion[2]), 
                                    rule = rule, tau = tau,
                                    ert = ert, rs = rs, ntree = Ntree,
                                    nodesize = nodesize, mindeath = mindeath)
  if (cv==1) saveRDS(CSK.leuk.i, paste0(fnm, "dtr.csk.leuk_", nm, "_", cv, rds))
  if (cv==1) saveRDS(GKRF.leuk.i, paste0(fnm, "dtr.gkRF.leuk_", nm, "_", cv, rds))
  if (cv==1) saveRDS(GKLM.leuk.i, paste0(fnm, "dtr.gkLM.leuk_", nm, "_", cv, rds))
  if (cv==1) saveRDS(DW.leuk.i, paste0(fnm, "dtr.dw.leuk_", nm, "_", cv, rds))
  rm(CSK.leuk.i, GKRF.leuk.i, GKLM.leuk.i, DW.leuk.i, zom.leuk.i); gc()
  saveRDS(values, paste0(fnm, "values_", nm, rds))
}
saveRDS(values, paste0(fnm, "values_", nm, rds))