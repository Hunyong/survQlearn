library(caret); library(purrr); library(dplyr); library(survival)
library(DTRreg); library(randomForestSRC); library(dtrSurv)
source("LF.helper.R")    # All the library and source files
source("F00.generic.R")
source("F21.GK1.R")
source("F21.GK2.R")
source("F22.DwSurv.R")


### 0. parameters
  Tx.nm = "ACAS"
  Tx.bin.nm = "AC" # For DW which does not admit more than two Tx arms.
  Tx.nm.list = paste(Tx.nm, 1:2, sep = ".")
  imp = 1
  
  
  ## time points
  tau = 2700
  timepoints = seq(0, sqrt(tau), length.out = 1000)^2
  
  ## other parameters
  nodesize = 50
  mindeath = round(sqrt(c(nodesize)), 0)
  Ntree = 300
  ert = TRUE; rs = 0.2 # randomSplit = 0.2
  

### 0.2 all_cause / criterion
  for (all_cause in c(FALSE, TRUE)) {
    all_cause_nm = if (all_cause) "_allcause" else ""
    for (value.criterion in c("mean", "surv.mean")) {
      if (value.criterion[1] != "mean") {
        value.criterion[2] = 2200 # six-year survival
        rule = "logrank"
      } else {
        rule = "mean"
      }
      
      cat("all_cause = ", all_cause,
          " value.criterion =",value.criterion[1],
          " value.s =", value.criterion[2], "tau =", tau, "\n")
      
      
      ## PM models
      modelPM = "Surv(V.%d, d.%d) ~ ACAS.%d + AA.%d + PREVHF01 + PRVCHD05 + hf + CENTER + gender + race + age + cig_yrs + bmi.%d + wth.%d + drink.%d + hypert.%d + glucose.%d + smoke.%d + hdl.%d + AAprev.%d + ACprev.%d + ASprev.%d"
      form.CSK <- lapply(1:2, function(i) {
        modelPM %>% gsub("%d", i, .) %>% 
          # Removing the xxprev.1, as they are only available from stage 2.
          {if (i == 1) {gsub("\\+ [A-z0-9]*prev\\.1", "", .)} else .} %>% as.formula
      })
      form.GK <- modelPM %>% gsub("\\.%d", "", .) %>% gsub("\\(V, d\\)", "(V, delta)", .) %>% as.formula
      print(form.CSK)
      
      ## For DW, due to the highdimensionality issue, only a few covariates are included.
      modelPM.DW = "~ AA.%d + AS.%d + PREVHF01 + PRVCHD05 + hf"
      form.DW <- lapply(1:2, function(i) gsub("%d", i, modelPM.DW) %>% as.formula)
      form.DW.tx <- lapply(1:2, function(i) paste0("factor(AC.%d) ", modelPM.DW) %>% gsub("%d", i, .) %>% as.formula)
      form.DW.cens <- lapply(1:2, function(i) paste0("delta ", modelPM.DW) %>% gsub("%d", i, .) %>% as.formula)
      print(form.DW)
      
      ## propensity models
      modelPr = "factor(ACAS.%d) ~ AA.%d + PREVHF01 + PRVCHD05 + hf + CENTER + gender + race + age + cig_yrs + bmi.%d + wth.%d + drink.%d + hypert.%d + glucose.%d + smoke.%d + hdl.%d"
      form.weight <- list(modelPr %>% gsub("%d", 1, .) %>% as.formula, modelPr %>% gsub("%d", 2, .) %>% as.formula)
      
      modelPrbin = "factor(AC.%d) ~ AA.%d + AS.%d + PREVHF01 + PRVCHD05 + hf + CENTER + gender + race + age + cig_yrs + bmi.%d + wth.%d + drink.%d + hypert.%d + glucose.%d + smoke.%d + hdl.%d"
      form.weight.bin <- list(modelPrbin %>% gsub("%d", 1, .) %>% as.formula, modelPrbin %>% gsub("%d", 2, .) %>% as.formula)
      print(form.weight)  
      
      
      ### 1. data preprocessing
      dat.aric.list <- readRDS(sprintf("../data_ARIC/03.aric.comp%s.rds", all_cause_nm)) ###TBD
      dat.aric = dat.aric.list[[imp]] %>%   # First imputed dataset.
        mutate(T.0 = pmin(T.0, tau), V.1 = pmin(V.1, tau), V.2 = T.0 - V.1, # Treat those administratively censored as failures for RMST.
               last.fu = ifelse(V.1 == tau, 1, last.fu),
               d.1 = ifelse(T.0 == tau & last.fu == 1, 1, d.1), d.2 = ifelse(T.0 == tau, ifelse(last.fu == 1, NA, 1), d.2),
               delta = ifelse(T.0 == tau, 1, delta)) %>% 
        mutate(ACAS.1 = ifelse(is.na(AC.1), NA, paste0(AC.1, AS.1)) %>% factor,
               ACAS.2 = ifelse(is.na(AC.2), NA, paste0(AC.2, AS.2)) %>% factor)
      
      mean(dat.aric$d.2==0, na.rm = TRUE)
      
      lvls = lapply(Tx.nm.list, function(x) levels(dat.aric[, x]))
      for (x in 1:length(lvls)) {if (is.null(lvls[[x]])) lvls[[x]] = 0:1}
      
      ## cross-validation
      K = 300
      set.seed(100, kind="Mersenne-Twister")
      cv.insample  = createDataPartition(dat.aric[, Tx.nm.list[1]], p = .8, list = FALSE, times = K)
      cv.outsample = map_dfc(1:K, ~which(!(1:dim(dat.aric)[1] %in% cv.insample[, .]))) %>% as.matrix
      
      
      # skeleton
      values <- matrix(NA, K, 11, 
                       dimnames = list(1:K, c("CSK", "GKRF", "GKLM", "DW", "ZOM", "observed", "cens1", "cens2",
                                              "ns.CSK", "ns.GK", "ns.GKLM")))
      
      # common parameters
      stg = 1:2  # stages
      
      # common file name
      nm = paste0(value.criterion[1], "_", 
                  round(as.numeric(value.criterion[2]), 1), "_", 
                  rule, "Split_", "tau_", tau, all_cause_nm, "_imp", imp)
      fnm = paste0("output/", Sys.Date(), "/") # folder name
      if (!dir.exists("output")) dir.create("output")
      if (!dir.exists(fnm)) dir.create(fnm)
      rds = paste0("_", Sys.Date(), ".rds")
      
      tab = data.frame(s1 = rep(NA, 300), s2 = NA, total = NA)
      
      for (cv in 1:300) {
        cat(cv, "th cv.\n")
        set.seed(cv)
        in.cv = cv.insample[, cv]   # insample index
        out.cv = cv.outsample[, cv] # outsample index
        train = dat.aric[in.cv, ]
        test  = dat.aric[out.cv,]
        values[cv, "cens1"] = mean(train$d.1==0); values[cv, "cens2"] = mean(train$d.2==0, na.rm = T)
        
        
        ############################################################################################################
        #### A. rule estimation ####################################################################################
        ############################################################################################################
        
        ### A1. The proposed method
        
        args.CSK <- list(data = train, 
                         txName = Tx.nm.list,
                         models = form.CSK,
                         usePrevTime = FALSE, tau = tau, timePoints = timepoints,
                         criticalValue = value.criterion[1], evalTime = as.numeric(value.criterion[2]), 
                         splitRule = ifelse(value.criterion[1] == "mean", "mean", "logrank"),
                         ERT = ert, uniformSplit = ert, replace = !ert,
                         randomSplit = rs, nTree = Ntree, mTry = c(6, 6),
                         pooled = FALSE, stratifiedSplit = FALSE)
        
        # actual fitting
        values[cv, "ns.CSK"] = nodesize
        set.seed(cv)
        
        CSK.aric.i <- 
          try(do.call(dtrSurv, c(args.CSK, list(nodeSize = nodesize, minEvent = mindeath))))
        err.CSK = class(CSK.aric.i)[1] == "try-error"
        
        
        ### A2. Goldberg-Kosorok RF
        args.GK = list(common.formula = form.GK,
                       common.Tx.label = Tx.nm, stage.label = 1:2, tau = tau,
                       data = train, stage.sep = ".", regress.prev.time = T)
        
        # actual fitting
        values[cv, "ns.GK"] = nodesize
        set.seed(cv)
        GKRF.aric.i <-
          try(do.call(gk.separate, c(args.GK, list(nodesize = nodesize, method = "rf"))))
        err.GKRF = class(GKRF.aric.i)[1] == "try-error"
        
        ### A3. Goldberg-Kosorok linear
        set.seed(cv)
        GKLM.aric.i <-
          try(gk.separate(common.formula = form.GK,
                          common.Tx.label = Tx.nm, stage.label = 1:2, tau = tau,
                          data = train, stage.sep = ".", method = "lm", regress.prev.time = T))
        err.GKLM = class(GKLM.aric.i)[1] == "try-error"
        
        
        ### A4. Simoneau et al.
        set.seed(cv)
        DW.aric.i <-
          try(DWSurv(time = list(~ V.1, ~ V.2),
                     blip.mod = form.DW,
                     treat.mod = form.DW.tx,
                     tf.mod = form.DW,
                     cens.mod = form.DW.cens,
                     data = train))
        err.DW = class(DW.aric.i)[1] == "try-error"
        
        ### A5. zero-order model
        zom.aric.i <-
          try(do.call(dtrSurv, c(args.CSK, list(nodeSize = 1e+4, minEvent = 1e+4))))
        err.zom = class(zom.aric.i)[1] == "try-error"
        
        if (!err.zom) {
          zom.pred =
            data.frame(
              zom.aric.i@stageResults[[1]]@optimal@optimalTx[1],
              zom.aric.i@stageResults[[2]]@optimal@optimalTx[1])
          tab[cv, 1] = zom.aric.i@stageResults[[1]]@optimal@optimalTx[1]
          tab[cv, 2] = zom.aric.i@stageResults[[2]]@optimal@optimalTx[1]
          tab[cv, 3] = paste0(tab[cv, 1], "-", tab[cv, 2])
          cat("ZOM freq table\n"); tab[, 3] %>% table %>% print
        }
        
        ############################################################################################################
        #### B. Rules applied to a test set ########################################################################
        ############################################################################################################
        
        ### B0. skeletons for the predicted A's
        # opt.rule.pred = matrix(NA, dim(cv.outsample)[1], 2, dimnames = list(NULL, stg))
        opt.rule.pred = data.frame(A.1 = rep(NA, dim(cv.outsample)[1]), A.2 = NA)
        for (q in 1:2) {
          opt.rule.pred[, q] = factor(NA, levels = lvls[[q]])
        }
        opt.rule.CSK <- opt.rule.GKRF <- opt.rule.GKLM <- opt.rule.zom <- opt.rule.pred
        opt.rule.DW =  data.frame(A.1 = rep(NA, dim(cv.outsample)[1]), A.2 = NA)
        
        for (q in seq_along(stg)) {
          elig = !test[, Tx.nm.list[q]] %>% is.na
          
          ### B1. the proposed method
          
          if (!err.CSK) {
            opt.CSK =
              predict(CSK.aric.i, 
                      # newdata = test[elig, ],
                      newdata = test[elig, ], 
                      stage = q)
            opt.rule.CSK[elig, q] = factor(opt.CSK$optimal@optimalTx, levels = lvls[[q]]) # %>% as.numeric()
          }
          ## B2. Goldberg-Kosorok RF
          if (!err.GKRF) {
            opt.GKRF = predict.opt.sep(GKRF.aric.i$survRF[[q]],
                                       newdata = test[elig, ] %>% dplyr::select(-V.2, -d.2) %>%
                                         mutate(prevTime.2 = V.1),
                                       Tx.label = paste0(Tx.nm, ".", q))
            opt.rule.GKRF[elig, q] = lvls[[q]][opt.GKRF$optimal.Tx]
          }
          ### B3. Goldberg-Kosorok linear
          if (!err.GKLM) {
            opt.GKLM = predict.opt.sep(GKLM.aric.i$survRF[[q]],
                                       newdata = test[elig, ] %>% dplyr::select(-V.2, -d.2) %>%
                                         mutate(prevTime.2 = V.1),
                                       Tx.label = paste0(Tx.nm, ".", q))
            opt.rule.GKLM[elig, q] = lvls[[q]][opt.GKLM$optimal.Tx]
          }
          
          ### B4. Simoneau et al.
          # DW (effect coded -> dummy)
          if (!err.DW) {
            opt.DW = decision.rule.dw(DW.aric.i, newdata = test[elig, ] %>% dplyr::select(-V.2, -d.2) %>%
                                        mutate(lfast.2 = ifelse(is.na(lfast.2), -1, lfast.2)), stage = q)$rule.fix
            opt.rule.DW[elig, q] = opt.DW*.5 + .5
          }
          
          ### B5. zero-order model
          if (!err.zom) {
            opt.rule.zom[elig, q] = zom.pred[1, q]
            opt.rule.zom[elig, q] = factor(zom.pred[1, q], levels = lvls[[q]]) # %>% as.numeric()
          }
        }
        print(table(opt.rule.CSK, useNA = "always"))
        print(table(opt.rule.zom, useNA = "always"))
        
        
        ############################################################################################################
        #### C. Value estimation ###################################################################################
        ############################################################################################################
        # weights
        weight = weights_aric (data = test, weight.formula.list = form.weight, 
                               weight.formula.bin.list = form.weight.bin)
        
        
        test.tmp = test %>% transmute(A.1 = test[, Tx.nm.list[[1]]] , A.2 = test[, Tx.nm.list[[2]]])
        test.tmp.DW = test %>% transmute(A.1 = test[, paste0(Tx.bin.nm, ".1")] , A.2 = test[, paste0(Tx.bin.nm, ".2")])
        
        arg.val = list(test = test, actual = test.tmp, propensity = weight$propensity,
                       weight.censor = weight$weight.censor, criterion = value.criterion,
                       tau = tau)
        if (!err.CSK)  values[cv, "CSK"] = do.call(getValue, c(arg.val, list(estimated = opt.rule.CSK)))
        if (!err.GKRF) values[cv, "GKRF"] = do.call(getValue, c(arg.val, list(estimated = opt.rule.GKRF)))
        if (!err.GKLM) values[cv, "GKLM"] = do.call(getValue, c(arg.val, list(estimated = opt.rule.GKLM)))
        values[cv, "observed"] = do.call(getValue, c(arg.val[-which(names(arg.val) == "propensity")], list(estimated = test.tmp, propensity = 1)))
        if (!err.zom)  values[cv, "ZOM"] = do.call(getValue, c(arg.val, list(estimated = opt.rule.zom)))
        
        if (!err.DW) {
          arg.val.DW = list(test = test, actual = test.tmp.DW, 
                            propensity = weight$propensity.DW,
                            weight.censor = weight$weight.censor, 
                            criterion = value.criterion,
                            tau = tau, estimated = opt.rule.DW)
          values[cv, "DW"] = do.call(getValue, arg.val.DW)
        }
        
        print(values[cv, ])
        print(apply(values, 2, mean, na.rm = TRUE))
        
        ############################################################################################################
        #### C. saving the results #################################################################################
        ############################################################################################################
        
        
        attr(values, "spec") = data.frame(criterion = value.criterion[1], criterion.s = as.numeric(value.criterion[2]), 
                                          rule = rule, tau = tau,
                                          ert = ert, rs = rs, ntree = Ntree,
                                          nodesize = nodesize, mindeath = mindeath)
        if (cv==1) saveRDS(CSK.aric.i, paste0(fnm, "dtr.csk.aric_", nm, "_", cv, rds))
        if (cv==1) saveRDS(GKRF.aric.i, paste0(fnm, "dtr.gkRF.aric_", nm, "_", cv, rds))
        if (cv==1) saveRDS(GKLM.aric.i, paste0(fnm, "dtr.gkLM.aric_", nm, "_", cv, rds))
        if (!err.DW & cv < 10) saveRDS(DW.aric.i, paste0(fnm, "dtr.dw.aric_", nm, "_", cv, rds))
        saveRDS(values, paste0(fnm, "values_", nm, rds))
        
      }
      saveRDS(values, paste0(fnm, "values_", nm, rds))
      saveRDS(tab, paste0(fnm, "tab_ZOM_", nm, rds))
    }
  }

  
  
  
