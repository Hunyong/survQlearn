#### C21.simulation_body.R is to be run by source("C21.simulation_run.R")

#### library
# install.packages("survival")
# install.packages("randomForestSRC")
# install.packages("DTRreg")
# install.packages("dtrSurv")
library(dtrSurv)
source("F00.generic.R")
source("F02.multiStage.R")
source("F21.GK1.R")
source("F21.GK2.R")
source("F22.DwSurv.R")

    if (criterion == "mean") {
      val.fn <- function(x) {mean(x, na.rm = TRUE)}
    } else if (criterion %in% c("surv.prob", "surv.mean")) {
      val.fn <- function(x) {mean(x >= crit.value, na.rm = TRUE)}
    }
    filename.tmp <- gsub("\\.rds", "_tmp.rds", filename)
    stat.stage <- matrix(NA, nrow = n.sim, ncol = n.stages, 
                         dimnames = list(1:n.sim, 1:n.stages))
### simulation
    result <- data.frame(no = 1:n.sim, observed = NA, csk = NA, gkLM = NA, gkRF = NA, 
                         dw = NA, zom = NA,
                         time.obs = NA, time.csk = NA, time.gkLM = NA, time.gkRF = NA,
                         time.dw = NA, time.zom = NA, percent.censor = NA)  # time for each method (both policy est and eval)
    attr(result, "criterion") <- list(criterion = criterion, crit.value = crit.value)
    
    for (i in 1:n.stages) result[[paste0("n.", i)]] = NA
    arg.obs <- arg.csk <- arg.obs.no.censor <- arg.gk.lm <- arg.gk.rf <-  arg.dw <-  arg.zom <- 
      list(
        n.sample = n, n.stages = n.stages, tau = tau, tick = 0.01,  # structural parameters
        at.risk = 1,      # initial state vector
        p = p, corr = -0.5,                              # cor of two error processes
        predHazardFn = predHazardFn, predPropensityFn = predPropensityFn, 
        predCensorFn = predCensorFn,             # list of predictor functions
        hidden_data = TRUE, printFlag = FALSE
      )
    arg.obs.no.censor$predCensorFn <- arg.csk$predCensorFn <- arg.gk.lm$predCensorFn <- 
      arg.gk.rf$predCensorFn <- arg.dw$predCensorFn <- arg.zom$predCensorFn <-
      noCensorFn
    arg.obs.no.censor$n.sample <- arg.csk$n.sample <- arg.gk.lm$n.sample <- 
      arg.gk.rf$n.sample <- arg.dw$n.sample <- arg.zom$n.sample <-
      n.eval
    nodesize = 5
    mindeath = round(sqrt(c(nodesize)), 0)
    
    # flow: obs (1) -> optimal (n.mc), obs.no.censor (n.mc)
    print(Sys.time())
    for (sim in 1:n.sim) {
      cat("###########################  simulation ", sim, "########################### \n")
      cat("########################### (criterion ", criterion, crit.value, ")################## \n")
      
      cat ("1. Data \n")
        tt(1)
        # obs.data
        set.seed(sim*10000)
        obs.data <- do.call(multiStageDynamics, arg.obs)
        
        # observed policy value
        set.seed(sim*10000)
        obs.data.rep <- do.call(multiStageDynamics, arg.obs.no.censor)
        result[sim, "observed"] <- val.fn(obs.data.rep$summary$cumulative.event.time)
        result[sim, "percent.censor"] <- 1 - mean(obs.data$summary$censor.status, na.rm = TRUE)
        result[sim, paste0("n.", 1:n.stages)] <- 
          sapply(1:n.stages, function(s) mean(obs.data$summary$terminal.stage == s))
        result[sim, "time.obs"] <- tt(2, reset = TRUE)["elapsed"]
        print(flowchart(obs.data$output))
        
        # transforming data from an array format to a data.frame format
        data.df = output2observable(obs.data$output)
        rm(obs.data.rep); gc()
      
      # estimation
      cat ("2. csk \n")
      if (!skip.csk) {
        cat ("  2. csk - Policy estimation \n")
        # new package dtrSurv
        models = 
          lapply(1:n.stages, function(x) {
            paste0(sprintf("Surv(T_%s, delta_%s) ~ ", x, x), 
                   paste(paste0(c(if(x>1) "lB_", paste0("Z", 1:p, "_")), x), collapse = " + ")) %>% 
              as.formula
          })
        arg.csk2 = list(data = data.df, 
                       txName = paste("A", 1:n.stages, sep = "_"),
                       models = models,
                       usePrevTime = TRUE, tau = tau, timePoints = "uni", nTimes = 200,
                       criticalValue = criterion, evalTime = crit.value, 
                       splitRule = ifelse(criterion == "mean", "mean", "logrank"),
                       ERT = TRUE, uniformSplit = TRUE, replace = FALSE,
                       randomSplit = 0.2, nTree = 300, mTry = rep(sqrt(p), n.stages),
                       pooled = FALSE, stratifiedSplit = 0.1)
        
        set.seed(sim*10000 + 1)
        optimal.csk <- do.call(dtrSurv, c(arg.csk2, list(nodeSize = nodesize, minEvent = mindeath )))
        csk.error <- class(optimal.csk)[1] == "try-error"
        arg.csk$policy <- if (!csk.error) optimal.csk
        rm(optimal.csk); gc()
        
        cat ("  2. csk - Evaluation \n")
        set.seed(sim*10000 + 10)
        if (!csk.error) csk.data.rep <- do.call(multiStageDynamics, arg.csk)
        if (!csk.error) result[sim, "csk"] <- val.fn(csk.data.rep$summary$cumulative.event.time)
        result[sim, "time.csk"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
        arg.csk$policy <- NULL; gc()
        rm(csk.data.rep); gc()
      }
        
      cat ("3. Goldberg & Kosorok - lm \n")
      if (!skip.gk) {
        cat ("  3. Goldberg & Kosorok - lm - Policy estimation \n")
        set.seed(sim*10000 + 2)
        optimal.gk.lm <-
          try(gk.separate(common.formula = formula.lm(p),
                          common.Tx.label = "A", stage.label = 1:n.stages, tau = tau,
                          data = data.df, stage.sep = "_", method = "lm", regress.prev.time = F))
        # optimal.gk.lm <- try(gk.Q(data = obs.data$output, est.fn = gk.lm, formula = formula.lm(p), tau = tau))
        gklm.error <- class(optimal.gk.lm)[1] == "try-error"
        arg.gk.lm$policy <- if (!gklm.error) optimal.gk.lm$survRF
        attr(arg.gk.lm$policy, "class") = "GKLM"
        rm(optimal.gk.lm); gc()
        
        cat ("  3. Goldberg & Kosorok - lm - Evaluation \n")
        set.seed(sim*10000 + 10)
        if (!gklm.error) gklm.data.rep <- do.call(multiStageDynamics, arg.gk.lm)
        if (!gklm.error) result[sim, "gkLM"]    <- val.fn(gklm.data.rep$summary$cumulative.event.time)
        result[sim, "time.gkLM"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
        arg.gk.lm$policy <- NULL; gc()
        rm(gklm.data.rep); gc()
        
      cat ("4. Estimation - Goldberg & Kosorok  - rf \n")
      
        cat ("  4. Estimation - Goldberg & Kosorok  - rf - Policy estimation \n")
        set.seed(sim*10000 + 3)
        arg.gk.rf2 = list(common.formula = formula.rf(p),
                         common.Tx.label = "A", stage.label = 1:n.stages, tau = tau,
                         data = data.df, stage.sep = "_", regress.prev.time = F)
        optimal.gk.rf <-
          try(do.call(gk.separate, c(arg.gk.rf2, list(nodesize = nodesize, method = "rf"))))
        # optimal.gk.lm <- try(gk.Q(data = obs.data$output, est.fn = gk.lm, formula = formula.lm(p), tau = tau))
        gkrf.error <- class(optimal.gk.rf)[1] == "try-error"
        arg.gk.rf$policy <- if (!gkrf.error) optimal.gk.rf$survRF
        if (!gkrf.error) attr(arg.gk.rf$policy, "class") = "GKRF"
        rm(optimal.gk.rf); gc()
        
        cat ("  4. Estimation - Goldberg & Kosorok  - rf - Evaluation \n")
        set.seed(sim*10000 + 10)
        if (!gkrf.error) gkrf.data.rep <- do.call(multiStageDynamics, arg.gk.rf)
        if (!gkrf.error) result[sim, "gkRF"]    <- val.fn(gkrf.data.rep$summary$cumulative.event.time)
        result[sim, "time.gkRF"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
        arg.gk.rf$policy <- NULL; gc()
        rm(gkrf.data.rep); gc()
      }
      
      cat ("5. Estimation - Simoneau et al. \n")
      if (!skip.dw) {
        cat ("  5. Estimation - Simoneau et al. - Policy estimation \n")
        set.seed(sim*10000 + 4)
        optimal.dw <- try(dwSurv(data = obs.data$output, tau = tau))
        dw.error <- class(optimal.dw)[1] == "try-error"
        arg.dw$policy <- if (!dw.error) optimal.dw$dw.est
        rm(optimal.dw); gc()
        
        cat ("  5. Estimation - Simoneau et al. - Evaluation \n")
        set.seed(sim*10000 + 10)
        if (!dw.error) dw.data.rep <- do.call(multiStageDynamics, arg.dw)
        if (!dw.error) result[sim, "dw"]  <- val.fn(dw.data.rep$summary$cumulative.event.time)
        result[sim, "time.dw"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
        arg.dw$policy <- NULL; gc()
        rm(dw.data.rep); gc()
      }
      
      cat ("6. Estimation - zero-order model\n")
      if (!skip.zom) {
        cat ("  6. zero-order model - Policy estimation \n")
        set.seed(sim*10000 + 5)
        optimal.zom <- try(dtrSurv:::dtrSurv(data = data.df, 
                                             txName = paste("A", 1:n.stages, sep = "_"),
                                             models = models,
                                             usePrevTime = TRUE,
                                             tau = tau,
                                             timePoints = "uni",
                                             nTimes = 200,
                                             criticalValue = criterion, evalTime = crit.value, 
                                             splitRule = ifelse(criterion == "mean", "mean", "logrank"),
                                             ERT = TRUE, uniformSplit = TRUE, replace = FALSE,
                                             randomSplit = 0.2, 
                                             minEvent = 1, nodeSize = 1e+4, # zero-order model
                                             nTree = 300, mTry = rep(1, n.stages),
                                             pooled = FALSE, stratifiedSplit = 0))
        zom.error <- class(optimal.zom)[1] == "try-error"
        arg.zom$policy <- if (!zom.error) optimal.zom
        rm(optimal.zom); gc()
        
        cat ("  6. zero-order model - Evaluation \n")
        set.seed(sim*10000 + 10)
        if (!zom.error) zom.data.rep <- do.call(multiStageDynamics, arg.zom)
        if (!zom.error) result[sim, "zom"] <- val.fn(zom.data.rep$summary$cumulative.event.time)
        result[sim, "time.zom"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
        arg.zom$policy <- NULL; gc()
        rm(zom.data.rep); gc()
      }
      
      ### saving and cleaning 
        print(result[sim, ])
        print(apply(result, 2, mean, na.rm = TRUE))
        saveRDS(result, filename.tmp) # saving the temporary results
        gc()
      # }
    }
    result
    print(Sys.time())
    saveRDS(list(statistics = result, settings = setting), filename)
    