gk.separate <-
  function(common.formula, 
           common.Tx.label = "A",
           regress.prev.time = !any(grepl("prevTime", names(data))),
           stage.label = NULL,
           formula.list = NULL, 
           Tx.label.list = NULL,
           data = NULL, 
           stage.sep = ".",
           # timepoints = NULL,
           # timeticks = c("quadratic", "uniform", "exponential"),
           # ntime = 100L,
           method = c("rf", "lm"),
           tau,
           # value.crit = c("mean", "surv.prob", "surv.mean"),
           # value.s.time = if (value.crit %in% c("surv.prob", "surv.mean")) tau/2,
           tie.method = c("random", "first", "NA"), 
           # printStages = FALSE,
           # do.trace = printStages,
           ...
           # , na.action = na.fail) {
  ) {
    ### part of the code originally from svm.formula (package e1071).
    
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    n    = dim(data)[1]
    tie.method = match.arg(tie.method)
    method = match.arg(method)
    methodFun = switch(method, rf = gk.rf, lm = gk.lm)
    # value.crit = match.arg(value.crit)
    
    # getting Y labels # Y.1, Y.2, ...
    if (!is.null(common.formula)) {
      # Getting common.Y.label  # Y.1, Y.2, ...
      y.lab = paste0(deparse(common.formula[[2]][[2]]), stage.sep)
      # getting stage.label first
      if (is.null(stage.label)) { # If there is no stage.label given:
        stage.label = as.numeric(greps(y.lab, "", names(data), fixed = TRUE))
        if (any(is.na(stage.label))) warning("There is a non-numeric stage.label index. That stage index is ignored")
        stage.label <- sort(stage.label)
      } else {
        if (any(is.na(stage.label))) warning("There is a non-numeric stage index. That stage index is ignored")
        stage.label <- as.numeric(sort(stage.label))
      }
      Y.label = paste0(y.lab, stage.label)
    } else {
      Y.label <-
        sapply(formula.list, function(s) {
          deparse(s[[2]][[2]])
        })
      stage.label <- 1:length(Y.label) # artificially name the stages.
    }
    n.stages = length(stage.label)
  
    # Force 0 values in place of NA.
    data[, Y.label][is.na(data[, Y.label])] = 0 
    
    # adding the prevTime into data.
    if (n.stages < 2 || is.null(common.formula)) regress.prev.time = FALSE # Not applicable.
    if (regress.prev.time) {
      # prevTime = t(apply(data[, Y.label[-n.stages], drop = FALSE], 1, cumsum))
      prevTime = matrix(0, n, n.stages - 1)
      prevTime[, 1] = data[, Y.label[1]] 
      for (i in seq_along(stage.label)[-c(1, n.stages)]) {
        prevTime[, i] = prevTime[, i - 1] + data[, Y.label[i]]
      }
      colnames(prevTime) = paste0("prevTime", stage.sep, stage.label[-1])   #prevTime.2, prevTime.2, ...
      data = cbind(data, prevTime)
    }
    # update data in m
    m$data = if (is.matrix(data)) as.data.frame(data) else data
    
    # getting all the variable names in data
    vars = names(data)
    
    # getting formula.list
    if (!is.null(common.formula)) {
      if (!inherits(common.formula, "formula")) stop("method is only for formula objects")
      if (!is.null(formula.list)) stop("Provide only one of formula.list and common.formula")
      if (regress.prev.time) common.formula = update.formula(common.formula, ~ . + prevTime)
      xvars = attr(terms(common.formula), "term.labels")
      if (!any(grepl(common.Tx.label, xvars))) stop("common.Tx.label is not included in the common.formula")
      for (i in xvars) {
        if (!any(grepl(i, vars))) stop(paste0 (i, " in common.formula is not in data"))
      }
      
      # Checking common.Tx.label
      for (i in stage.label) {
        if (! paste0(common.Tx.label[1], stage.sep, i) %in% vars) 
          stop(paste0("There does not exist ", common.Tx.label[1], stage.sep, i, " in the data."))
      }
      Tx.label.list = paste0(common.Tx.label[1], stage.sep, stage.label)
      
      # completing formula.list
      formula.list <- lapply(stage.label, function(s) common.formula)
      names(formula.list) <- as.character(stage.label)
      
      for (i in seq_along(stage.label)) {
        # removing non-existing variables for each stage.
        xvars.i = paste0(xvars, stage.sep, i)
        xvars.index  = (xvars.i %in% vars)
        xvars.common = (xvars   %in% vars)  # to make sure the plain variable names (without the stage index) can be considered.
        xvars.i = ifelse(xvars.index, xvars.i, ifelse(xvars.common, xvars, NA))
        xvars.i = xvars.i[!is.na(xvars.i)]
        
        LHS = paste0(common.formula[[2]][[1]], "(", 
                     Y.label[i] , ", ", common.formula[[2]][[3]], ")") # Surv(Y.i, delta)
        RHS = paste0(xvars.i, collapse = " + ")
        formula.list[[i]] <- 
          as.formula(paste0(LHS, " ~ ", RHS))
      }
    } else { # if formula.list is provided
      if (any(sapply(formula.list, function(s) !inherits(s, "formula"))))
        stop("method is only for formula objects")
      names(formula.list) <- as.character(stage.label)
      # regress.prev.time not supported for formula.list
      
      for (i in seq_along(formula.list)) {
        if (!grepl(Tx.label.list, attr(terms(formula.list[[i]]), "term.labels"))) 
          stop("Tx.labels is not included in the elements of formula.list")
      }
      
      stop("TBD")
    }
    
    delta.label = deparse(formula.list[[1]][[2]][[3]]) #"delta"

    # stage eligibility index
    # n x n.stages eligibility matrix
    eligibility <- data[, Y.label] != 0
    for (i in seq_along(stage.label)) {
      xvars.i = attr(terms(formula.list[[i]]), "term.labels")
      xvars.i = c(delta.label, xvars.i)
      eligibility[, i] = eligibility[, i] & complete.cases(data[, xvars.i])
    }    

    # censoring KM (magrinal)
    Y.cum <- t(apply(data[, Y.label], 1, cumsum))
    terminal <- apply(data[, Y.label], 1, function(x) max(c(which(x > 0), 1)))
    delta.vec <- data[, delta.label]
    terminal.delta1 <- ifelse(delta.vec, terminal, terminal - 1)
    delta.mat <- t(sapply(terminal.delta1, function(x) c(rep(1, x), rep(0, n.stages - x))))
    Sc.hat <- lapply(1:n.stages, function(s) survival::survfit(survival::Surv(Y.cum[, s], 1 - delta.mat[, s]) ~ 1))
    weight = sapply(1:n.stages, function(s) S.t(Y.cum[, s], Sc.hat[[s]], rightcts = FALSE)) #S_C(sum(R_t))
    weight = delta.mat/weight # delta_t / S_C(sum(R_t))
    
    ### Transforming data to an auxilliary problem
    ## trimming Y only up to tau.
    Y.cum2 <- pmin(cbind(0, Y.cum), tau)
    R2 <- Y.cum2[, -1] - Y.cum2[, -(n.stages + 1)]
    data[, Y.label] <- R2
    
    ## Tx levels
    Tx.levels = lapply(Tx.label.list, function(x) sort(unique(data[, x])))
    # 1.2 Adding random draws to the forced-to-be-missing fields
    overflow.index <- is.na0(R2)
    for (q in seq_along(stage.label)) {
      data[overflow.index[, q], Tx.label.list[q]] <- 
        sample(Tx.levels[[q]], sum(overflow.index[, q]), replace = TRUE)
    }
    
    # return(list(formula.list = formula.list, data = data))    
    for (i in names(m)) {
      if (!i %in% c("", "data")) m[[i]] <- NULL
    }
    m[[1]] <- as.name("model.frame")
    m$na.action <- na.pass  # NA's are handled in eligibility matrix (they are omitted).
    
    # decision rule
    survRF.list <- structure(vector("list", n.stages), names = stage.label)
    res <- list(survRF = survRF.list,
                # optimal.survival = survRF.list,
                optimal.mean = survRF.list,
                optimal.Tx = matrix(NA, n, n.stages, dimnames = list(rownames(data), stage.label)),
                value.train = NULL,
                stages = stage.label,
                Terms  = list(time = Y.label,
                              delta = delta.label,
                              Tx   = Tx.label.list,
                              x    = lapply(formula.list, function(s) attr(terms(s), "term.labels"))),
                formula = formula.list,
                # timepoints = timepoints,
                eligibility = eligibility)
    
    # New formulae for lm and rfsrc (y only without delta)
    formula.list.new <- formula.list
    for (q in seq_along(stage.label)) {
      formula.list.new[[q]][[2]] <- formula.list.new[[q]][[2]][[2]]  # replace Surv(Y.q, delta) with Y.q
    }
    
    for (q in n.stages:1L) {
      
      m$formula <- formula.list[[q]]
      m.eval <- eval(m, parent.frame())
      #rn <- 1:nrow(m.eval)
      y <- model.response(m.eval)
      attr(y, "na.action") <- attr(m.eval, "na.action")
      Terms <- attr(m.eval, "terms")
      attr(Terms, "intercept") <- 0
      ## Drop any "negative" terms in the formula.
      x.q   <- model.frame(terms(reformulate(attributes(Terms)$term.labels)),
                           data.frame(m.eval), na.action = na.pass)
      x.q   <- x.q[eligibility[, q], ]
      
      y.q   <- y[eligibility[, q], 1L]
      delta <- y[eligibility[, q], 2L]
      
      for (i in seq_along(x.q)) {
        if (is.ordered(x.q[[i]])) x.q[[i]] <- as.numeric(x.q[[i]])
      }
      n.q <- sum(eligibility[, q])
      # if (printStages)
      #   cat("Stage ", q, ", sample size: ", n.q, ", model: ", deparse(formula.list[[q]]), "\n")
      # Augmenting
      po = y.q # pseudo-outcome
      if (q == n.stages) {
        # pr = y2prob(y.q, timepoints)
      } else {
        # initialize the surv.matrix: by default S(t) = 1(t <= 0)
        # survMatrix = matrix(0, n.q, ntime)
        # survMatrix[, 1] = 1
        
        # fill in surv.matrix
        elig.next = eligibility[, q + 1][eligibility[, q]] # eligible on q+1st stage among those eligible on qth stage.
        # survMatrix[elig.next, ] = res$optimal.survival[[q + 1]]    # next stage survival probability (estimated from the previous iteration)
        # pr = shiftMat (timepoints = timepoints, surv.matrix = survMatrix, 
        #                shift.vector = y.q, ntime = ntime, nsample = n.q,
        #                surv2prob = TRUE)
        po [elig.next] = po [elig.next] + res$optimal.mean[[q + 1]]
      }
      
      Tx.lab = Tx.label.list[q]
      Tx.levs = sort(unique(x.q[, Tx.lab]))
      Tx.group = lapply(Tx.levs, function(x) x.q[, Tx.lab] == x)
      # storing the surv object (To be returned as one of the final outputs for later use)
      # res$survRF[[q]] <- survRF.default(x = x.q, pr = pr, delta = delta, timepoints = timepoints, 
      #                                   do.trace = do.trace, ...)
      formula.list.new[[q]] <- update.formula(formula.list.new[[q]], as.formula(paste0("~. -", Tx.lab)))
      res$survRF[[q]] <- 
        lapply(Tx.levs, function(s) {
          cat("Tx :", s, "\n")
          group = data[eligibility[, q], ][, Tx.lab] == s
          methodFun(data = data[eligibility[, q], ][group, ], 
                    formula = formula.list.new[[q]], 
                    weight = weight[eligibility[, q], q][group], ...)
        })
      # class(res$survRF[[q]]) = method
      
      # storing the estimated optimal decision
      rule <- predict.opt.sep(res$survRF[[q]], newdata = data[eligibility[, q], ],
                              Tx.label = Tx.label.list[[q]], tie.method = tie.method)
      res$optimal.mean[[q]] <- rule$optimal.value
      res$optimal.Tx[eligibility[, q], q] <- rule$optimal.Tx
      gc()
    }
    
    cl <- match.call()
    cl[[1]] <- as.name("GK")
    res$call <- cl
    class(res) <- c("GK")
    return(res)
  }

predict.opt.sep <- function(object.list, newdata, Tx.label, tie.method = c("random", "first", "NA")) {
  tie.method = match.arg(tie.method)
  ntest = dim(newdata)[1]
  
  # Check treatment levels
  Tx.levels = unique(sort(newdata[ , Tx.label]))
  n.Tx      = length(Tx.levels)
  if (n.Tx <= 1) {
    n.Tx = length(object.list)
    Tx.levels = 1:n.Tx
  }
  
  # pseudo data
  newdata2 <- newdata
  output.value = matrix(NA, ntest, n.Tx, dimnames = list(NULL, Tx.levels))
  # tmp2 <<- list(object.list = object.list, newdata = newdata)
  predictor = switch(class(object.list[[1]])[1], 
                     "rf" = function(...) predict(...)$predicted,
                     "rfsrc" = function(...) predict(...)$predicted,
                     "lm" = function(...) predict(...))
  
  for (i in 1:n.Tx) {
    # newdata2[, Tx.label] <- Tx.levels[i] # not needed here.
    output.value [ ,i] = predictor(object.list[[i]], newdata = newdata2)
  }
# print(head(output.value))
  optimal.Tx = output.value %>% apply(1, whichMax, tie.method = tie.method)
  optimal.value = output.value[cbind(1:ntest, optimal.Tx)]
  return(list(optimal.Tx = optimal.Tx, optimal.value = optimal.value))
}
