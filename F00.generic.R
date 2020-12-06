library(dplyr)
search2 <- function(s, xvec, include = TRUE) {
  # Finding the index of "the largest x smaller or equal to s"
  # Fn is the cdf table (x (should be ordered) and y)
  # index = which (s >= c(-Inf, xvec) & s <= c(xvec, Inf)) - 1
  index = if (include) {which (s >= c(-Inf, xvec)) - 1} else {which (s > c(-Inf, xvec)) - 1}
  max(index)
}
search2.vec <- Vectorize(search2, vectorize.args = "s")

# search2(49.25, xvec =  48:52)
lin.interpolate <- function(s, xvec, yvec = NULL, return.y.only = FALSE, exponential.tail = TRUE) {
  x.len = length(xvec)
  low.index <- search2(s, xvec)
  lower <- xvec[low.index]
  if (low.index == x.len) {
    upper <- lower
    proportion <- 0
  } else {
    upper <- xvec[low.index + 1]
    proportion <- (s - lower) / (upper - lower)
  }
  result <- c(low.index = low.index, proportion = proportion)
  if (!is.null(yvec)) {
    lower.y <- yvec[low.index]
    upper.y <- yvec[min(x.len, low.index + 1)]
    if (is.infinite(upper.y) && exponential.tail) {
      y.interpolate <- lower.y * log (1 - s) / log (1 - lower)
    } else {
      y.interpolate <- lower.y + proportion * (upper.y - lower.y)
    }
    if (return.y.only) return(y.interpolate)
    result["y.interpolate"] <- y.interpolate
  }
  return(result)
}
lin.interpolate.vec <- Vectorize(lin.interpolate, vectorize.args = "s")

# time calculator
tt <- function(s, reset = FALSE, units = "auto"){
  if (s==1) {time.tmp <<- Sys.time() # record time
  } else if (s==2) { # calculate time
    result <- data.frame(begin = time.tmp, end = Sys.time(), elapsed = difftime(Sys.time(), time.tmp, units = units))
    if (reset) time.tmp <<- Sys.time()
    return(result)
  }
}

# either 0 or na gives TRUE
# used in G01.Goldberg.R::auxiliary()
# is.na0(c(1:3,NA, 0, 0, 3))
is.na0 <- function(x) {
  is.na(x) | (x == 0)
}

# complete the delta vector so that the final status is inherited to the next NA values
# complete.delta.vec(c(1,0,NA, 0))
complete.delta.vec <- function(vec, n.stages = length(vec)) {
  ref <- min(which(is.na(vec)), n.stages + 1) - 1
  delta.final <- vec[ref]
  if (ref < n.stages)
  vec[ref:n.stages] <- delta.final
  vec
}

complete.delta.mat <- function(mat, n.stages = dim(mat)[2]) {
  apply(mat, 1, complete.delta.vec) %>% t
}
# tmp[1:10,"delta",]
# complete.delta.mat(tmp[1:10,"delta",])

S.t <- function(t, S, rightcts = TRUE) {
  index = search2.vec(t, S$time, include = rightcts)
  S$surv[pmax(index, 1)]
}

flowchart <- function(outputData) {
  n.stages = dim(outputData)[3]
  n.sample = dim(outputData)[1]
  result <- data.frame(stage = c(1:n.stages, "Total", "Percent"),
                       total = rep(NA, n.stages + 2), percentage = NA,
                       censored = NA, died = NA, nextStage = NA)
  for (stage in 1:n.stages) {
    result[stage, "total"]  = dim(outputData)[1]
    result[stage, "percentage"]  = round(result[stage, "total"] / n.sample,2)
    cens <- outputData[, "delta", stage] == 0
    # drop those censored
    outputData <- outputData[!cens,,, drop = FALSE]
    died <- outputData[, "gamma", stage] == 1
    # drop those died
    outputData <- outputData[!died,,, drop = FALSE]
    result[stage, "censored"]  = sum(cens, na.rm = TRUE) 
    result[stage, "died"]      = sum(died, na.rm = TRUE)
    result[stage, "nextStage"] = sum(!died, na.rm = TRUE)
  }
  result[n.stages + 1, c("censored", "died")] <- apply(result[, c("censored", "died")], 2, sum, na.rm = TRUE)
  result[n.stages + 2, c("censored", "died")] <- round(result[n.stages + 1, c("censored", "died")]/n.sample, 2)
  # result[c("Total", "Percent"), 1:2]
  result
}

output2observable <- function(output, stage = NULL, cumulative = FALSE) {
  n.stages = dim(output)[3]
  if (is.null(stage)) {
    stage = 1:n.stages
  } else if (cumulative) {
    stage = 1:max(stage)
  }
  
  n = dim(output)[1]
  nm = dimnames(output)[[2]]
  nm.covar = c("lB", grep("Z[0-9]+", nm, value = TRUE))
  nm.stage =  c("event.time", "delta", "action", nm.covar)
  
  df <- data.frame(subject.id = output[, "subj.id", 1])
  for (i in stage) {
    df.i <- data.frame(output[, nm.stage, i])
    names(df.i) <- paste(c("T", "delta", "A", nm.covar), i, sep = "_")
    df <- cbind(df, df.i)
  }
  df$delta = 
    df %>% dplyr::select(starts_with("delta_")) %>% as.matrix %>% 
    apply(1, function(s) 1 - any(s == 0, na.rm = TRUE))
  df
}

# output2observable(obs.data$output)

whichMax <- function (x, tie.method = "random") {
  ind = which(x == max(x))
  if (length(ind) == 1L) return(ind)
  if (tie.method == "random") {
    sample(ind, 1)
  } else if (tie.method == "first") {
    ind[1]
  } else {
    NA
  }
}
