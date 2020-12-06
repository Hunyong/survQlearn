
### 1. settings
date = "2020-12-5"

library(dplyr); library(ggplot2); library(cowplot)

settings <- 
  expand.grid(# model = c("full", "base"),
    crit  = c("mean", "surv.mean"),
    tau   = 450) %>% 
  mutate(crit.s = ifelse(crit == "mean", NA, 180),
         rule  = ifelse(crit == "mean", "mean", "logrank"),
         tau = tau,
         nm = paste0(crit, "_", as.numeric(crit.s), "_", 
                     rule, "Split_", "tau_", tau))
rds = paste0("_", date, ".rds")

### 2. some functions
meanSD <- function(mat) {
  cbind(mean = apply(mat, 2, mean, na.rm = T),
        sd   = apply(mat, 2, sd, na.rm = T),
        na   = apply(mat, 2, function(s) mean(is.na(s)))
  ) %>% as.data.frame %>% 
    tibble::rownames_to_column("method")
}
gathr <- function(mat) tidyr::gather(as.data.frame(mat), key = "method", value = "value")

### 3. summarize
value_summary <-
  lapply(1:dim(settings)[1], function(i) {
    readRDS(paste0("output_leuk/",date,"/values_", settings[i, "nm"], rds)) %>% 
      gathr %>% 
      # meanSD %>% 
      data.frame(scenario = i, settings[i, ], .)
  }) %>% do.call(rbind, .)
(value_summary %>% dplyr::filter(method == "CSK", crit == "mean"))$value %>% is.na %>% mean # 72% not done!

p <- list()
for (crit1 in c("mean", "surv.mean")) {
  binary1 = 0
  
  max.y = if (crit1 == "mean") 400 else 0.9 # y-coordinate for the mean labels
  min.y = if (crit1 == "mean") 100 else 0.3
  round.y = if (crit1 == "mean") 0 else 2
  lab.y = if (crit1 == "mean") "Truncated mean survival time" else "Six month survival probability"
  
  val.tmp <-
    value_summary %>% 
    dplyr::filter(crit == crit1)
 
  p[[crit1]] <-
    val.tmp %>% 
    dplyr::filter(!method %in% c("cens1", "cens2", "ns.CSK", "ns.GK", "ns.GKLM")) %>% 
    mutate(method = factor(method, levels = c("CSK", "GKRF", "GKLM","DW", "ZOM", "observed"),
                           labels = c("the proposed method", "Goldberg & Kosorok (2012), RF", 
                                      "Goldberg & Kosorok (2012), linear", "Simoneau et al. (2019)", "zero-order model", "observed"))) %>% 
    ggplot(aes(method, value, col = method)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, height = 0) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    stat_summary(aes(x = as.numeric(method), y = value),
                 fun.y = mean, geom = 'point',  col = "black", shape = "square", size = 2) +
    stat_summary(aes(x = as.numeric(method) + 0.3, y = value, label = round(..y.., round.y)),
                 fun.y = mean, geom = 'text',  col = "black", position = 'dodge', vjust = -3) +
    theme_bw() + 
    ylab(lab.y) + xlab("") +
    guides(col = FALSE)
}

p2 <- plot_grid(p[[1]], p[[2]], align = "h", nrow = 1)
save_plot(paste0("figure_leuk/leuk_", date,".eps"), p2, base_height = 7, base_width = 13)
save_plot(paste0("figure_leuk/leuk.eps"), p2, base_height = 7, base_width = 13)
