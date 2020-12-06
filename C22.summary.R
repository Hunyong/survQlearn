library(dplyr);library(ggplot2);library(cowplot); library(tidyr)

crit.no = 1

lab.date = "2020-12-05"
files <- list.files(path = "output", pattern = paste0(lab.date, ".*\\.rds"), full.names = TRUE)

method.nm.abc = 
  c("A", "B", "C", "D", "E", "F")
method.nm.simple = 
  c("csk", "gkRF", "gkLM", "dw", "zom", "observed")
method.nm.formal = 
  c("the proposed method", "Goldberg & Kosorok (2012), RF", 
    "Goldberg & Kosorok (2012), linear", "Simoneau et al. (2019)", "zero-order model", 
    "observed policy")

  
  filename = function(lab.date, crit.no) 
    paste0("figure/C22_summary_", gsub("-", "", lab.date), "_crit", crit.no, ".eps")
  fn <- 
    expand.grid(beta = 1:6, prop = 1:2, n = 1:2, crit.no = 1:2) %>% 
    mutate(nm = paste0(beta, "-", prop, "-", n),
           fn = paste0("output/simResult_", lab.date, "_beta", beta, 
                       "_prop", prop, "_n", n, "_crit", crit.no, ".rds"))
    
  result.comb <- 
      lapply(1:dim(fn)[1], function(i) {
        fn.i = fn$fn[i]
        if (file.exists(fn.i)) {
          a <- readRDS(fn.i)$statistics
        } else if (file.exists(gsub("\\.rds", "_tmp.rds", fn.i))) {
          a <- readRDS(gsub("\\.rds", "_tmp.rds", fn.i))
        } else {
          a <- NULL
        }
        if (is.null(a)) {
          NULL
        } else {
          as.data.frame(a) %>% 
            dplyr::select(rep = no, observed, 
                          csk, gkLM, gkRF, dw, zom, percent.censor, starts_with("n.")) %>% 
            mutate(beta = fn$beta[i],
                   prop = fn$prop[i],
                   n    = fn$n[i]   ,
                   crit = fn$crit.no[i]) %>% 
            pivot_longer(names_to = "method", values_to = "value", 
                         c("observed", "csk", "gkLM", "gkRF", "dw", "zom"))
        }
      }) %>% 
      do.call(rbind, .) %>% 
      mutate(
        method = factor(method, levels = method.nm.simple, labels = method.nm.abc),
        setting = factor(beta, levels = c(1,3,4,2), # reordering so that high censoring comes last.
                          labels = c("1" = "moderate censoring rate\n d = 5",
                                     # "2" = "Pr(C)\u2191", 
                                     "3" = "moderate censoring rate\n d = 2", 
                                     "4" = "moderate censoring rate\n d = 10",
                                     "2" = "high censoring rate\n d = 5")), 
        n = factor (n, levels = 1:2, labels = c("1" = "n=300", "2" = "n = 1000")),
        design = factor(prop, levels = 1:2, labels = c("1" = "observational", "2" = "RCT")),
        crit.label = factor(crit, levels = 1:3, labels = c("1" = "Truncated mean, E[T]", 
                                                           "2" = "Survival probability at t=5, S(5)",
                                                           "3" = "Survival probability at t=7, S(7)")))
    
    result.stat <- 
      result.comb %>% 
      aggregate(cbind(value, percent.censor, n.1, n.2, n.3) ~ method + setting + n + design + 
                  crit + crit.label, data = ., 
                FUN = function(x) round(mean(x, na.rm = TRUE), 2)) %>% 
      mutate(progress = paste(n.1, n.2, n.3, sep = " / "))
    result.stat.sd <- 
      result.comb %>% 
      aggregate(value ~ method + setting + n + design + 
                  crit + crit.label, data = ., 
                FUN = function(x) round(sd(x, na.rm = TRUE), 3)) %>% 
      rename(sd = value)
    result.stat <- left_join(result.stat, result.stat.sd)
    rm(result.stat.sd)
    
    p.list <- list()
    for (crit.no in 1:2) {
      file.name.crit = filename(lab.date, crit.no)
      design.filter = c("observational", "RCT")
      ylabs = if (crit.no == 1) {
        "Truncated mean survival time"
      } else {
        paste0("Survival probability at t = ", crit[[crit.no]]$crit.value)
      }
      result.stat.i = 
        result.stat %>% filter(design %in% design.filter, crit %in% crit.no)
      
      p <- 
        result.comb %>% 
        dplyr::filter(design %in% design.filter, crit == crit.no) %>% 
        ggplot(aes(method, value, group = method, color = method)) +
        facet_grid(setting ~ design + n) +
        geom_boxplot() +
        # geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
        geom_jitter(width = 0.1, height = 0) + # EPS does not support alpha.
        scale_color_discrete(labels = paste0(method.nm.abc, ": ", method.nm.formal)) +
        ylab(ylabs) + theme_bw() + theme(legend.position = "bottom")
        # guides(col = guide_legend(nrow = 2, byrow = TRUE, title = NULL))
      rng = suppressWarnings(layer_scales(p)$y$range$range)
      rng[3] = rng[2] - rng[1]
      rng[4] = rng[3] * 0.4 + rng[2] # y coordiate for censoring %
      rng[5] = rng[3] * 0.2 + rng[2] # y coordiate for flowchart
      p.list[[crit.no]] <-
        p + 
        stat_summary(aes(x = as.numeric(method), y = value),
                     fun.y = mean, geom = 'point',  col = "black", shape = "square", size = 2) +
        stat_summary(aes(x = as.numeric(method) + 0.3, y = value*1.1, label = round(..y.., 2)),
                     fun.y = mean, geom = 'text',  col = "black", position = 'dodge')
        
      ggsave(file.name.crit, p.list[[crit.no]], device="eps", width = 10, height = 10)
    } 
    p.grid <- plot_grid(p.list[[1]] + guides(col = FALSE), 
                        p.list[[2]] + guides(col = FALSE), 
                        align = "h", nrow = 1, ncol = 2, common.legend = TRUE)
    
    p.grid <- plot_grid(p.grid, get_legend(p.list[[1]]), align = "v", nrow = 2, rel_heights = c(1, 0.1))
    save_plot("figure/Fig1.eps", p.grid, base_height = 10, base_width = 20)
    