library(dplyr); library(tibble)
library(pracma) # for linking NA values in geom_lines
library(ggplot2); library(cowplot)

  # Specify date of the outputs and all_cause (T/F)
  date = "2022-04-06"; all_cause = FALSE
  date = "2022-04-09"; all_cause = TRUE
  
  all_cause_nm = if (all_cause) "_allcause" else ""
  mthd = c("CSK", "GKRF", "GKLM", "DW", "ZOM", "observed")
  labs = c("the proposed method", "Goldberg & Kosorok (2012), RF",
            "Goldberg & Kosorok (2012), linear", "Simoneau et al. (2019)", "zero-order model", "observed")
  cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
  names(cols) = labs
  crit = "mean"
  nodesize = 50
  Tx = "ACAS"
  Tx.full = "Anti-coagulant and Statin"
  clipping = "5%"  # clipping at 5%

  p = list()
  for (crit in c("mean", "surv")) {
    nm.tmp =
      sprintf("output/%s/values_%s_tau_2700%s_imp1_%s.rds", 
              date, if (crit == "mean") "mean_NA_meanSplit" else "surv.mean_2200_logrankSplit", 
              all_cause_nm, date) 
    data.long = 
      nm.tmp %>%
      readRDS() %>% as.data.frame %>% 
      dplyr::select(-cens1, -cens2, -ns.CSK, -ns.GK, -ns.GKLM) %>% 
      mutate(rep = 1:n()) %>% 
      tidyr::gather(key = "method", value = "value", -rep) %>% 
      mutate(method = factor(method, levels = c("CSK", "GKRF", "GKLM","DW", "ZOM", "observed"),
                             labels = labs),
             crit = crit)
    
    p[[crit]] =
      ggplot(data.long) +
      geom_boxplot(aes(method, value, col = method)) +
      geom_jitter(aes(method, value, col = method, group = rep), width = 0.2, alpha = 1) +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      scale_color_manual(values = cols) +
      stat_summary(aes(x = as.numeric(method), y = value),
                   fun.y = mean, geom = 'point',  col = "black", shape = "square", size = 2) +
      stat_summary(aes(x = as.numeric(method) + 0.3, y = value, label = sprintf( paste0("%3.", if (crit == "mean") 0 else 3, "f"), ..y..)),
                   fun.y = mean, geom = 'text',  col = "black", position = 'dodge', vjust = 2) +
      theme_bw() + theme(axis.title.x = element_blank()) +
      ylim(c(if (crit == "mean") { if (all_cause) 2000 else 2200} else { if (all_cause) 0.45 else 0.6}, NA)) +
      ylab(if (crit == "mean") "Truncated mean survival time" else "Six year survival probability") +
      guides(col = "none", group = "none")
    ggsave(nm.tmp %>% gsub("output/.*/", "figure/", .) %>% gsub("\\.rds", "\\.png", .),
           width = 6, height = 4)
  }
  
  p2 <- plot_grid(p[[1]], p[[2]], align = "h", nrow = 1)
  save_plot(paste0("figure/aric", all_cause_nm, "_", date,".eps"), p2, base_height = 3, base_width = 13)
  save_plot(paste0("figure/aric", all_cause_nm, ".eps"), p2, base_height = 3, base_width = 13)
              