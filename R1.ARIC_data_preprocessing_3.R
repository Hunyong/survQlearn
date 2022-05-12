### STEP 3. 
### Handling the missing values.
library(dplyr)
library(mice)  # multiple imputation
all_cause = TRUE
dat.56 = readRDS(sprintf("../data_ARIC/02.aric.56%s.rds", if (all_cause) "_allcause" else "" ))
m = 5 # Five imputation sets

### Composition of the missing values
  ##  type A. Intrinsic missing values -- For those who were censored or failed at a previous stage, the current stage values are inherently missing.
  ##  type B. Undue missing values -- Lack of information due to any other reasons than type A.
### Approches
  ##  Only type B missings are imputed.
  ##  For the remaining, multiple imputation beginging from baseline, stage 1, 2, and 3.
  ##  The multiple imputation is done for the first stage (baseline), and for each imputed set, the next stage imputation is done once.

### 1. Staging the variables
  names(dat.56)
  var.y = c("ID", "T.0", "delta", "last.fu", "UCOD", "d.1", "V.1", "d.2", "V.2")
  var.0 = c("PREVHF01", "PRVCHD05", "TIAB01", "CENTER", "gender", "race", "age", "cig_yrs")
  var.1 = c("AA.1", "AS.1", "AC.1", "bmi.1", "wth.1", "drink.1", "hypert.1", "glucose.1", "smoke.1", "hdl.1", "fast.1", "lfast.1", "hf")
  var.2 = c("AA.2", "AS.2", "AC.2", "bmi.2", "wth.2", "drink.2", "hypert.2", "glucose.2", "smoke.2", "hdl.2", "fast.2", "lfast.2")
  # var.3 = c("AA.3", "AS.3", "AC.3", "bmi.3", "wth.3", "drink.3", "hypert.3", "glucose.3", "smoke.3", "hdl.3", "fast.3", "lfast.3")
  var.0a = c("ID", var.0); var.1a = c("ID", var.1); var.2a = c("ID", var.2); #var.3a = c("ID", var.3)
  
### 2. Imputation
  ## 2.0 Stage 0 -- outcomes and the baseline variables
  dat.56[!complete.cases(dat.56[, var.y]), var.y] # complete
  
  # All less than 3% missing. PREVHF01 (1.3%), PRVCHD05 (1.6%), TIAB01 (2.2%), cig_yrs (1.0%)
  #  missing fractions
  apply(dat.56[, var.0], 2, function(x) mean(is.na(x)))
  
  set.seed(0)
  imp.0 = mice(dat.56[, var.0a], m = m, block = list("ID", var.0))  # Using block=, "ID" is not used to guess the other variables.
  dat.comp.list0 = 
    lapply(1:m, function(i) complete(imp.0, action = i))
  saveRDS(dat.comp.list0, sprintf("../data_ARIC/03.tmp0%s.rds", if (all_cause) "_allcause" else "" ))
  
  ## 2.1 Stage 1
  #  missing fractions
  apply(dat.56[, var.1], 2, function(x) mean(is.na(x)))
  
  set.seed(1)
  dat.comp.list1 = 
    lapply(1:m, function(i) {
      augmented.i = left_join(dat.comp.list0[[i]], dat.56[, var.1a], by = "ID")
      imp.i = mice(augmented.i, m = 1, block = list("ID", c(var.0, var.1))) # For each previously imputed set, imputation is done once.
      comp.i = complete(imp.i, action = 1)
      comp.i
    })
  saveRDS(dat.comp.list1, sprintf("../data_ARIC/03.tmp1%s.rds", if (all_cause) "_allcause" else "" ))
  
  ## 2.2 Stage 2
  #  missing fractions
  avail2 = !dat.56$d.2 %>% is.na      # Those who are available (no failure or dropout so far)
  apply(dat.56[avail2, var.2], 2, function(x) mean(is.na(x)))
  
  set.seed(2)
  dat.comp.list2 = 
    lapply(1:m, function(i) {
      augmented.i = 
        left_join(dat.comp.list1[[i]], dat.56[, var.2a], by = "ID") %>% 
        dplyr::filter(avail2)    ### Consider only those who were observed!
      imp.i = mice(augmented.i, m = 1, block = list("ID", c(var.0, var.1, var.2))) # For each previously imputed set, imputation is done once.
      comp.i = complete(imp.i, action = 1)
      left_join(dat.comp.list1[[i]], comp.i[, var.2a], by = "ID")
    })
  saveRDS(dat.comp.list2, sprintf("../data_ARIC/03.tmp2%s.rds", if (all_cause) "_allcause" else "" ))
  
  
  
  dat.comp.list =
    dat.comp.list2 %>% 
    lapply(function(i) {
      i %>%
        left_join(dat.56[, var.y], ., by = "ID") %>% 
        mutate(#ACFS.1 = combineTx(AC.1 == 1, lfast.1 == 1) %>% factor, 
               #ACFS.2 = combineTx(AC.2 == 1, lfast.2 == 1) %>% factor, 
               #ACFS8.1 = combineTx(AC.1 == 1, fast.1 == 1) %>% factor, 
               #ACFS8.2 = combineTx(AC.2 == 1, fast.2 == 1) %>% factor,
               AAprev.2 = AA.1,
               ACprev.2 = AC.1,
               ASprev.2 = AS.1,
               fastprev.2 = fast.1,
               lfastprev.2 = lfast.1,
               fast3prev.2 = ifelse(is.na(lfast.1), NA, ifelse(lfast.1, 2, fast.1)) %>% factor,
               ACFS.1 = paste0(AC.1, lfast.1) %>% {ifelse(. == "NANA", NA, .)} %>% factor,
               ACFS.2 = paste0(AC.2, lfast.2) %>% {ifelse(. == "NANA", NA, .)} %>% factor)
      }) 
  saveRDS(dat.comp.list, sprintf("../data_ARIC/dat56/03.aric.comp%s.rds",  if (all_cause) "_allcause" else "" ))
  saveRDS(dat.comp.list, sprintf("../data_ARIC/03.aric.comp%s.rds",  if (all_cause) "_allcause" else "" ))
  
  