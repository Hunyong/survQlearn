library(dplyr)

### STEP 2.
### Subsetting the latest visits only Visit 5 (2011 ~ 2012) and after. (Visits 4 were in 1996-1997, visits 1 were 1986-1988)
dat.clean = readRDS("../data_ARIC/01.aric.clean.rds")
all_cause = TRUE

### 1. trimmer(): stages 1:7 to 5:7
trimmer = function(data, full = 1:7, reduced = 4:6, 
                   stage.length = "V", status = "d",
                   sep = ".",
                   total.time = "T.0", last.fu = "last.fu") {
  require(dplyr)
  # reduced should be continuous (no intermittent missing).
  
  # 0. defining the variables
  vars = names(data)
  # reduced = reduced         # 4, 5, 6
  # to.be = 1:length(reduced) # 1, 2, 3
  
  # V.full = paste0(stage.length, sep, full)
  # d.full = paste0(status, sep, full)
  V.new = paste0(stage.length, sep, reduced)
  # d.new = paste0(status, sep, reduced)
  
  # if reduced does not end with the final stage, the last stage length is augmented.
  if (max(reduced) < max(full)) {
    V.residual = paste0(stage.length, sep, max(reduced) : max(full))     # if reduced = 2:4, obtain c("V.4", "V.5", "V.6", "V.7").
    data[, V.residual[1]] = apply(data[, V.residual], 1, sum, na.rm = TRUE)
  }
  
  # 1. updating the total time
  data[, total.time] = apply(data[, V.new], 1, sum, na.rm = TRUE)
  
  new.i = 0
  for (i in full) {
    if (i %in% reduced) {
      new.i = new.i + 1
      vars = gsub(paste0("\\.", i), paste0("\\.NEW", new.i), vars)
    } else {
      vars = gsub(paste0("\\.", i), paste0("\\.DEPRECATED", i), vars)
    }
  }
  names(data) = vars
  
  ## 2. Drop the unnecessary variables, and update the last.fu
  # drop the unused visits
  data = data %>% 
    dplyr::select(-contains("DEPRECATED"))
  # fix the names by dropping "NEW"
  names(data) = gsub("\\.NEW([0-9]+)", ".\\1", names(data))
  
  # updating the last.fu
  data[, last.fu] = pmin(data[, last.fu], max(reduced))
  data[, last.fu] = match(data[, last.fu], reduced)
  
  # updating the statuses
  for (i in 1:length(reduced)) {
    data[, sprintf("%s%s%s", status, sep, i)] = 
      ifelse(data[, last.fu] > i, 1, ifelse(data[, last.fu] == i, data[, "DEAD19"], NA))
  }
  
  # Removing those NA's (who failed or dropped out)
  data = data %>% 
    dplyr::filter(!is.na(last.fu))
  
  return(data)
}

### 2. specific.failure: Treat failure from other causes as censoring.
##     See https://health.mo.gov/data/documentation/death/death-icd10.php for the list of the codes
specific.failure = function(data, cause = "^I.*") {
  fail = grepl(cause, data$UCOD)
  cens = data$UCOD == ""
  cens.other.cause = !(fail | cens)
  data[cens.other.cause, "DEAD19"] = 0
  for (i in 1:max(data$last.fu, na.rm = TRUE)) {
    matched = data[cens.other.cause, "last.fu"] == i
    after.matched = data[cens.other.cause, "last.fu"] < i
    data[cens.other.cause, paste0("d.", i)][matched] = 0         # Change death = 1 to = 0
    data[cens.other.cause, paste0("d.", i)][after.matched] = NA  # Change delta = NA afterwards
  }
  return(data)
}

### dat.567 has only seven observed deaths. So we move the window to 4:6.
### dat.456 still has only seven observed deaths. So we move the window to 5:6.
# dat.567 = dat.clean %>% trimmer(reduced = 5:7) %>% specific.failure %>% filter(PRVCHD05 == 1| PREVHF01 == 1) # 92% censoring
# dat.567 = dat.clean %>% trimmer(reduced = 5:7) %>% specific.failure %>% dplyr::filter(hf.1 == 1) # 81% censoring
# dat.567 = dat.clean %>% trimmer(reduced = 5:7) %>% specific.failure %>% dplyr::filter(hf.1 == 1 | PRVCHD05 == 1| PREVHF01 == 1) # 82% censoring
# dat.456 = dat.clean %>% trimmer(reduced = 4:6) %>% specific.failure %>% dplyr::filter( PRVCHD05 == 1| PREVHF01 == 1) # 76% censoring
dat.56 = dat.clean %>% trimmer(reduced = 5:6) %>% {if (all_cause) . else specific.failure(.)} %>% dplyr::filter(hf.1 == 1 | PRVCHD05 == 1| PREVHF01 == 1) # 82% censoring
# dat.567 %>% filter(DEAD19 == 1) %>% {table(.$last.fu)} # 147, 20, 7 (2011, 2017, 2019) # Too few deaths at the end
# dat.456 %>% filter(DEAD19 == 1) %>% {table(.$last.fu)} # 187, 37, 8 (1996, 2011, 2017) # Too few deaths at the end
dat.56 %>% filter(DEAD19 == 1) %>% {table(.$last.fu)} # 147, 27 (2011, 2017) / 324, 65 for all cause deaths.
# dat.45 %>% filter(DEAD19 == 1) %>% {table(.$last.fu)} # Too old (1996, 2011)

# combineTx = function(vec1, vec2) {ifelse(is.na(vec1) | is.na(vec2), NA, paste0(vec1, vec2))}
combineTx = function(vec1, vec2) {vec1 + vec2}
# dat.456 = dat.456 %>% 
dat.56 = dat.56 %>% 
  dplyr::filter(!race %in% c("A", "I")) %>% mutate(race = race %>% droplevels) %>%   # Removing the only one "A" and no "I"s from analysis.
  rename(hf = hf.1) %>%               # Baseline feature, only available at stage 5.  #hf is not available at visit 4.
  rename(delta = DEAD19)              # Replace DEAD19 with delta
saveRDS(dat.56, sprintf("../data_ARIC/02.aric.56%s.rds", if (all_cause) "_allcause" else "" ))

# 1. raw data:  n=15,760, % non-censoring = 55%
dim(dat.clean) ; mean(dat.clean$DEAD19)

# 2. visits 5-7 only:  n=5,890, % non-censoring = 23%
dim(trimmer(dat.clean, reduced = 5:6)) ; mean(trimmer(dat.clean, reduced = 5:6)$DEAD19)

# 3. circular death only:  n=5,890, % non-censoring = 7.6%
dim(specific.failure(trimmer(dat.clean, reduced = 5:6))) ; mean(specific.failure(trimmer(dat.clean, reduced = 5:6))$DEAD19)

# 3. the cases at visit 1: n = 945,  % non-censoring = 18% / 41% for all cause mortality
dim(dat.56) ; mean(dat.56$delta)
