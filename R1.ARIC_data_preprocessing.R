### 0. library
  library(dplyr)

### 1. Raw data read-in
  visit1 = read.delim("../data_ARIC/csv/visit_1.csv", sep = "\t")
  visit2 = read.delim("../data_ARIC/csv/visit_2.csv", sep = "\t")
  visit3 = read.delim("../data_ARIC/csv/visit_3.csv", sep = "\t")
  visit4 = read.delim("../data_ARIC/csv/visit_4.csv", sep = "\t")
  visit5 = read.delim("../data_ARIC/csv/visit_5.csv", sep = "\t")
  visit6 = read.delim("../data_ARIC/csv/visit_6.csv", sep = "\t")
  visit7 = read.delim("../data_ARIC/csv/visit_7.csv", sep = "\t")
  sta = read.delim("../data_ARIC/csv/status71.csv", sep = "\t")   # status
  inc = read.delim("../data_ARIC/csv/inc_by19.csv", sep = "\t") # incidence

  # Dimension of the raw data
  dim(visit1) # 15,760 x 216
  dim(visit2) # 14,316 x 185
  dim(visit3) # 12,858 x 194
  dim(visit4) # 11,629 x 213
  dim(visit5) #  6,527 x 415
  dim(visit6) #  3,996 x 368
  dim(visit7) #  3,583 x 547   # Missing subjects removed.
  dim(sta)    # 15,760 x 67
  dim(inc)    # 15,760 x 67


### 2. Screening variables
  if (0) {
    # key variables
    key = readxl::read_xlsx("../data_ARIC/key_variables.xlsx")
    key2 = key %>% dplyr::filter(category != "Tx2") # removing the detailed Tx information
    membership = function(x, ref) x[x %in% ref]
    
    # Key variables only. # For a screening purpose only. Can be skipped.
    visit1 = visit1[, membership(names(visit1), key2$variable)]
    visit2 = visit2[, membership(names(visit2), key2$variable)]
    visit3 = visit3[, membership(names(visit3), key2$variable)]
    visit4 = visit4[, membership(names(visit4), key2$variable)]
    visit5 = visit5[, membership(names(visit5), key2$variable)]
    visit6 = visit6[, membership(names(visit6), key2$variable)]
    visit7 = visit7[, membership(names(visit7), key2$variable)]
    sta = sta[, membership(names(sta), key2$variable)]
    inc = inc[, membership(names(inc), key2$variable)]
    
    # Dimension of the raw data
    dim(visit1) # 15,760 x 41
    dim(visit2) # 14,316 x 23
    dim(visit3) # 12,858 x 20
    dim(visit4) # 11,629 x 21
    dim(visit5) #  6,527 x 21
    dim(visit6) #  3,006 x 20
    dim(visit7) #  3,583 x 20   # Missing subjects removed.
    dim(sta) # 15,760 x 21
    dim(inc) # 15,760 x 15
  }
  
### 2B. description of the key variables
  ## (1) inc
  # ID
  # C7_FUTIME : follow-up time for MI or fatal CHD
  # FUTIMED   : follow-up time for death (longer than or equal to C7_FUTIME)
  # DEAD19    : status of death as of 2019-12-31
  # UCOD      : underlying cause of death (coded)
  # C7_SOURCINC : the source of event when C7_INC_BY19=1  ‘MI’ if the event is definite or probable MI
  # C7_SOURCIP : the source of event when C7_IN_BY19P=1   ‘MI’ (‘FATCHD’) if there is no cardiac procedure or cardiac procedure is after the MI.
  # C7_SOURCISP : the source of event when C7_IN_19SP=1   ‘MI’ (‘FATCHD’) if the event is definite or probable MI (definite fatal CHD) 
  # PREVHF01  : prevalent heart failure at visit 1.
  # PRVCHD05  : prevalent coronary heart disease at visit 1.
  # TIAB01    : prevalent stroke at visit 1.
  
  ## (2) sta
  # CENTER    : center four centers
  # KNWNDEADBYVISIT21 : known to be dead at visit 2 => delta_1
  # KNWNDEADBYVISIT31 : known to be dead at visit 3 => delta_2
  # KNWNDEADBYVISIT41 : known to be dead at visit 4 ...
  # KNWNDEADBYVISIT51 : known to be dead at visit 5
  # KNWNDEADBYVISIT61 : known to be dead at visit 6
  # KNWNDEADBYVISIT71 : known to be dead at visit 7 => delta_6

  
  ## (3) visit1
  # ID, CENTER,
  # ASPIRINCODE01   : Aspirin use in the past 2 weeks based on 2004 medication codes
  # STATINCODE01    : Statin use in the past 2 weeks based on 2004 medication codes
  # ANTICOAGCODE01  : Used Anticoagulates (At Visit 1) Last 2 Weeks (0=no, 1=yes) Based On 2004 Med Code
  # PREVHF01, PRVCHD05 (overlaps with inc)
  # RACEGRP, GENDER, V1AGE01
  # V1DATE01,
  # CIGTYR01    : Cigarette years of smoking
  # BMI01       : 
  # WSTHPR01    : Waist-To-Hip Ratio
  # DRNKR01     : Drinker Status (1: current drinker, 2: former, 3: never, 4: unknown)
  # HYPERT04    : Hypertension, definition 4; replaces HYPERT01
  # GLUCOS01
  # 
  # GLEFH01 FAST1202 TCHSIU01 WORK_I02 SPRT_I02 ELEVEL01 ELEVEL02 HYPTMD01 ...
  # OCCUPN01 MENOPS01 INTPLQ01 HYPTMDCODE01 MOMHISTORYSTR DADHISTORYSTR MOMHISTORYCHD DADHISTORYCHD MOMHISTORYDIA DADHISTORYDIA 

  
### 3. Screening variables - Phase II.
  #### Keeping necessary variables only.
  inc =
    inc %>% 
    transmute(ID, T.0 = FUTIMED, DEAD19, UCOD, PREVHF01, PRVCHD05, TIAB01 = factor(TIAB01, levels = c("N", "Y")))
    # "C7_FUTIME", "DEAD19", "C7_SOURCINC"
  
  sta = 
    sta %>% 
    transmute(ID, CENTER = CENTER %>% as.factor, kd.1 = 0, # known not dead at the beginning of stage 1.
              kd.2 = KNWNDEADBYVISIT21, kd.3 = KNWNDEADBYVISIT31, kd.4 = KNWNDEADBYVISIT41, 
              kd.5 = KNWNDEADBYVISIT51, kd.6 = KNWNDEADBYVISIT61, kd.7 = KNWNDEADBYVISIT71)

  visit1 =
    visit1 %>% 
    transmute(ID, gender = GENDER %>% as.factor, race = RACEGRP %>% as.factor, age = V1AGE01,
              cig_yrs = CIGTYR01, 
              date.1 = V1DATE01,
              AA.1 = ASPIRINCODE01 %>% as.factor,
              AS.1 = STATINCODE01 %>% as.factor,
              AC.1 = ANTICOAGCODE01 %>% as.factor,
              bmi.1 = BMI01, wth.1 = WSTHPR01,
              drink.1 = ordered(DRNKR01, levels = c(3, 2, 1), labels = c("never", "former drinker", "current drinker")),
              hypert.1 = HYPERT04,
              glucose.1 = GLUCOS01,
              smoke.1 = CURSMK01,
              hdl.1 = HDL01,
              fast.1 = FAST0802,
              lfast.1 = FAST1202)
  
  visit2 =
    visit2 %>% 
    transmute(ID, 
              date.2 = V2DATE21,
              AA.2 = ASPIRINCODE21 %>% as.factor,
              AS.2 = STATINCODE21 %>% as.factor,
              AC.2 = ANTICOAGCODE21 %>% as.factor,
              bmi.2 = BMI21, wth.2 = WSTHPR21,
              drink.2 = ordered(DRNKR21, levels = c(3, 2, 1), labels = c("never", "former", "current")),
              hypert.2 = HYPERT24,
              glucose.2 = GLUSIU21,
              smoke.2 = CURSMK21,
              hdl.2 = HDL221,
              fast.2 = FAST0823,
              lfast.2 = FAST1223)
  
  visit3 =
    visit3 %>% 
    transmute(ID, 
              date.3 = V3DATE31,
              AA.3 = ASPIRINCODE31 %>% as.factor,
              AS.3 = STATINCODE31 %>% as.factor,
              AC.3 = ANTICOAGCODE31 %>% as.factor,
              bmi.3 = BMI32,
              wth.3 = WSTHPR31,
              drink.3 = ordered(DRNKR31, levels = c(3, 2, 1), labels = c("never", "former", "current")),
              hypert.3 = HYPERT34,
              # glucose.3 = GLUSIU31, # not available
              smoke.3 = CURSMK31,
              hdl.3 = HDLSIU31,
              fast.3 = FAST0834,
              lfast.3 = FAST1234)
  
  visit4 =
    visit4 %>% 
    transmute(ID, 
              date.4 = V4DATE41,
              AA.4 = ASPIRINCODE41 %>% as.factor,
              AS.4 = STATINCODE41 %>% as.factor,
              AC.4 = ANTICOAGCODE41 %>% as.factor,
              bmi.4 = BMI41,
              wth.4 = WSTHPR41,
              drink.4 = ordered(DRNKR41, levels = c(3, 2, 1), labels = c("never", "former", "current")),
              hypert.4 = HYPERT44,
              glucose.4 = GLUSIU41,
              smoke.4 = CURSMK41,
              hdl.4 = HDLSIU41,
              fast.4 = FAST0841,
              lfast.4 = FAST1241)
  
  visit5 =
    visit5 %>% 
    transmute(ID, 
              date.5 = V5DATE51,
              AA.5 = ASPIRINCODE51 %>% as.factor,
              AS.5 = STATINCODE51 %>% as.factor,
              AC.5 = ANTICOAGCODE51 %>% as.factor,
              bmi.5 = BMI51,
              wth.5 = WSTHPR51,
              drink.5 = ordered(DRNKR51, levels = c(3, 2, 1), labels = c("never", "former", "current")),
              hypert.5 = HYPERT54,
              glucose.5 = GLUSIU51,
              smoke.5 = CURSMK52,
              hdl.5 = HDLSIU51,
              fast.5 = FAST0851,
              hf.5 = PREVHF52,
              lfast.5 = FAST1251)
  
  visit6 =
    visit6 %>% 
    transmute(ID, 
              date.6 = V6DATE61,
              AA.6 = ASPIRINCODE61 %>% as.factor,
              AS.6 = STATINCODE61 %>% as.factor,
              AC.6 = ANTICOAGCODE61 %>% as.factor,
              bmi.6 = BMI61,
              wth.6 = WSTHPR61,
              drink.6 = ordered(DRNKR61, levels = c(3, 2, 1), labels = c("never", "former", "current")),
              hypert.6 = HYPERT64,
              glucose.6 = GLUSIU61,
              smoke.6 = CURSMK62,
              hdl.6 = HDLSIU61,
              fast.6 = FAST0861,
              lfast.6 = FAST1261)
  
  visit7 =
    visit7 %>% 
    transmute(ID, 
              date.7 = V7DATE71,
              AA.7 = ASPIRINCODE71 %>% as.factor,
              AS.7 = STATINCODE71 %>% as.factor,
              AC.7 = ANTICOAGCODE71 %>% as.factor,
              bmi.7 = BMI71,
              wth.7 = WSTHPR71,
              drink.7 = ordered(DRNKR71, levels = c(3, 2, 1), labels = c("never", "former", "current")),
              hypert.7 = HYPERT74,
              glucose.7 = GLUSIU71,
              smoke.7 = CURSMK72,
              hdl.7 = HDLSIU71,
              fast.7 = FAST0871,
              lfast.7 = FAST1271)

  
### 4. Stage length calculation
  outcome =
    left_join(inc %>% select(ID, T.0, DEAD19, UCOD), 
              sta %>% select(ID, kd.1, kd.2, kd.3, kd.4, kd.5, kd.6, kd.7), by = "ID") %>% 
    left_join(visit1 %>% select(ID, date.1), by = "ID") %>% 
    left_join(visit2 %>% select(ID, date.2), by = "ID") %>% 
    left_join(visit3 %>% select(ID, date.3), by = "ID") %>% 
    left_join(visit4 %>% select(ID, date.4), by = "ID") %>% 
    left_join(visit5 %>% select(ID, date.5), by = "ID") %>% 
    left_join(visit6 %>% select(ID, date.6), by = "ID") %>% 
    left_join(visit7 %>% select(ID, date.7), by = "ID") 
   
  
  ## 4.1 Obtain the last follow up visits
  visits = function(x) which(!is.na(x))
  # This returns 
  #  1) the last continuous followup visit (no intermittent missing) and 
  #  2) the next follow-up time if there is intermittent missings (for calculation of the censoring time).
  last.continuous.followup = function(x) {
    visits = visits(x)
    max.visit = max(visits)
    if (max.visit == 1) 
      return(c(lastfu = 1, nextfu = 1))
    for (i in 1:max.visit) {
      if (i != visits[i]) 
        return(c(lastfu = i-1, nextfu = visits[i]))
    }
    return(c(lastfu = max.visit, nextfu = max.visit))
  }
  
  last.visits = 
    outcome %>% 
    dplyr::select(date.1, date.2, date.3, date.4, date.5, date.6, date.7) %>% 
    apply(1, last.continuous.followup) %>% 
    t
  outcome$last.fu = last.visits[, 1]
  outcome$next.fu = last.visits[, 2]
  outcome$date.next.fu = 
    outcome %>% 
    dplyr::select(date.1, date.2, date.3, date.4, date.5, date.6, date.7) %>% 
    "["(cbind(1:dim(outcome)[1], outcome$next.fu))
  
  ## 4.2 New dates for intermittent missing subjects to be treated as right censored.
  #      The next earliest available visit is the observed new censoring time.
  #      e.g., If visit dates are 11, 15, NA, NA, 30, 31, NA, NA, NA, 52,
  #            then the final data are recorded as 11, 15, 30+, NA, NA, NA, NA, NA, NA, NA.
  outcome =
    outcome %>% 
    # STEP1: replacing the missing intermediate date with the next follow up date.
    mutate(new.date.1 = date.1,
           new.date.2 = ifelse(last.fu == 1 & next.fu > last.fu, date.next.fu, date.2),
           new.date.3 = ifelse(last.fu == 2 & next.fu > last.fu, date.next.fu, date.3), 
           new.date.4 = ifelse(last.fu == 3 & next.fu > last.fu, date.next.fu, date.4),
           new.date.5 = ifelse(last.fu == 4 & next.fu > last.fu, date.next.fu, date.5),
           new.date.6 = ifelse(last.fu == 5 & next.fu > last.fu, date.next.fu, date.6),
           new.date.7 = ifelse(last.fu == 6 & next.fu > last.fu, date.next.fu, date.7)) %>% 
    # STEP2: removing the later than the last followup dates (for the intermediate missings)
    mutate(new.date.1 = ifelse(1 > last.fu, NA, new.date.1),
           new.date.2 = ifelse(2 > last.fu, NA, new.date.2),
           new.date.3 = ifelse(3 > last.fu, NA, new.date.3),
           new.date.4 = ifelse(4 > last.fu, NA, new.date.4),
           new.date.5 = ifelse(5 > last.fu, NA, new.date.5),
           new.date.6 = ifelse(6 > last.fu, NA, new.date.6),
           new.date.7 = ifelse(7 > last.fu, NA, new.date.7))
  
  
  ## 4.3 The cumulative time calculation
  # na.futime(): Filling the residual follow up times.
  na.futime = function(x, fu) {ifelse(is.na(x), fu, x)}
  outcome =
    outcome %>% 
    mutate(cum.1 = difftime(new.date.2, new.date.1, units = "day") %>% as.numeric,
           cum.2 = difftime(new.date.3, new.date.1, units = "day") %>% as.numeric,
           cum.3 = difftime(new.date.4, new.date.1, units = "day") %>% as.numeric,
           cum.4 = difftime(new.date.5, new.date.1, units = "day") %>% as.numeric,
           cum.5 = difftime(new.date.6, new.date.1, units = "day") %>% as.numeric,
           cum.6 = difftime(new.date.7, new.date.1, units = "day") %>% as.numeric,
           cum.7 = NA) %>% 
    # Cumulative time is completed by filling in the residual times (until event or censoring)
    mutate(cum.1 = na.futime(cum.1, T.0),
           cum.2 = na.futime(cum.2, T.0),
           cum.3 = na.futime(cum.3, T.0),
           cum.4 = na.futime(cum.4, T.0),
           cum.5 = na.futime(cum.5, T.0),
           cum.6 = na.futime(cum.6, T.0),
           cum.7 = na.futime(cum.7, T.0))

  
  ## 4.4 delta (status) and V (stage lengths) calculation for each stage
  
  # lastfu2delta: If censoring or failure happened at stage 3 (=last.fu),
  #    all previous stages should be given 1, all later stages NA.
  #    and the stage should be given by the "dead" status.
  lastfu2delta = function(i, last.fu, dead) {ifelse(i < last.fu, 1, ifelse(i == last.fu, dead, NA))}
  
  outcome =
    outcome %>% 
    # Stagewise death status
    mutate(d.1 = lastfu2delta(1, last.fu, DEAD19),
           d.2 = lastfu2delta(2, last.fu, DEAD19),
           d.3 = lastfu2delta(3, last.fu, DEAD19),
           d.4 = lastfu2delta(4, last.fu, DEAD19),
           d.5 = lastfu2delta(5, last.fu, DEAD19),
           d.6 = lastfu2delta(6, last.fu, DEAD19),
           d.7 = lastfu2delta(7, last.fu, DEAD19)) %>% 
    # stage lengths
    mutate(V.1 = cum.1,
           V.2 = cum.2 - cum.1,
           V.3 = cum.3 - cum.2,
           V.4 = cum.4 - cum.3,
           V.5 = cum.5 - cum.4,
           V.6 = cum.6 - cum.5,
           V.7 = cum.7 - cum.6)
    
  ### Some followup times (FUTIMED)  are shorter than the longest intervals (last visit date - first visit date)
  ##  In this case, the last stage length is negative (FUTIME - longest interval < 0).
  ##  We adjust the second to the last visit time so that it matches the followup time, and remove the latest statuses (d.7, e.g.)
  outcome = 
    outcome %>% 
    # V.7 < 0  n = 583 / 16K
    mutate(V.6 = ifelse(!is.na(V.7) & V.7 < 0, V.6 + V.7, V.6),
           d.6 = ifelse(!is.na(V.7) & V.7 < 0, d.7, d.6),
           last.fu = ifelse(!is.na(V.7) & V.7 < 0, 6, last.fu),
           next.fu = ifelse(!is.na(V.7) & V.7 < 0, 6, next.fu),
           d.7 = ifelse(!is.na(V.7) & V.7 < 0, NA, d.7),
           V.7 = ifelse(!is.na(V.7) & V.7 < 0, 0, V.7)) %>% 
    # For V.6 < 0, n = 1
    mutate(V.5 = ifelse(!is.na(V.6) & V.6 < 0, V.5 + V.6, V.5),
           d.5 = ifelse(!is.na(V.6) & V.6 < 0, d.6, d.5),
           last.fu = ifelse(!is.na(V.6) & V.6 < 0, 5, last.fu),
           next.fu = ifelse(!is.na(V.6) & V.6 < 0, 5, next.fu),
           d.6 = ifelse(!is.na(V.6) & V.6 < 0, NA, d.6),
           V.6 = ifelse(!is.na(V.6) & V.6 < 0, 0, V.6)) %>% 
    # For V.4 < 0, n = 1
    mutate(V.3 = ifelse(!is.na(V.4) & V.4 < 0, V.3 + V.4, V.3),
           d.3 = ifelse(!is.na(V.4) & V.4 < 0, d.4, d.3),
           last.fu = ifelse(!is.na(V.4) & V.4 < 0, 3, last.fu),
           next.fu = ifelse(!is.na(V.4) & V.4 < 0, 3, next.fu),
           d.4 = ifelse(!is.na(V.4) & V.4 < 0, NA, d.4),
           V.4 = ifelse(!is.na(V.4) & V.4 < 0, 0, V.4)) %>% 
    # For V.2 < 0, n = 1
    mutate(V.1 = ifelse(!is.na(V.2) & V.2 < 0, V.1 + V.2, V.1),
           d.1 = ifelse(!is.na(V.2) & V.2 < 0, d.2, d.1),
           last.fu = ifelse(!is.na(V.2) & V.2 < 0, 1, last.fu),
           next.fu = ifelse(!is.na(V.2) & V.2 < 0, 1, next.fu),
           d.2 = ifelse(!is.na(V.2) & V.2 < 0, NA, d.2),
           V.2 = ifelse(!is.na(V.2) & V.2 < 0, 0, V.2))
  
  
    ### (Code to be erased.)
    # # stage-wide deltas: 1 for death and 0 for censoring.
    # mutate(d.1 = kd2delta(kd.1, DEAD19), d.2 = kd2delta(kd.2, DEAD19), 
    #        d.3 = kd2delta(kd.3, DEAD19), d.4 = kd2delta(kd.4, DEAD19), 
    #        d.5 = kd2delta(kd.5, DEAD19), d.6 = kd2delta(kd.6, DEAD19), 
    #        d.7 = kd2delta(kd.7, DEAD19)) %>% 
    # # If censoring happened at stage 3, all previous stages should be given 1 and all later stages NA.
    # mutate(d.1 = ifelse(DEAD19 == 0 & last.fu > 1, 1, d.1),
    #        d.2 = ifelse(DEAD19 == 0 & last.fu > 2, 1, ifelse(DEAD19 == 0 & last.fu < 2, NA, d.2)),
    #        d.3 = ifelse(DEAD19 == 0 & last.fu > 3, 1, ifelse(DEAD19 == 0 & last.fu < 3, NA, d.3)),
    #        d.4 = ifelse(DEAD19 == 0 & last.fu > 4, 1, ifelse(DEAD19 == 0 & last.fu < 4, NA, d.4)),
    #        d.5 = ifelse(DEAD19 == 0 & last.fu > 5, 1, ifelse(DEAD19 == 0 & last.fu < 5, NA, d.5)),
    #        d.6 = ifelse(DEAD19 == 0 & last.fu > 6, 1, ifelse(DEAD19 == 0 & last.fu < 6, NA, d.6)),
    #        d.7 = ifelse(DEAD19 == 0 & last.fu > 7, 1, ifelse(DEAD19 == 0 & last.fu < 7, NA, d.7)))
    # 
  
  ## 9. Final set.
  dat.clean =
    outcome %>% 
    ##    9.1 Removes kd.?, date.?, new.date.?, cum.?, next.fu, date.next.fu.
    dplyr::select(ID, T.0, DEAD19, UCOD, last.fu, d.1:d.7, V.1:V.7) %>% 
    ##    9.2 Adding necessary information back
    left_join(inc %>% select(ID, PREVHF01, PRVCHD05, TIAB01), by = "ID") %>% 
    left_join(sta %>% select(ID, CENTER), by = "ID") %>% 
    left_join(visit1 %>% select(-date.1), by = "ID") %>% 
    left_join(visit2 %>% select(-date.2), by = "ID") %>% 
    left_join(visit3 %>% select(-date.3), by = "ID") %>% 
    left_join(visit4 %>% select(-date.4), by = "ID") %>% 
    left_join(visit5 %>% select(-date.5), by = "ID") %>% 
    left_join(visit6 %>% select(-date.6), by = "ID") %>% 
    left_join(visit7 %>% select(-date.7), by = "ID")
  
  
  saveRDS(dat.clean, "../data_ARIC/01.aric.clean.rds")
