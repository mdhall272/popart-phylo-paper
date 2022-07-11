library(tidyverse)
library(glue)
library(lubridate)
library(anesrake)
library(spatstat)
library(janitor)
library(stargazer)
library(dtplyr)
library(magrittr)
library(wesanderson)
library(devtools)
library(grImport)
library(rsvg)
library(scales)
library(showtext)

# Input files used here are:
# 1) IBM outputs (labelled 1 to 6 and 9 to 12 for communities), in folders IBM_results_1 to IBM_results_12
# 2) Recency estimates (PopART_recency.csv)
# 3) Table of participant marital statuses (PopART_martial_status.csv)
# 4) Table of drug resistance results (PopART_resistance.csv))
# 4) Table of other participant demographics (PopART_demographic_data.csv)
# 5) Summary of phyloscanner pair results with phyloscanner scores (PopART_phyloscanner_pairs.csv)
# 6) SVG images from fontawesome.com: person-solid, person-dress-solid, stopwatch-solid, arrow-right-arrow-left-solid), recoloured

# One figure requires a custom ggupset installation to allow weighted bar graphs. This is available on requests.


# Load of IBM output and setup of raking targets ----

## Load of IBM files ----

ibm.results.raw <- map(1:12, function(number){
  
  folder <- glue("IBM_results_{number}")
  print(folder)
  txtno <- number %>% as.character()
  txtno <- ifelse(nchar(txtno) == 1, glue("0{txtno}"), txtno)
  
  
  infile <- list.files(path=folder, recursive = T, full.names = T, pattern = glue("CL{txtno}"))[1]
  
  read_csv(infile)
  
}) %>% bind_rows() 

## Wrangling for raking targets. Warnings here are fine and NAs are overwritten. ----

ibm.results.processed <- ibm.results.raw %>%
  select(Year,CommunityID, any_of(starts_with("Ninc"))) %>% # new infections; people infected during the trial are the target
  filter(Year >= 2013 & Year <=2017) %>% # trial period
  filter(!(CommunityID %in% 7:9)) %>% # not Lusaka
  pivot_longer(3:32, values_to = "count") %>% 
  mutate(gender = map_chr(name, function(nm) str_match(nm, "Ninc_([MF])_.*")[,2]  )) %>%
  mutate(agegp.lower = map_dbl(name, function(nm) str_match(nm, "Ninc_[MF]_(.*)_.*")[,2] %>% as.numeric)) %>%
  mutate(agegp.upper = map_dbl(name, function(nm) str_match(nm, "Ninc_[MF]_.*_(.*)")[,2] %>% as.numeric)) %>%
  select(-name) %>%
  mutate(agegp.upper = replace(agegp.upper, is.na(agegp.upper), 85)) %>%
  mutate(birth.year.lower = Year - agegp.upper - 1, birth.year.upper = Year - agegp.lower - 1) %>%
  mutate(birth.year.estimate = (birth.year.lower+birth.year.upper)/2) %>%
  mutate(birth.year.estimate.group = cut(birth.year.estimate, breaks = c(1925, seq(1960, 1995, by=5), 2020), right = FALSE)) # nonempty bins in the real pairs

## Wrangling for overall HIV+ population. This takes a long time in dplyr but fine in dtplyr (proof that can work!). Warnings OK as above. ----

ibm.results.processed.prevalent <- ibm.results.raw

ibm.results.processed.prevalent <- ibm.results.processed.prevalent %>%
  lazy_dt(immutable = T) %>%
  select(Year,CommunityID, any_of(starts_with("Ntot"))) %>%
  filter(Year >= 2013 & Year <=2017) %>%
  filter(!(CommunityID %in% 7:9)) %>%
  pivot_longer(3:5252, values_to = "count")  %>%
  mutate(gender = map_chr(name, function(nm) str_match(nm, "Ntot_([MF])_.*")[,2]  )) %>%
  mutate(agegp.lower = map_dbl(name, function(nm) str_match(nm, "Ntot_[MF]_(.*)_.*_C.*")[,2] %>% as.numeric)) %>%
  mutate(agegp.upper = map_dbl(name, function(nm) str_match(nm, "Ntot_[MF]_.*_(.*)_C.*")[,2] %>% as.numeric)) %>%
  mutate(agegp.upper = replace(agegp.upper, is.na(agegp.upper), 85)) %>%
  filter(!endsWith(name, "_CD40_SPVL0_ART0")) %>%
  mutate(birth.year.lower = Year - agegp.upper - 1, birth.year.upper = Year - agegp.lower - 1) %>%
  mutate(birth.year.estimate = (birth.year.lower+birth.year.upper)/2) %>%
  mutate(birth.year.estimate.group = cut(birth.year.estimate, breaks = c(1925, seq(1960, 1995, by=5), 2020), right = FALSE)) %>%
  as_tibble()

## Raking targets for YEAR OF BIRTH ----

yob.props<- ibm.results.processed %>%
  group_by(Year, birth.year.estimate.group) %>%
  summarise(count = sum(count)) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(fraction = count/total)

yob.props.mean <- yob.props %>% # Averaged over the years
  group_by(birth.year.estimate.group) %>%
  summarise(fraction = mean(fraction))

## Raking targets for SEX ----

sex.props <- ibm.results.processed %>% 
  group_by(Year, gender) %>%
  summarise(count = sum(count)) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(fraction = count/total)

sex.props.mean <- sex.props %>% 
  group_by(gender) %>%
  summarise(fraction = mean(fraction))

## Raking targets for COMMUNITY ----

community.props <- ibm.results.processed %>% 
  group_by(Year, CommunityID) %>%
  summarise(count = sum(count)) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(fraction = count/total)

community.props.mean <- community.props %>% 
  group_by(CommunityID) %>%
  summarise(fraction = mean(fraction)) %>%
  mutate(community = glue("com{CommunityID}"))

## Raking targets for MARITAL STATUS (directly derived from Zambia Health and Demographic Survey) ----

f.status.all <- list(never.married = 31.2, married.or.cohabiting = 55.9, divorced.or.separated = 10.1, widowed = 2.9)
m.status.all <- list(never.married = 46.0, married.or.cohabiting = 49.8, divorced.or.separated = 3.7, widowed = 0.4)

f.prop.hiv = list(never.married = 0.07, married.or.cohabiting = 0.141, divorced.or.separated = 0.291, widowed = 0.435)
m.prop.hiv = list(never.married = 0.028, married.or.cohabiting = 0.106, divorced.or.separated = 0.226, widowed = 0.372)

f.hiv.pos <- unlist(f.status.all)*unlist(f.prop.hiv)
m.hiv.pos <- unlist(m.status.all)*unlist(m.prop.hiv)

f.ms.hiv.pos <- 0.583*f.hiv.pos/sum(f.hiv.pos)
m.ms.hiv.pos <- 0.417*m.hiv.pos/sum(m.hiv.pos)

ms.hiv.pos = f.ms.hiv.pos + m.ms.hiv.pos

# Pairs import ----

## Participants marital status ----

all.marriage <- read_csv("PopART_marital_status.csv") %>%
  mutate(CONID = as.character(CONID))

## Participants drug resistance ----

drpairs <- read_csv("PopART_resistance.csv") %>%
  mutate(fl.res = map2_chr(firstline, hl_firstline, function(fl, hlfl){
    if(is.na(fl)){
      NA
    } else if(hlfl == "TRUE"){
      "high"
    } else if(fl == "TRUE"){
      "low"
    } else {
      "WT"
    }
    
  }))  %>%
  select(GLOBALID, fl.res) %>%
  mutate(GLOBALID = as.character(GLOBALID)) %>%
  group_by(GLOBALID) %>%
  dplyr::slice(1)

## Recency estimates ----

recency.estimates <- read_csv("PopART_recency.csv") %>%
  mutate(host.id = as.character(host.id))

## Other metadata ----

participant.reference <- read_csv("PopART_demographic_data.csv") %>%
  mutate(CONID = as.character(CONID))

## PAIRS ----

### Preliminary step comparing directionality estimates ----

pairs.for.tsi.calcs <- read_csv("PopART_phyloscanner_pairs.csv") %>% 
  mutate(host.1 = as.character(host.1), host.2 = as.character(host.2)) %>%
  left_join(recency.estimates, by=c("host.1" = "host.id"), suffix = c(".1", ".2")) %>%
  left_join(recency.estimates, by=c("host.2" = "host.id"), suffix = c(".1", ".2")) %>%
  mutate(tsi.norm.sd.1 = (RF_pred_max_linear.1 - RF_pred_min_linear.1)/(2*1.96),
         tsi.norm.sd.2 = (RF_pred_max_linear.2 - RF_pred_min_linear.2)/(2*1.96)) %>%
  left_join(participant.reference, by= c("host.1" = "CONID"), suffix=c(".1", ".2")) %>%
  left_join(participant.reference, by= c("host.2" = "CONID"), suffix=c(".1", ".2")) %>% 
  mutate(estimated.infection.date.1 = decimal_date(sample.date.1) - (RF_pred_linear.1),
         estimated.infection.date.2 = decimal_date(sample.date.2) - (RF_pred_linear.2)) 

### Calculate generation times ----

calc.tsi.direction <- function(x){
  if(is.na(x)) return(NA)
  if(x>0.8) return("12")
  else if (x<0.2) return("21")
  else return("None")
}

pairs.for.tsi.calcs <- pairs.for.tsi.calcs %>% 
  filter(!is.na(community.1) & !is.na(community.2)) %>%
  mutate(generation.time = estimated.infection.date.2 - estimated.infection.date.1) %>%
  mutate(tsi.score = pmap_dbl(list(generation.time, tsi.norm.sd.1, tsi.norm.sd.2), 
                              function(gt, tsn1, tsn2) pnorm(gt, mean = 0, sd = sqrt((tsn1^2)+(tsn2^2))))) %>%
  mutate(tsi.dir = map_chr(tsi.score, calc.tsi.direction)) %>% 
  mutate(tsi.arrow = tsi.dir != "None")

### Consensus direction of transmission ----

consensus.direction.2 <- function(x,y){
  if(x){
    if(is.na(y)) return(NA)
    if(y == "12") return("12")
    if(y == "21") return("Conflict")
    if(y == "None") return("12")
  }
  else{
    if(is.na(y)) return(NA)
    if(y == "12") return("12")
    if(y == "21") return("21")
    if(y == "None") return("None")
  }
}

pairs.for.tsi.calcs <- pairs.for.tsi.calcs %>% mutate(consensus.direction = map2_chr(arrow, tsi.dir, consensus.direction.2)) %>%
  mutate(direction.combination = case_when(arrow & tsi.dir == "12" ~ "Both - Agreement",
                                           arrow & tsi.dir == "21" ~ "Both - Disagreement",
                                           arrow & tsi.dir == "None" ~ "Phylo only",
                                           !arrow & tsi.dir != "None" ~ "Recency only",
                                           TRUE ~ "Neither"))

### Table of all identified pairs, directed or no, opposite-sex or no ----

final.pairs <- pairs.for.tsi.calcs %>% select(host.1, host.2,  cluster, label, consensus.direction, direction.combination) %>%
  mutate(temp.host.1 = ifelse(consensus.direction == "21", host.2, host.1)) %>%
  mutate(temp.host.2 = ifelse(consensus.direction == "21", host.1, host.2)) %>%
  select(-host.1, -host.2) %>%
  rename(host.1 = temp.host.1, host.2 = temp.host.2) %>%
  left_join(recency.estimates, by=c("host.1" = "host.id"), suffix = c(".1", ".2")) %>%
  left_join(recency.estimates, by=c("host.2" = "host.id"), suffix = c(".1", ".2")) %>%
  left_join(participant.reference, by= c("host.1" = "CONID"), suffix=c(".1", ".2")) %>%
  left_join(participant.reference, by= c("host.2" = "CONID"), suffix=c(".1", ".2")) %>%
  mutate(arrow = ifelse(consensus.direction != "None" & consensus.direction != "Conflict", "Directed", "Undirected" )) %>%
  mutate(genders = map2_chr(gender.1, gender.2, function(g1, g2){
    if(is.na(g1) | is.na(g2)){
      NA
    } else if(g1 != g2){
      "Opposite-sex"
    } else {
      "Same-sex"
    }
  }))

### Restrict to directed opposite-sex pairs ----

directed.mf.pairs <- final.pairs %>% filter(arrow == "Directed" & (gender.1 != gender.2)) %>%
  mutate(source.age.at.infection = pmap_dbl(list(dob.1, sample.date.2, RF_pred_linear.2), function(x,y,z){
    decimal_date(y) - z - decimal_date(x)
  } )) %>%
  mutate(recipient.age.at.infection = pmap_dbl(list(dob.2, sample.date.2, RF_pred_linear.2), function(x,y,z){
    decimal_date(y) - z - decimal_date(x)
  } ))  %>%
  mutate(source.age.at.sampling = map2_dbl(dob.1, sample.date.1, function(x,y){
    decimal_date(y) - decimal_date(x)
  } )) %>%
  mutate(recipient.age.at.sampling = map2_dbl(dob.2, sample.date.2, function(x,y){
    decimal_date(y) - decimal_date(x)
  } ))  %>%
  mutate(sampling.age.gap = pmap_dbl(list(gender.1, source.age.at.sampling, recipient.age.at.sampling), function(x,y,z){
    if(x=="M"){
      y - z
    } else {
      z - y
    }
  })) %>%
  mutate(fancy.age.gap = pmap_dbl(list(gender.1, source.age.at.infection, recipient.age.at.infection), function(x,y,z){
    if(x=="M"){
      y - z
    } else {
      z - y
    }
  })) %>%
  mutate(origin = map2_chr(origin.1, origin.2, function(x,y){
    if(x==y){
      x
    } else {
      "Mixed"
    }
  })) %>%
  mutate(male.age.at.infection = pmap_dbl(list(gender.1, source.age.at.infection, recipient.age.at.infection), function(x,y,z){
    ifelse(x=="M", y, z)
  })) %>%
  mutate(female.age.at.infection = pmap_dbl(list(gender.1, source.age.at.infection, recipient.age.at.infection), function(x,y,z){
    ifelse(x=="F", y, z)
  })) %>%
  mutate(source.agegp = cut_interval(source.age.at.infection, length = 5, right = F)) %>%
  mutate(source.sampling.agegp = cut_interval(source.age.at.sampling, length = 5, right = F)) %>%
  mutate(score.pair = map_dbl(label, function(sc){
    if(!is.na(as.numeric(sc))){
      as.numeric(sc)
    } else {
      as.numeric(str_split_fixed(sc, "/", 2)[2])
    }
    
  })) %>%
  group_by(host.2) %>% arrange(desc(score.pair)) %>% dplyr::slice(1) %>% ungroup()


arm.mapper <- function(community){
  case_when(community %in% c("02", "05", "10") ~ "A",
            community %in% c("01", "06", "11") ~ "B",
            TRUE ~ "C")
}

### Make master pairs table ----

pairs <- directed.mf.pairs %>%
  mutate(fact.com.1 = as.character(fact.com.1),  fact.com.2 =  as.character(fact.com.2)) %>%
  mutate(host.1 = as.character(host.1), host.2 = as.character(host.2)) %>%
  mutate(source.infection.date = pmap_dbl(list(sample.date.1, RF_pred_linear.1), function(y,z){
    decimal_date(y) - z
  } )) %>%
  mutate(recipient.infection.date = pmap_dbl(list(sample.date.2, RF_pred_linear.2), function(y,z){
    decimal_date(y) - z
  } )) %>%
  mutate(source.age.at.infection = pmap_dbl(list(dob.1, sample.date.2, RF_pred_linear.2), function(x,y,z){
    decimal_date(y) - z - decimal_date(x)
  } )) %>%
  mutate(recipient.age.at.infection = pmap_dbl(list(dob.2, sample.date.2, RF_pred_linear.2), function(x,y,z){
    decimal_date(y) - z - decimal_date(x)
  } ))  %>%
  mutate(infection.date = map2_dbl(sample.date.2, RF_pred_linear.2, function(x,y){
    decimal_date(x) - y 
  } )) %>%
  mutate(est.tsi.source = floor(RF_pred_linear.1)) %>%
  mutate(est.tsi.recipient = floor(RF_pred_linear.2)) %>%
  mutate(recent = pmap_lgl(list(infection.date, RF_pred_linear.1,sample.date.1), function(id, rfpl, sd2){
    decimal_date(sd2) - rfpl >= id - 1
  })) %>%
  mutate(age.gap = source.age.at.infection- recipient.age.at.infection) %>%
  mutate(gender.age.gap = map2_dbl(age.gap, gender.2, function(ag, rs){
    ifelse(rs == "F", ag, -ag)
  })) %>%
  mutate(source.age.at.infection.hr = source.age.at.infection >= 25 & source.age.at.infection <= 40) %>%
  mutate(source.age.at.infection.hr.2 = source.age.at.infection >= 20 & source.age.at.infection <= 35) %>%
  mutate(source.age.at.infection.hr.3 = source.age.at.infection <25) %>%
  filter(infection.date >= 2013.833) %>%
  mutate(recipient.id = host.2) %>%
  mutate(source.id = host.1) %>%
  mutate(source.yob = year(dob.1)) %>%
  mutate(source.yob.group = cut(source.yob, breaks = c(1925, seq(1960, 1995, by=5), 2020), right = FALSE)) %>%
  mutate(recipient.yob = year(dob.2)) %>%
  mutate(recipient.yob.group = cut(recipient.yob, breaks = c(1925, seq(1960, 1995, by=5), 2020), right = FALSE)) %>%
  mutate(age.gap.group = cut(gender.age.gap, breaks = seq(-50,50, by=5))) %>%
  mutate(recipient.yob = year(dob.2)) %>%
  mutate(recipient.sex = gender.2) %>%
  mutate(recipient.community = glue("com{as.character(as.numeric(fact.com.2))}")) %>%
  mutate(recipient.age.at.infection.group = cut(recipient.age.at.infection, breaks = c(0, seq(15, 80, by = 5)), right = FALSE)) %>%
  mutate(source.age.at.infection.group = cut(source.age.at.infection, breaks = c(0, seq(15, 80, by = 5)), right = FALSE)) %>%
  mutate(source.community = glue("com{as.character(as.numeric(fact.com.1))}")) %>%
  mutate(source.different.community = source.community != recipient.community) %>%
  left_join(all.marriage, by= c("host.2" = "CONID")) %>%
  rename(recipient.marital.status = married) %>%
  filter(!is.na(recipient.marital.status)) %>%
  left_join(all.marriage, by= c("host.1" = "CONID")) %>%
  rename(source.marital.status = married) %>%
  left_join(drpairs, by= c("source.id" = "GLOBALID"), suffix = c(".source", ".recipient")) %>%
  left_join(drpairs, by= c("recipient.id" = "GLOBALID"), suffix = c(".source", ".recipient")) %>%
  mutate(male.infector = recipient.sex == "F") %>%
  mutate(high.res = fl.res.source == "high") %>%
  mutate(recipient.arm = map_chr(fact.com.2, arm.mapper)) %>%
  mutate(source.arm = map_chr(fact.com.1, arm.mapper)) %>%
  mutate(source.sample.year = year(sample.date.1)) %>%
  mutate(recipient.sample.year = year(sample.date.2)) %>%
  mutate(direction = case_when(recipient.sex == "M" ~ "Female to male",
                               TRUE ~ "Male to female"))


pairs <- pairs %>% 
  mutate(recipient.yob.group = fct_drop(recipient.yob.group)) %>%
  mutate(source.yob.group = fct_drop(source.yob.group)) %>%
  mutate(recipient.marital.status = case_when(recipient.marital.status == "Never married"~"never.married",
                                              recipient.marital.status == "Married/living as married"~"married.or.cohabiting",
                                              recipient.marital.status == "Divorced/separated"~"divorced.or.separated",
                                              recipient.marital.status == "Widowed"~"widowed")) %>%
  mutate(source.marital.status = case_when(source.marital.status == "Never married"~"never.married",
                                           source.marital.status == "Married/living as married"~"married.or.cohabiting",
                                           source.marital.status == "Divorced/separated"~"divorced.or.separated",
                                           source.marital.status == "Widowed"~"widowed")) %>%
  mutate(dec.dob.1 = decimal_date(dob.1)) %>%
  mutate(age.2017.1 = 2017.5 - dec.dob.1) %>%
  mutate(agegp.2017 = cut(age.2017.1, c(13, seq(20, 50, by=5), 90), right = F )) %>%
  mutate(agegp.2017 = as.character(agegp.2017)) %>%
  mutate(agegp.2017.midpoint = map_dbl(agegp.2017, function(ag2017){
    tmp <- as.character(ag2017)
    tmp2 <- substr(tmp, 2, nchar(tmp)-1)
    str_split_fixed(tmp2, ",", 2)[1,] %>% as.numeric() %>% mean()
  })) %>%
  mutate(agegp.2017 = map_chr(agegp.2017, function(ag17){
    tmp <- substr(ag17, 2, nchar(ag17) - 1)
    tmp <- str_split_fixed(tmp, ",", 2)[1,] %>% as.numeric()
    glue("{tmp[1]}-{tmp[2] - 1}")
  })) %>%
  mutate(agegp.2017 = replace(agegp.2017, agegp.2017 == "50-89", "50+"))

pairs <- pairs %>%
  mutate(high.risk.male.infector = male.infector & source.age.at.infection.hr) %>%
  mutate(high.risk.female.infector = !male.infector & source.age.at.infection.hr.2) %>%
  mutate(cat.1 = case_when(high.risk.female.infector ~ "hr.f",
                           high.risk.male.infector ~ "hr.m",
                           TRUE ~ "lr")) %>%
  mutate(fl.res.source.2 = case_when(is.na(fl.res.source) ~ "Unknown",
                                     fl.res.source=="high" ~ "yes",
                                     fl.res.source=="low" ~ "yes",
                                     TRUE ~ "no"))

## Set up the rake ----

pairs$recipient.sex <- pairs$recipient.sex %>% factor
pairs$source.marital.status <- pairs$source.marital.status %>% factor
pairs$recipient.marital.status <- pairs$recipient.marital.status %>% factor
pairs$recipient.community <- pairs$recipient.community %>% factor

sex.targets <- sex.props.mean %>% pull(fraction)
names(sex.targets) <- c("F", "M")

community.targets <- community.props.mean %>% pull(fraction)
names(community.targets) <- as.character(community.props.mean %>% pull(community))

yob.props.mean.2 <- yob.props.mean %>% filter( birth.year.estimate.group %in% unique(pairs$recipient.yob.group)) %>% mutate( birth.year.estimate.group = fct_drop( birth.year.estimate.group))
yob.targets <- yob.props.mean.2 %>% pull(fraction) 
names(yob.targets) <- yob.props.mean.2 %>% pull( birth.year.estimate.group) 

married.targets <- ms.hiv.pos

targets <- list(sex.targets, community.targets, yob.targets, married.targets)
names(targets) <- c("recipient.sex", "recipient.community", "recipient.yob.group", "recipient.marital.status")

## Rake ----

raked <- anesrake(targets, pairs %>% as.data.frame, caseid = pairs$recipient.id, verbose = T)
pairs$weight <- raked$weightvec

# Bootstrapping analyses (figs 3B/3C) ----

## Point estimates ----

ibm.results.processed.prevalent.2017 <- ibm.results.processed.prevalent %>% 
  mutate(art.code = map_dbl(name,function(nm){
    str_match(nm, ".*_ART(\\-?[0-9]+)")[,2] %>% as.numeric()
  })) %>%
  filter(Year == 2017) %>%
  mutate(agegp.midpoint = (agegp.lower + agegp.upper + 1)/2) %>%
  mutate(agegp.2017 = cut(agegp.midpoint, c(13, seq(20, 50, by=5), 90), right = F )) %>%
  mutate(agegp.2017.midpoint = map_dbl(agegp.2017, function(ag2017){
    tmp <- as.character(ag2017)
    tmp2 <- substr(tmp, 2, nchar(tmp)-1)
    str_split_fixed(tmp2, ",", 2)[1,] %>% as.numeric() %>% mean()
  })) %>%
  mutate(agegp.2017 = as.character(agegp.2017)) %>%
  mutate(agegp.2017 = map_chr(agegp.2017, function(ag17){
    tmp <- substr(ag17, 2, nchar(ag17) - 1)
    tmp <- str_split_fixed(tmp, ",", 2)[1,] %>% as.numeric()
    glue("{tmp[1]}-{tmp[2] - 1}")
  })) %>% 
  mutate(agegp.2017 = replace(agegp.2017, agegp.2017 == "50-89", "50+"))


ibm.results.processed.prevalent.2017.noart <- ibm.results.processed.prevalent.2017 %>% 
  filter(art.code %in% c(0, 4, 5, 6))


ibm.results.sex.2017.summarised <- ibm.results.processed.prevalent.2017 %>%
  filter(count > 0) %>%
  mutate(total = sum(count)) %>%
  mutate(prop = count/total) %>%
  ungroup() %>%
  group_by(gender, total) %>% 
  dplyr::summarise(pop = sum(count), prop = sum(prop))

ibm.results.sex.2017.noart.summarised <- ibm.results.processed.prevalent.2017.noart %>%
  filter(count > 0) %>%
  mutate(total = sum(count)) %>%
  mutate(prop = count/total) %>%
  ungroup() %>%
  group_by(gender, total) %>% 
  dplyr::summarise(pop = sum(count), prop = sum(prop))

ibm.results.agesex.2017.summarised <- ibm.results.processed.prevalent.2017 %>%
  filter(count > 0) %>%
  mutate(total = sum(count)) %>%
  mutate(prop = count/total) %>%
  ungroup() %>%
  group_by(gender, agegp.2017, total) %>% 
  dplyr::summarise(pop = sum(count), prop = sum(prop))

ibm.results.agesex.2017.noart.summarised <- ibm.results.processed.prevalent.2017.noart %>%
  filter(count > 0) %>%
  mutate(total = sum(count)) %>%
  mutate(prop = count/total) %>%
  ungroup() %>%
  group_by(gender, agegp.2017, total) %>% 
  dplyr::summarise(pop = sum(count), prop = sum(prop))

pairs.2017.summarised.agesex <- pairs %>%
  group_by(agegp.2017, gender.1) %>%
  summarise(pop = n(), weightpop = sum(weight))

pairs.2017.summarised.agesex <- pairs.2017.summarised.agesex %>%
  mutate(prop = pop/sum(pairs.2017.summarised.agesex$pop), weightprop = weightpop/sum(pairs.2017.summarised.agesex$weightpop)) %>%
  rename(gender = gender.1) 

pairs.2017.summarised.sex <- pairs %>%
  group_by(gender.1) %>%
  summarise(pop = n(), weightpop = sum(weight))

pairs.2017.summarised.sex <- pairs.2017.summarised.sex %>%
  mutate(prop = pop/sum(pairs.2017.summarised.sex$pop), weightprop = weightpop/sum(pairs.2017.summarised.sex$weightpop)) %>%
  rename(gender = gender.1)

### Sex only, ART included ----

rho.tbl.sex <- ibm.results.sex.2017.summarised %>%
  left_join(pairs.2017.summarised.sex, by= c("gender"), suffix = c(".all", ".infectors")) %>%
  mutate(pop.infectors = replace(pop.infectors, is.na(pop.infectors), 0)) %>%
  mutate(prop.infectors = replace(prop.infectors, is.na(prop.infectors), 0)) %>%
  mutate(weightpop = replace(weightpop, is.na(weightpop), 0)) %>%
  mutate(weightprop = replace(weightprop, is.na(weightprop), 0)) %>%
  mutate(rho = weightprop/prop.all)

### Age and sex, ART included ----

rho.tbl.agesex <- ibm.results.agesex.2017.summarised %>%
  left_join(pairs.2017.summarised.agesex, by= c("agegp.2017", "gender"), suffix = c(".all", ".infectors")) %>%
  mutate(pop.infectors = replace(pop.infectors, is.na(pop.infectors), 0)) %>%
  mutate(prop.infectors = replace(prop.infectors, is.na(prop.infectors), 0)) %>%
  mutate(weightpop = replace(weightpop, is.na(weightpop), 0)) %>%
  mutate(weightprop = replace(weightprop, is.na(weightprop), 0)) %>%
  mutate(rounded.weightpop = round(weightpop)) %>%
  mutate(rho = weightprop/prop.all)

### Sex only, no ART ----

rho.tbl.sex.noart <- ibm.results.sex.2017.noart.summarised %>%
  left_join(pairs.2017.summarised.sex, by= c("gender"), suffix = c(".all", ".infectors")) %>%
  mutate(pop.infectors = replace(pop.infectors, is.na(pop.infectors), 0)) %>%
  mutate(prop.infectors = replace(prop.infectors, is.na(prop.infectors), 0)) %>%
  mutate(weightpop = replace(weightpop, is.na(weightpop), 0)) %>%
  mutate(weightprop = replace(weightprop, is.na(weightprop), 0)) %>%
  mutate(rho = weightprop/prop.all)

### Age and sex, no ART ----

rho.table.agesex.noart <- ibm.results.agesex.2017.noart.summarised %>%
  ungroup() %>%
  left_join(pairs.2017.summarised.agesex, by= c("agegp.2017", "gender"), suffix = c(".all", ".infectors")) %>%
  mutate(pop.infectors = replace(pop.infectors, is.na(pop.infectors), 0)) %>%
  mutate(prop.infectors = replace(prop.infectors, is.na(prop.infectors), 0)) %>%
  mutate(weightpop = replace(weightpop, is.na(weightpop), 0)) %>%
  mutate(weightprop = replace(weightprop, is.na(weightprop), 0)) %>%
  mutate(rounded.weightpop = round(weightpop)) %>%
  mutate(rho = weightprop/prop.all)

## Bootstrapping ----

### Sex only, ART included ----

bootstrapped.results.sex <- map(1:200, function(i){
  bs.pairs <- pairs[sample(1:nrow(pairs), replace = T),] %>%
    mutate(rn = as.character(row_number()))
  
  bs.targets <- targets
  
  if(any(bs.pairs$recipient.yob.group %>%table() == 0)){
    empty.categories <- names(which(bs.pairs$recipient.yob.group %>% table() == 0))
    bs.pairs <- bs.pairs %>% mutate(recipient.yob.group = map_chr(recipient.yob.group, function(ryg){
      tmp <- as.character(ryg)
      
      case_when(tmp == "[1925,1960)"~ "[1925,1965)",
                tmp == "[1960,1965)"~ "[1925,1965)",
                TRUE ~ tmp)
    }))
    
    
    bs.targets$recipient.yob.group[10] <- bs.targets$recipient.yob.group[1] + bs.targets$recipient.yob.group[2]
    bs.targets$recipient.yob.group <- bs.targets$recipient.yob.group[3:10]
    names(bs.targets$recipient.yob.group)[8] <- "[1925,1965)"
    bs.pairs <- bs.pairs %>% mutate(recipient.yob.group = factor(recipient.yob.group))
  }
  
  raked <- anesrake(bs.targets, bs.pairs %>% as.data.frame, caseid = bs.pairs$rn, verbose = F)
  
  caseweights <- data.frame(cases=raked$caseid, weights=raked$weightvec)
  
  bs.pairs$weight <- raked$weightvec
  bs.pairs <- bs.pairs %>%
    mutate(direction = case_when(recipient.sex == "M" ~ "Female to male",
                                 TRUE ~ "Male to female"))
  
  bs.pairs.2017.summarised <- bs.pairs %>%
    group_by(gender.1) %>%
    summarise(pop = n(), weightpop = sum(weight))
  
  bs.pairs.2017.summarised <- bs.pairs.2017.summarised %>%
    mutate(prop = pop/sum(bs.pairs.2017.summarised$pop), weightprop = weightpop/sum(bs.pairs.2017.summarised$weightpop)) %>%
    rename(gender = gender.1) 
  
  rho.tbl <- ibm.results.sex.2017.summarised %>%
    left_join(bs.pairs.2017.summarised, by= c("gender"), suffix = c(".all", ".infectors")) %>%
    mutate(pop.infectors = replace(pop.infectors, is.na(pop.infectors), 0)) %>%
    mutate(prop.infectors = replace(prop.infectors, is.na(prop.infectors), 0)) %>%
    mutate(weightpop = replace(weightpop, is.na(weightpop), 0)) %>%
    mutate(weightprop = replace(weightprop, is.na(weightprop), 0)) %>%
    mutate(rounded.weightpop = round(weightpop)) %>%
    mutate(rho = weightprop/prop.all)
  
  rho.tbl %>% ungroup %>%
    select(gender, rho) %>%
    pivot_wider(names_from = "gender", values_from = "rho") %>%
    rename(Male = M, Female = "F") %>%
    mutate(ratio = Male/Female) %>%
    mutate(bs.number = i)
  
  
}) %>% bind_rows()

### Age and sex, ART included ----

bootstrapped.results.agesex <- map(1:200, function(i){
  
  bs.pairs <- pairs[sample(1:nrow(pairs), replace = T),] %>%
    mutate(rn = as.character(row_number()))
  
  bs.targets <- targets
  
  if(any(bs.pairs$recipient.yob.group %>%table() == 0)){
    empty.categories <- names(which(bs.pairs$recipient.yob.group %>% table() == 0))
    bs.pairs <- bs.pairs %>% mutate(recipient.yob.group = map_chr(recipient.yob.group, function(ryg){
      tmp <- as.character(ryg)
      
      case_when(tmp == "[1925,1960)"~ "[1925,1965)",
                tmp == "[1960,1965)"~ "[1925,1965)",
                TRUE ~ tmp)
    }))
    
    
    bs.targets$recipient.yob.group[10] <- bs.targets$recipient.yob.group[1] + bs.targets$recipient.yob.group[2]
    bs.targets$recipient.yob.group <- bs.targets$recipient.yob.group[3:10]
    names(bs.targets$recipient.yob.group)[8] <- "[1925,1965)"
    bs.pairs <- bs.pairs %>% mutate(recipient.yob.group = factor(recipient.yob.group))
  }

  
  raked <- anesrake(bs.targets, bs.pairs %>% as.data.frame, caseid = bs.pairs$rn, 
                    verbose = F
  )
  
  caseweights <- data.frame(cases=raked$caseid, weights=raked$weightvec)
  
  bs.pairs$weight <- raked$weightvec
  bs.pairs <- bs.pairs %>%
    mutate(direction = case_when(recipient.sex == "M" ~ "Female to male",
                                 TRUE ~ "Male to female"))
  
  bs.pairs.2017.summarised <- bs.pairs %>%
    group_by(agegp.2017, gender.1) %>%
    summarise(pop = n(), weightpop = sum(weight))
  
  bs.pairs.2017.summarised <- bs.pairs.2017.summarised %>%
    mutate(prop = pop/sum(bs.pairs.2017.summarised$pop), weightprop = weightpop/sum(bs.pairs.2017.summarised$weightpop)) %>%
    rename(gender = gender.1) 
  
  rho.tbl.agesex <- ibm.results.agesex.2017.summarised %>%
    left_join(bs.pairs.2017.summarised, by= c("agegp.2017", "gender"), suffix = c(".all", ".infectors")) %>%
    mutate(pop.infectors = replace(pop.infectors, is.na(pop.infectors), 0)) %>%
    mutate(prop.infectors = replace(prop.infectors, is.na(prop.infectors), 0)) %>%
    mutate(weightpop = replace(weightpop, is.na(weightpop), 0)) %>%
    mutate(weightprop = replace(weightprop, is.na(weightprop), 0)) %>%
    mutate(rounded.weightpop = round(weightpop)) %>%
    mutate(rho = weightprop/prop.all)
  
  rho.tbl.agesex %>% ungroup %>%
    select(agegp.2017, gender, rho) %>%
    pivot_wider(names_from = "gender", values_from = "rho") %>%
    rename(Male = M, Female = "F") %>%
    mutate(ratio = Male/Female) %>%
    mutate(bs.number = i)
  
}) %>% bind_rows()


rho.ratio.agesex <- rho.tbl.agesex %>% ungroup() %>%
  select(agegp.2017, gender, rho) %>%
  pivot_wider(names_from = "gender", values_from = "rho") %>%
  rename(Male = M, Female = "F") %>%
  mutate(ratio = Male/Female) %>%
  mutate(ratio.bs.ci = map(agegp.2017, function(agmp){
    (bootstrapped.results.agesex %>% filter(agegp.2017 == agmp) %>% pull(ratio) %>% sort())[c(5,195)]
  })) %>%
  mutate(ratio.lower.ci = map_dbl(ratio.bs.ci, extract(1)), ratio.upper.ci = map_dbl(ratio.bs.ci, extract(2))) %>%
  mutate(male.bs.ci = map(agegp.2017, function(agmp){
    (bootstrapped.results.agesex %>% filter(agegp.2017 == agmp) %>% pull(Male) %>% sort())[c(5,195)]
  })) %>%
  mutate(male.lower.ci = map_dbl(male.bs.ci, extract(1)), male.upper.ci = map_dbl(male.bs.ci, extract(2))) %>%
  mutate(female.bs.ci = map(agegp.2017, function(agmp){
    (bootstrapped.results.agesex %>% filter(agegp.2017 == agmp) %>% pull(Female) %>% sort())[c(5,195)]
  })) %>%
  mutate(female.lower.ci = map_dbl(female.bs.ci, extract(1)), female.upper.ci = map_dbl(female.bs.ci, extract(2))) %>%
  mutate(denominator = "As proportion of HIV+ individuals")

rho.ratio.pt1 <- rho.ratio.agesex %>% 
  select(agegp.2017, Female, female.lower.ci, female.upper.ci, Male, male.lower.ci, male.upper.ci) %>%
  rename(female.estimate = Female, male.estimate = Male) %>%
  pivot_longer(2:7) %>%
  mutate(sex = map_chr(name, function(nm) str_split_fixed(nm, "\\.", 2)[,1])) %>%
  mutate(val = map_chr(name, function(nm) str_split_fixed(nm, "\\.", 2)[,2])) %>%
  select(-name) %>%
  pivot_wider(names_from = val, values_from = value) %>%
  mutate(denominator = "As proportion of HIV+ individuals")

### Sex only, no ART ----

bootstrapped.results.sex.art <- map(1:200, function(i){
  bs.pairs <- pairs[sample(1:nrow(pairs), replace = T),] %>%
    mutate(rn = as.character(row_number()))
  
  bs.targets <- targets
  
  if(any(bs.pairs$recipient.yob.group %>% table() == 0)){
    empty.categories <- names(which(bs.pairs$recipient.yob.group %>% table() == 0))
    bs.pairs <- bs.pairs %>% mutate(recipient.yob.group = map_chr(recipient.yob.group, function(ryg){
      tmp <- as.character(ryg)
      
      case_when(tmp == "[1925,1960)"~ "[1925,1965)",
                tmp == "[1960,1965)"~ "[1925,1965)",
                TRUE ~ tmp)
    }))
    
    
    bs.targets$recipient.yob.group[10] <- bs.targets$recipient.yob.group[1] + bs.targets$recipient.yob.group[2]
    bs.targets$recipient.yob.group <- bs.targets$recipient.yob.group[3:10]
    names(bs.targets$recipient.yob.group)[8] <- "[1925,1965)"
    bs.pairs <- bs.pairs %>% mutate(recipient.yob.group = factor(recipient.yob.group))
  }
  
  raked <- anesrake(bs.targets, bs.pairs %>% as.data.frame, caseid = bs.pairs$rn, verbose = F)
  
  caseweights <- data.frame(cases=raked$caseid, weights=raked$weightvec)
  
  bs.pairs$weight <- raked$weightvec
  bs.pairs <- bs.pairs %>%
    mutate(direction = case_when(recipient.sex == "M" ~ "Female to male",
                                 TRUE ~ "Male to female"))
  
  bs.pairs.2017.summarised <- bs.pairs %>%
    group_by(gender.1) %>%
    summarise(pop = n(), weightpop = sum(weight))
  
  bs.pairs.2017.summarised <- bs.pairs.2017.summarised %>%
    mutate(prop = pop/sum(bs.pairs.2017.summarised$pop), weightprop = weightpop/sum(bs.pairs.2017.summarised$weightpop)) %>%
    rename(gender = gender.1) 
  
  rho.tbl <- ibm.results.sex.2017.noart.summarised %>%
    left_join(bs.pairs.2017.summarised, by= c("gender"), suffix = c(".all", ".infectors")) %>%
    mutate(pop.infectors = replace(pop.infectors, is.na(pop.infectors), 0)) %>%
    mutate(prop.infectors = replace(prop.infectors, is.na(prop.infectors), 0)) %>%
    mutate(weightpop = replace(weightpop, is.na(weightpop), 0)) %>%
    mutate(weightprop = replace(weightprop, is.na(weightprop), 0)) %>%
    mutate(rounded.weightpop = round(weightpop)) %>%
    mutate(rho = weightprop/prop.all)
  
  rho.tbl %>% ungroup %>%
    select(gender, rho) %>%
    pivot_wider(names_from = "gender", values_from = "rho") %>%
    rename(Male = M, Female = "F") %>%
    mutate(ratio = Male/Female) %>%
    mutate(bs.number = i)
  
  
}) %>% bind_rows()

### Age and sex, no ART ----

bootstrapped.results.agesex.noart <- map(1:200, function(i){
  bs.pairs <- pairs[sample(1:nrow(pairs), replace = T),] %>%
    mutate(rn = as.character(row_number()))
  
  bs.targets <- targets
  
  if(any(bs.pairs$recipient.yob.group %>%table() == 0)){
    empty.categories <- names(which(bs.pairs$recipient.yob.group %>% table() == 0))
    bs.pairs <- bs.pairs %>% mutate(recipient.yob.group = map_chr(recipient.yob.group, function(ryg){
      tmp <- as.character(ryg)
      
      case_when(tmp == "[1925,1960)"~ "[1925,1965)",
                tmp == "[1960,1965)"~ "[1925,1965)",
                TRUE ~ tmp)
    }))
    
    
    bs.targets$recipient.yob.group[10] <- bs.targets$recipient.yob.group[1] + bs.targets$recipient.yob.group[2]
    bs.targets$recipient.yob.group <- bs.targets$recipient.yob.group[3:10]
    names(bs.targets$recipient.yob.group)[8] <- "[1925,1965)"
    bs.pairs <- bs.pairs %>% mutate(recipient.yob.group = factor(recipient.yob.group))
  }
  
  raked <- anesrake(bs.targets, bs.pairs %>% as.data.frame, caseid = bs.pairs$rn, 
                    verbose = F
  )
  
  caseweights <- data.frame(cases=raked$caseid, weights=raked$weightvec)
  
  bs.pairs$weight <- raked$weightvec
  bs.pairs <- bs.pairs %>%
    mutate(direction = case_when(recipient.sex == "M" ~ "Female to male",
                                 TRUE ~ "Male to female"))
  
  bs.pairs.2017.summarised <- bs.pairs %>%
    group_by(agegp.2017, gender.1) %>%
    summarise(pop = n(), weightpop = sum(weight))
  
  bs.pairs.2017.summarised <- bs.pairs.2017.summarised %>%
    mutate(prop = pop/sum(bs.pairs.2017.summarised$pop), weightprop = weightpop/sum(bs.pairs.2017.summarised$weightpop)) %>%
    rename(gender = gender.1) 
  
  rho.tbl.agesex <- ibm.results.agesex.2017.noart.summarised %>%
    left_join(bs.pairs.2017.summarised, by= c("agegp.2017", "gender"), suffix = c(".all", ".infectors")) %>%
    mutate(pop.infectors = replace(pop.infectors, is.na(pop.infectors), 0)) %>%
    mutate(prop.infectors = replace(prop.infectors, is.na(prop.infectors), 0)) %>%
    mutate(weightpop = replace(weightpop, is.na(weightpop), 0)) %>%
    mutate(weightprop = replace(weightprop, is.na(weightprop), 0)) %>%
    mutate(rounded.weightpop = round(weightpop)) %>%
    mutate(rho = weightprop/prop.all)
  
  rho.tbl.agesex %>% ungroup %>%
    select(agegp.2017, gender, rho) %>%
    pivot_wider(names_from = "gender", values_from = "rho") %>%
    rename(Male = M, Female = "F") %>%
    mutate(ratio = Male/Female) %>%
    mutate(bs.number = i)
  
}) %>% bind_rows()

## Figures 3B and 3C ----

rho.ratio.agesex.noart <- rho.table.agesex.noart %>% ungroup() %>%
  select(agegp.2017, gender, rho) %>%
  pivot_wider(names_from = "gender", values_from = "rho") %>%
  rename(Male = M, Female = "F") %>%
  mutate(ratio = Male/Female) %>%
  mutate(ratio.bs.ci = map(agegp.2017, function(agmp){
    (bootstrapped.results.agesex.noart %>% filter(agegp.2017 == agmp) %>% pull(ratio) %>% sort())[c(5,195)]
  })) %>%
  mutate(ratio.lower.ci = map_dbl(ratio.bs.ci, extract(1)), ratio.upper.ci = map_dbl(ratio.bs.ci, extract(2))) %>%
  mutate(male.bs.ci = map(agegp.2017, function(agmp){
    (bootstrapped.results.agesex.noart %>% filter(agegp.2017 == agmp) %>% pull(Male) %>% sort())[c(5,195)]
  })) %>%
  mutate(male.lower.ci = map_dbl(male.bs.ci, extract(1)), male.upper.ci = map_dbl(male.bs.ci, extract(2))) %>%
  mutate(female.bs.ci = map(agegp.2017, function(agmp){
    (bootstrapped.results.agesex.noart %>% filter(agegp.2017 == agmp) %>% pull(Female) %>% sort())[c(5,195)]
  })) %>%
  mutate(female.lower.ci = map_dbl(female.bs.ci, extract(1)), female.upper.ci = map_dbl(female.bs.ci, extract(2))) %>%
  mutate(denominator = "As proportion of HIV+ individuals not on ART")

rho.ratio.pt2 <- rho.ratio.agesex.noart %>% 
  select(agegp.2017, Female, female.lower.ci, female.upper.ci, Male, male.lower.ci, male.upper.ci) %>%
  rename(female.estimate = Female, male.estimate = Male) %>%
  pivot_longer(2:7) %>%
  mutate(sex = map_chr(name, function(nm) str_split_fixed(nm, "\\.", 2)[,1])) %>%
  mutate(val = map_chr(name, function(nm) str_split_fixed(nm, "\\.", 2)[,2])) %>%
  select(-name) %>%
  pivot_wider(names_from = val, values_from = value) %>%
  mutate(denominator = "As proportion of HIV+ individuals not on ART")

rho.ratios <- bind_rows(rho.ratio.agesex, rho.ratio.agesex.noart)


mf.rhos <- bind_rows(rho.ratio.pt1, rho.ratio.pt2)

ggplot(mf.rhos %>% filter(estimate > 0)) + 
  geom_linerange(aes(x = agegp.2017, y=estimate, ymin = lower.ci, ymax = upper.ci, col = sex), size = 2, alpha = 0.33, position = position_dodge2(width = 0.25)) +
  facet_wrap(~denominator, ncol = 1) +
  geom_line(aes(x = agegp.2017, y = estimate, col = sex, group = sex), position = position_dodge2(width = 0.25)) +
  geom_point(aes(x = agegp.2017, y = estimate, col = sex), position = position_dodge2(width = 0.25)) +
  scale_colour_manual(values = c("#33a02c", "#1f78b4"), name = "Sex", labels = c("Female", "Male")) +
  ylab("Relative transmission rate") +
  xlab("Age group") +
  theme_bw() + 
  theme(legend.position="top")

ggsave("rho_by_sex.pdf", height = 4, width = 6)

ggplot(rho.ratios %>% filter(ratio > 0) %>% mutate(ratio.upper.ci = replace(ratio.upper.ci, ratio.upper.ci == Inf, 1000))) + 
  geom_linerange(aes(x = agegp.2017, ymin = ratio.lower.ci, ymax = ratio.upper.ci, col = denominator, group  =denominator), size = 2, alpha = 0.25, position = position_dodge2(width = 0.25)) +
  geom_line(aes(x = agegp.2017, y = ratio, col = denominator, group  =denominator), position = position_dodge2(width = 0.25)) +
  geom_point(aes(x = agegp.2017, y = ratio, col = denominator, group  =denominator), position = position_dodge2(width = 0.25)) +
  scale_color_manual(values = wes_palette("Darjeeling1")[1:2], name = "Denominator", labels = c("All HIV+", "HIV+, not on ART")) +
  scale_fill_manual(values = wes_palette("Darjeeling1")[1:2], name = "Denominator", labels = c("All HIV+", "HIV+, not on ART")) +
  ylab("Relative transmission rate\nratio (male/female)") +
  xlab("Age group") +
  theme_bw() +
  coord_cartesian(y = c(0.1, 100)) +
  scale_y_log10() + 
  theme(legend.position="top")


ggsave("rho_ratio.pdf", height = 2.5, width = 6)

# Histograms ----

## Function to calculate weighted CIs ----

ci_clt <- function(df, ci_level, colname) {
  ### The following calculations are based on the following setup for each category given in 'df':
  ### 1) Let X(i) be the measurement of variable X for sampled case i, i = 1 ... n (n=500 in this case)
  ### where X is a 0/1 variable indicating absence or presence of a selected category.
  ### From the X(i) samples we would like to estimate the
  ### true proportion p of the presence of the category in the population.
  ### Therefore X(i) are iid random variables with Binomial(1,p) distribution
  ###
  ### 2) Let Y(i) = w(i)*X(i)
  ### where w(i) is the sampling weight applied to variable X(i).
  ###
  ### We apply the CLT to the sum of the Y(i)'s, using:
  ### - E(Y(i)) = mu(i) = w(i) * E(X(i)) = w(i) * p (since w(i) is a constant and the X(i) are identically distributed)
  ### - Var(Y(i)) = sigma2(i) = w(i)^2 * Var(X(i)) = w(i)^2 * p*(1-p) (since the X(i) iid)
  ###
  ### Hence, by CLT:
  ###   Sum{Y(i) - mu(i)} / sigma -> N(0,1)
  ### where:
  ###   sigma = sqrt( Sum{ sigma2(i) } ) = sqrt( Sum{ w(i)^2 } ) * sqrt( p*(1-p) )
  ### and note that:
  ###   Sum{ mu(i) } = Sum{ w(i) } * p = n*p
  ### since the sampling weights are assumed to sum up to the sample size.
  ###
  ### Note: all the Sums are from i = 1, ..., n
  ###
  ### 3) Compute the approximate confidence interval for p based on the N(0,1) distribution
  ### in the usual way, by first estimating sigma replacing p for the estimated p.
  ###
  
  alpha = 1 - ci_level                                         # area outside the confidence band
  z = qnorm(1 - alpha/2)                                       # critical z-quantile from Normal(0,1)
  
  n = nrow(df)                                                 # Sample size (assuming no missing values)
  df$ws = df$weight / sum(df$weight) * n                          # Weights scaled to sum the sample size (assumed for sampling weights)
  Sdf <-  df %>% group_by_at(colname) %>% summarise(ws = sum(ws))   # Weighted-base estimate of the total by category (Sum{ Y(i) })
  sigma2 = sum( df$ws^2 )                                         # Sum of squared weights (note that we must NOT sum by category)
  Sdf[,"p"] = Sdf[,"ws"] / n                                       # Estimated proportion by category
  Sdf[,"Delta"] =  z * sqrt( sigma2 ) *
    sqrt( Sdf$p * (1 - Sdf$p) ) / n               # Semi-size of the CI by category
  LB_name = paste(formatC(alpha/2*100, format="g"), "%")       # Name for the CI's Lower Bound column
  UB_name = paste(formatC((1 - alpha/2)*100, format="g"), "%") # Name for the CI's Upper Bound column
  Sdf[,LB_name] = Sdf[,"p"] - Sdf[,"Delta"]                          # CI's Lower Bound
  Sdf[,UB_name] = Sdf[,"p"] + Sdf[,"Delta"]                          # CI's Upper Bound
  
  return(Sdf)
}

## Figure 3A ----

source.ages.and.sexes <- pairs %>%
  group_by(direction, source.age.at.infection.group) %>%
  summarise(total.weight = sum(weight)) %>%
  ungroup() %>%
  mutate(prop = total.weight/sum(total.weight))

nice.x.labels <- map_chr(levels(source.ages.and.sexes$source.age.at.infection.group), function(x){
  temp <- substr(x, 2, nchar(x)-1)
  comps <- str_split_fixed(temp,",", 2)[1,]
  glue("{comps[1]}-{as.numeric(comps[2])-1}")
})

ggplot(source.ages.and.sexes) + geom_col(aes(x=source.age.at.infection.group, y= prop*100, fill = direction), size = 0.5, col = "black") +
  facet_grid(~direction) +
  theme_minimal() +
  scale_fill_manual(values = c("#33a02c", "#1f78b4"), guide = "none") +
  scale_x_discrete(drop = F, labels = nice.x.labels) +
  xlab("Source age group at transmission") +
  ylab("Weighted percentage") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(panel.spacing = unit(3, "lines"))

ggsave("Source_age_adjusted_trialarm.pdf", width = 7, height =4)

## ANOVA for age ----

aov(source.age.at.infection ~source.arm ,  data = pairs %>% filter(gender.1 == "M"), weight = weight) %>% summary()
aov(source.age.at.infection ~source.arm ,  data = pairs %>% filter(gender.1 == "F"), weight = weight) %>% summary()

## Figure 4A ----

recent.tbl <- pairs %>% 
  ungroup() %>% 
  group_by(recent) %>%
  summarise(count = n(), weighted.count = sum(weight)) %>%
  ungroup() %>%
  mutate(total = sum(weighted.count)) %>%
  mutate(prop = 100*weighted.count/total)

for.chisq <- recent.tbl %>% filter(name == "weighted.count") %>%
  select(recipient.arm, count) %>%
  mutate(count = round(count)) %>%
  select(source.different.community, count, recipient.arm) %>%
  pivot_wider(names_from = "recipient.arm", values_from = "count") %>%
  ungroup() %>%
  select(-source.different.community) %>%
  as.matrix()

chisq.test(for.chisq)

ggplot(recent.tbl) + 
  geom_col(aes(x = recent, y=prop), position = "dodge", size = 0.5, col = "black",  fill = "red3") +
  theme_minimal() +
  scale_x_discrete(name = "Source infected\n< 1 year prior") +
  ylab("Weighted percentage") +
  ylim(c(0, 100))

ggsave("outside.png", width = 2, height = 4)

## chi-squared for recency ----

for.chisq <- pairs %>% 
  ungroup() %>% 
  group_by(recent, recipient.arm) %>%
  summarise(count = n(), weighted.count = sum(weight)) %>%
  ungroup() %>%
  group_by(recipient.arm) %>%
  mutate(count = round(weighted.count)) %>%
  select(recent, count, recipient.arm) %>%
  pivot_wider(names_from = "recipient.arm", values_from = "count") %>%
  ungroup() %>%
  select(-recent) %>%
  as.matrix()

chisq.test(for.chisq)

## Figure 4B ----

fl.res.tbl <- pairs %>% 
  ungroup() %>% 
  mutate(fl.res.source = factor(fl.res.source, levels = c("WT", "low", "high"))) %>% 
  group_by(fl.res.source) %>% 
  filter(!is.na(fl.res.source)) %>%
  summarise(count = n(), weighted.count = sum(weight)) %>%
  mutate(prop = 100*weighted.count/289) 

ggplot(fl.res.tbl)  + 
  geom_col(aes(x = fl.res.source, y=prop), position = "dodge", size = 0.5, col = "black", fill = "violetred3") +
  theme_minimal() + 
  # facet_wrap(~recipient.arm, ncol = 3) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c"), name = "Version", labels = c("Raw", "Weighted")) +
  scale_x_discrete(name = "First-line drug\nresistance (source)", labels = c("Wild-type", "Low", "High"))   +
  ylab("Weighted percentage") +
  ylim(c(0, 100))

ggsave("resistance.png", width = 3, height = 4)

## chi-squared for drug resistance ----

for.chisq <- pairs %>% 
  filter(!is.na(fl.res.source)) %>%
  ungroup() %>% 
  group_by(fl.res.source, recipient.arm) %>%
  summarise(count = n(), weighted.count = sum(weight)) %>%
  ungroup() %>%
  group_by(recipient.arm) %>%
  mutate(count = round(weighted.count)) %>%
  select(fl.res.source, count, recipient.arm) %>%
  pivot_wider(names_from = "recipient.arm", values_from = "count") %>%
  ungroup() %>%
  select(-fl.res.source) %>%
  as.matrix()

chisq.test(for.chisq)

## Figure 4C ----

outside.tbl <- pairs %>% 
  ungroup() %>% 
  group_by(source.different.community) %>% 
  summarise(count = n(), weighted.count = sum(weight)) %>%
  mutate(prop = 100*weighted.count/300) 

ggplot(outside.tbl) + 
  geom_col(aes(x = source.different.community, y=prop), position = "dodge", size = 0.5, col = "black", fill = "darkorange1") +
  theme_minimal() +
  # facet_wrap(~recipient.arm, ncol = 3) +
  scale_fill_brewer(palette = "Paired", name = "Version", labels = c("Raw", "Weighted"), guide = "none") +
  scale_x_discrete(name = "Source from\noutside community") +
  ylab("Weighted percentage") +
  ylim(c(0, 100))

ggsave("outside.png", width = 2, height = 4)

## chi-squared for import ----

for.chisq <- pairs %>% 
  # filter(!is.na(fl.res.source)) %>%
  ungroup() %>% 
  group_by(source.different.community, recipient.arm) %>%
  summarise(count = n(), weighted.count = sum(weight)) %>%
  ungroup() %>%
  group_by(recipient.arm) %>%
  mutate(count = round(weighted.count)) %>%
  select(source.different.community, count, recipient.arm) %>%
  pivot_wider(names_from = "recipient.arm", values_from = "count") %>%
  ungroup() %>%
  select(-source.different.community) %>%
  as.matrix()

chisq.test(for.chisq)

# Upset plots ----

## Setup ----

### Upset plot table ----

pairs.for.upset <- pairs %>%
  mutate(rn = row_number()) %>%
  select(rn, weight, recent, high.risk.male.infector, high.risk.female.infector, high.res, source.different.community) %>%
  pivot_longer(3:7) %>%
  mutate(nice.name = case_when(name == "recent" ~ "Infected < 1 year prior",
                               name == "high.risk.male.infector" ~ "Male aged 25-40",
                               name == "high.risk.female.infector" ~ "Female aged 20-35",
                               name == "high.res" ~ "High level first line DR",
                               name == "source.different.community" ~ "From different community"
  )) %>%
  group_by(rn, weight) %>%
  summarise(riskfactors = list(nice.name[which(value)])) %>%
  mutate(origin = "real")

### Null distribution table ----

cat.1.proportions <- pairs %>% group_by(cat.1) %>% summarise(pop = sum(weight)/300)
recent.proportions <- pairs %>% group_by(recent) %>% summarise(pop = sum(weight)/300)
high.res.proportions <- pairs %>% filter(!is.na(high.res)) %>% group_by(high.res) %>% summarise(pop = sum(weight)/300)
diff.com.proportions <- pairs  %>% group_by(source.different.community) %>% summarise(pop = sum(weight)/300)

cat.1.vector <- unique(cat.1.proportions$cat.1)
recent.vector <- glue("recent_{c(TRUE, FALSE)}")
high.res.vector <- glue("high.res_{c(TRUE, FALSE)}")
diff.com.vector <- glue("diff.com_{c(TRUE, FALSE)}")

expected.props.table <- expand.grid(cat.1.vector, c(TRUE, FALSE), c(TRUE, FALSE), c(TRUE, FALSE)) %>%
  as_tibble() %>%
  rename(cat.1 = Var1, recent = Var2, high.res = Var3, source.different.community = Var4) %>%
  left_join(cat.1.proportions) %>%
  left_join(recent.proportions, by="recent") %>%
  left_join(high.res.proportions, by="high.res") %>%
  left_join(diff.com.proportions, by="source.different.community") %>%
  mutate(expected.prop = pop.x*pop.y*pop.x.x*pop.y.y) %>%
  mutate(weight = expected.prop*302) %>%
  mutate(high.risk.male.infector = cat.1 == "hr.m") %>%
  mutate(high.risk.female.infector = cat.1 == "hr.f")  %>%
  mutate(rn = row_number()) %>%
  select(rn, weight, high.risk.male.infector, high.risk.female.infector, recent, high.res, source.different.community) %>%
  pivot_longer(3:7) %>%
  mutate(nice.name = case_when(name == "recent" ~ "Infected < 1 year prior",
                               name == "high.risk.male.infector" ~ "Male aged 25-40",
                               name == "high.risk.female.infector" ~ "Female aged 20-35",
                               name == "high.res" ~ "High level first line DR",
                               name == "source.different.community" ~ "From different community"
  )) %>%
  group_by(rn, weight) %>%
  summarise(riskfactors = list(nice.name[which(value)])) %>%
  mutate(origin  = "expected")

pairs.duo <- bind_rows(pairs.for.upset, expected.props.table) %>%
  ungroup() %>%
  mutate(tempref = map_chr(riskfactors, function(rf) paste0(sort(rf), collapse = "-"))) %>%
  group_by(tempref) %>%
  mutate(count=length(unique(origin))) %>%
  ungroup() %>%
  filter(count > 1) %>%
  mutate(origin = factor(origin, levels = c("real", "expected")))

# This requires a custom ggupset, which is available on request.

detach(package:ggupset, unload=TRUE)
install("~/ggupset/")
library(ggupset)

## Figure 5A ----

pdf(file = "upset_weighted_May2022.pdf", width = 7, height = 4)

ggplot(pairs.for.upset, aes(x=riskfactors, y=100*(..count..)/290, weight = weight)) +
  geom_bar(fill = "indianred3", col = "black") +
  xlab("Source risk factors") +
  ylab("Count") +
  ylab("Weighted percentage") +
  scale_x_upset(order_by = "weighted_freq") +
  theme_minimal() +
  axis_combmatrix(clip = "off") +
  theme_combmatrix(
    combmatrix.label.text = element_text(size = 7))
dev.off()

## Figure 5B ----

pdf(file = "upset_weighted_wExpected.pdf", width = 7, height = 4)

ggplot(pairs.duo, aes(x=riskfactors, y=100*(..count..)/290, fill = origin, group = origin, weight = weight)) +
  geom_bar(col = "black", position = "dodge", size = 0.25 ) +
  xlab("Source risk factors") +
  ylab("Count") +
  ylab("Weighted percentage") +
  scale_x_upset(order_by = "weighted_freq") +
  theme_minimal() +
  axis_combmatrix(clip = "off") +
  scale_fill_manual(values = wes_palette("Chevalier1"), name = "Value", labels = c("Actual", "Expected")) +
  theme_combmatrix(
    combmatrix.label.text = element_text(size = 7))
dev.off()

# Figure 5C ----

pairs.for.fake.upset <- pairs %>%
  mutate(male.source = gender.1 == "M") %>%
  mutate(rn = row_number()) %>%
  select(rn,weight, recent, male.source, high.risk.male.infector, high.risk.female.infector, high.res, source.different.community) %>%
  pivot_longer(3:8) %>%
  mutate(nice.name = case_when(name == "male.source" ~ "Male",
                               name == "recent" ~ "Infected < 1 year prior",
                               name == "high.risk.male.infector" ~ "Male aged 25-40",
                               name == "high.risk.female.infector" ~ "Female aged 20-35",
                               name == "high.res" ~ "High level first line DR",
                               name == "source.different.community" ~ "From different community"
  )) %>%
  group_by(rn, weight) %>%
  summarise(riskfactors = list(nice.name[which(value)])) 

## Rounding so there are 100 individuals and determining icon order ----

typical.sources <- pairs.for.fake.upset %>% group_by(riskfactors) %>% summarise(count = sum(weight)/3) %>%
  mutate(stringrf = map_chr(riskfactors, function(rfl){
    paste(rfl, collapse = " - ") 
    
  })) %>%
  mutate(whole.count = floor(count)) %>%
  mutate(remainder = count - whole.count) %>%
  arrange(-remainder) %>%
  mutate(rn = row_number()) %>%
  mutate(converted = map2_dbl(rn, whole.count, function(r, rc){
    ifelse(r <=9, rc +1 , rc)
  })) %>%
  filter(converted > 0) %>%
  mutate(factrf = factor(stringrf)) %>%
  mutate(factrf = fct_drop(factrf)) %>%
  mutate(factrf = fct_relevel(factrf, 
                              c("Male - Male aged 25-40",
                                "Male - Male aged 25-40 - High level first line DR",
                                "Infected < 1 year prior - Male - Male aged 25-40",
                                "Male - Male aged 25-40 - From different community",
                                "Infected < 1 year prior - Male - Male aged 25-40 - High level first line DR",
                                "Infected < 1 year prior - Male - Male aged 25-40 - From different community",
                                "Male - Male aged 25-40 - High level first line DR - From different community",
                                "Male",
                                "Infected < 1 year prior - Male",
                                "Male - From different community",
                                "Female aged 20-35",
                                "Female aged 20-35 - High level first line DR",
                                "Infected < 1 year prior - Female aged 20-35",
                                "Female aged 20-35 - From different community",
                                "Infected < 1 year prior - Female aged 20-35 - High level first line DR",
                                "Infected < 1 year prior - Female aged 20-35 - From different community",
                                "",
                                "Infected < 1 year prior"
                              ))) 

gridify <- function(element, wrap.threshold, rowwise = T){
  
  dim1 <- ifelse(element %% wrap.threshold !=0, element %% wrap.threshold, wrap.threshold)
  dim2 <- ifelse(element %% wrap.threshold !=0, (element %/% wrap.threshold) + 1, (element %/% wrap.threshold))
  if(rowwise){
    list(row = wrap.threshold-dim2 + 1, col = dim1)
  } else {
    list(row = wrap.threshold-dim1 + 1, col = dim2)
  }
}


## Expand to 100 individuals and set up the grid ----

typical.sources.grid <- typical.sources %>%
  uncount(converted) %>%
  arrange(factrf) %>%
  mutate(rn = row_number()) %>%
  mutate(grid.position = map(rn, function(n) gridify(n, 5, F))) %>%
  mutate(row = map_dbl(grid.position, function(gp) extract2(gp, 1)),
         col = map_dbl(grid.position, function(gp) extract2(gp, 2))) 

## Determine which icons should appear ----

typical.sources.grid <- typical.sources.grid %>%
  mutate(agesex = case_when(str_detect(stringrf, "Male") & !str_detect(stringrf, "aged 25-40") ~ "OM",
                            str_detect(stringrf, "Male")  ~ "YM",
                            str_detect(stringrf, "Female")  ~ "YW",
                            TRUE ~ "OW")) %>%
  mutate(outside = str_detect(stringrf, "community")) %>%
  mutate(recent = str_detect(stringrf, "prior")) %>%
  mutate(resistant = str_detect(stringrf, "level"))

## Load in the images ----

# .svg files here are from font-awesome and recoloured

src.file <- file("person-solid_dblue.svg")

rsvg_eps(charToRaw(read_lines(src.file)), file = "person-solid_dblue.eps")
PostScriptTrace("person-solid_dblue.eps")

yman.image <- readPicture("person-solid_dblue.eps.xml")

src.file <- file("person-solid_lblue.svg")

rsvg_eps(charToRaw(read_lines(src.file)), file = "person-solid_lblue.eps")
PostScriptTrace("person-solid_lblue.eps")

oman.image <- readPicture("person-solid_lblue.eps.xml")

src.file <- file("person-dress-solid_dgreen.svg")

rsvg_eps(charToRaw(read_lines(src.file)), file = "person-dress-solid_dgreen.eps")
PostScriptTrace("person-dress-solid_dgreen.eps")

ywoman.image <- readPicture("person-dress-solid_dgreen.eps.xml")



src.file <- file("person-dress-solid_lgreen.svg")

rsvg_eps(charToRaw(read_lines(src.file)), file = "person-dress-solid_lgreen.eps")
PostScriptTrace("person-dress-solid_lgreen.eps")

owoman.image <- readPicture("person-dress-solid_lgreen.eps.xml")

src.file <- file("stopwatch-solid.svg")

rsvg_eps(charToRaw(read_lines(src.file)), file = "stopwatch-solid.eps")
PostScriptTrace("stopwatch-solid.eps")

clock.image <- readPicture("stopwatch-solid.eps.xml")

src.file <- file("arrow-right-arrow-left-solid.svg")

rsvg_eps(charToRaw(read_lines(src.file)), file = "arrow-right-arrow-left-solid.eps")
PostScriptTrace("arrow-right-arrow-left-solid.eps")

arrows.image <- readPicture("arrow-right-arrow-left-solid.eps.xml")

xlims <- c(0,21)
ylims <- c(0,6)

ymen.grob <- symbolsGrob(yman.image,
                         x=rescale(typical.sources.grid %>% filter(agesex == "YM") %>% pull(col), from = xlims),
                         y=rescale(typical.sources.grid %>% filter(agesex == "YM") %>% pull(row), from = ylims),
                         size=0.15)


omen.grob <- symbolsGrob(oman.image,
                         x=rescale(typical.sources.grid %>% filter(agesex == "OM") %>% pull(col), from = xlims),
                         y=rescale(typical.sources.grid %>% filter(agesex == "OM") %>% pull(row), from = ylims),
                         size=0.15)

ywomen.grob <- symbolsGrob(ywoman.image,
                           x=rescale(typical.sources.grid %>% filter(agesex == "YW") %>% pull(col), from = xlims),
                           y=rescale(typical.sources.grid %>% filter(agesex == "YW") %>% pull(row), from = ylims),
                           size=0.15)


owomen.grob <- symbolsGrob(owoman.image,
                           x=rescale(typical.sources.grid %>% filter(agesex == "OW") %>% pull(col), from = xlims),
                           y=rescale(typical.sources.grid %>% filter(agesex == "OW") %>% pull(row), from = ylims),
                           size=0.15)

clocks.grob <- symbolsGrob(clock.image,
                           x=rescale(typical.sources.grid %>% filter(recent) %>% pull(col) + 0.25, from = xlims),
                           y=rescale(typical.sources.grid %>% filter(recent) %>% pull(row) - 0.3, from = ylims),
                           size=0.05)

arrows.grob <- symbolsGrob(arrows.image,
                           x=rescale(typical.sources.grid %>% filter(outside) %>% pull(col) - 0.25, from = xlims),
                           y=rescale(typical.sources.grid %>% filter(outside) %>% pull(row) - 0.3, from = ylims),
                           size=0.04)

dr.grob <- textGrob("DR",
                    x=rescale(typical.sources.grid %>% filter(resistant) %>% pull(col), from = xlims),
                    y=rescale(typical.sources.grid %>% filter(resistant) %>% pull(row) +0.06, from = ylims),
                    gp = gpar(col= "white", fontsize = 3)
)


ggplot() + 
  annotation_custom(ymen.grob) +
  annotation_custom(omen.grob) +
  annotation_custom(ywomen.grob) +
  annotation_custom(owomen.grob) +
  annotation_custom(clocks.grob) +
  annotation_custom(arrows.grob) +
  annotation_custom(dr.grob) +
  xlim(c(0, 21)) + ylim(c(0, 6)) +
  theme_void()

ggsave("hundred.pdf", width = 3.5, height = 1.5)

# Sensitivity analysis ----

pairs.hcf.only <- pairs %>% filter(cohort.1 == "HCF" & cohort.2 == "HCF")

raked <- anesrake(targets, pairs.hcf.only %>% as.data.frame, caseid = pairs.hcf.only$recipient.id, verbose = T)
pairs.hcf.only$weight <- raked$weightvec


