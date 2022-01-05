library(tidyverse)

load("Data/w1network_dvn.RData")
dat_net = x
load("Data/w1inhome_dvn.RData")
dat_home = x
remove(x)

#load("Data/w1context_dvn.RData")
#load("Data/w1weight.RData")

# Network data?
# https://networks.skewed.de/net/add_health

# Sex
# Section 1, pg 5
# BIO_SEX: Interviewer, please confirm that R's sex is (male) female. (Ask if necessary.)

# Race (make binary -- {1=white,0=non-white})
# Section 1, pg 11
# H1GI4: Are you of Hispanic or Latino origin?
# What is your race?
# H1GI6A: White
# H1GI6B: Black or African American
# H1GI6C: American Indian or Native American
# H1GI6D: Asian or Pacific Islander
# H1GI6E: Other

# Grade (make discrete -- {6,7,...,12}?)
# Section 5, pg 60
# H1ED11: what was your grade in English or language arts?
# H1ED12: And what was your grade in mathematics?
# H1ED13: And what was your grade in history or social studies? 
# H1ED14: And what was your grade in science?
# GPA is the mean grade across four core subjects from the in-school questionnaire (items S10a - S10d). Grades are weighted as follows: A = 4, B = 3, C = 2, D or F = 1. GPA was calculated using only valid responses.

# Network (S&R = union of send and receive network (includes ego))
# NESR: Size of ego send- and receive-network
# AXGPA: S&R alter mean: gpa 
# NAGPA: Ego Net Denominator axgpa
# AXNUMACT: S&R alter mean: numact (# extracurricular activities)
# AXS2: S&R alter mean: s2 (sex)
# NAS2: Ego Net Denominator axs2

# Saliency (tendency to nominate others with similar characteristics)

clean_sex_f = function(x)
{
  # Assuming AXS2=0 is all male and AXS2=1 is all female
  sex_map = c("(1) Male"=0, "(2) Female"=1)
  return(sex_map[as.character(x)])
}
clean_grade_f = function(x)
{
  grade_map = c("(1) A"=4, "(2) B"=3, "(3) C"=2, "(4) D or lower"=1)
  return(grade_map[as.character(x)])
}
dat_net = dat_net %>%
  select(AID, NESR, AXGPA, AXS2, NAGPA, NAS2)
# Note: denominators are not always the same
dat_home = dat_home %>%
  select(AID, BIO_SEX,
         #H1GI4, H1GI6A, H1GI6B, H1GI6C, H1GI6D, H1GI6E,
         H1ED11, H1ED12, H1ED13, H1ED14) %>%
  mutate_at(vars(H1ED11, H1ED12, H1ED13, H1ED14), clean_grade_f) %>%
  rowwise() %>%
  mutate(BIO_SEX = clean_sex_f(BIO_SEX),
         GPA = mean(c(H1ED11, H1ED12, H1ED13, H1ED14), na.rm=T)) %>%
  filter(!is.na(BIO_SEX) & !is.nan(GPA)) %>%
  select(AID, BIO_SEX, GPA)
dat = inner_join(dat_home, dat_net) %>%
  # Scale GPA from [1,4] to [6,12]
  mutate(GPA = 2*GPA+4, AXGPA = 2*AXGPA+4) %>%
  na.omit()


# Simulation

sim = dat

# Scenario 1
expit = function(p) {1 / (1+exp(-p))}
sim$Z = rbernoulli(nrow(sim), expit(-18 + 2*sim$GPA + 3*sim$BIO_SEX))*1
