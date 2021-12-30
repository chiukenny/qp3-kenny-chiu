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

# Network (S&R = union of send and receive network)
# NESR: Size of ego send- and receive-network
# AXGPA: S&R alter mean: gpa 
# AXNUMACT: S&R alter mean: numact (# extracurricular activities)
# AXS2: S&R alter mean: s2 (sex)
# Saliency (tendency to nominate others with similar characteristics)

dat_net = dat_net %>%
  select(AID, AXGPA, AXS2)
dat_home = dat_home %>%
  select(AID, H1GI4, H1GI6A, H1GI6B, H1GI6C, H1GI6D, H1GI6E,
         H1ED11, H1ED12, H1ED13, H1ED14)
dat = dat_home %>%
  inner_join(dat_net)
