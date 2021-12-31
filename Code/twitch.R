# Games -> Self-promotion -> Views?

# https://snap.stanford.edu/data/twitch_gamers.html
# https://arxiv.org/pdf/1909.13021.pdf
# Nodes: 7126
# Features: 3169

# Top features:
#   feats Freq
#1   920 6742
#2   224 3816
#3   569 3710
#4   861 3242
#5  3152 3180
#6  2645 2461

library(rjson)
library(tidyverse)
library(Matrix)

dat_json = fromJSON(file = "musae_ENGB_features.json")
dat_mat = matrix(0, 7126, 3)

for (user_str in names(dat_json))
{
  id = as.numeric(user_str) + 1
  dat_mat[id,1] = id
  #feat_i = unname(unlist(dat_json[user]))
  for (g in unlist(dat_json[user_str]))
  {
    if (g == 224) {dat_mat[id,2] = 1}
    else if (g == 569) {dat_mat[id,3] = 1}
  }
  #dat = bind_rows(dat, c(id=id, game1=game1, game2=game2))
  #feats = c(feats, feat_i)
  #users = c(users, rep(user_i, length(feat_i)))
}
dat = as.data.frame(dat_mat)
colnames(dat) = c("id", "game1", "game2")

remove(dat_json)
remove(dat_mat)

edge_csv = read.csv("musae_ENGB_edges.csv", header=T) + 1
dat_edMat = sparseMatrix(i=edge_csv$from, j=edge_csv$to, symmetric=T)
remove(edge_csv)
#degrees = as.data.frame(table(c(dat_csv$from,dat_csv$to))) %>%
#  arrange(desc(Freq))
#summary(degrees$Freq)


# Simulation

sim = dat
sim$Ngame1 = 0
sim$Ngame2 = 0
sim$N = 0

for (i in 1:nrow(sim))
{
  id = sim$id[i]
  ii = which(dat_edMat[id,])
  sim$Ngame1[i] = mean(sim$game1[ii])
  sim$Ngame2[i] = mean(sim$game2[ii])
  sim$N[i] = length(ii)
}


# Scenario 1
expit = function(p) {1 / (1+exp(-p))}
sim$Z = rbernoulli(nrow(sim), expit(-3 + 3*sim$game1 + 5*sim$game2))*1

# Function to compute neighbourhood covariates
nbrhd_trt = function(x)
{
  return(sum(sim$Z[which(dat_edMat[x,])]))
}
sim$G = sapply(sim$id, nbrhd_trt)
