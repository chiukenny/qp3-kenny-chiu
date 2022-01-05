# Helper functions
# ---------------------------------------------------------

# Compute expit
expit = function(x) {1 / (1+exp(-x))}

# Compute individual propensity score (Scenario 1)
ind_prop = function(x1, x2)
{
  return(expit(-3 + 3*x1 + 4*x2))
}

# Compute outcome mean (model for main effects simulations)
outcome_mean = function(z, g, x1, x2, beta)
{
  prop = ind_prop(x1, x2)
  return(5 + 6*(prop>=0.7) + 10*z - 3*z*(prop>=0.7) + beta*g)
}



# Bias functions
# ---------------------------------------------------------

# Variable names used in bias functions:
#   - i_X: indices of nodes that satisfy X
#   - n_X: number of nodes that satisfy X
#   - PY_X: (proportional) probability of Y given X

# Compute bias (Theorem 2B, eqn.12)
bias_T2B = function(sim, beta)
{
  n = nrow(sim)
  i_1 = which(sim$Z==1)
  n_1 = length(i_1)
  i_0 = sim$id[-i_1]
  n_0 = n - n_1
  
  bias = 0
  for (g in unique(sim$G))
  {
    i_g = which(sim$G==g)
    n_g = length(i_g)
    i_1g = intersect(i_1, i_g)
    n_1g = length(i_1g)
    i_0g = intersect(i_0, i_g)
    n_0g = n_g - n_1g
    
    # Subset by unmeasured confounders (game1, game2)
    for (u1 in 0:1)
    {
      for (u2 in 0:1)
      {
        i_uu = which(sim$game1==u1 & sim$game2==u2)
        
        # Compute empirical probabilities
        Pug_1 = length(intersect(i_1g, i_uu)) / n_1
        Pug_0 = length(intersect(i_0g, i_uu)) / n_0
        Pug = length(intersect(i_g, i_uu)) / n
        
        # Take g'=0, u'=(0,0)
        bias = bias +
          (outcome_mean(1,g,u1,u2,beta)-outcome_mean(1,0,0,0,beta))*(Pug_1-Pug) -
          (outcome_mean(0,g,u1,u2,beta)-outcome_mean(0,0,0,0,beta))*(Pug_0-Pug)
      }
    }
  }
  return(bias)
}

# Compute bias given individual covariates (Corollary 2, eqn.11)
bias_C2I = function(sim, beta)
{
  i_1 = which(sim$Z==1)
  i_0 = sim$id[-i_1]
  
  bias = 0
  # Subset by individual covariates (game1, game2)
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      i_xx = which(sim$game1==x1 & sim$game2==x2)
      n_xx = length(i_xx)
      i_1xx = intersect(i_1, i_xx)
      n_1xx = length(i_1xx)
      i_0xx = intersect(i_0, i_xx)
      n_0xx = n_xx - n_1xx
      
      bias_xx = 0
      for (g in unique(sim$G))
      {
        # Compute empirical probabilities
        Pg_1xx = length(which(sim$G[i_1xx]==g)) / n_1xx
        Pg_0xx = length(which(sim$G[i_0xx]==g)) / n_0xx
        
        # Take g'=0
        bias_xx = bias_xx +
          (outcome_mean(0,g,x1,x2,beta)-outcome_mean(0,0,x1,x2,beta))*(Pg_1xx-Pg_0xx)
      }
      bias = bias + bias_xx*n_xx
    }
  }
  return(bias/nrow(sim))
}


# Compute bias given individual and neighbourhood covariates (Corollary 2, eqn.11)
bias_C2N = function(sim, beta)
{
  bias = 0
  # Subset by (game1, game2, N, Ngame1, Ngame2)
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      for (N in unique(sim$N[which(sim$game1==x1 & sim$game2==x2)]))
      {
        i_xxN = which(sim$game1==x1 & sim$game2==x2 & sim$N==N)
        sim_xxN = sim[i_xxN,]
        i_1xxN = which(sim_xxN$Z==1)
        i_0xxN = which(sim_xxN$Z==0)
        
        for (n1 in unique(sim_xxN$Ngame1))
        {
          for (n2 in unique(sim_xxN$Ngame2))
          {
            i_xxnnN = which(sim_xxN$Ngame1==n1 & sim_xxN$Ngame2==n2)
            n_xxnnN = length(i_xxnnN)
            if (n_xxnnN==0) {next}
            i_1xxnnN = intersect(i_xxnnN,i_1xxN)
            n_1xxnnN = length(i_1xxnnN)
            i_0xxnnN = intersect(i_xxnnN,i_0xxN)
            n_0xxnnN = length(i_0xxnnN)
            
            bias_xxnnN = 0
            for (g in unique(sim_xxN$G[i_xxnnN]))
            {
              i_gxxN = which(sim_xxN$G==g)
              
              # Compute empirical probabilities
              p = 0
              if (n_1xxnnN>0) {p = p+length(intersect(i_1xxnnN,i_gxxN))/n_1xxnnN}
              if (n_0xxnnN>0) {p = p-length(intersect(i_0xxnnN,i_gxxN))/n_0xxnnN}
              
              # Take g'=0
              bias_xxnnN = bias_xxnnN +
                (outcome_mean(0,g,x1,x2,beta)-outcome_mean(0,0,x1,x2,beta))*p
            }
            bias = bias + bias_xxnnN*n_xxnnN
          }
        }
      }
    }
  }
  return(bias/nrow(sim))
}