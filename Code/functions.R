# Helper functions
# ---------------------------------------------------------

# Compute expit
expit = function(x) {1 / (1+exp(-x))}

# Compute individual propensity score (Scenario 1)
ind_prop = function(x1, x2)
{
  return(expit(-3 + 3*x1 + 4*x2))
}

# Compute outcome mean (model for main effects datulations)
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
bias_T2B = function(dat, beta)
{
  n = nrow(dat)
  i_1 = dat$Z==1
  n_1 = sum(i_1)
  i_0 = !i_1
  n_0 = n - n_1
  
  bias = 0
  for (g in unique(dat$G))
  {
    i_g = dat$G==g
    n_g = sum(i_g)
    i_1g = i_1 & i_g
    n_1g = sum(i_1g)
    i_0g = i_0 & i_g
    n_0g = n_g - n_1g
    
    # Subset by unmeasured confounders (game1, game2)
    for (u1 in 0:1)
    {
      for (u2 in 0:1)
      {
        i_uu = dat$game1==u1 & dat$game2==u2
        
        # Compute empirical probabilities
        Pug_1 = sum(i_1g & i_uu) / n_1
        Pug_0 = sum(i_0g & i_uu) / n_0
        Pug = sum(i_g & i_uu) / n
        
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
bias_C2I = function(dat, beta)
{
  i_1 = dat$Z==1
  i_0 = !i_1
  
  bias = 0
  # Subset by individual covariates (game1, game2)
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      i_xx = dat$game1==x1 & dat$game2==x2
      n_xx = sum(i_xx)
      i_1xx = i_1 & i_xx
      n_1xx = sum(i_1xx)
      i_0xx = i_0 & i_xx
      n_0xx = n_xx - n_1xx
      
      bias_xx = 0
      for (g in unique(dat$G))
      {
        # Compute empirical probabilities
        Pg_1xx = sum(dat$G[i_1xx]==g) / n_1xx
        Pg_0xx = sum(dat$G[i_0xx]==g) / n_0xx
        
        # Take g'=0
        bias_xx = bias_xx +
          (outcome_mean(0,g,x1,x2,beta)-outcome_mean(0,0,x1,x2,beta))*(Pg_1xx-Pg_0xx)
      }
      bias = bias + bias_xx*n_xx
    }
  }
  return(bias/nrow(dat))
}

# Compute bias given individual and neighbourhood covariates (Corollary 2, eqn.11)
bias_C2N = function(dat, beta)
{
  bias = 0
  # Subset by (game1, game2, N, Ngame1, Ngame2)
  for (x1 in 0:1)
  {
    for (x2 in 0:1)
    {
      for (N in unique(dat$N[dat$game1==x1 & dat$game2==x2]))
      {
        i_xxN = dat$game1==x1 & dat$game2==x2 & dat$N==N
        dat_xxN = dat[i_xxN,] # Save subsetted data for performance
        i_1xxN = dat_xxN$Z==1

        for (n1 in unique(dat_xxN$Ngame1))
        {
          for (n2 in unique(dat_xxN$Ngame2))
          {
            i_xxnnN = dat_xxN$Ngame1==n1 & dat_xxN$Ngame2==n2
            n_xxnnN = sum(i_xxnnN)
            if (n_xxnnN==0) {next}
            i_1xxnnN = i_xxnnN & i_1xxN
            n_1xxnnN = sum(i_1xxnnN)
            i_0xxnnN = i_xxnnN & !i_1xxN
            n_0xxnnN = n_xxnnN - n_1xxnnN

            bias_xxnnN = 0
            for (g in unique(dat_xxN$G[i_xxnnN]))
            {
              i_gxxN = dat_xxN$G==g

              # Compute empirical probabilities
              p = 0
              if (n_1xxnnN>0) {p = p+sum(i_1xxnnN & i_gxxN)/n_1xxnnN}
              if (n_0xxnnN>0) {p = p-sum(i_0xxnnN & i_gxxN)/n_0xxnnN}

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
  return(bias/nrow(dat))
}



# Subclassification functions
# ---------------------------------------------------------

# Estimates treatment effect by classifying units into j roughly
# equal sized subclasses based on individual covariates
subcl_ind = function(dat, j)
{
  if (j > 1)
  {
    classes = cut(glm(Z~game1+game2,data=dat,family=binomial)$fitted.values,
                  j, labels=F)
  } else {
    classes = rep(1, nrow(dat))
  }
  return(classes)
}

# Estimates treatment effect by classifying units into j roughly
# equal sized subclasses based on all covariates
subcl_all = function(dat, j)
{
  if (j > 1)
  {
    classes = cut(glm(Z~game1+game2+Ngame1+Ngame2+N,
                      data=dat,family=binomial)$fitted.values,
                  j, labels=F)
  } else {
    classes = rep(1, nrow(dat))
  }
  return(classes)
}

subcl_est = function(dat, classes)
{
  i_z = dat$Z==1
  tau_est = 0
  for (c in unique(classes))
  {
    i_c = classes==c
    tau_c = 0
    Y_c1 = mean(dat$Y[i_c & i_z])
    Y_c0 = mean(dat$Y[i_c & !i_z])
    if (!is.na(Y_c1)) {tau_c = tau_c + Y_c1}
    if (!is.na(Y_c0)) {tau_c = tau_c - Y_c0}
    tau_est = tau_est + tau_c*sum(i_c)
  }
  return(tau_est/nrow(dat))
}

subcl_gps = function(dat, j, sum_exposure=F)
{
  # Subclassification based on individual propensity
  ind_classes = subcl_all(dat, j)
  
  n = nrow(dat)
  
  obs_g = sort(unique(dat$G))
  n_g = length(obs_g)
  tau_g = rep(0, n_g)
  
  for (c in unique(ind_classes))
  {
    dat_c = dat[ind_classes==c,]
    n_c = nrow(dat_c)
    
    if (sum_exposure)
    {
      # TODO
    } else {
      # Note: assumes V_g = all units in subclass
      prop_fit = glm(G~Z+Ngame1+Ngame2+N, data=dat_c, family=binomial, weights=N)
      dat_c$nbr_prop = prop_fit$fitted.values
      lin_fit = lm(Y~Z+G+nbr_prop, data=dat_c)
      
      for (i in 1:n_g)
      {
        g = obs_g[i]
        dat_c$G = g
        
        dat_c$Z = 1
        dat_c$nbr_prop = predict(prop_fit, newdata=dat_c)
        Y_c1g = mean(predict(lin_fit,newdata=dat_c))
        
        dat_c$Z = 0
        dat_c$nbr_prop = predict(prop_fit, newdata=dat_c)
        Y_c0g = mean(predict(lin_fit,newdata=dat_c))
        
        tau_g[i] = tau_g[i] + (Y_c1g-Y_c0g)*n_c/n
      }
    }
  }
  return(sum(tau_g * table(dat$G)/n))
}
