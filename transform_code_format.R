mdat_mixingTSP <- function(ctdata) {
    ## make data for stan fit
    head(ctdata)
    list(NRNA=nrow(ctdata),
         tot_obs=as.integer(round(ctdata$Tot)),
         sup_obs=as.integer(round(ctdata$Sup)),
         p100_obs=as.integer(round(ctdata$P100))
         )
    
}

make_mixingTSP <- function(ctdata) {
    stan_dat <- mdat_mixingTSP(ctdata)
    stan(model_code='// -*- mode: C -*-
data {
  // Number of RNAs
  int<lower=1> NRNA;     
  
  // Note: These are all integers
  // columns t, s, p
  int<lower=0> tot_obs[NRNA];
  int<lower=0> sup_obs[NRNA];
  int<lower=0> p100_obs[NRNA];
}
parameters {
  // Unnormalized mixing proportions
  // real<lower=0> mixing_t;
  real<lower=0> mixing_sup;
  real<lower=0> mixing_p100;
  
  // dispersion parameter for counts
  real phi;
}
model{
  // mixing ratios
  mixing_sup ~ gamma(1,1);
  mixing_p100 ~ gamma(1,1);
  // Cauchy prior for negbin dispersion parameter
  phi ~ cauchy(0,3);
  
  for(idx in 1:NRNA){ 
    // count distn negative binomial with specified means
    // Total
    tot_obs[idx] ~ neg_binomial_2(mixing_sup * sup_obs[idx] + mixing_p100 * p100_obs[idx], phi);
  }

}
generated quantities{
  // print("Mixing pars (sup,p100) = (", mixing_sup,",",mixing_p100,")");
  // print("dispersion phi = ", phi);
  // print("------------------");
}
',
data=stan_dat,chains = 1,iter = 10)
}

fit_mixingTSP <- function(ctdata,stan_mixing=NULL,...) {
    stan_dat <- mdat_mixingTSP(ctdata)
    if (is.null(stan_mixing)) {
        stan_mixing <- make_mixingTSP(ctdata)
    }
    stan_mixing_fit <- stan(fit=stan_mixing,data=stan_dat,chains = 4,...)
    return(stan_mixing_fit)
}

getmixingratiosTSP <- function(ctdata,iter=1000,
                                control=list(adapt_delta=0.85),...) {
    # head(ctdata) %>% print()
    stansummary <- fit_mixingTSP(ctdata=ctdata,iter=iter,control=control,...) %>%
        summary()
    # return medians
    data.frame(
        mixing.Sup=stansummary$summary["mixing_sup","50%"],
        mixing.P100=stansummary$summary["mixing_p100","50%"],
        lp.n_eff  =stansummary$summary["lp__","n_eff"],
        lp.Rhat   =stansummary$summary["lp__","Rhat"])
}
