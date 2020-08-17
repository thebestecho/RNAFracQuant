# run stan mixing ratio inference, reusing same model
mixing_ratios_TSP <- TSP_Count_FracBySample %>%
    ddply(~Condition,getmixingratiosTSP,
          stan_mixing=stan_mixingTSP_wd5,iter=1000,seed=myseed)
