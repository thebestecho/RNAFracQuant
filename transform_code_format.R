TSP_Count_FracBySample_pSup <- 
    TSP_Count_FracBySample %>%
    left_join(select(mixing_ratios_TSP,Condition,mixing.Sup,mixing.P100)) %>%
    mutate(mixSum=mixing.Sup*Sup+mixing.P100*P100,
           pSup=mixing.Sup*Sup/mixSum)
