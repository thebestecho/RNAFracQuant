# Filter for only selected ORFs, and calculate TPMs
TSP_Count_TPM <- TSP_Count %>%
    inner_join(ORFselect) %>%
    group_by(Sample,Condition,Fraction) %>%
    mutate(Density=Count/Length, 
           TPM=Density/sum(Density)*1e6  ) %>%
    select(Sample,Condition,Fraction,ORF,Count,TPM)

# Transform the data frame into the wider one
TSP_TPM_wide <- 
    TSP_Count_TPM %>%
    ungroup() %>%
    select(Condition,Fraction,ORF,TPM) %>%
    unite(Var,Condition,Fraction,sep="_") %>%
    spread(Var,TPM)

