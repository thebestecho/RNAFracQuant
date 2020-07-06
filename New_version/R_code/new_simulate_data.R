# poisson mean = total * scaling factor * psup, or *(1-psup) for pellet.
# Tot = scaling_factor_1 * Sup + scaling_factor_2 * Pellet
# We assume there is only one condition ( one replicate)
# And there are two fractions: Sup & Pellet
# Create a folder to store those simulated data
dir.create(here::here("new_simulated_data_in"))
ORF <- paste0("G",1:20) # The first column (mRNA-name)
set.seed(5)
Tot <- floor(rlnorm(20, log(8000), log(1.5)))
# pSup = scaling_factor_1 * Sup / (scaling_factor_1 * Sup + scaling_factor_2 * Pellet)
pSup <- runif(20, min=0, max=1) # simulated pSup, the proportion of mRNA in supernatant fraction1
scaling_factor_1 <- runif(1, min=0, max=1) #scaling factor for Sup
scaling_factor_2 <- runif(1, min=0, max=1) #scaling factor for Pellet
Sup <- rpois(20, lambda = Tot * scaling_factor_1 * pSup)
Pellet <- rpois(20, lambda = Tot * scaling_factor_2 * (1-pSup))
sample1 <- data.frame(ORF, Tot)
write_tsv(sample1,
          here::here("new_simulated_data_in","sample1.txt"),col_names = FALSE) # total mRNA counts
sample2 <- data.frame(ORF, Sup)
write_tsv(sample2,
          here::here("new_simulated_data_in","sample2.txt"),col_names = FALSE) # mRNA counts of Sup
sample3 <- data.frame(ORF, Pellet)
write_tsv(sample3,
          here::here("new_simulated_data_in","sample3.txt"),col_names = FALSE) # mRNA counts of Pellet
