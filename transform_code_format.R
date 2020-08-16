# Create a new directory in the current working directory
dir.create(here::here("simulated_data"))
ORF <- paste0("G",1:20)
set.seed(5) # For reproducing the same data
Total <- floor (rlnorm (20, log (8000), log (1.5)))
pSup <- runif (20, min=0, max=1) # pSup for condition 30C
# scaling factor for Sup in condition 30C
scaling_factor1S <- runif (1, min=0, max=1)
# scaling factor for Pellet in condition 30C
scaling_factor1P <- runif (1, min=0, max=1) 
Tot1 <- rpois (20, lambda = 1000) # N total for condition 30C
Sup1 <- rpois (20, lambda = Total*scaling_factor1S*pSup) # N sup
Pellet1 <- rpois (20, lambda = Total*scaling_factor1P*(1-pSup)) # N pellet
# Repeat for condition 42C
set.seed(10)
scaling_factor2S <- runif (1, min=0, max=1)
scaling_factor2P <- runif (1, min=0, max=1)
Tot2 <- rpois (20, lambda = 1000) # N total for condition 42C
Sup2 <- rpois (20, lambda = Total*scaling_factor2S*pSup) # N sup
Pellet2 <- rpois (20, lambda = Total*scaling_factor2P*(1-pSup)) # N pellet
# Write data into sample files (for example)
sample1 <- data.frame(ORF, Tot1) 
write_tsv (sample1,
          here::here("simulated_data","sample1.txt"),col_names = FALSE)
