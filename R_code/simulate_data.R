# poisson mean = total * scaling factor * psup, or *(1-psup) for pellet.
# We assume there are two conditions, 30C and 42C
# Create a folder to store those simulated data
dir.create(here::here("simulated_data_in"))
ORF <- paste0("G",1:20)
set.seed(5)
Tot1 <- floor(rlnorm(20, log(5000), log(1.5)))
pSup1 <- runif(20, min=0, max=1) #pSup for condition 30C
scaling_factor1S <- runif(1, min=0, max=1) #scaling factor for Sup in condition 30C
scaling_factor1P <- runif(1, min=0, max=1) #scaling factor for Pellet in condition 30C
Sup1 <- rpois(20, lambda = Tot1*scaling_factor1S*pSup1)
Pellet1 <- rpois(20, lambda = Tot1*scaling_factor1P*(1-pSup1))
sample1 <- data.frame(ORF, Tot1)
write_tsv(sample1,
          here::here("simulated_data_in","sample1.txt"),col_names = FALSE)
sample2 <- data.frame(ORF, Sup1)
write_tsv(sample2,
          here::here("simulated_data_in","sample2.txt"),col_names = FALSE)
sample3 <- data.frame(ORF, Pellet1)
write_tsv(sample3,
          here::here("simulated_data_in","sample3.txt"),col_names = FALSE)

Length <- c(floor(rlnorm(20, log(1000), log(1.5))))
fake_ORF_select <- data.frame(ORF,Length)
write_tsv(fake_ORF_select,
          here::here("simulated_data_in","Scer_ORF_length.txt"))

# The same for condition 42C
set.seed(10)
Tot2 <- floor(rlnorm(20, log(5000), log(1.5)))
pSup2 <- runif(20, min=0, max=1) #pSup for condition 42C
scaling_factor2S <- runif(1, min=0, max=1) #scaling factor for Sup in condition 42C
scaling_factor2P <- runif(1, min=0, max=1) #scaling factor for Pellet in condition 42C
Sup2 <- rpois(20, lambda = Tot2*scaling_factor2S*pSup2)
Pellet2 <- rpois(20, lambda = Tot2*scaling_factor2P*(1-pSup2))
sample4 <- data.frame(ORF, Tot2)
write_tsv(sample4,
          here::here("simulated_data_in","sample4.txt"),col_names = FALSE)
sample5 <- data.frame(ORF, Sup2)
write_tsv(sample5,
          here::here("simulated_data_in","sample5.txt"),col_names = FALSE)
sample6 <- data.frame(ORF, Pellet1)
write_tsv(sample6,
          here::here("simulated_data_in","sample6.txt"),col_names = FALSE)
