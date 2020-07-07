# To prepare the packages
.First <- function() {
  # In case of the users have loaded dplyr befor plyr
  if("dplyr" %in% (.packages())){
    detach("package:dplyr", unload=TRUE) 
    detach("package:plyr", unload=TRUE) 
  } 
  library(plyr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tidyverse)
  library(rstan)
  library(here)
}


# Read the Samplesheet file that is provided by users
# Return a single data frame
# There are two arguments, dir_in: the input directory & file: Samplesheet file name
# Users need to make sure that this input directory is under the working directory
read_samplesheet <- function(dir_in,file)
{
  load_samplesheet = read_tsv(here::here(dir_in, file), comment = "#")
      return(load_samplesheet)
}


# Use the file names that is read from Samplesheet to read all the count files
# All the count files should be within the same directory as the Samplesheet file
# There are two arguments, 
# dir_in: the input directory & file: Samplesheet file name
# This function is dependent on function read_samplesheet
# Return a single data frame
read_count_files <- function(dir_in, file)
{
  count_data = read_samplesheet(dir_in,file) %>%
    dplyr::group_by_all() %>%
    do(read_tsv(here::here(dir_in, .$File[1]),col_names = c("ORF", 
                                                          "Count")))
}


# Get tidy data
# Change the data format
# There are two arguments, 
# dir_in: the input directory & file: Samplesheet file name
# This function is dependent on function read_count_files & read_samplesheet
get_tidy_data <- function(dir_in, file)
{
  data <- read_count_files(dir_in, file) %>%
    ungroup() %>%
    # Select the columns that we need
    select(Condition,Fraction,ORF,Count) %>%
    # Generate a wider data format to get the variables of fractions
    pivot_wider(names_from = Fraction,values_from = Count)
}




# Model fit
# There is only one argument, tidydata: the input data
# This function is dependent on function read_count_files & read_samplesheet & get_tidy_data
tidy_fit <- function(tidydata)
{
  # Use the function compose_data from package tidybayes to get the right format of data for model fitting
  # tidydata = get_tidy_data(read_count_files, dir_in, file)
  reformated_tidydata <- tidybayes::compose_data(tidydata)
  tidy_sampling <- sampling(compile_model, data = reformated_tidydata,
                            chains = 4, iter = 1000, 
                            control = list(adapt_delta = 0.85))
}



# Get the statistical results of parameters
# There is only one argument, tidydata: the input data
# This function is dependent on function read_count_files & read_samplesheet & get_tidy_data & tidy_fit
get_para_sta <- function(tidydata)
{
  para_sta <- tidy_fit(tidydata) %>% summary()
  # return medians
  data.frame(scaling.factor.Sup = para_sta$summary["scaling_factor_sup", 
                                                   "50%"], 
             scaling.factor.Pellet = para_sta$summary["scaling_factor_pellet",
                                                      "50%"], 
             lp.n_eff = para_sta$summary["lp__", 
                                         "n_eff"], 
             lp.Rhat = para_sta$summary["lp__", 
                                        "Rhat"])
}


# For different conditions
# There is only one argument, tidydata: the input data
# This function is dependent on all the five functions above 
paras_mean <- function(tidydata)
{
  plyr::ddply(tidydata, ~Condition, get_para_sta) # Apply function in each condition
}




# Calculate pSup for each transcript
# There is only one argument, tidydata: the input data
# This function is dependent on all the functions above 
calculate_pSup <- function(tidydata)
{
  paras_mean = paras_mean(tidydata)
 # tidydata <- get_tidy_data(read_count_files, dir_in, file)
  tidydata_pSup <- tidydata %>%
    dplyr::left_join(select(paras_mean,Condition,scaling.factor.Sup,scaling.factor.Pellet)) %>%
    mutate(pSup=scaling.factor.Sup*Sup/(scaling.factor.Sup*Sup+scaling.factor.Pellet*Pellet))
}


# Output pSup for each transcript
# Clear format
# There is only one argument, data: the input data
# This function is dependent on all the functions above 
each_mRNA_pSup <- function(tidydata)
{
  calculate_pSup(tidydata) %>%
    select(Condition,ORF,pSup) %>%
    pivot_wider(names_from = Condition,values_from = pSup)
}


# The last function to save history automatically
.Last <- function() {
  history.name <- paste(paste("AutoSave", Sys.Date(), strsplit(date(), " ")[[1]][4], sep="-"), ".Rhistory", sep="")
  savehistory(paste(getwd(), history.name))
}


