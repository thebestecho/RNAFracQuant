# This script includes the definitions of all the functions
# The starter function
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

# Read the original mRNA count data from the files provided by users
# The file type is restricted to CSV only
# Return a single data frame
read_samplesheet <- function(dir_in)
{
  # load data from the samplesheet file
  # I dont know whether we need to change the ORFselect file later if deal with other experiments
  # Ask users if their experiments are about yeasts
  # Maybe I can provide other choices later, e.g. prepare different reference files as ORFselect
  while (TRUE) {
    input = readline('Are your data files from yeast?\n Please type Y for yes and N for no. ')
    if (input == "Y") {
      samplesheet = paste0(dir_in, "/Samplesheet.txt")
      tryCatch(read_tsv(here::here(dir_in, "Samplesheet.txt"), comment = "#"), error = function(c) {
        c$message <- paste0(" Error ",c$message, " (error in ", samplesheet, ")") # point out where the error is from
        stop(c)
      })
      load_samplesheet = read_tsv(here::here(dir_in, "Samplesheet.txt"), comment = "#")
      return(load_samplesheet)
    } else if (input == "N"){
      return(message("Sorry, you data are not eligible for using this model."))
    } else {
      message("Warning: ")
    }
  }
}
  


# Use the file names that is read from Samplesheet.txt to read all the count files
read_count_files <- function(dir_in)
{# samplesheet <- read_samplesheet(dir_in=input)
  group = samplesheet %>%
                     group_by(File, Sample, Condition, Fraction)
  if (sum(is.na(group)) > 0) {
    stop("NA values found! Require complete values of File, Sample, Condition and Fraction in file Samplesheet.txt.")
  } else if (sum(complete.cases(group$File)) < 2) {
    stop("Require at least two values of File in file Samplesheet.txt.")
  } else if (sum(complete.cases(group$Sample)) < 2) {
    stop("Require at least two values of Sample in file Samplesheet.txt.")
  } else if (sum(complete.cases(group$Condition)) < 2) {
    stop("Require at least two values of Condition in file Samplesheet.txt.")
  } else if (sum(complete.cases(group$Fraction)) < 2) {
    stop("Require at least two values of Fraction in file Samplesheet.txt.")
  } else {
    count_data = group %>%
      do(read_tsv(here::here(dir_in, .$File[1]),col_names = c("ORF", 
                                                        "Count"), comment = "__"))
    if (ncol(count_data) != 6) {
      stop("Require ORF and its Count value in all data files.")
    } else if (typeof(count_data$ORF) != "character"){
      stop("Require ORF as character input in all data files.")
    } else if (typeof(count_data$Count) != "double"){
      stop("Require Count values for each ORF in all data files.")
    } else {
      count_data = mutate(count_data,ORF = str_sub(ORF)) # Used be end = -5, need to consider this later
      return(count_data)
    }
  }
}


# Calculate the TPM value and add it to the data frame
get_TPM <- function(dir_in){
  # load data from the ORFselect file
  # I dont know whether we need to change the file here if deal with other experiments
  # Maybe I can provide other choices later
  while(TRUE){
    choice = readline('Please choose the organism your data are from: \n Please type A for yeast and B for others. ')
    if (choice == "A") {
      ORFselect <- read_tsv(here::here(dir_in,"Scer_ORF_length.txt")) 
      if (sum(is.na(ORFselect)) > 0) {
        stop("NA values found! Require complete values of ORF and Length in file orf_coding_EW_select.txt.")
      } else {
        getcount = get_count %>% inner_join(ORFselect)
        if (ncol(getcount) != 7) {
          stop("OOps! Failed to join ORFselect data to the raw count data.")
        } else {
          data_with_TPM = getcount %>%
            group_by(Sample,Condition,Fraction) %>%
            dplyr::mutate(TPM=(Count/Length)/sum(Count/Length)*1e6) 
          if (sum(is.na(data_with_TPM)) > 0) {
            warning("NA values calculated in TPM!")
          } else {
            data_TPM = select(data_with_TPM,Sample,Condition,Fraction,ORF,Count,TPM)
            return(data_TPM)
          }
        }
      }
    } else if (choice == "B"){
      return(message("Dont worry! I will load other ORFselect later to match the choice later\n (not completed yet)."))
    } else {
      message("Warning: ")
    }
  }
  }


# Get the ORF list with high TPM from the data 
admissible_ORF <- function(dir_in){
# Get the fraction of Total count
# count_with_TPM <- get_TPM(dir_in=dir_in)
  get_Tot <- count_with_TPM %>%
    dplyr::filter(Fraction=="Tot")
  if (nrow(get_Tot) == 0){
    stop("Fraction Tot is not found. Please make sure your columns are listed as Tot,Sup,P100 in file Samplesheet.txt.")
  } else {
# samplesheet <- read_samplesheet(dir_in=dir_in)
    n_condition = length(unique(samplesheet$Condition)) # Get the number of conditions from samplesheet
    tryCatch(samplesheet$Condition, warning = function(a) {
      a$message <- paste0(" Warning! ",a$message, " (error in file Samplesheet.txt)") # point out where the error is from
      stop(a)
    })
    True_ORF_list = get_Tot %>% 
      group_by(ORF) %>%
# Filter by median tpm > 5, min tpm > 0, across all conditions.
# This is to remove low-abundance transcripts
      dplyr::summarise(n_good_condition = sum(TPM > 5),ORF_list=(n_good_condition == n_condition)) %>%
      dplyr::filter(ORF_list) %>%
      .$ORF
    return(True_ORF_list)
  }
}



# Filter data using the ORF list with high TPM
# Get tidy data
get_tidy_data <- function(dir_in){
# count_with_TPM <- get_TPM(dir_in=dir_in)
  tidy_data <- count_with_TPM %>%
    ungroup() %>%
# Select the columns that we need
    select(Condition,Fraction,ORF,Count) %>%
# Generate a wider data format to get the variables of fractions
    pivot_wider(names_from = Fraction,values_from = Count) %>%
# Remove low-abundance transcripts to reduce the time to fit
# admissible_ORF_list <- admissible_ORF(dir_in=dir_in)
    filter(ORF %in% admissible_ORF_list)
  return(tidy_data)
}


# Model fit
tidy_fit <- function(tidydata, chains = 4, iter = 1000, 
                     control = list(adapt_delta = 0.85), seed = 39, 
                     ...)
{
# Use the function compose_data from package tidybayes to get the right format of data for model fitting
  reformated_tidydata <- tidybayes::compose_data(tidydata)
  tidy_sampling <- sampling(compile_model, data = reformated_tidydata, 
                            iter = iter, control = control, seed = seed, chains = chains,
                            ...)
}


# Get the statistical results of parameters
get_para_sta <- function(tidydata)
{
  para_sta <- tidy_fit(tidydata = tidydata) %>% summary()
  # return medians
  data.frame(mixing.Sup = para_sta$summary["mixing_sup", 
                                                "50%"], 
             mixing.P100 = para_sta$summary["mixing_p100",
                                                 "50%"], 
             lp.n_eff = para_sta$summary["lp__", 
                                              "n_eff"], 
             lp.Rhat = para_sta$summary["lp__", 
                                             "Rhat"])
}


# The last function to save history automatically
.Last <- function() {
  history.name <- paste(paste("AutoSave", Sys.Date(), strsplit(date(), " ")[[1]][4], sep="-"), ".Rhistory", sep="")
  savehistory(paste(getwd(), history.name))
}
