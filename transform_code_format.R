# Read the raw data
raw_data = readr::read_tsv("./GSE141029_counts_matrix.txt")
# Replace the column names and covert it into longer format
wide_published_data <- tidyr::pivot_longer ( raw_data,
cols = "IDR_30C_1_pellet_262":"WT_42C_2_total_219",
names_to = c (“Condition","Replicate","Fraction","File"), 
names_pattern = "(. * _.*)_(.*)_(.*)_(.*)",
values_to = "Count") %>% 
  dplyr::select(-File) %>% 
  tidyr::pivot_wider(names_from = "Fraction",values_from = "Count")

colnames(wide_published_data) [c (1,4:6)] <- c("ORF","Pellet","Sup","Tot")
# Filter data
ORFselect = readr::read_tsv("./original_Scer_ORF_length.txt")
wide_published_filtered <- wide_published_data %>%
  tidyr::unite(Condition,Condition,Replicate,sep = "_") %>%
  dplyr::filter(ORF %in% ORFselect$ORF, Pellet > 0, Sup > 0, Tot > 0)
# Get pSup values
published_data_pSup <- each_mRNA_pSup (wide_published_filtered)
