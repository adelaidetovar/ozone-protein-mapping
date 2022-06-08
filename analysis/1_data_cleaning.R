source("0_functions.R")

setwd("../data")

#################################
#### Data Input and Cleaning ####
#################################

dir.create("../output")

# BAL data
bal_data <- read.xlsx('raw_data.xlsx',sheet=1)
bal_data$batch <- as.factor(bal_data$batch)
colnames(bal_data)[5] <- "harvest"
colnames(bal_data)[4] <- "rx"
write.table(bal_data, "../output/bal_data.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#---------------

# Protein data
protein_data <- read.xlsx('raw_data.xlsx',sheet=2)
protein_data$harvest <- as.factor(protein_data$harvest)
protein_data$rx <- as.factor(protein_data$rx)
protein_data$protein.concentration <- as.numeric(protein_data$protein.concentration)

# Manually exclude value for 5202, missing its match
protein_data$protein.concentration[protein_data$mouse_id == "5202"] <- NA

# Manually recode harvest IDs for strains included more than once in a batch
## CC021 in batch 17 twice
protein_data$harvest_update <- as.character(protein_data$harvest)

protein_data$harvest_update[protein_data$mouse_id %in% c("5590","5591","5592")] <- "17A"
protein_data$harvest_update[protein_data$mouse_id %in% c("5631","5632","5633")] <- "17B"

## CC072 in batch 15 twice
protein_data$harvest_update[protein_data$mouse_id %in% c("5566","5567")] <- "15A"
protein_data$harvest_update[protein_data$mouse_id %in% c("5568","5569")] <- "15B"
write.table(protein_data, "../output/protein_data.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#---------------

# Delta for mapping
proteinPairs <- protein_data %>%
  filter(!is.na(protein.concentration)) %>%
  group_by(strain, harvest_update, sex, rx) %>%
  mutate(protein.mean = mean(protein.concentration)) %>%
  ungroup() %>%
  group_by(strain, harvest_update, sex) %>%
  mutate(protein.delta = first(protein.mean[rx == "O3"]) - first(protein.mean[rx == "FA"])) %>%
  ungroup() %>%
  group_by(strain, harvest, harvest_update, sex) %>%
  dplyr::select(strain, harvest, harvest_update, sex, protein.delta) %>%
  filter(!is.na(protein.delta)) %>%
  distinct(protein.delta)

# Recode multiple harvests to merged
proteinPairs <- as.data.frame(proteinPairs)
proteinPairs$harvest <- proteinPairs$harvest_update
proteinPairs$harvest[proteinPairs$harvest %in% c("17A","17B")] <- "17"
proteinPairs$harvest[proteinPairs$harvest %in% c("15A","15B")] <- "15"
proteinPairs$harvest <- as.factor(proteinPairs$harvest)

# Round up difference, create pseudo-mouse IDs for pairs
proteinPairs <- proteinPairs[complete.cases(proteinPairs),]
proteinPairs$protein.delta <- round(proteinPairs$protein.delta, digits = 3)
proteinPairs$mouse_id <- rep("M",length(proteinPairs$protein.delta))
proteinPairs$mouse_id <- paste0(proteinPairs$mouse_id,seq(1,length(proteinPairs$strain)))

write.table(proteinPairs, "../output/protein_pairs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#---------------

## FA only data for heritability
fa_protein <- protein_data[protein_data$rx == "FA",]

write.table(fa_protein, "../output/fa_protein.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#---------------

## O3 only data for heritability
O3_protein <- protein_data[protein_data$rx == "O3",]

write.table(O3_protein, "../output/O3_protein.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
