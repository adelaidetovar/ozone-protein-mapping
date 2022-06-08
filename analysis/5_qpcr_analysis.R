source("0_functions.R")

dir.create("./qPCR")
setwd("../data")

###############################
#### qPCR Sample Selection ####
###############################

# protein_mean <- read.table("../output/protein_data.tsv")
# protein_mean <- merge(protein_mean, chr15_founder, by = "strain")
# 
# strain_haplo <- protein_mean[,c("strain","haplotype")]
# strain_haplo <- strain_haplo[!duplicated(strain_haplo),]
# 
# set.seed(32)
# random_strain <- strain_haplo %>%
#   # select from non-B6 or CAST haplotypes, ie not B or F
#   filter(haplotype %in% c("A", "C", "D", "E", "G", "H")) %>%
#   group_by(haplotype) %>%
#   do(sample_n(.,2)) 
# 
# # select all B6 strains
# pcr_strain <- strain_haplo %>%
#   filter(haplotype %in% c("B"))
# 
# # combine randomly selected and B6 strains, then add CC071 (highest responder)
# pcr_strain <- bind_rows(pcr_strain, random_strain)
# pcr_strain[nrow(pcr_strain) + 1,] = c("CC071", "F")
# 
# #### With strains, select random samples ####
# 
# protein_data <- read.table("../output/protein_pairs.tsv")
# protein_pairs <- proteinPairs %>%
#   group_by(strain) %>%
#   count()
# 
# # figure out which strain has fewer than 2 pairs
# protein_pairs <- protein_pairs %>%
#   filter(strain %in% pcr_strain$strain) %>%
#   mutate(size = if_else(n >= 2, 2, 1))
# ## it's CC010, exclude from sampling b/c we'll take all samples
# 
# pcr_strain <- pcr_strain$strain[pcr_strain$strain!="CC010"]
# 
# pcr_sample <- protein_data %>%
#   filter(strain %in% pcr_strain) %>%
#   group_by(strain, rx, harvest) %>%
#   do(sample_n(., 1)) %>%
#   group_by(strain, harvest) %>%
#   add_count(harvest) %>%
#   filter(n == 2) %>%
#   group_by(strain) %>%
#   arrange(strain, harvest) %>%
#   filter(row_number() <= 2) %>%
#   bind_rows(., protein_data[protein_data$strain=="CC010",])
# 
# write.xlsx(pcr_sample[,c(1:5)], "pcr_samples.xlsx", overwrite = TRUE)

############################
#### qPCR Data Analysis ####
############################

qpcr_samples <- read.xlsx("qpcr_raw-data.xlsx", sheet = 1)
qpcr_a <- read.xlsx("qpcr_raw-data.xlsx", sheet = 2) # Angpt1
qpcr_b <- read.xlsx("qpcr_raw-data.xlsx", sheet = 3) # Oxr1, both sets
qpcr_c <- read.xlsx("qpcr_raw-data.xlsx", sheet = 4) # Rspo2

# keep select columns
qpcr_a <- qpcr_a[,c(3, 5, 7)]
qpcr_b <- qpcr_b[,c(3, 5, 7)]
qpcr_c <- qpcr_c[,c(3, 5, 7)]

# fill NAs with 0
qpcr_a$Cq[is.na(qpcr_a$Cq)] <- 0
qpcr_b$Cq[is.na(qpcr_b$Cq)] <- 0
qpcr_c$Cq[is.na(qpcr_c$Cq)] <- 0

# rename columns
colnames(qpcr_a) <- colnames(qpcr_b) <- colnames(qpcr_c) <- c("gene", "sample_id", "Cq")
colnames(qpcr_samples)[1] <- "sample_id"
qpcr_samples$sample_id <- as.character(qpcr_samples$sample_id)

qpcr_a$sample_id <- as.character(qpcr_a$sample_id)
qpcr_b$sample_id <- as.character(qpcr_b$sample_id)
qpcr_c$sample_id <- as.character(qpcr_c$sample_id)

qpcr_a$plate <- rep("plate_1", length(qpcr_a$gene))
qpcr_b$plate <- rep("plate_2", length(qpcr_b$gene))
qpcr_c$plate <- rep("plate_3", length(qpcr_c$gene))

# calculate mean cq by sample
cq_a <- mean_cq(cq_a)
cq_a <- cq_a[,c(1:4)]
colnames(cq_a)[4] <- c("Rps20_Plate1")

cq_b <- mean_cq(cq_b)
cq_b <- cq_b[,c(1:4)]
colnames(cq_b)[4] <- c("Rps20_Plate2")

cq_c <- mean_cq(cq_c)
cq_c <- cq_c[,c(1:3)]
colnames(cq_c)[3] <- c("Rps20_Plate3")

all_cq <- merge(cq_a, cq_b, by = "sample_id")
all_cq <- merge(all_cq, cq_c, by = "sample_id")
all_cq <- merge(all_cq, qpcr_samples, by = "sample_id")
all_cq <- all_cq[order(as.numeric(all_cq$sample_id)),]

# Check between plate correlation for housekeeping genes

# Plate 1 and Plate 2
cor(all_cq$Rps20_Plate1, all_cq$Rps20_Plate2)
# [1] 0.8706175
cor(all_cq$Rps20_Plate1, all_cq$Rps20_Plate2, method = "spearman")
# [1] 0.8665062

# Plate 2 and Plate 3
cor(all_cq$Rps20_Plate2, all_cq$Rps20_Plate3)
# [1] 0.9286525
cor(all_cq$Rps20_Plate2, all_cq$Rps20_Plate3, method = "spearman")
# [1] 0.8846701

# Plate 1 and Plate 3
cor(all_cq$Rps20_Plate1, all_cq$Rps20_Plate3)
# [1] 0.913412
cor(all_cq$Rps20_Plate1, all_cq$Rps20_Plate3, method = "spearman")
# [1] 0.9002079

# calculate dCq by gene
cq_a[cq_a == 0] <- NA
cq_b[cq_b == 0] <- NA
cq_c[cq_c == 0] <- NA

dcq_a <- cq_a %>%
  mutate(dCq_Angpt1_Rps20 = Angpt1 - Rps20_Plate1)

dcq_b <- cq_b %>%
  mutate(dCq_Rspo2_Rps20 = Rspo2 - Rps20_Plate2,
         dCq_Oxr1_Set1_Rps20 = Oxr1_Set1 - Rps20_Plate2) %>%
  filter(sample_id %in% seq(1,38,by=1))

dcq_c <- cq_cc %>%
  mutate(dCq_Oxr1_Set2_Rps20 = Oxr1_Set2 - Rps20_Plate3) %>%
  filter(sample_id %in% seq(1,38,by=1))

# created merged dCq across all
dcq_full <- merge(dcq_a[,c(1,5,6)], dcq_b[,c(1,5,6)], by = "sample_id")
dcq_full <- merge(dcq_full, dcq_c[,c(1,4)], by = "sample_id")
dcq_full <- merge(dcq_full, qpcr_samples, by = "sample_id")

# merge to add haplotype
dcq_full <- merge(dcq_full, chr15_founder, by = "strain")

# add dichotomous grouping
dcq_full <- dcq_full %>%
  mutate(dichot = case_when(haplotype == "B"| haplotype == "F" ~ "group_a",
                            TRUE ~ "group_b"))

#---------------

# Perform regression analyses
dir.create("../output/qPCR_ANOVA_tables")

# Rspo2 with strain term and interaction
sink(file = "../output/qPCR_ANOVA_tables/rspo2_dichot-rx-strain-interact.txt")
anova(lm(dCq_Rspo2_Rps20 ~ dichot*rx + strain, data = dcq_full))
sink()

# Angpt1 with strain term and interaction
sink(file = "../output/qPCR_ANOVA_tables/Angpt1_dichot-rx-strain-interact.txt")
anova(lm(dCq_Angpt1_Rps20 ~ dichot*rx + strain, data = dcq_full))
sink()

# Oxr1_Set1 with strain term and interaction
sink(file = "../output/qPCR_ANOVA_tables/Oxr1_Set1_dichot-rx-strain-interact.txt")
anova(lm(dCq_Oxr1_Set1_Rps20 ~ dichot*rx + strain, data = dcq_full))
sink()

# Oxr1_Set2 with strain term and interaction
sink(file = "../output/qPCR_ANOVA_tables/Oxr1_Set2_dichot-rx-strain-interact.txt")
anova(lm(dCq_Oxr1_Set2_Rps20 ~ dichot*rx + strain, data = dcq_full))
sink()

#---------------

# Plots
## Rspo2
rspo2_interact_plot <- dcq_interact_plot(dcq_full, dCq_Rspo2_Rps20, ylab = "Rspo2 dCq (Rps20) * -1", 3.8, 4.02)
rspo2_rx_plot <- dcq_rx_plot(dcq_full, dCq_Rspo2_Rps20, ylab = "Rspo2 dCq (Rps20) * -1")
rspo2_dichot_plot <- dcq_dichot_plot(dcq_full, dCq_Rspo2_Rps20, ylab = "Rspo2 dCq (Rps20) * -1")

rspo2_arranged <- grid.arrange(arrangeGrob(rspo2_dichot_plot, rspo2_rx_plot, nrow=2), rspo2_interact_plot, ncol = 2,
                               widths = c(1,1.5))

ggsave(rspo2_arranged, filename="../output/plots/rspo2_arranged.png", dpi = 300,
       units = "in", width = 10, height = 6)

## Angpt1
angpt1_interact_plot <- dcq_interact_plot(dcq_full, dCq_Angpt1_Rps20, ylab = "Angpt1 dCq (Rps20) * -1", 4, 4.2)
angpt1_rx_plot <- dcq_rx_plot(dcq_full, dCq_Angpt1_Rps20, ylab = "Angpt1 dCq (Rps20) * -1")
angpt1_dichot_plot <- dcq_dichot_plot(dcq_full, dCq_Angpt1_Rps20, ylab = "Angpt1 dCq (Rps20) * -1")

angpt1_arranged <- grid.arrange(arrangeGrob(angpt1_dichot_plot, angpt1_rx_plot, nrow=2), angpt1_interact_plot, ncol = 2,
                               widths = c(1,1.5))

ggsave(angpt1_arranged, filename="../output/plots/angpt1_arranged.png", dpi = 300,
       units = "in", width = 10, height = 6)

## Oxr1_Set1
oxr1_set1_interact_plot <- dcq_interact_plot(dcq_full, dCq_Oxr1_Set1_Rps20, ylab = "Oxr1-Set1 dCq (Rps20) * -1", 4.9, 5.13)
oxr1_set1_rx_plot <- dcq_rx_plot(dcq_full, dCq_Oxr1_Set1_Rps20, ylab = "Oxr1-Set1 dCq (Rps20) * -1")
oxr1_set1_dichot_plot <- dcq_dichot_plot(dcq_full, dCq_Oxr1_Set1_Rps20, ylab = "Oxr1-Set1 dCq (Rps20) * -1")

oxr1_set1_arranged <- grid.arrange(arrangeGrob(oxr1_set1_dichot_plot, oxr1_set1_rx_plot, nrow=2), oxr1_set1_interact_plot, ncol = 2,
                                widths = c(1,1.5))

ggsave(oxr1_set1_arranged, filename="../output/plots/oxr1_set1_arranged.png", dpi = 300,
       units = "in", width = 10, height = 6)


## Oxr1_Set2
oxr1_set2_interact_plot <- dcq_interact_plot(dcq_full, dCq_Oxr1_Set2_Rps20, ylab = "Oxr1-Set2 dCq (Rps20) * -1", 4.03, 4.23)
oxr1_set2_rx_plot <- dcq_rx_plot(dcq_full, dCq_Oxr1_Set2_Rps20, ylab = "Oxr1-Set2 dCq (Rps20) * -1")
oxr1_set2_dichot_plot <- dcq_dichot_plot(dcq_full, dCq_Oxr1_Set2_Rps20, ylab = "Oxr1-Set2 dCq (Rps20) * -1")

oxr1_set2_arranged <- grid.arrange(arrangeGrob(oxr1_set2_dichot_plot, oxr1_set2_rx_plot, nrow=2), oxr1_set2_interact_plot, ncol = 2,
                                   widths = c(1,1.5))

ggsave(oxr1_set2_arranged, filename="../output/plots/oxr1_set2_arranged.png", dpi = 300,
       units = "in", width = 10, height = 6)