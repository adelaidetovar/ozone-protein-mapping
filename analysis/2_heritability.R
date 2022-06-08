source("0_functions.R")
dir.create("../output/pheno_anova_tables")
dir.create("../output/heritability/")

# Normalization factors determined via Box-Cox

###############################
#### ANOVA for phenotypes #####
###############################

# BAL neutrophils (%)
sink(file = "../output/pheno_anova_tables/percent_neut_anova.txt")
anova(lm((p.pmn)^(1/3) ~ strain + rx + sex, data = bal_data))
sink()

sink(file = "../output/pheno_anova_tables/percent_neut_anova_interact.txt")
anova(lm((p.pmn)^(1/3) ~ strain + rx + sex + strain*rx, data = bal_data))
sink()

#---------------

# BAL neutrophils (#)
sink(file = "../output/pheno_anova_tables/number_neut_anova.txt")
anova(lm((log10((p.pmn*total/10)+1)) ~ strain + rx + sex, data = bal_data))
sink()

sink(file = "../output/pheno_anova_tables/number_neut_anova_interact.txt")
anova(lm((log10((p.pmn*total/10)+1)) ~ strain + rx + sex + strain*rx, data = bal_data))
sink()

#---------------

# BAL protein
sink(file = "../output/pheno_anova_tables/protein_anova.txt")
anova(lm((protein.concentration)^(1/3) ~ strain + rx + sex, data = protein_data))
sink()

sink(file = "../output/pheno_anova_tables/protein_anova_interact.txt")
anova(lm((protein.concentration)^(1/3) ~ strain + rx + sex + strain*rx, data = protein_data))
sink()

###############################
#### Heritability via INLA ####
###############################

bal_data = read.table("../output/bal_data.tsv", header = TRUE)

# BAL O3 values (% neutrophils)
K <- read.delim('../data/neutrophil_ozone_kinship.txt', header=T, sep="\t", check.names=FALSE, row.names=1)
K[is.na(K)] <- 0

traits <- c("p.pmn.norm", "t.pmn.norm")
bal_data_O3 <- bal_data[bal_data$rx=="O3",]
bal_data_O3 <- bal_data_O3[complete.cases(bal_data_O3),]
bal_data_O3$p.pmn.norm <- bal_data_O3$p.pmn^(1/3)
bal_data_O3$t.pmn <- bal_data_O3$p.pmn*bal_data_O3$total/10
bal_data_O3$t.pmn.norm <- log10(bal_data_O3$t.pmn+1)

new.K <- K + diag(rep(0.001, nrow(K)))

heritability(bal_data_O3, traits = c("p.pmn.norm", "t.pmn.norm"))

## Return summary of results
sink(file = "../output/heritability/percent_neut_heritability.txt")
print(summary(inla.results[[1]]))
sink()

sink(file = "../output/heritability/number_neut_heritability.txt")
print(summary(inla.results[[2]]))
sink()

#---------------

protein_data = read.table("../output/protein_data.tsv", header = TRUE)

######## BAL protein (FA only)

## BAL FA protein values
protein_data_FA <- protein_data[protein_data$rx=="FA",]
protein_data_FA$protein.concentration.norm <- protein_data_FA$protein.concentration^(1/3)
K <- read.delim('../data/protein_FA_kinship.txt', header=T, sep="\t", check.names=FALSE, row.names=1)
K[is.na(K)] <- 0

new.K <- K + diag(rep(0.001, nrow(K)))

heritability(protein_data_FA, "protein.concentration.norm")

# Return summary of results
# BAL protein for FA
sink(file = "../output/heritability/protein_FA_heritability.txt")
print(summary(inla.results))
sink()

#---------------

######## BAL protein (O3)
protein_data_O3 <- protein_data[protein_data$rx=="O3",]
K <- read.delim('../data/protein_ozone_kinship.txt', header=T, sep="\t", check.names=FALSE, row.names=1)
K[is.na(K)] <- 0

new.K <- K + diag(rep(0.001, nrow(K)))

heritability(protein_data_O3, "protein.concentration")

# Return summary of results
# BAL protein for O3
sink(file = "../output/heritability/protein_O3_heritability.txt")
print(summary(inla.results))
sink()

#---------------

proteinPairs = read.table("../output/protein_pairs.tsv", header = TRUE)

######## BAL protein (difference)
K <- read.delim('../data/protein_pairs_kinship.txt', header=T, sep="\t", check.names=FALSE, row.names=1)
K[is.na(K)] <- 0

new.K <- K + diag(rep(0.001, nrow(K)))

heritability(proteinPairs, "protein.delta")

# Return summary of results
# BAL protein delta/difference
sink(file = "../output/heritability/protein_difference_heritability.txt")
print(summary(inla.results))
sink()

#----------------

# Misc calc

summ_bal <- bal_data %>%
  group_by(strain, rx) %>%
  count()

mean(summ_bal$n[summ_bal$rx=="FA"])
mean(summ_bal$n[summ_bal$rx=="O3"])
