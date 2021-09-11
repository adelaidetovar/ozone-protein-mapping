source("0_functions.R")

#############################
#### About running miqtl ####
#############################

# Initial genome scans are quick to run in miqtl, and thus can be run on a desktop
# Running threshold scans and bootstrap scans for locus intervals require more memory and thus take more time
# I usually perform these processes with cluster computing

#####################################
#### ROP with protein difference ####
#####################################

dir.create("../output/mapping")
proteinPairs = read.table("../output/protein_pairs.tsv", header = TRUE)

# ROP scan with harvest/batch as covariate (explicit batch correction), using replicates (reflected with pheno.id argument)
difference_rop <- scan.h2lmm(genomecache=cache_path, data=proteinPairs, formula=protein.delta ~ 1 + harvest + sex, 
                             model="additive", use.multi.impute = FALSE, locus.as.fixed=TRUE, return.allele.effects=TRUE, 
                             pheno.id = "mouse_id", geno.id="strain")

saveRDS(difference_rop,"../output/mapping/difference_rop.rds")

# View scan to see where top peaks are - output to viewer
genome.plotter.whole(list(difference_rop))

# Create permutation matrix for threshold scans
difference_rop_perms <- generate.sample.outcomes.matrix(scan.object=difference_rop,
                                                        method="permutation", num.samples = 1000)

saveRDS(difference_rop_perms,"../output/mapping/difference_rop_perms.rds")

# Threshold scans
difference_rop_perm_scans <- run.threshold.scans(sim.threshold.object = difference_rop_perms, keep.full.scans = TRUE,
                                                 use.multi.impute = FALSE, data = proteinPairs, genomecache = cache_path)

saveRDS(difference_rop_perm_scans,"../output/mapping/difference_rop_perm_scans.rds")

# Grab top marker for chr 10, 15 locus
difference_ten <- grab.locus.from.scan(difference_rop, chr=10, return.value="marker")
difference_fift <- grab.locus.from.scan(difference_rop, chr=15, return.value="marker")

#---------------

# Chr 15 locus

# Subsample for top locus (Chr 15)
difference_rop_fift <- scan.h2lmm(genomecache=cache_path, data=difference, formula=delta ~ 1 + harvest + sex, model="additive",
                                  use.multi.impute = FALSE, locus.as.fixed=TRUE, return.allele.effects=TRUE,
                                  just.these.loci = difference_fift, pheno.id = "mouse_id", geno.id="strain")

saveRDS(difference_rop_fift, "../output/mapping/difference_rop_fift.rds")

# Generate bootstraps from subsampled scan
difference_bootstraps_fift <- generate.sample.outcomes.matrix(scan.object=difference_rop_fift,model.type="alt",
                                                              method="bootstrap",subsample.chr=15,num.samples=1000)

saveRDS(difference_bootstraps_fift, "../output/mapping/difference_bootstraps_fift.rds")

# Positional scans to define confidence interval
difference_positional_fift <- run.positional.scans(difference_bootstraps_fift,keep.full.scans=TRUE,use.multi.impute=FALSE,
                                                   data=difference,genomecache=cache_path)

saveRDS(difference_positional_fift, "../output/mapping/difference_positional_fift.rds")

#---------------

# Chr 10 Locus

# Subsample for top locus (Chr 10)
difference_rop_ten <- scan.h2lmm(genomecache=cache_path, data=difference, formula=delta ~ 1 + harvest + sex, model="additive",
                                 use.multi.impute = FALSE, locus.as.fixed=TRUE, return.allele.effects=TRUE,
                                 just.these.loci = difference_ten, pheno.id = "mouse_id", geno.id="strain")

saveRDS(difference_rop_ten,"../output/mapping/difference_rop_ten.rds")

# Generate bootstraps from subsampled scan
difference_bootstraps_ten <- generate.sample.outcomes.matrix(scan.object=difference_rop_tenScan,model.type="alt",
                                                             method="bootstrap",subsample.chr=10,num.samples=1000)

saveRDS(difference_bootstraps_ten,"../output/mapping/difference_bootstraps_ten.rds")

# Positional scans to define confidence interval
difference_positional_ten <- run.positional.scans(difference_bootstraps_ten,keep.full.scans=TRUE,use.multi.impute=FALSE,
                                                  data=difference,genomecache=cache_path)

saveRDS(difference_positional_ten,"../output/mapping/difference_positional_ten.rds")

###########################################
#### Effect size for Chr 10 and 15 QTL ####
###########################################

# Get founder haplotypes at Chr 10 marker
x = load(paste0(cache_path,'additive/chr10/data/',difference_ten,".Rdata"))
chr10_founder <- as.data.frame(get(x))
chr10_founder <- chr10_founder %>%
  rownames_to_column('strain') %>%
  gather(haplotype, cnt, A:H) %>%
  group_by(strain) %>%
  filter(cnt == max(cnt)) %>%
  arrange(strain) %>%
  rename(chr10_haplo = haplotype)

# Get founder haplotypes at Chr 15 marker
x = load(paste0(cache_path,'additive/chr15/data/',difference_fift,".Rdata"))
chr15_founder <- as.data.frame(get(x))
chr15_founder <- chr15_founder %>%
  rownames_to_column('strain') %>%
  gather(haplotype, cnt, A:H) %>%
  group_by(strain) %>%
  filter(cnt == max(cnt)) %>%
  arrange(strain) %>%
  rename(chr15_haplo = haplotype)

protein_QTL <- left_join(chr15_founder, proteinPairs)

## Effect size for Chr 10 locus

af <- anova(lm(protein.delta ~ harvest + sex + chr10_haplo, data = protein_QTL))
afss <- af$"Sum Sq"
sink("../output/mapping/chr10_effectsize.txt")
print(cbind(af, PctExp = afss/sum(afss)*100))
sink()

## Effect size for Chr 15 locus

af <- anova(lm(protein.delta ~ harvest + sex + chr15_haplo, data = protein_QTL))
afss <- af$"Sum Sq"
sink("../output/mapping/chr15_effectsize.txt")
print(cbind(af, PctExpt = afss/sum(afss)*100))
sink()

###############
#### TIMBR ####
###############

# Eliminate CC082 -- no info in locus matrix
timbr_df <- proteinPairs[proteinPairs$strain != "CC082",]

# Pull base data from TIMBR for model settings
data("mcv.data")
## TIMBR O3/difference (same marker, diff y)
prior.M <- list(model.type = "crp", # crp - Chinese Restaurant Process
                prior.alpha.type = "gamma", 
                prior.alpha.shape = 1, 
                prior.alpha.rate = 2.333415)

# Order based on strain
timbr_df <- timbr_df %>%
  arrange(strain)
timbr_df$strain <- as.character(timbr_df$strain)

# Gather number of replicates for each strain
W <- as.numeric(as.vector(unname(table(timbr_df$strain))))

timbr_df$strain <- as.factor(timbr_df$strain)
timbr_strains <- levels(timbr_df$strain)

# Get haplotype probs for peak marker (Chr 15)
prior.D <- list()
prior.D$P <- load(paste0(cache_path,'full/chr15/data/',difference_fift,".Rdata"))
prior.D$P <- get(prior.D$P)
prior.D$P <- prior.D$P[timbr_strains,]
prior.D$P <- as.data.frame(prior.D$P) %>% uncount(W)
rownames(prior.D$P) <- colnames(prior.D$P) <- c()
prior.D$P <- as.matrix(prior.D$P)
dimnames(prior.D$P) <- NULL
prior.D$A <- mcv.data$prior.D$A
prior.D$fixed.diplo <- FALSE

# Make dataframe for fitting
df <- timbr_df
# Regress out batch and sex effects
df$delta.resids <- resid(lm(protein.delta ~ 1 + harvest + sex, data=df))

timbr_chr15 <- list(prior.D = prior.D,
                    y = df$delta.resids)

str(timbr_chr15$prior.D)

results <- TIMBR(timbr_chr15$y, timbr_chr15$prior.D, prior.M, verbose = TRUE)

# Report posterior probabilities for the top allelic series models
head(results$p.M.given.y)

# Posterior density of number of alleles
results$p.K.given.y

save(results,file='../output/mapping/chr15_timbr.rds')

####################
#### Diploffect ####
####################

## Chr 15
# Load marker locus matrix
locus_matrix = load(paste0(cache_path,'full/chr15/data/',difference_fift,".Rdata"))
locus_matrix <- get(locus_matrix)

# Repurpose elements of TIMBR

locus_matrix <- locus_matrix[timbr_strains,]
locus_matrix <- as.data.frame(locus_matrix) %>% uncount(W)

diplo_df <- df
colnames(diplo_df)[7] <- "SUBJECT.NAME"
diplo_df$SUBJECT.NAME <- rownames(locus_matrix)

inla.diploffect <- run.diploffect.inla(formula=delta.resids~1+(1|strain), add.on=FALSE,                     
                                       data=diplo_df, K=NULL, prob.matrix=as.matrix(locus_matrix),
                                       num.draws=10, use.dip.lincomb=TRUE, seed=1, 
                                       gamma.rate=1, impute.on="strain")

inla.diploffect.summary <- run.diploffect.inla.summary.stats(inla.diploffect)

save(inla.diploffect.summary, file='../output/mapping/chr15_diploffect.rds')
