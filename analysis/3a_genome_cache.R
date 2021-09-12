source("0_functions.R")

######################################
#### Genome Cache (ONLY RUN ONCE) ####
######################################

# Download MegaMUGA SNPs and MRCAs from csbio.unc.edu

# Load in CC MRCAs
CC.paths <- list.files("/path/to/strain_mrcas/", full.names=TRUE)
CC.names <- gsub(x=list.files("/path/to/strain_mrcas/", full.names=FALSE), 
                 pattern="\\_.+$", replacement="", perl=TRUE)

CC.1 <- read.csv(CC.paths[1], header=TRUE)

chr <- CC.1$chromosome

chr.keep <- c(1:19, "X")
## Keep loci in the standard mouse chromosomes and that have probabilities
keep <- chr %in% chr.keep & rowSums(CC.1[,-(1:3)]) == 1

final.chr <- chr[keep]
final.loci <- CC.1$marker[keep]
final.Mb <- CC.1$position.B37.[keep]

result.array <- array(NA, dim=c(length(CC.paths), ncol(CC.1) - 3, length(final.loci)),
                      dimnames=list(CC.names, names(CC.1)[-(1:3)], final.loci))
for(i in 1:length(CC.paths)){
  this.csv <- read.csv(CC.paths[i], header=TRUE)
  result.array[i,,] <- t(this.csv[keep, -(1:3)])
}

# Load Megamuga SNP chip info
load("/path/to/megamuga.Rdata")

convert.full.DOQTL.array.to.HAPPY(full.array=result.array, map=snps,
                                  map.locus_name.colname="marker", map.chr.colname="chr", 
                                  map.physical_dist.colname="pos", map.genetic_dist.colname="cM",
                                  HAPPY.output.path="/full_genome_cache/", remove.chr.from.chr=TRUE,
                                  physical_dist.is.Mb=FALSE,
                                  allele.labels=LETTERS[1:8],
                                  chr=c(1:19, "X"),
                                  diplotype.order=c("CC"))

# Use L2 norm to collapse genome cache
collapse.genomecache(original.cache="/full_genome_cache/", 
                     new.cache="/reduced_genome_cache/", subjects=NULL,
                     criterion="l2.norm",
                     model="additive",
                     proportion.tol=0.1)

cache_path <- file.path("/reduced_genome_cache/")