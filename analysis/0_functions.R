#####################
#### Environment ####
#####################

#### USING R VERSION 4.0.3

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
devtools::install_github('gkeele/miqtl')
devtools::install_github('gkeele/Diploffect.INLA')
devtools::install_github('wesleycrouse/TIMBR')
BiocManager::install('Gviz')


tovar_packages <- c("openxlsx", "tidyverse", "miqtl", "TIMBR", "Diploffect.INLA", "reshape2",
                    "MASS","Gviz","biomaRt", "INLA", "patchwork", "scales", "RColorBrewer")
lapply(tovar_packages, require, character.only=TRUE)

cc_palette <- c("#ffff00","#888888","#ff8888","#1111ff",
                "#00ccff","#00aa00","#ff0000","#9900ee")

############################
#### Plotting Functions ####
############################

# Strain means plot

balPlot <- function(df, var, ylab){
  ggplot(df,aes(x=strain,y={{var}},fill=rx,group=interaction(rx,strain))) + 
    geom_point(aes(shape=rx,stroke=1),size=4) + 
    geom_errorbar(aes(ymin={{var}}-se,ymax={{var}}+se)) + 
    scale_shape_manual(values=c(21,22)) +
    scale_fill_manual(values=c("#FFFFFF","#000000"),
                      name="Treatment",
                      labels=c("FA",expression("O"[3]))) +
    theme_linedraw(base_size=18) +
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1),
          plot.margin=unit(c(0.4,0.2,0.4,0.4),"cm"),
          legend.justification=c(0,1),
          legend.position=c(0.02,0.96),
          legend.box.background=element_rect(fill="#FFFFFF",color="#000000",size=1),
          legend.margin=margin(0.15,0.2,0.15,0.2,"cm"),
          legend.title=element_text(size=15),
          legend.text=element_text(size=15)) + 
    labs(x="Strain",
         y=ylab) +
    scale_x_discrete(expand=c(0,0.1)) +
    guides(shape="none",
           fill=guide_legend(override.aes=list(shape=c(21,22),stroke=1,size=4)),
           color="none")
}

# BAL histogram

balHisto <- function(df, var, xlab) {
  ggplot(df,aes(x={{var}})) + 
  geom_histogram(bins=25,fill="#0049c3",color="black",alpha=0.75) + 
  theme_linedraw(base_size=16) +
  theme(plot.margin=unit(c(0.4,1.2,0.4,1.2),"cm")) +
  labs(x=xlab,y="Count")
}

################################
#### Heritability with INLA ####
################################

heritability <- function(df, traits) {
  inla.results <- list()
  for (i in 1:length(traits)) {
    y <- df[, traits[i]]
    pheno.K <- new.K[as.character(df$mouse_id), as.character(df$mouse_id)]
    
    remove.na <- is.na(y)
    y <- y[!(remove.na)]

    pheno.K <- pheno.K[!(remove.na), !(remove.na)]

    y <- as.vector(scale(y))
  
    inla.data <- list(
      prec.R = list(prior="loggamma", param=c(1,1)),
      prec.e = list(prior="loggamma", param=c(1,1)),
      inv.R = solve(pheno.K),
      y = y,
      x = rep(1, length(y)),
      idx = 1:length(y))  
  
    formula <- y ~ 1 + f(idx, model="generic2", constr=TRUE,
                         Cmatrix=inv.R, hyper=list(theta1=prec.R, 
                                                   theta2=prec.e))
    
    inla.fit <- inla(formula=formula, data=inla.data)
  
    inla.results[[i]] = inla.fit
  }
  names(inla.results) <- traits
}

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