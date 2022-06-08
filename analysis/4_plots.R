source("0_functions.R")

dir.create("../output/plots")

#########################
#### Phenotype plots ####
#########################

protein_data = read.table("../output/protein_data.tsv", header = TRUE)
proteinPairs = read.table("../output/protein_pairs.tsv", header = TRUE)

# Summarize protein data
protein_data <- protein_data %>%
  group_by(strain, rx) %>%
  summarize(meanConc = mean(protein.concentration),
            sd = sd(protein.concentration),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by FA values
strain_order <- protein_data %>%
  filter(rx == "FA", .preserve = TRUE) %>%
  arrange(meanConc) %>%
  dplyr::select(strain)
# Reorder
protein_data <- protein_data %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

# Protein plot
protein_plot <- balPlot(protein_data, meanConc,
                        expression(atop("BAL Total Protein", "Concentration (ug/mL)")))

# Protein histogram
protein_histo <- balHisto(proteinPairs, protein.delta,
                          xlab = expression("Total Protein Concentration Difference (" * O[3] * " - FA)"))

# Fig with arranged plots
protein_arranged <- protein_plot + protein_histo +
  plot_layout(widths = c(6/9, 3/9))

ggsave(protein_arranged, filename= "../plots/protein_arranged.png", dpi = 300, units = "in",
       height = 6.5, width = 18.5)

#---------------

bal_data = read.table("../output/bal_data.tsv", header = TRUE)

# Summarize BAL data - neutrophil percentage
neutrophil_per <- bal_data %>%
  filter(!is.na(p.pmn)) %>%
  filter(mouse_id != c("5168","5386")) %>%
  group_by(strain, rx) %>%
  summarize(meanNeutPer = mean(p.pmn),
            sd = sd(p.pmn),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by O3 values
strain_order <- neutrophil_per %>%
  filter(rx == "O3", .preserve = TRUE) %>%
  arrange(meanNeutPer) %>%
  dplyr::select(strain)
# Reorder
neutrophil_per <- neutrophil_per %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

# Neutrophil percentage plot
neup_plot <- balPlot(neutrophil_per, meanNeutPer, ylab = "% neutrophils in BAL") + 
  scale_y_continuous(limits=c(0,90),breaks=c(0,25,50,75))

# Summarize BAL data - total neutrophils
neutrophil_num <- bal_data %>%
  filter(!is.na(p.pmn)) %>%
  filter(mouse_id != c("5168","5386")) %>%
  group_by(strain, rx) %>%
  mutate(total.pmn = total * p.pmn/10) %>%
  summarize(meanNeutNum = mean(total.pmn),
            sd = sd(total.pmn),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by O3 values
strain_order <- neutrophil_num %>%
  filter(rx == "O3", .preserve = TRUE) %>%
  arrange(meanNeutNum) %>%
  dplyr::select(strain)
# Reorder
neutrophil_num <- neutrophil_num %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

neut_plot <- balPlot(neutrophil_num, meanNeutNum,
                     ylab = parse(text=paste0('"# of neutrophils"','~ (10^6)'))) +
  scale_y_continuous(limits = c(0, 8e6), breaks = seq(0,8e6,by=1.5e6), 
                     labels = c(0, 1.5, 3, 4.5, 6, 7.5))
  
# Arranged BAL data plot
neutrophils_arranged <- neup_plot / neut_plot

ggsave(neutrophils_arranged, filename="../output/plots/neutrophils_arranged.png", dpi = 300,
       units = "in", width = 13, height = 8)

#---------------

# Macrophage plots

# Summarize BAL data - macrophage percentage
macrophage_per <- bal_data %>%
  filter(!is.na(p.mac)) %>%
  filter(mouse_id != c("5168","5386")) %>%
  group_by(strain, rx) %>%
  summarize(meanMacPer = mean(p.mac),
            sd = sd(p.mac),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by FA values
strain_order <- macrophage_per %>%
  filter(rx == "O3", .preserve = TRUE) %>%
  arrange(meanMacPer) %>%
  dplyr::select(strain)
# Reorder
macrophage_per <- macrophage_per %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

# Macrophage percentage plot
macp_plot <- balPlot(macrophage_per, meanMacPer, ylab = "% macrophages in BAL") + 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
  theme(legend.position="none")

# Summarize BAL data - total macrophages
macrophage_num <- bal_data %>%
  filter(!is.na(p.mac)) %>%
  filter(mouse_id != c("5168","5386")) %>%
  group_by(strain, rx) %>%
  mutate(total.mac = total * p.mac/10) %>%
  summarize(meanMacNum = mean(total.mac),
            sd = sd(total.mac),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by FA values
strain_order <- macrophage_num %>%
  filter(rx == "O3", .preserve = TRUE) %>%
  arrange(meanMacNum) %>%
  dplyr::select(strain)
# Reorder
macrophage_num <- macrophage_num %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

mact_plot <- balPlot(macrophage_num, meanMacNum,
                     ylab = parse(text=paste0('"# of macrophages"','~ (10^7)'))) +
  scale_y_continuous(limits = c(0, 3.5e7), breaks = seq(0,3.5e7,by=5e6), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
  theme(legend.position = "none")

#---------------

# Total cell number

# Summarize BAL data - total cells
total_cells <- bal_data %>%
  filter(!is.na(total)) %>%
  filter(mouse_id != c("5168","5386")) %>%
  group_by(strain, rx) %>%
  mutate(total.cells = total * 10) %>%
  summarize(meanNum = mean(total.cells),
            sd = sd(total.cells),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by FA values
strain_order <- total_cells %>%
  filter(rx == "O3", .preserve = TRUE) %>%
  arrange(meanNum) %>%
  dplyr::select(strain)
# Reorder
total_cells <- total_cells %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

total_plot <- balPlot(total_cells, meanNum,
                     ylab = parse(text=paste0('"Total # of BAL cells"','~ (10^7)'))) +
  scale_y_continuous(limits = c(0, 3.5e7), breaks = seq(0,3.5e7,by=5e6), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5))

bal_arranged <- total_plot / macp_plot / mact_plot

ggsave(bal_arranged, filename="../output/plots/bal_arranged.png", dpi = 300,
       units = "in", width = 13, height = 12)


#---------------

# Main data plotted by sex

# Summarize BAL data w sex - neutrophil percentage
neutrophil_per_sex <- bal_data %>%
  filter(!is.na(p.pmn)) %>%
  filter(mouse_id != c("5168","5386")) %>%
  group_by(strain, rx, sex) %>%
  summarize(meanNeutPer = mean(p.pmn),
            sd = sd(p.pmn),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by FA values (ignoring sex)
strain_order <- neutrophil_per %>%
  filter(rx == "O3", .preserve = TRUE) %>%
  arrange(meanNeutPer) %>%
  dplyr::select(strain)
# Reorder
neutrophil_per_sex <- neutrophil_per_sex %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

neup_sex <- balSexDiffPlot(df = neutrophil_per_sex, var = meanNeutPer,
                          ylab = "% neutrophils in BAL")

# Summarize BAL data w sex - total neutrophils
neutrophil_num_sex <- bal_data %>%
  filter(!is.na(p.pmn)) %>%
  filter(mouse_id != c("5168","5386")) %>%
  group_by(strain, rx, sex) %>%
  mutate(total.pmn = total * p.pmn/10) %>%
  summarize(meanNeutNum = mean(total.pmn),
            sd = sd(total.pmn),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by O3 values
strain_order <- neutrophil_num %>%
  filter(rx == "O3", .preserve = TRUE) %>%
  arrange(meanNeutNum) %>%
  dplyr::select(strain)
# Reorder
neutrophil_num_sex <- neutrophil_num_sex %>%
  mutate(strain = factor(strain, levels = strain_order$strain))


neut_sex <- balSexDiffPlot(df = neutrophil_num_sex, var = meanNeutNum,
                           ylab = parse(text=paste0('"# of neutrophils"','~ (10^6)'))) +
  scale_y_continuous(limits = c(0, 8e6), breaks = seq(0,8e6,by=1.5e6), 
                     labels = c(0, 1.5, 3, 4.5, 6, 7.5)) +
  theme(legend.position="none")

# Summarize protein data w sex
protein_data_sex <- protein_data %>%
  group_by(strain, rx, sex) %>%
  summarize(meanConc = mean(protein.concentration),
            sd = sd(protein.concentration),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by FA values
strain_order <- protein_data %>%
  filter(rx == "FA", .preserve = TRUE) %>%
  arrange(meanConc) %>%
  dplyr::select(strain)
# Reorder
protein_data_sex <- protein_data_sex %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

# Protein plot
protein_sex <- balSexDiffPlot(df = protein_data_sex, var = meanConc,
                                   ylab = expression(atop("BAL Total Protein", "Concentration (ug/mL)"))) +
  theme(legend.position="none")

# Arranged protein, neu data by sex plot
neu_prot_sex_arranged <- neup_sex / neut_sex / protein_sex

ggsave(neu_prot_sex_arranged, filename="../output/plots/neu_prot_sex_arranged.png", dpi = 300,
       units = "in", width = 13, height = 12)

#---------------

# Other BAL data by sex

# Total cells by sex
total_cells_sex <- bal_data %>%
  filter(!is.na(total)) %>%
  filter(mouse_id != c("5168","5386")) %>%
  group_by(strain, rx, sex) %>%
  mutate(total.cells = total * 10) %>%
  summarize(meanNum = mean(total.cells),
            sd = sd(total.cells),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by FA values
strain_order <- total_cells %>%
  filter(rx == "O3", .preserve = TRUE) %>%
  arrange(meanNum) %>%
  dplyr::select(strain)
# Reorder
total_cells_sex <- total_cells_sex %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

# Plot
total_sex <- balSexDiffPlot(df = total_cells_sex, var = meanNum,
                                   ylab = parse(text=paste0('"Total # of BAL cells"','~ (10^7)'))) +
  scale_y_continuous(limits = c(0, 3.5e7), breaks = seq(0,3.5e7,by=5e6), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5))

# Summarize BAL data - macrophage percentage by sex
macrophage_per_sex <- bal_data %>%
  filter(!is.na(p.mac)) %>%
  filter(mouse_id != c("5168","5386")) %>%
  group_by(strain, rx, sex) %>%
  summarize(meanMacPer = mean(p.mac),
            sd = sd(p.mac),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by FA values
strain_order <- macrophage_per %>%
  filter(rx == "O3", .preserve = TRUE) %>%
  arrange(meanMacPer) %>%
  dplyr::select(strain)
# Reorder
macrophage_per_sex <- macrophage_per_sex %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

# Macrophage percentage plot
macp_sex <- balSexDiffPlot(df = macrophage_per_sex, var = meanMacPer,
                           ylab = "% macrophages in BAL") + 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
theme(legend.position="none")

# Summarize BAL data - total macrophages by sex
macrophage_num_sex <- bal_data %>%
  filter(!is.na(p.mac)) %>%
  filter(mouse_id != c("5168","5386")) %>%
  group_by(strain, rx, sex) %>%
  mutate(total.mac = total * p.mac/10) %>%
  summarize(meanMacNum = mean(total.mac),
            sd = sd(total.mac),
            n = n(),
            se = sd/sqrt(n))
# Extract strain order by FA values
strain_order <- macrophage_num %>%
  filter(rx == "O3", .preserve = TRUE) %>%
  arrange(meanMacNum) %>%
  dplyr::select(strain)
# Reorder
macrophage_num_sex <- macrophage_num_sex %>%
  mutate(strain = factor(strain, levels = strain_order$strain))

mact_sex <- balSexDiffPlot(df = macrophage_num_sex, var = meanMacNum,
                     ylab = parse(text=paste0('"# of macrophages"','~ (10^7)'))) +
  scale_y_continuous(limits = c(0, 3.5e7), breaks = seq(0,3.5e7,by=5e6), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
  theme(legend.position="none")

bal_sex_arranged <- total_sex / macp_sex / mact_sex

ggsave(bal_sex_arranged, filename="../output/plots/bal_sex_arranged.png", dpi = 300,
       units = "in", width = 13, height = 12)

###########################
#### QTL Mapping Plots ####
###########################

load(file = "../output/mapping/difference_rop.rds")
load(file = "../output/mapping/difference_rop_perm_scans.rds")
load(file = "../output/mapping/difference_positional_ten.rds")
load(file = "../output/mapping/difference_positional_fift.rds")

# Define permutation-based significance threshold for difference
difference_90 <- get.gev.thresholds(difference_rop_perm_scans,use.lod=FALSE,percentile=0.9,type="min")
difference_95 <- get.gev.thresholds(difference_rop_perm_scans,use.lod=FALSE,percentile=0.95,type="min")
difference_99 <- get.gev.thresholds(difference_rop_perm_scans,use.lod=FALSE,percentile=0.99,type="min")

# Grab markers for top loci
difference_ten <- grab.locus.from.scan(difference_rop, chr = 10)
difference_fift <- grab.locus.from.scan(difference_rop, chr = 15)

# Difference scan
png('../plots/genome_scan.png',units="in",res=300,height=3.5,width=12)
par(cex=1.5)
par(mar=c(1.75,3,1.5,0.75))
genome.plotter.whole(list(difference_rop),
                     hard.thresholds=c(difference_90,difference_95),
                     thresholds.col=c("#0049c3","#f45d48"),
                     main.colors="#4b4b4b",
                     distinguish.chr.type = "box",
                     distinguish.box.col="#d7d7d7",
                     my.y.lab.cex = 1,
                     y.max.manual=5,
                     thresholds.lwd=c(2,2),
                     my.legend.lwd=1.25,
                     override.title=expression(bold("Total Protein Concentration Difference (" * O[3] * " - FA)")))
dev.off()

# Haplotype dosage plot (Chr 10)
merge.by <- "strain"
phenotype <- "protein.delta"
png('../output/plots/chr10_haplotypedosage.png',units="in",height=6,width=6.45,res=300)
par(oma=c(0.25,3.5,0.25,0.5))
prob.heatmap(difference_ten,
             genomecache = cache_path,
             phenotype=phenotype, phenotype.data=proteinPairs,
             merge.by=merge.by,
             founder.labels = c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ",
                                "NZO/HILtJ","CAST/EiJ","PWK/PhJ","WSB/EiJ"),
             founder.cex=1.25, phenotype.lab.cex=1.25,phenotype.num.cex=1.25,prob.axis.cex=1.1)
dev.off()

# 80% credible interval plot (Chr 10)
png('../output/plots/chr10_credibleinterval.png',units="in",height=6,width=6.45,res=300)
par(cex=1.5)
single.chr.plotter.w.ci(difference_rop,difference_positional_ten,
                        ci.type="Parametric Bootstrap",
                        scan.type="ROP",alpha = 0.2,
                        these.col=c("#0049ce","#f45d48"))
dev.off()

# Straw plot (Chr 10)
png('../output/plots/chr10_strawplot.png',res=300,width=6.45,height = 6,units="in")
par(cex=1.5)
plot_locus.effect.from.scan(difference_rop,
                            difference_ten,
                            names=c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ",
                                    "NZO/HILtJ","CAST/EiJ","PWK/PhJ","WSB/EiJ"),
                            main = "Allele Effects for Chr. 10 Peak Marker (Difference)")
dev.off()

# Haplotype dosage plot (Chr 15)
merge.by <- "strain"
phenotype <- "protein.delta"
png('../output/plots/chr15_haplotypedosage.png',units="in",height=6.1,width=6,res=300)
par(mar=c(0.1,0.8,0.1,3.2))
prob.heatmap(difference_fift,
             genomecache = cache_path,
             phenotype=phenotype,phenotype.data=proteinPairs,
             merge.by=merge.by,alternative.phenotype.label = "Protein Difference",
             founder.labels = NULL,
             founder.cex=1.25, phenotype.lab.cex=1.25,phenotype.num.cex=1.25,prob.axis.cex=1.1)
dev.off()

# 80% credible interval plot (Chr 15)
png('../output/plots/chr15_credibleinterval_wide.png',units="in",height=5,width=5,res=300)
par(cex=1.25)
single.chr.plotter.w.ci(difference_rop, difference_positional_fift,
                        ci.type="Parametric Bootstrap",
                        scan.type="ROP",alpha = 0.2,
                        these.col=c("#0049ce","#f45d48"))
dev.off()

# Straw plot (Chr 15)
png('../output/plots/chr15_strawplot.png',res=300,width=8,height = 6,units="in")
par(cex=1.5)
plot_locus.effect.from.scan(difference_rop,
                            difference_fift,
                            names=c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ",
                                    "NZO/HILtJ","CAST/EiJ","PWK/PhJ","WSB/EiJ"),
                            main = "Allele Effects for Chr. 15 Peak Marker (Difference)")
dev.off()


# Plot posterior haplotype effects from TIMBR
load(file='../output/mapping/chr15_timbr.rds')

png('../output/plots/TIMBR_chr15.png',width=6,height=4,units="in",res=300)
TIMBR.plot.haplotypes(results)
title(main="TIMBR Results for Chr. 15 QTL")
dev.off()

# Diploffect plot
load(file='../output/mapping/chr15_diploffect.rds')

png('../output/plots/Diploffect_chr15.png', width = 6.25, height = 6, units = "in", res = 300)
par(mar=c(0.1,4,0.1,0.5))
par(cex=1.5)
plot.straineff.ci(inla.diploffect.summary,main =NULL,
                  flip=FALSE,sn = c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ",
                                    "NZO/HILtJ","CAST/EiJ","PWK/PhJ","WSB/EiJ"))
dev.off()

################################
#### Subspecies origin plot ####
################################

# Requires input from Mouse Phylogeny Viewer (msub.csbio.unc.edu)
## Coordinates reported as Build 37, use liftOver to make Build 38

subspecies <- read.csv('../data/subspecies_origin_coordinates.csv',header=TRUE)

subspecies <- subspecies %>%
  mutate(length = end_positionb38 - start_positionb38,
         subspecies = as.factor(subspecies),
         strain = as.factor(strain))

subspecies <- subspecies %>%
  mutate(subspecies = factor(subspecies, levels = c("Dom", "Cast", "Mus", "Undef", "Marker")),
         strain = factor(strain, levels = c("WSB/EiJ", "PWK/PhJ", "CAST/EiJ", "NZO/HlLtJ",
                                            "NOD/ShiLtJ", "129S1SvlmJ", "C57BL/6J", "A/J")))

levels(subspecies$subspecies) <- c("Domesticus","Castaneus","Musculus","Undef","Marker")

subspecies_plot <- ggplot() + 
  geom_bar(aes(y=length,x=strain,fill=subspecies,group=strain),data=subspecies,stat="identity",width=0.95) + 
  coord_flip() + 
  scale_fill_manual(breaks=c("Domesticus","Castaneus","Musculus"),
                    values=c("#1111ff","#00aa00","#ff0000","#FFFFFF","#000000"),
                    name="Subspecies Origin") + 
  labs(y="Chr. 15 (Mb)",x=NULL) + 
  theme_minimal(base_size=16) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.5,color="black"),
        legend.position="top",axis.ticks.x = element_line(size=0.5),
        axis.text.y=element_text(color="black")) +
  scale_y_continuous(breaks=seq(833546,13833546,by=1e6),
                     labels=c("41","42","43","44","45","46","47","48","49","50","51","52","53","54"),
                     expand=c(0,0),limits=c(-0.01e6,14.72e6)) + 
  scale_x_discrete(expand=c(0,0)) +
  guides(fill=guide_legend(label.theme=element_text(face="italic",size=15))) + 
  geom_hline(yintercept=6952000,linetype="dashed") +
  geom_text(aes(4.5,6900000,label="peak marker",angle=90,vjust=-0.5),size=5.5)

ggsave(subspecies_plot, filename = "../output/plots/subspecies_origin_chr15qtl.png", dpi = 300,
       width = 8, height = 4)

############################
#### Edited Merge Plots ####
############################

# Altered original merge code (found at 'yanweicai/MergeAnalysisCC')
ensembl <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
filters = listFilters(ensembl)

ph='y';peakI <- c(40.2,54.9);CHR <- '15';Pcutoff = 4

diplo <- read.table(file="../data/merge_top-level.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE, colClasses=c("sdp"="character"))

diplo <- diplo[which(diplo$gene_name!='unknown_intergenic'),c('pos','sdp','logP.merge','logP.founder','gene_name','cq','transcript_name')]

# Add annotations to gene names
for (i in 1:dim(diplo)[1]){ diplo$gene_name[i] <-  strsplit(diplo$gene_name[i],split=',')[[1]][1] }

geneanno <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol","description"),values=unique(diplo$gene_name),mart= ensembl)

diplo <- merge(diplo,geneanno,by.x='gene_name',by.y='ensembl_gene_id')
diplo <- diplo[order(-diplo$logP.merge),c(2,3,4,5,1,8,9,6,7)]
diplo$posm<-diplo$pos/10^6

# Select genotype order from table
geno <- diplo[which(diplo$logP.merge>=Pcutoff),]
geno <- geno[which(geno$posm>=peakI[1] & geno$posm<=peakI[2]),]

diplo<- read.table(file="../data/merge_diplo-res-sample.txt",header=TRUE,stringsAsFactors=FALSE)
diplo$posm<-diplo$pos/10^6
diplo <- diplo[order(diplo$posm),c('posm','efs')]
colnames(diplo)[2] <- 'logP.founder'

# Order SDP by number of variants
StrDis_df <- aggregate(logP.merge ~ sdp, data=geno,length )
StrDis_df <- StrDis_df[order(StrDis_df$logP.merge),]
StrDis <- StrDis_df$sdp

rincol <- brewer.pal(8, 'Dark2')[1:8]; # select eight-color base palette
rincol <- c(rep('darkblue',(length(StrDis)-8)),rincol)
names(rincol) <- StrDis;
geno$col <- rincol[geno$sdp]

#---------------

# Create top plot
png('../output/plots/merge_panel_a_b.png', res = 300, units = "in", height = 12, width = 16)
par(mfrow=c(2,1),mar=c(2.5,5,0.5,1.1))
plot(diplo$posm,diplo$logP.founder,xlim=peakI,ylim=c(0,max(geno$logP.merge,diplo$logP.founder)),
     col='red',ylab="-log10(P-value)",xlab=NA,type='l',lwd=2.5)
points(diplo$posm,diplo$logP.merge,pch=20,col='grey',cex=0.8)
points(geno$posm,geno$logP.merge,pch=20,col=geno$col)
#abline(v=peakP,lty=3,lwd=2)

blcknm <- length(StrDis)
plot(geno$posm,geno$posm,xlim=peakI,ylim=c(0,(blcknm-0.5)),type='n',yaxt = "n",xlab=paste0('Mb on Chromosome',CHR),ylab='')
axis(2, at = seq(0.5, (blcknm-0.5) , by = 1), label = StrDis, las=1,cex.axis=1)
for (ii in c(1:floor(blcknm/2))){
  rect(peakI[1], (2*ii-1), peakI[2], (2*ii), col='gray', angle = 0,lwd=0)
}

for ( i in c(1:length(geno$posm))){
  if (geno$sdp[i] %in% StrDis){
    lineid=which(StrDis==geno$sdp[i])
    segments(geno$posm[i], (lineid-1) ,geno$posm[i],(lineid),col=alpha(geno$col[i]),lwd=2)
  }
}

dev.off()

#---------------

# Create bottom plot (variant consequence names added in Photoshop)

# Create data frame where each variant is separated out

data_4b <- geno %>% 
  separate(cq, c("Consequence_A","Consequence_B","Consequence_C"), ",")

data_csq <- melt(data_4b, id.vars = c(1:3), measure.vars = c(8:10))
data_csq <- data_csq[,c(1:3,5)] # 10 levels
data_csq$posm<-data_csq$pos/10^6
data_csq <- data_csq[!is.na(data_csq$value),]

StrDis_df <- aggregate(logP.merge ~ value, data=data_csq,length )
StrDis_df <- StrDis_df[order(StrDis_df$logP.merge),]
StrDis <- StrDis_df$value

data_csq <- merge(data_csq, geno[,c(1,11)], by = "pos")
data_csq <- data_csq %>%
  add_count(col) %>%
  arrange(desc(n)) %>%
  dplyr::select(-n)

blcknm <- length(StrDis)

png('../output/plots/merge_panel_c.png', res = 300, units = "in", height = 3.5, width = 16)
par(mar=c(2.5,5,0.5,1.1))
plot(data_csq$posm,data_csq$posm,xlim=peakI,ylim=c(0,(blcknm-0.5)),type='n',yaxt = "n",xlab=paste0('Mb on Chromosome',CHR),ylab='')
axis(2, at = seq(0.5, (blcknm-0.5) , by = 1), label = StrDis, las=1,cex.axis=1)
for (ii in c(1:floor(blcknm/2))){
  rect(peakI[1], (2*ii-1), peakI[2], (2*ii), col='gray', angle = 0,lwd=0)
}

for ( i in c(1:length(data_csq$posm))){
  if (data_csq$value[i] %in% StrDis){
    lineid=which(StrDis==data_csq$value[i])
    segments(data_csq$posm[i], (lineid-1) ,data_csq$posm[i],(lineid),col=alpha(data_csq$col[i]),lwd=2)
  }
}

dev.off()

####################
#### Gene Track ####
####################

# This plot was created, printed to a PDF, and edited subsequently in Illustrator to remove pseudogenes

axTrack <- GenomeAxisTrack()

biomTrack <- BiomartGeneRegionTrack(genome="mm10",chromosome=15,
                                    start=40000000,end=4500000,name=NULL,
                                    collapseTranscripts="meta",
                                    biomart=useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl"))

pdf('../output/plots/merge_panel_d.pdf', width=9, height=3, pointsize=300)
plotTracks(list(axTrack,biomTrack),from=40000000,to=45000000,
           transcriptAnnotation="symbol",
           background.panel="transparent",showTitle=FALSE)
dev.off()
