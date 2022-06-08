#####################
#### Environment ####
#####################

#### USING R VERSION 4.1.1

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

haplo_labels <- c("FA.A"="FA",
                  "O3.A"=expression("O"[3]),
                  "FA.B"="FA",
                  "O3.B"=expression("O"[3]),
                  "FA.C"="FA",
                  "O3.C"=expression("O"[3]),
                  "FA.D"="FA",
                  "O3.D"=expression("O"[3]),
                  "FA.E"="FA",
                  "O3.E"=expression("O"[3]),
                  "FA.F"="FA",
                  "O3.F"=expression("O"[3]),
                  "FA.G"="FA",
                  "O3.G"=expression("O"[3]),
                  "FA.H"="FA",
                  "O3.H"=expression("O"[3]))

dichot_interact_labels <- c("FA.group_a"="FA",
                   "O3.group_a"=expression("O"[3]),
                   "FA.group_b"="FA",
                   "O3.group_b"=expression("O"[3]))

haplotypes <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "WSB", "PWK")

exposure_labels <- c("FA" ="FA",
                   "O3" = expression("O"[3]))

haplo_groups <- c("group_a"="High protein response",
                  "group_b"="Low protein response")


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
          legend.background=element_blank(),
          legend.box.background=element_rect(color="black"),
          legend.margin=margin(0.15,0.2,0.15,0.2,"cm"),
          legend.title=element_text(size=15),
          legend.text=element_text(size=15),
          legend.spacing.y=unit(0,"mm")) + 
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

# BAL sex differences
balSexDiffPlot <- function(df, var, ylab){
  ggplot(df,aes(x=strain,y={{var}},group=interaction(rx,strain))) + 
    geom_point(size=4,aes(shape=factor(rx),
                   color=factor(sex)),stroke=1) +
    geom_errorbar(aes(ymin={{var}}-se,ymax={{var}}+se,
                      color=factor(sex),linetype=factor(sex))) + 
    scale_shape_manual(values=c(1,15),
                       name="Treatment",
                       labels=c("FA",expression("O"[3]))) +
    scale_color_manual(values=c("red", "blue"),
                       name = "Sex") +
    scale_linetype_manual(values=c("solid","dashed")) +
    theme_linedraw(base_size=18) +
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1),
          plot.margin=unit(c(0.4,0.2,0.4,0.4),"cm"),
          legend.justification=c(0,1),
          legend.position=c(0.02,0.96),
          legend.background=element_blank(),
          legend.box.background=element_rect(color="black"),
          legend.margin=margin(0.15,0.2,0.15,0.2,"cm"),
          legend.title=element_text(size=15),
          legend.text=element_text(size=15),
          legend.spacing.y=unit(0,"mm")) + 
    labs(x="Strain",
         y=ylab) +
    scale_x_discrete(expand=c(0,0.1)) +
    guides(color=guide_legend(override.aes=list(shape=17,size=4,linetype=0)),
           linetype="none")
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

#################################
#### qPCR Analysis Functions ####
#################################

# Calculate mean Cq per sample
mean_cq <- function(df) {
  df %>%
    group_by(gene, sample_id) %>%
    mutate(Cq_sd = sd(Cq),
           Cq_mean = if_else(Cq_sd < 1, mean(Cq), max(Cq))) %>%
    distinct(Cq_mean, .keep_all = TRUE) %>%
    pivot_wider(id_cols = sample_id,
                names_from = gene,
                values_from = c(Cq_mean, Cq_sd)) %>%
    filter(sample_id %in% seq(1,38,by=1))
}

# dCq interaction plot
dcq_interact_plot <- function(df, gene, ylab, bot.val, top.val){
  y_axis_max = -as.vector(df %>% select({{gene}}) %>% arrange(desc({{gene}})) %>% filter(!is.na({{gene}})) %>% slice(1))[[1]]
  y_axis_min = -as.vector(df %>% select({{gene}}) %>% arrange({{gene}}) %>% filter(!is.na({{gene}})) %>% slice(1))[[1]]
  df %>% 
    ggplot(aes(x = interaction(rx, dichot),
               y = -{{gene}}, fill = dichot,
               group = interaction(rx, dichot))) +
    geom_boxplot(aes(group = interaction(rx, dichot)),outlier.shape = NA) +
    geom_jitter(aes(shape=rx, stroke=1),height=0.0,width=0.2,size=3) +
    scale_fill_manual(values = c("gray80", "gray30")) +
    scale_shape_manual(values = c(21, 22)) +
    labs(y = ylab) +
    theme_linedraw(base_size = 16) +
    theme(legend.position = "none",
          plot.margin=unit(c(1,1,4.5,1),"lines"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    coord_cartesian(ylim=c(y_axis_max+(y_axis_max*0.025),
                           y_axis_min-(y_axis_min*0.025)), clip = "off", expand = FALSE) +
    annotate(geom = "text", x = seq(1,4,by=1), y = y_axis_min+(y_axis_max+y_axis_min)/{{top.val}}, label = dichot_interact_labels, size = 4) +
    annotate(geom = "text", x = 1.5 + 2 * (0:1),y = y_axis_min+(y_axis_max+y_axis_min)/{{bot.val}}, label = haplo_groups, size= 5, fontface=2)
}

# dCq exposure effect plot
dcq_rx_plot <- function(df, gene, ylab){
  y_axis_max = as.vector(df %>% select({{gene}}) %>% arrange(desc({{gene}})) %>% filter(!is.na({{gene}})) %>% slice(1))[[1]]
  y_axis_min = as.vector(df %>% select({{gene}}) %>% arrange({{gene}}) %>% filter(!is.na({{gene}})) %>% slice(1))[[1]]
  df %>% 
    ggplot(aes(x = rx,
               y = -{{gene}}, shape = rx,
               fill = dichot, group = rx)) +
    geom_boxplot(fill = "white", outlier.shape = NA) +
    geom_jitter(aes(stroke=1),height=0.0,width=0.2,size=3) +
    scale_fill_manual(values = c("gray80", "gray30")) +
    scale_shape_manual(values = c(21, 22)) +
    labs(y = ylab, x = NULL) +
    theme_linedraw(base_size = 14) +
    theme(plot.title=element_text(hjust=0.5,face="bold"),
          legend.position = "none",
          axis.text.x =element_text(size = 14),
          plot.margin=unit(c(1,1,1,1),"lines")) +
    scale_x_discrete(labels = exposure_labels) +
    ggtitle("Exposure")
}

# dCq haplotype grouping effect plot
dcq_dichot_plot <- function(df, gene, ylab){
  y_axis_max = as.vector(df %>% select({{gene}}) %>% arrange(desc({{gene}})) %>% filter(!is.na({{gene}})) %>% slice(1))[[1]]
  y_axis_min = as.vector(df %>% select({{gene}}) %>% arrange({{gene}}) %>% filter(!is.na({{gene}})) %>% slice(1))[[1]]
  df %>% 
    ggplot(aes(x = dichot,
               y = -{{gene}}, shape = rx,
               fill = dichot, group = dichot)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(stroke=1),height=0.0,width=0.2,size=3) +
    scale_fill_manual(values = c("gray80", "gray30")) +
    scale_shape_manual(values = c(21, 22)) +
    labs(y = ylab, x = NULL) +
    theme_linedraw(base_size = 14) +
    theme(plot.title=element_text(hjust=0.5,face="bold"),
          legend.position = "none",
          axis.text.x =element_text(size = 10),
          plot.margin=unit(c(1,1,1,1),"lines")) +
    scale_x_discrete(labels = c("High protein resp.", "Low protein resp.")) +
    ggtitle("Haplotype grouping")
}
