## illustration with figures
## compositional analysis
compositional_analysis <- function(PhyloseqObj, target_rank='order', rel_abundance=0.03){
  ## PhyloseqObj: phyloseq object after filtering, decontamination and normalization
  ## target_rank: one of : order, family, genus, species
  ## rel_abundance: a threshold which less than that to be considered as 'others'.
  physeq.fam <- transform_sample_counts(PhyloseqObj, function(x){x/sum(x)})
  df.fam.melt <- psmelt(physeq.fam)
  if(tolower(target_rank)=='order'){
    df.fam.melt$target_rank <- ifelse(df.fam.melt$Abundance>rel_abundance,as.character(df.fam.melt$order), "others")
  } else if(tolower(target_rank)=='family'){
    df.fam.melt$target_rank <- ifelse(df.fam.melt$Abundance>rel_abundance,as.character(df.fam.melt$family), "others")
  } else if(tolower(target_rank)=='genus'){
    df.fam.melt$target_rank <- ifelse(df.fam.melt$Abundance>rel_abundance,as.character(df.fam.melt$genus), "others")
  } else if(tolower(target_rank)== 'species'){
    df.fam.melt$target_rank <- ifelse(df.fam.melt$Abundance>rel_abundance,as.character(df.fam.melt$species), "others")
  } else{
    stop("Unknown target rank. Must be one of the list: order, family, genus, species.")
  }
  
  df.fam.melt$target_rank <- factor(df.fam.melt$target_rank, levels=rev(unique(df.fam.melt$target_rank)))
  cols <- c("#9d547c","#56ca63","#a357d6","#419d2a","#525fd6","#8cbe3a","#c944aa","#5ba557","#9e66cb","#c1b735","#6d82ec","#e69728","#6654b0","#799330","#da7fdf","#3c782c","#e44586","#63c996","#dc3f53","#49cbc8","#cf3f29","#4fabda","#da6c2b","#598bd1","#b78c24","#8d4191","#a0b971","#b2386a","#479d71","#ae4341","#2ba198","#e07557","#5361a3","#dda353","#aa98df","#5b6114","#dc89bf","#327243","#e57b94","#277257","#9b62a0","#bbab59","#98495a","#526229","#d8827d","#857624","#9a4a22","#7c7d46","#e3a073","#9e6b33")
  
  
  ggplot(df.fam.melt, aes(x=id, y=Abundance, fill=target_rank)) +
    geom_bar(stat="identity") +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1.02)) +
    scale_fill_manual(values = cols) +
    xlab ("") +
    ylab("Relative abundance") +
    theme_bw() +
    ggtitle(paste0("Microbial composition at ",target_rank, " level"))+
    theme( axis.text.x = element_text(angle=45,  vjust = .5),
           axis.text.y = element_text (size=12),
           axis.title = element_text(size=14, face="bold"))+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
          legend.key.size = unit(4,"mm"))
}
##------------------------------------------------------
## Alpha analysis
alpha_analysis <- function(PhyloseqObj, choosen_metric='shannon'){
  rich = estimate_richness(PhyloseqObj)
  rich$id <- gsub("X","", rownames(rich))
  rownames(rich) <- rich$id
  rich_meta <- merge(sample_data(PhyloseqObj),rich, by="row.names")
  rich_meta %>% group_by(sample_side) %>% count()
  
  sh.cm <- NULL # #Shapiro-Wilk
  l.cm <- NULL #Levene Test
  cm <- NULL # stat-test
  title_name <- NULL
  group.colors <- c(normal = "dodgerblue3", tumor = "firebrick2")
  my.labels <- c("Normal pancreas", "Tumor tissue")
  if(tolower(choosen_metric)=="observed"){
    title_name <- "observed species"
    sh.cm <- rich_meta %>% shapiro_test(Observed)
    l.cm <- rich_meta %>% levene_test(Observed ~ sample_side)
    cm <- rich_meta %>% wilcox_test(Observed ~ sample_side) %>% adjust_pvalue(method = "BH") %>%  add_significance("p.adj") %>% add_xy_position()
  } else if(tolower(choosen_metric)== "shannon"){
    title_name <- "Shannon index"
    sh.cm <- rich_meta %>% shapiro_test(Shannon)
    l.cm <- rich_meta %>% levene_test(Shannon ~ sample_side)
    cm <- rich_meta %>% t_test(Shannon ~ sample_side) %>% adjust_pvalue(method = "BH") %>%  add_significance("p.adj") %>% add_xy_position()
  } else if(tolower(choosen_metric)=="invsimpson"){
    title_name <- "Inverse Simpson index"
    sh.cm <- rich_meta %>% shapiro_test(InvSimpson)
    l.cm <- rich_meta %>% levene_test(InvSimpson ~ sample_side)
    cm <- rich_meta %>% wilcox_test(InvSimpson ~ sample_side) %>% adjust_pvalue(method = "BH") %>%  add_significance("p.adj") %>% add_xy_position()
  } else{
    stop("choosen_metric must be in: Observed, Shannon, InvSimpson")
  }
  
  #plot
  ggplot (rich_meta, aes (x=sample_side, y=choosen_metric, fill=sample_side))+ 
    geom_boxplot()+
    geom_point (position=position_jitterdodge( jitter.width = 0.05))+
    theme_gray() + 
    xlab("")+
    # facet_grid(.~tp)+
    scale_x_discrete(labels= my.labels)+
    scale_fill_manual(values=group.colors, labels = c("\nNormal pancreas\nn = 15\n", "\nTumor tissue\nn = 25\n"))+
    ggtitle(paste0("FFPE human tumor tissue vs. normal pancreas - ", title_name)) +
    theme(axis.text.y = element_text (size=12),
          axis.title = element_text(size=14, face="bold"))+
    theme(legend.text = element_text(size = 12),
          legend.title = element_blank(),
          # legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
          legend.key.size = unit(4,"mm"),
          axis.text.x = element_text(size=12, hjust = .5, angle=0, vjust = .5),
          plot.title = element_text(size = 14))+
    stat_pvalue_manual(cm, label = "p.adj.signif", inherit.aes = FALSE, tip.length = 0.01)
}
##-----------------------------------------------------------------------------
## Beta analysis
beta_analysis <- function(PhyloseqObj, distance_metric="unweighted UF"){
  panc_count <- as_tibble(sample_data(PhyloseqObj)) %>% filter(sample_side == "normal") %>% count(sample_side) %>% pull(n)
  tum_count <- as_tibble(sample_data(PhyloseqObj)) %>% filter(sample_side == "tumor") %>% count(sample_side) %>% pull(n)
  tmp_dist <- NULL
  tmp_ordination <- NULL
  dist_name <- NULL
  if(tolower(distance_metric)=="unweighted_uf"){
    dist_name <- "unweighted UniFrac"
    tmp_dist <- phyloseq::distance(PhyloseqObj, method="unifrac", weighted=F)
  }else if(tolower(distance_metric)=="weighted_uf"){
    #dist_name <- "weighted UniFrac"
    tmp_dist <- phyloseq::distance(PhyloseqObj, method="unifrac", weighted=T)
  }else if(tolower(distance_metric)=="bray"){
    dist_name <- 'Bray curtis'
    tmp_dist <- phyloseq::distance(PhyloseqObj, method="bray")
  }else{
    stop("Unknown distance metric. It must be in: unweighted_uf, weighted_uf, bray")
  }
  tmp_ordination <- phyloseq::ordinate(PhyloseqObj, method="PCoA", distance=tmp_dist)
  
  pcoa1 <- paste("PCoA 1", " [", round(tmp_ordination[[3]]$Relative_eig[1], digits = 3)*100, "%]", sep = "")
  pcoa2 <- paste("PCoA 2", " [", round(tmp_ordination[[3]]$Relative_eig[2], digits = 3)*100, "%]", sep = "")
  
  set.seed(40)
  p.adonis <- adonis2(tmp_dist ~ sample_data(PhyloseqObj)$sample_side)
  p <- case_when(
    p.adonis$`Pr(>F)`[1] > 0.05 ~ "n.s.",
    p.adonis$`Pr(>F)`[1] < 0.05 &  p.adonis$`Pr(>F)`[1] > 0.01 ~ paste(p.adonis$`Pr(>F)`[1], "*", sep = " "),
    p.adonis$`Pr(>F)`[1] <= 0.01 & p.adonis$`Pr(>F)`[1] > 0.001  ~ paste(p.adonis$`Pr(>F)`[1], "**", sep = " "),
    p.adonis$`Pr(>F)`[1] <= 0.001 ~ paste(p.adonis$`Pr(>F)`[1], "***", sep = " "),
  )
  annotations <- data.frame(
    xpos = c(-Inf),
    ypos =  c(Inf),
    annotateText = p,
    hjustvar = c(-0.2) ,
    vjustvar = c(1.5))
  
  plot_ordination(PhyloseqObj, tmp_ordination, color="sample_side") + 
    geom_point(aes(colour=sample_side), size=3) +
    theme(aspect.ratio=1) +
    theme_bw()+
    scale_color_manual(values=group.colors, labels=c(paste("Normal pancreas\nn =", panc_count, sep = " "), paste("Tumor tissue\nn =", tum_count, sep = " ")))+
    stat_ellipse() +
    xlab(pcoa1)+
    ylab(pcoa2)+
    theme(panel.grid =  element_blank())+
    ggtitle(paste0("FFPE human tumor tissue vs. normal pancreas - ", dist_name)) +
    theme (axis.text=element_text(size=14),
           axis.title=element_text(size=16,face="bold"))+
    theme(
      # legend.justification=c(1.15,-0.15), legend.position=c(1,0),
      legend.text = element_text(size = 12),
      legend.title = element_blank())+
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), inherit.aes = FALSE)
  
}
