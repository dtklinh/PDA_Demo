#Installing all packages.

bio_pkgs <- c('PERFect', 'microbiome', 'phyloseq', 'decontam', 'tidyverse', 'pairwiseAdonis', 'kableExtra', 'vegan',
              'tidytree', 'ape', 'ggrepel', 'ggpubr', 'rstatix')
BiocManager::install(bio_pkgs)
## install pairwiseAdonis
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(PERFect)
library(microbiome)
library(phyloseq)
library(decontam)
library(tidyverse)
library(pairwiseAdonis)
library(kableExtra)
library(vegan)
library(tidytree)
library(ape)
library(ggrepel)
library(ggpubr) # for manual p-value assigning
library(rstatix) # for wilcox test

###------------------------------------

physeq.comb.filt <- merge_phyloseq(p.true.RmLowAbun, physeq.nct.filt.s.prune)
sample_data(physeq.comb.filt)$is.neg <- sample_data(physeq.comb.filt)$TRUE_control == "control"

prev.nct <- tibble(
  prev = prevalence(physeq.nct.filt.s.prune, detection  = 0, sort = TRUE, count = FALSE),
  OTU =as.factor(names(prevalence(physeq.nct.filt.s.prune , detection = 0, sort = TRUE, count = FALSE)))
)

prev.trueSample <- tibble(
  prev = prevalence(p.true.RmLowAbun, detection  = 0, sort = TRUE, count = FALSE),
  OTU =as.factor(names(prevalence(p.true.RmLowAbun , detection = 0, sort = TRUE, count = FALSE)))
)

x <- prev.nct$prev
factorx <- factor(cut(x, breaks=nclass.Sturges(x)))
xout <- as.data.frame(table(factorx)) %>% map_df(rev)
xout <- mutate(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))
xout <- xout %>% mutate(xxx=cumFreq/685)

##---------------------------

sample_idx <- metadata[metadata$DNAex_round==1,"uniqueID"]
physeq.comb.DNAex.1 <- subset_samples(physeq.comb.filt, uniqueID %in% sample_idx)
physeq.comb.DNAex.1 <- prune_taxa(taxa_sums(physeq.comb.DNAex.1)>0, physeq.comb.DNAex.1)
metadata[metadata$DNAex_round==1,"sample_side"] %>% table()



sub_metadata <- metadata[metadata$uniqueID %in% sample_idx,]
sample.true.idx <- sub_metadata[sub_metadata$sample_type=="true_sample","uniqueID"]
sample.paraffin.idx <- sub_metadata[sub_metadata$sample_type=="paraffin","uniqueID"]

p.true.DNAex.1 <- subset_samples(p.true.RmLowAbun.RmHighPrev, uniqueID %in% sample.true.idx)
p.true.DNAex.1 <- prune_taxa(taxa_sums(p.true.DNAex.1)>0, p.true.DNAex.1)

physeq.nct.DNAex.1 <- subset_samples(physeq.nct.filt.s.prune, uniqueID %in% sample.paraffin.idx)
physeq.nct.DNAex.1 <- prune_taxa(taxa_sums(physeq.nct.DNAex.1)>0, physeq.nct.DNAex.1)

prev.nct.DNAex.1 <- tibble(
  my_p = prevalence(physeq.nct.DNAex.1, detection  = 0, sort = TRUE, count = FALSE),
  OTU =as.factor(names(prevalence(physeq.nct.DNAex.1 , detection = 0, sort = TRUE, count = FALSE)))
)

prev.true.DNAex.1 <- tibble(
  my_x = prevalence(p.true.DNAex.1, detection  = 0, sort = TRUE, count = TRUE),
  OTU =as.factor(names(prevalence(p.true.DNAex.1 , detection = 0, sort = TRUE, count = TRUE)))
)

Tab.merge <- merge(x=prev.nct.DNAex.1, y=prev.true.DNAex.1, by.x = "OTU", by.y = "OTU",
                   all.x = TRUE, all.y = FALSE) %>% tibble()
Tab.merge[is.na(Tab.merge)] <- 0
Tab.merge$my_x <- as.integer(Tab.merge$my_x)
n <- nrow(sample_data(p.true.DNAex.1))

Tab.merge$pval <- mapply(MyBinomTest, n, Tab.merge %>% pull(my_x), Tab.merge %>% pull(my_p))
binom.test(Tab.merge$my_x, 10, Tab.merge$my_p, alternative = "greater")

BinomTest_Wrapper(1, "DNA", "all", metadata, p.true.RmLowAbun.RmHighPrev, physeq.nct.filt.s.prune)

##--------------- ConQuR----------------------------



otu = as(otu_table(p.true.RmLowAbun.RmHLPrev), "matrix")
t_otu <- NULL
if(taxa_are_rows(p.true.RmLowAbun.RmHLPrev)){
  t_otu <- t(otu); 
  rownames(t_otu) <- colnames(otu)
  otu <- t_otu
}
otu <- as.data.frame(otu)
otu$uniqueID <- lapply(otu %>% rownames(), function(x){strsplit(x,split = ".", fixed = TRUE)[[1]][1]})
otu$uniqueID <- as.integer(otu$uniqueID)
metadata <- as.data.frame(metadata)
otu.merge <- merge(x=otu, y=metadata[, c("uniqueID", "DNAex_round", "PCR_round")],  
                   by.x = "uniqueID", by.y = "uniqueID", all.x = TRUE, all.y = FALSE)
otu.merge$DNAex_round <- as.factor(otu.merge$DNAex_round)
otu.merge$PCR_round <- as.factor(otu.merge$PCR_round)

batchid <- otu.merge[,"DNAex_round"]
covar <- otu.merge[, "PCR_round"]
batchid2 <- otu.merge[,"PCR_round"]
covar2 <- otu.merge[, "DNAex_round"]
taxa_tab <- otu.merge[,c(2:512)] 

options(warn=-1)
taxa_corrected1 = ConQuR(tax_tab=taxa_tab, batchid=batchid, covariates=covar, batch_ref="1")
taxa_corrected2 = ConQuR(tax_tab=taxa_corrected1, batchid=batchid2, covariates=covar2, batch_ref="1")

taxa_corrected2 = ConQuR(tax_tab=taxa_corrected1, batchid=batchid, covariates=covar, batch_ref="1",
                         logistic_lasso=T, quantile_type="lasso", interplt=T)


par(mfrow=c(2, 3))

Plot_PCoA(TAX=taxa_tab, factor=batchid, main="Before Correction, Bray-Curtis")
Plot_PCoA(TAX=taxa_corrected1, factor=batchid, main="ConQuR (Default), Bray-Curtis")
Plot_PCoA(TAX=taxa_corrected2, factor=batchid, main="ConQuR (Penalized), Bray-Curtis")

Plot_PCoA(TAX=taxa_tab, factor=batchid, dissimilarity="Aitch", main="Before Correction, Aitchison")
Plot_PCoA(TAX=taxa_corrected1, factor=batchid, dissimilarity="Aitch", main="ConQuR (Default), Aitchison")
Plot_PCoA(TAX=taxa_corrected2, factor=batchid, dissimilarity="Aitch", main="ConQuR (Penalized), Aitchison")

tax_table(p.true.RmLowAbun.RmHLPrev) <- t(taxa_corrected2)

p.true.RmLowAbun.RmHLPrev.rar <- rarefy_even_depth(p.true.RmLowAbun.RmHLPrev,sample.size = 1000,rngseed = 711)

compositional_analysis(p.true.RmLowAbun.RmHLPrev.rar,target_rank = "species", rel_abundance = 0.04)

p.true.2<- Subtract_Species(p.true.RmLowAbun.RmHLPrev, physeq.nct.filt.s.prune)
p.true.2.rar <- rarefy_even_depth(p.true.2,sample.size = 1000,rngseed = 911)
compositional_analysis(p.true.2.rar,target_rank = "species", rel_abundance = 0.04)

tmp <- psmelt(p.true.2.rar)


