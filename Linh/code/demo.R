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
