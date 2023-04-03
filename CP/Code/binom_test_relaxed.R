phy <- p.comb.filt.s.1500.80.r
phy
Rm_taxa <- list()
# for DNAex
IDv <- sample_data(phy)$DNAex_round %>% unique()
for (i in IDv) {
  b <- subset_samples(phy, DNAex_round == i)
  check.true <- length(which(meta(b)$TRUE_control == "TRUE"))
  check.nct <- length(which(meta(b)$TRUE_control == "control"))
  if (check.true == 0) {
    print(paste0("No TRUE samples found in DNAex batch: ","",i,""))    
    next}
  if (check.nct == 0) {
    print(paste0("No NCT samples found in DNAex batch: ","",i,""))
    next}
  b <- prune_taxa(taxa = taxa_sums(b)>0,x = b) #this removes taxa with taxa_sums == 0 from the overall batch
  b.true <- subset_samples(b,TRUE_control == "TRUE") #subset of true samples in batch
  b.nct <- subset_samples(b,TRUE_control == "control") #subset of nct samples in batch
  b.sub <- merge_phyloseq(b.true,b.nct) # we take the subsets we want, merge them to allow for thorough removal of bacteria
  if (any(taxa_sums(b.sub) < 1)) {
    b.sub <- prune_taxa(taxa = taxa_sums(b.sub)>0,x = b.sub) #this removes taxa with taxa_sums == 0 from our batch subset
  }
  b.true <- subset_samples(b.sub,TRUE_control == "TRUE") #create new subset of true samples in batch
  b.nct <- subset_samples(b.sub,TRUE_control == "control") #create new subset of nct samples in batch  
  # Im honest: this is a bit redundant
  df.filt <- data.frame(OTU = taxa_names(b.true),
                        x = prevalence(b.true, detection = 0, count = TRUE),
                        n = nsamples(b.true),
                        p = prevalence(b.nct, detection = 0, count = FALSE)) #df stores OTU handle - x,n,p for binom.test
  vec <- seq(1,nrow(df.filt))
  for (j in vec) {
    vec[j] <- binom.test(x = df.filt[j,"x"], n = df.filt[j,"n"], p = df.filt[j,"p"],alternative = "greater")$p.value
  }
  df.filt$pval <- vec
  Rm_taxa <- list.append(Rm_taxa,df.filt[df.filt$pval > 0.05,"OTU"])
}

Rm_taxa <- unique(unlist(Rm_taxa))
length(Rm_taxa)

# for PCR
IDv <- sample_data(phy)$PCR_round %>% unique()
for (i in IDv) {
  b <- subset_samples(phy, PCR_round == i)
  check.true <- length(which(meta(b)$TRUE_control == "TRUE"))
  check.nct <- length(which(meta(b)$TRUE_control == "control"))
  if (check.true == 0) {
    print(paste0("No TRUE samples found in PCR batch: ","",i,""))    
    next}
  if (check.nct == 0) {
    print(paste0("No NCT samples found in PCR batch: ","",i,""))
    next}
  b <- prune_taxa(taxa = taxa_sums(b)>0,x = b) #this removes taxa with taxa_sums == 0 from the overall batch
  b.true <- subset_samples(b,TRUE_control == "TRUE") #subset of true samples in batch
  b.nct <- subset_samples(b,TRUE_control == "control") #subset of nct samples in batch
  b.sub <- merge_phyloseq(b.true,b.nct) # we take the subsets we want, merge them to allow for thorough removal of bacteria
  if (any(taxa_sums(b.sub) < 1)) {
    b.sub <- prune_taxa(taxa = taxa_sums(b.sub)>0,x = b.sub) #this removes taxa with taxa_sums == 0 from our batch subset
  }
  b.true <- subset_samples(b.sub,TRUE_control == "TRUE") #create new subset of true samples in batch
  b.nct <- subset_samples(b.sub,TRUE_control == "control") #create new subset of nct samples in batch  
  # Im honest: this is a bit redundant
  df.filt <- data.frame(OTU = taxa_names(b.true),
                        x = prevalence(b.true, detection = 0, count = TRUE),
                        n = nsamples(b.true),
                        p = prevalence(b.nct, detection = 0, count = FALSE)) #df stores OTU handle - x,n,p for binom.test
  vec <- seq(1,nrow(df.filt))
  for (j in vec) {
    vec[j] <- binom.test(x = df.filt[j,"x"], n = df.filt[j,"n"], p = df.filt[j,"p"],alternative = "greater")$p.value
  }
  df.filt$pval <- vec
  Rm_taxa <- list.append(Rm_taxa,df.filt[df.filt$pval > 0.05,"OTU"])
}

Rm_taxa <- unique(unlist(Rm_taxa))
length(Rm_taxa)