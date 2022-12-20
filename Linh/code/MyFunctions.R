# functions 
## filter by low abundance
filter_by_low_abundance <- function(phyloseqObj, thres=1e-4){
  physeq.abun.filt <- phyloseq::transform_sample_counts(phyloseqObj, function(x){x/sum(x)})
  physeq.abun.filt <- phyloseq::filter_taxa(physeq.abun.filt, function(x) mean(x) > 1e-4, TRUE)
  return(physeq.abun.filt)
}
###--------------------------------------------------------------------------------
## filtering sample from low read count
filter_samples_low_reads <- function(PhyObj, theshold=750){
  # use to remove sample from NTC whose number of reads are less than a  certain threshold
  bac.count.ntc <- tibble(bac.count = sample_sums(physeq.filt.nct.s), 
                          sample = sample_data(physeq.filt.nct.s)$id) %>% 
    arrange (bac.count)
}
##------------------------------------------

## wrapper of low abundance filtering by PERfect package
Wrapper_FERfect <- function(PhyObj){
  otu = as(otu_table(PhyObj), "matrix")
  if(taxa_are_rows(PhyObj)){otu <- t(otu)}
  otu = as_tibble(otu)
  res_sim <- PERFect_sim(X = otu)
  # dim(res_sim$filtX) 
  ids.sim<- colnames(res_sim$filtX) 
  return(prune_taxa(ids.sim, p.true.filt.s))
}

##--------------------------------------

## determine a threshold for high prevalence filtering.
## return an array, how many percent of taxa in true sample I remove if I choose that threshold 
HighPrevalence_Data <- function(True.Sample, NCT.Sample){
  num_taxa <- nrow(otu_table(True.Sample))
  prev.nct <- tibble(
    prev = prevalence(NCT.Sample, detection  = 0, sort = TRUE, count = FALSE),
    OTU =as.factor(names(prevalence(NCT.Sample , detection = 0, sort = TRUE, count = FALSE)))
  )
  x <- prev.nct$prev
  factorx <- factor(cut(x, breaks=nclass.Sturges(x)))
  xout <- as.data.frame(table(factorx)) %>% map_df(rev)
  xout <- mutate(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))
  xout <- xout %>% mutate(rel_tax=cumFreq/num_taxa)
  return(xout[,c(1,5)])
}

##--------------------------------------
## filter by high prevalence in NTC

