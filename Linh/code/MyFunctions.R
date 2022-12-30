# functions 
inspect_physeqObj <- function(Obj){
  print(nrow(sample_data(Obj)))
  print(nrow(otu_table(Obj)))
}
##---------------

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
  # subset of taxa which in in overlap
  idx_overlap <- intersect(rownames(otu_table(True.Sample)) %>% unique(),
                           rownames(otu_table(NCT.Sample)) %>% unique())
  NCT.Sample.Subset <- prune_taxa(idx_overlap, NCT.Sample)
  prev.nct <- tibble(
    prev = prevalence(NCT.Sample.Subset, detection  = 0, sort = TRUE, count = FALSE),
    OTU =as.factor(names(prevalence(NCT.Sample.Subset , detection = 0, sort = TRUE, count = FALSE)))
  )
  x <- prev.nct$prev
  factorx <- factor(cut(x, breaks=nclass.Sturges(x)))
  xout <- as.data.frame(table(factorx)) %>% map_df(rev)
  xout <- mutate(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))
  xout <- xout %>% mutate(rel_tax=cumFreq/num_taxa)
  return(xout[,c(1,5)])
}

##--------------------------------------
## check if a batch of condition is a mixture of true sample and NCT
is.valid.batch <- function(batch_num, batch_type, meta_data ){
  
  tmp <- NULL
  if(grepl("dna", tolower(batch_type), fixed = TRUE)){
    tmp <- meta_data[meta_data$DNAex_round==batch_num, "TRUE_control"] %>% table() %>% length()
  }else if(grepl("pcr", tolower(batch_type), fixed = TRUE)){
    tmp <- meta_data[meta_data$PCR_round==batch_num, "TRUE_control"] %>% table() %>% length()
  }else{
    stop("Batch type must contain either dna or pcr.")
  }
  return(tmp==2)
}
###----------------------------------------------
##------------------------------------------------
## Binomial testing
MyBinomTest <- function(n,x,p){ 
  return(binom.test(x,n,p,alternative = "greater")$p.value)
}
#---------------------------------------------
## when the batch_type is fixed, extract true samples and all different kind of controls.
Extract_SampleIdx <- function(batch_num, batch_type, control_type="all", meta_data){
  sub_meta_data <- NULL
  if(grepl("dna", tolower(batch_type), fixed = TRUE)){
    sub_meta_data <- meta_data[meta_data$DNAex_round == batch_num,]
  }else if(grepl("pcr", tolower(batch_type), fixed = TRUE)){
    sub_meta_data <- meta_data[meta_data$PCR_round == batch_num,]
  }else{
    stop("Batch type must contain either dna or pcr.")
  }
  Sample.TRUE.idx <- sub_meta_data[sub_meta_data$sample_type=="true_sample","uniqueID"]
  Sample.NCT.idx <- NULL
  if(grepl("paraffin", tolower(control_type), fixed = TRUE)){

  }else if(grepl("buffer", tolower(control_type), fixed = TRUE)){

  }else if(grepl("pcr", tolower(control_type), fixed = TRUE)){

  }else if(grepl("all", tolower(control_type), fixed = TRUE)){
    Sample.NCT.idx <- sub_meta_data[sub_meta_data$sample_type != "true_sample", "uniqueID"]
  }
  else{stop("control type must be one of those: paraffin, buffer, pcr, all")}
  return(list("TrueSam"=Sample.TRUE.idx, "NCTSam"=Sample.NCT.idx))
}
##-----------------------
BinomTest_True.vs.NCT <- function(True.Sample, NCT.Sample){
  prev.NCT <- tibble(
    my_p = prevalence(NCT.Sample, detection  = 0, sort = TRUE, count = FALSE),
    OTU =as.factor(names(prevalence(NCT.Sample , detection = 0, sort = TRUE, count = FALSE)))
  )
  
  prev.TRUE <- tibble(
    my_x = prevalence(True.Sample, detection  = 0, sort = TRUE, count = TRUE),
    OTU =as.factor(names(prevalence(True.Sample , detection = 0, sort = TRUE, count = TRUE)))
  )
  
  Tab.merge <- merge(x=prev.NCT, y=prev.TRUE, by.x = "OTU", by.y = "OTU",
                     all.x = TRUE, all.y = FALSE) %>% tibble()
  Tab.merge[is.na(Tab.merge)] <- 0
  Tab.merge$my_x <- as.integer(Tab.merge$my_x)
  n <- nrow(sample_data(True.Sample))
  
  Tab.merge$pval <- mapply(MyBinomTest, n, Tab.merge %>% pull(my_x), Tab.merge %>% pull(my_p))
  return(list("Contaminant" = Tab.merge[Tab.merge$pval > 0.05,] %>% pull("OTU"), 
              "TrueSpecies" = Tab.merge[Tab.merge$pval <= 0.05,] %>% pull("OTU")))
}
##----------------------------------------------------------------
## when the condition is fixed
BinomTest_Wrapper <- function(batch_num, batch_type, control_type="all", meta_data, True.Sample, NCT.Sample){
  if(!is.valid.batch(batch_num, batch_type, meta_data)){
    print("In this batch with this condition, we have either all true sample or all control.")
    return(NULL)
  }
  Lst_idx <- Extract_SampleIdx(batch_num, batch_type,"all", meta_data)
  
  tmp_idx_true <- Lst_idx$TrueSam
  sub.True.Sample <- prune_samples(sample_data(True.Sample) %>% pull(uniqueID) %in% tmp_idx_true ,True.Sample)
  sub.True.Sample <- prune_taxa(taxa_sums(sub.True.Sample)>0, sub.True.Sample)
  
  tmp_idx_nct <- Lst_idx$NCTSam
  sub.NCT.Sample <- prune_samples(sample_data(NCT.Sample) %>% pull(uniqueID) %in% tmp_idx_nct ,NCT.Sample)
  sub.NCT.Sample <- prune_taxa(taxa_sums(sub.NCT.Sample)>0, sub.NCT.Sample)
  
  return(BinomTest_True.vs.NCT(sub.True.Sample, sub.NCT.Sample))
}
##----------------------------------------------------------
Wrapper_filter_by_low_prevalence <- function(True.Sample, NCT.Sample, metadata){
  ## filter taxa in true sample when their prevalence is low in NCT samples.
  tmp <- metadata$DNAex_round %>% as.numeric()
  DNAex_lst <- tmp[!is.na(tmp)] %>% unique()
  tmp <- metadata$PCR_round %>% as.numeric()
  PCR_lst <- tmp[!is.na(tmp)] %>% unique()
  
  Lst_contaminant = c() # list of contaminant taxa, empty at the beginning.
  tmp_lst <- lapply(DNAex_lst, BinomTest_Wrapper, batch_type="DNA", control_type="all", meta_data=metadata, 
                    True.Sample=True.Sample, NCT.Sample=NCT.Sample)
  for(item in tmp_lst){
    if(!is.null(item)){
      Lst_contaminant <- c(Lst_contaminant, item$Contaminant %>% levels())
    }
  }
  Lst_contaminant <- Lst_contaminant %>% unique()
  ## Second, we fix the condition PCR
  tmp_lst <- lapply(PCR_lst, BinomTest_Wrapper, batch_type="PCR", control_type="all", meta_data=metadata, 
                    True.Sample=True.Sample, NCT.Sample=NCT.Sample)
  for(item in tmp_lst){
    if(!is.null(item)){
      Lst_contaminant <- c(Lst_contaminant, item$Contaminant %>% levels())
    }
  }
  Lst_contaminant <- Lst_contaminant %>% unique()
  keep_taxa <- setdiff(rownames(otu_table(p.true.RmLowAbun.RmHighPrev)), Lst_contaminant)
  return(prune_taxa(keep_taxa, p.true.RmLowAbun.RmHighPrev))
}
##-------------------------------------------------------------
Wrapper_BatchEffectCorrection_ConQuR <- function(True.Sample, metadata, batchid, covar){
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
  
  taxa_corrected1 = ConQuR(tax_tab=taxa_tab, batchid=batchid, covariates=covar, batch_ref="1")
}
##----------------------------------------
Subtract_Species <- function(Sample1, Sample2){
  set_diff <- setdiff(Sample1 %>% otu_table() %>% rownames() %>% unique(),
    Sample2 %>% otu_table() %>% rownames() %>% unique()
  )
  return(prune_taxa(set_diff, Sample1))
}


