#combining dataframe 
library("readxl")

metaphyseq1 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2021_02_18_16s_kpc_tumor_invitro_est1_NEW.xlsx") %>% as.data.frame()
metaphyseq1 <- metaphyseq1 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq1 <- metaphyseq1 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

metaphyseq2 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2021_02_18_16s_kpc_tumor_invitro_est2_NEW.xlsx") %>% as.data.frame()
metaphyseq2 <- metaphyseq2 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq2 <- metaphyseq2 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

metaphyseq3 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2021_02_19_16s_kpc_panc_invitro_est1_NEW.xlsx") %>% as.data.frame()
metaphyseq3 <- metaphyseq3 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq3 <- metaphyseq3 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

metaphyseq4 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2021_02_19_16s_kpc_panc_invitro_est2_NEW.xlsx") %>% as.data.frame()
metaphyseq4 <- metaphyseq4 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq4 <- metaphyseq4 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

metaphyseq5 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2021_04_30_16s_kpc_tum-panc_invitro_NEW.xlsx") %>% as.data.frame()
metaphyseq5 <- metaphyseq5 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq5 <- metaphyseq5 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

metaphyseq6 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2021_06_28_16s_kpc_tum-panc_invitro_rep_NEW.xlsx") %>% as.data.frame()
metaphyseq6 <- metaphyseq6 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq6 <- metaphyseq6 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

metaphyseq7 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2022_02_22_16s_kpc_tum-panc_invitro_rep2_NEW.xlsx") %>% as.data.frame()
metaphyseq7 <- metaphyseq7 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq7 <- metaphyseq7 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

metaphyseq8 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2022_02_22_16s_kpc_tum-panc_invitro_rep3_NEW.xlsx") %>% as.data.frame()
metaphyseq8 <- metaphyseq8 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq8 <- metaphyseq8 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

metaphyseq9 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2022_02_23_16s_kpc_tum-panc_invitro_NEW.xlsx") %>%  as.data.frame()
metaphyseq9 <- metaphyseq9 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq9 <- metaphyseq9 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

tmp_phy_merge <- merge_samples(merge_phyloseq(physeq1,physeq2,physeq3,physeq4,physeq5,physeq6,physeq7,physeq8,physeq9),group = "uniqueID")

head(otu_table(tmp_phy_merge)) # transposed
head(tax_table(tmp_phy_merge))
head(sample_data(tmp_phy_merge)) # NA-filled

#prep sample_data
tmp_sample_data <- rbind(metaphyseq1,metaphyseq2,metaphyseq3,metaphyseq4,metaphyseq5,metaphyseq6,metaphyseq7,metaphyseq8,metaphyseq9) %>% distinct_at(.vars = "uniqueID",.keep_all = T)
rownames(tmp_sample_data) <- tmp_sample_data$id #When the argument is a data.frame, sample_data will create a sample_data-class object. In this case, the rows should be named to match the sample_names of the other objects to which it will ultimately be paired.
tmp_phy_sample <- sample_data(tmp_sample_data) #construct sample_data_class object

#prep otu_table
otu = as(otu_table(tmp_phy_merge), "matrix")
if(!taxa_are_rows(tmp_phy_merge)){otu <- t(otu)}
colnames(otu) <- tmp_sample_data$id
tmp_phy_otu <- otu_table(otu,taxa_are_rows = TRUE) # construct otu_table_class object

#prep tax_table
tmp_phy_taxa <- tax_table(tmp_phy_merge)

# combine elements and form new phyloseq
tmp_phy <- phyloseq(tmp_phy_otu,tmp_phy_taxa,tmp_phy_sample)
---
metaphyseq_ntc1 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2021_02_22_16s_kpc_tum-panc_invitro_ntc_NEW.xlsx") %>% as.data.frame()
metaphyseq_ntc1 <- metaphyseq_ntc1 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq_ntc1 <- metaphyseq_ntc1 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

metaphyseq_nct2 <- read_excel("bulk_tumor/batch_data_phyloseqs/xslx_file/2022_03_09_16s_kpc_tum-panc_invitro_ntc_NEW.xlsx") %>% as.data.frame()
metaphyseq_nct2 <- metaphyseq_nct2 %>% mutate(id = paste0(uniqueID,".",sample_side))
metaphyseq_ntc2 <- metaphyseq_ntc2 %>% mutate_at(.vars = c("DNAex_round","PCR_round"),function(x) as.character(x))

tmp_sample_data_nct <- rbind(metaphyseq_ntc1,metaphyseq_nct2)
rownames(tmp_sample_data_nct) <- tmp_sample_data_nct$id #When the argument is a data.frame, sample_data will create a sample_data-class object. In this case, the rows should be named to match the sample_names of the other objects to which it will ultimately be paired.
tmp_phy_sample <- sample_data(tmp_sample_data_nct)

