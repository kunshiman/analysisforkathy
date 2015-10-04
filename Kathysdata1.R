library(WriteXLS) # uses Perl
library(readxl)
library(biomaRt) # bioconductor
DPKD <- read_excel("DPKD_Ctr.xlsx") # reads the xlsx file
# make a list with corresponding ensembl genes
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- DPKD$GeneSymbol
Gene_list <- getBM(filters= "external_gene_name", 
                   attributes= c("external_gene_name", "ensembl_gene_id"), 
                   values=genes,mart= mart)
detach("package:biomaRt", unload=TRUE) # needs to be detached 
# as dplyr select function would not work

library(dplyr)
DPKD2 <- filter(DPKD, Status == "OK") # removes all genes with low read quality
# renames the following column names for easier handling in R
dpkd0 <- tolower(gsub(" ", "", names(DPKD2)))
names(DPKD2) = dpkd0
names(DPKD2)[names(DPKD2) == "fdradjustedpvalue"] <- "FDR"
# transforms the table and removes unwanted columns
DPKD2 <- select(DPKD2, genesymbol, log2foldchange, pvalue, FDR, description)
DPKD2 <- mutate(DPKD2, foldchange = 2^log2foldchange)
DPKD2$foldchange <- ifelse(DPKD2$foldchange < 1, -(1/DPKD2$foldchange), DPKD2$foldchange)
DPKD1 <- arrange(DPKD2, desc(foldchange), genesymbol)
DPKD2 <- merge(DPKD1, Gene_list, by.x = "genesymbol", by.y = "external_gene_name", all.x = TRUE)
DPKD2 <- select(DPKD2, genesymbol, ensembl_gene_id, foldchange, pvalue, FDR, description)

# makes 2 more sheets with different fold change and pvalue:
DPKD3 <- filter(DPKD2, foldchange >= 2 | foldchange <= -2, FDR <= 0.05)
DPKD4 <- filter(DPKD2, foldchange >= 5 | foldchange <= -5, FDR <= 0.05)
# writes everything to the xlsx file:
WriteXLS(c("DPKD2", "DPKD3", "DPKD4"), "DPKD_Ctr_simple.xlsx", SheetNames = c("FC>=1;FDR<=1", "FC>=2;FDR<=0.05", "FC>=5;FDR<=0.05"), row.names = FALSE, AdjWidth = TRUE)
rm(DPKD3); rm(DPKD4)

# same protocol for PGKD_Ctr:
PGKD <- read_excel("PGKD_Ctr.xlsx")
PGKD2 <- filter(PGKD, Status == "OK")
names(PGKD2) = dpkd0
names(PGKD2)[names(PGKD2) == "fdradjustedpvalue"] <- "FDR"

PGKD2 <- select(PGKD2, genesymbol, log2foldchange, pvalue, FDR, description)
PGKD2 <- mutate(PGKD2, foldchange = 2^log2foldchange)
PGKD2$foldchange <- ifelse(PGKD2$foldchange < 1, -(1/PGKD2$foldchange), PGKD2$foldchange)
PGKD1 <- arrange(PGKD2, desc(foldchange), genesymbol)
PGKD2 <- merge(PGKD1, Gene_list, by.x = "genesymbol", by.y = "external_gene_name", all.x = TRUE)
PGKD2 <- select(PGKD2, genesymbol, ensembl_gene_id, foldchange, pvalue, FDR, description)


PGKD3 <- filter(PGKD2, foldchange >= 2 | foldchange <= -2, FDR <= 0.05)
PGKD4 <- filter(PGKD2, foldchange >= 5 | foldchange <= -5, FDR <= 0.05)

WriteXLS(c("PGKD2", "PGKD3", "PGKD4"), "PGKD_Ctr_simple.xlsx", SheetNames = c("FC>=1;FDR<=1", "FC>=2;FDR<=0.05", "FC>=5;FDR<=0.05"), row.names = FALSE, AdjWidth = TRUE)
rm(PGKD3); rm(PGKD4)

# same for PKP2_Ctr
PKP2KD <- read_excel("PKP2KD_Ctr.xlsx")
PKP2KD2 <- filter(PKP2KD, Status == "OK")
names(PKP2KD2) = dpkd0
names(PKP2KD2)[names(PKP2KD2) == "fdradjustedpvalue"] <- "FDR"

PKP2KD2 <- select(PKP2KD2, genesymbol, log2foldchange, pvalue, FDR, description)
PKP2KD2 <- mutate(PKP2KD2, foldchange = 2^log2foldchange)
PKP2KD2$foldchange <- ifelse(PKP2KD2$foldchange < 1, -(1/PKP2KD2$foldchange), PKP2KD2$foldchange)
PKP2KD1 <- arrange(PKP2KD2, desc(foldchange), genesymbol)
PKP2KD2 <- merge(PKP2KD1, Gene_list, by.x = "genesymbol", by.y = "external_gene_name", all.x = TRUE)
PKP2KD2 <- select(PKP2KD2, genesymbol, ensembl_gene_id, foldchange, pvalue, FDR, description)

PKP2KD3 <- filter(PKP2KD2, foldchange >= 2 | foldchange <= -2, FDR <= 0.05)
PKP2KD4 <- filter(PKP2KD2, foldchange >= 5 | foldchange <= -5, FDR <= 0.05)

WriteXLS(c("PKP2KD2", "PKP2KD3", "PKP2KD4"), "PKP2KD_Ctr_simple.xlsx", 
         SheetNames = c("FC>=1;FDR<=1", "FC>=2;FDR<=0.05", "FC>=5;FDR<=0.05"), 
         row.names = FALSE, AdjWidth = TRUE)
rm(PKP2KD3); rm(PKP2KD4)

# merge the tables and renames the columns
DP_PG_PKP2_merge <- merge(DPKD1, merge(PGKD1, PKP2KD1, by = "genesymbol"), by = "genesymbol")
DP_PG_PKP2_merge <- select(DP_PG_PKP2_merge, genesymbol, foldchange, FDR, foldchange.x, FDR.x, 
                           foldchange.y, FDR.y, description)
DP_PG_PKP2_merge <- merge(DP_PG_PKP2_merge, Gene_list, by.x = "genesymbol", by.y = "external_gene_name", all.x = TRUE)
DP_PG_PKP2_merge <- select(DP_PG_PKP2_merge, genesymbol, ensembl_gene_id, foldchange, FDR, foldchange.x, FDR.x, 
                           foldchange.y, FDR.y, description)
# renames the columns and merges the data
names(DP_PG_PKP2_merge)[names(DP_PG_PKP2_merge) == "foldchange"] <- "FC_DPKD"
names(DP_PG_PKP2_merge)[names(DP_PG_PKP2_merge) == "foldchange.x"] <- "FC_PGKD"
names(DP_PG_PKP2_merge)[names(DP_PG_PKP2_merge) == "foldchange.y"] <- "FC_PKP2KD"
names(DP_PG_PKP2_merge)[names(DP_PG_PKP2_merge) == "FDR"] <- "FDR_DPKD"
names(DP_PG_PKP2_merge)[names(DP_PG_PKP2_merge) == "FDR.x"] <- "FDR_PGKD"
names(DP_PG_PKP2_merge)[names(DP_PG_PKP2_merge) == "FDR.y"] <- "FDR_PKP2KD"
DP_PG_PKP2_merge1 <- filter(DP_PG_PKP2_merge,
                           FC_DPKD >= 2 | FC_PGKD >= 2 | FC_PKP2KD >= 2 | FC_DPKD <= (-2) | FC_PGKD <= (-2) | FC_PKP2KD <= (-2))

WriteXLS(c("DP_PG_PKP2_merge1"), "DP_PG_PKP2vsCtr.xlsx", 
         SheetNames = "oneFC>=2", 
         row.names = FALSE, AdjWidth = TRUE)
rm(PKP2KD); rm(DPKD); rm(PGKD); rm(DPKD1); rm(PGKD1); rm(PKP2KD1)
####################################################################
# for PKP2KD_DPKD
PKP2KD_DPKD <- read_excel("PKP2KD_DPKD.xlsx") # reads the xlsx file
PKP2KD_DPKD2 <- filter(PKP2KD_DPKD, Status == "OK") # removes all genes with low read quality
# renames the following column names for easier handling in R
names(PKP2KD_DPKD2) = dpkd0
names(PKP2KD_DPKD2)[names(PKP2KD_DPKD2) == "fdradjustedpvalue"] <- "FDR"
# transforms the table and removes unwanted columns
PKP2KD_DPKD2 <- select(PKP2KD_DPKD2, genesymbol, log2foldchange, pvalue, FDR, description)
PKP2KD_DPKD2 <- mutate(PKP2KD_DPKD2, foldchange = 2^log2foldchange)
PKP2KD_DPKD2$foldchange <- ifelse(PKP2KD_DPKD2$foldchange < 1, -(1/PKP2KD_DPKD2$foldchange), PKP2KD_DPKD2$foldchange)
PKP2KD_DPKD1 <- arrange(PKP2KD_DPKD2, desc(foldchange), genesymbol)
PKP2KD_DPKD2 <- merge(PKP2KD_DPKD1, Gene_list, by.x = "genesymbol", by.y = "external_gene_name", all.x = TRUE)
PKP2KD_DPKD2 <- select(PKP2KD_DPKD2, genesymbol, ensembl_gene_id, foldchange, pvalue, FDR, description)

# makes 2 more sheets with different fold change and pvalue:
PKP2KD_DPKD3 <- filter(PKP2KD_DPKD2, foldchange >= 2 | foldchange <= -2, FDR <= 0.05)
PKP2KD_DPKD4 <- filter(PKP2KD_DPKD2, foldchange >= 5 | foldchange <= -5, FDR <= 0.05)
# writes everything to the xlsx file:
WriteXLS(c("PKP2KD_DPKD2", "PKP2KD_DPKD3", "PKP2KD_DPKD4"), "PKP2KD_DPKD_simple.xlsx", SheetNames = c("FC>=1;FDR<=1", "FC>=2;FDR<=0.05", "FC>=5;FDR<=0.05"), row.names = FALSE, AdjWidth = TRUE)
rm(PKP2KD_DPKD3); rm(PKP2KD_DPKD4)

PKP2KD_PGKD <- read_excel("PKP2KD_PGKD.xlsx") # reads the xlsx file
PKP2KD_PGKD2 <- filter(PKP2KD_PGKD, Status == "OK") # removes all genes with low read quality
# renames the following column names for easier handling in R
names(PKP2KD_PGKD2) = dpkd0
names(PKP2KD_PGKD2)[names(PKP2KD_PGKD2) == "fdradjustedpvalue"] <- "FDR"
# transforms the table and removes unwanted columns
PKP2KD_PGKD2 <- select(PKP2KD_PGKD2, genesymbol, log2foldchange, pvalue, FDR, description)
PKP2KD_PGKD2 <- mutate(PKP2KD_PGKD2, foldchange = 2^log2foldchange)
PKP2KD_PGKD2$foldchange <- ifelse(PKP2KD_PGKD2$foldchange < 1, -(1/PKP2KD_PGKD2$foldchange), PKP2KD_PGKD2$foldchange)
PKP2KD_PGKD1 <- arrange(PKP2KD_PGKD2, desc(foldchange), genesymbol)
PKP2KD_PGKD2 <- merge(PKP2KD_PGKD1, Gene_list, by.x = "genesymbol", by.y = "external_gene_name", all.x = TRUE)
PKP2KD_PGKD2 <- select(PKP2KD_PGKD2, genesymbol, ensembl_gene_id, foldchange, pvalue, FDR, description)

# makes 2 more sheets with different fold change and pvalue:
PKP2KD_PGKD3 <- filter(PKP2KD_PGKD2, foldchange >= 2 | foldchange <= -2, FDR <= 0.05)
PKP2KD_PGKD4 <- filter(PKP2KD_PGKD2, foldchange >= 5 | foldchange <= -5, FDR <= 0.05)
# writes everything to the xlsx file:
WriteXLS(c("PKP2KD_PGKD2", "PKP2KD_PGKD3", "PKP2KD_PGKD4"), "PKP2KD_PGKD_simple.xlsx", SheetNames = c("FC>=1;FDR<=1", "FC>=2;FDR<=0.05", "FC>=5;FDR<=0.05"), row.names = FALSE, AdjWidth = TRUE)
rm(PKP2KD_PGKD3); rm(PKP2KD_PGKD4)

PGKD_DPKD <- read_excel("PGKD_DPKD.xlsx") # reads the xlsx file
PGKD_DPKD2 <- filter(PGKD_DPKD, Status == "OK") # removes all genes with low read quality
# renames the following column names for easier handling in R
names(PGKD_DPKD2) = dpkd0
names(PGKD_DPKD2)[names(PGKD_DPKD2) == "fdradjustedpvalue"] <- "FDR"
# transforms the table and removes unwanted columns
PGKD_DPKD2 <- select(PGKD_DPKD2, genesymbol, log2foldchange, pvalue, FDR, description)
PGKD_DPKD2 <- mutate(PGKD_DPKD2, foldchange = 2^log2foldchange)
PGKD_DPKD2$foldchange <- ifelse(PGKD_DPKD2$foldchange < 1, -(1/PGKD_DPKD2$foldchange), PGKD_DPKD2$foldchange)
PGKD_DPKD1 <- arrange(PGKD_DPKD2, desc(foldchange), genesymbol)
PGKD_DPKD2 <- merge(PGKD_DPKD1, Gene_list, by.x = "genesymbol", by.y = "external_gene_name", all.x = TRUE)
PGKD_DPKD2 <- select(PGKD_DPKD2, genesymbol, ensembl_gene_id, foldchange, pvalue, FDR, description)

# makes 2 more sheets with different fold change and pvalue:
PGKD_DPKD3 <- filter(PGKD_DPKD2, foldchange >= 2 | foldchange <= -2, FDR <= 0.05)
PGKD_DPKD4 <- filter(PGKD_DPKD2, foldchange >= 5 | foldchange <= -5, FDR <= 0.05)
# writes everything to the xlsx file:
WriteXLS(c("PGKD_DPKD2", "PGKD_DPKD3", "PGKD_DPKD4"), "PGKD_DPKD_simple.xlsx", SheetNames = c("FC>=1;FDR<=1", "FC>=2;FDR<=0.05", "FC>=5;FDR<=0.05"), row.names = FALSE, AdjWidth = TRUE)
rm(PGKD_DPKD3); rm(PGKD_DPKD4);
#####
# merge the tables and renames the columns
merge_2 <- merge(PKP2KD_DPKD1, merge(PKP2KD_PGKD1, PGKD_DPKD1, by = "genesymbol"), by = "genesymbol")
merge_2 <- select(merge_2, genesymbol, foldchange, FDR, foldchange.x, FDR.x, 
                           foldchange.y, FDR.y, description)
merge_2 <- merge(merge_2, Gene_list, by.x = "genesymbol", by.y = "external_gene_name", all.x = TRUE)
merge_2 <- select(merge_2, genesymbol, ensembl_gene_id, foldchange, FDR, foldchange.x, FDR.x, 
                           foldchange.y, FDR.y, description)
# renames the columns and merges the data
names(merge_2)[names(merge_2) == "foldchange"] <- "FC_PKP2KD_DPKD"
names(merge_2)[names(merge_2) == "foldchange.x"] <- "FC_PKP2KD_PGKD"
names(merge_2)[names(merge_2) == "foldchange.y"] <- "FC_PGKD_DPKD"
names(merge_2)[names(merge_2) == "FDR"] <- "FDR_PKP2KD_DPKD"
names(merge_2)[names(merge_2) == "FDR.x"] <- "FDR_PKP2KD_PGKD"
names(merge_2)[names(merge_2) == "FDR.y"] <- "FDR_PGKD_DPKD"

all <- merge(DP_PG_PKP2_merge, merge_2)
all <- filter(all, FC_PKP2KD_DPKD >= 2 | FC_PKP2KD_PGKD >= 2 |
                FC_PGKD_DPKD >= 2 | FC_PKP2KD_DPKD <= (-2) | FC_PKP2KD_PGKD <= (-2) | FC_PGKD_DPKD <= (-2) |
                FC_DPKD >= 2 | FC_PGKD >= 2 | FC_PKP2KD >= 2 | FC_DPKD <= (-2) | FC_PGKD <= (-2) | FC_PKP2KD <= (-2))
all <- select(all, genesymbol, ensembl_gene_id, FC_DPKD, FC_PGKD, FC_PKP2KD, FC_PKP2KD_DPKD, FC_PKP2KD_PGKD, FC_PGKD_DPKD, description)
all <- arrange(all, desc(FC_DPKD))

WriteXLS(c("all"), "Allmerged.xlsx", 
         SheetNames = "oneFC>=2", 
         row.names = FALSE, AdjWidth = TRUE)

