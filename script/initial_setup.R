# Load the libraries
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(xlsx); packageVersion("xlsx")
library(tidyverse); packageVersion("tidyverse")
library(naniar); packageVersion("naniar")
library(dplyr); packageVersion("dplyr")
library(decontam); packageVersion("decontam")
"%!in%" <- Negate("%in%")

# path for data input
data_path <- "data/directory"

# load data and tables
taxa_table <- readRDS(file.path(data_path, "taxTables.rds"))
asv_table <- readRDS(file.path(data_path, "asvs.rds"))
# metadata
metadata <- read.xlsx(file.path(data_path, "metadata.xlsx"), sheetIndex = 1)

## Prepare the data to be analysed
taxa_table <- taxa_table$Species 
asv_table <- asv_table$Species

## Adjust the metadata
# set as rownames the samples IDs
rownames(metadata) <- metadata$Microomics.ID
# remove ID_sample column
metadata <- metadata[,-1]   

# modify the colnames by selecting only the middle part of the names
colnames(asv_table) <- str_replace(colnames(asv_table), pattern=".*-([^-]*)-.*", replacement = "\\1")
# prepare regular expression for digits
regexp <- "[[:digit:]]+"
# extract digits from specific column names and convert to numeric
num <- as.numeric(str_extract(colnames(asv_table)[startsWith(colnames(asv_table), "FLU")], regexp))
# adjust "FLU" samples names by adding leading zeros
colnames(asv_table)[startsWith(colnames(asv_table), "FLU")] <- paste0("FLU", sprintf("%03d", num))
rownames(metadata)[rownames(metadata) %!in% colnames(asv_table)]

# add control samples to metadata
cntrls <- colnames(asv_table)[colnames(asv_table) %!in% rownames(metadata)]
# create a dataframe with rownames = control samples and add 125 columns of NAs
cntrls_df <- as.data.frame(matrix(ncol = 125, nrow = 4), row.names = cntrls)
colnames(cntrls_df) <- colnames(metadata)
metadata <- rbind(metadata, cntrls_df)
# add to metadata a column with the four subgroups of interest
subgroups <- vector()
subgroups[startsWith(rownames(metadata), "CAP")] = "CAP"
subgroups[startsWith(rownames(metadata), "COV")] = "COV"
subgroups[startsWith(rownames(metadata), "FLU")] = "FLU"
subgroups[startsWith(rownames(metadata), "IAP")] = "IAP"
subgroups[startsWith(rownames(metadata), "CN")] = "CTRL"
subgroups[startsWith(rownames(metadata), "CM")] = "CTRL"
metadata <- cbind(metadata, subgroups)

# remove not relevant columns from metadata
metadata_columns_to_remove <- c("NA.","d")
metadata <- metadata[,-which(names(metadata) %in% metadata_columns_to_remove)]

## Calculate the sequencing depth
# divide samples from controls
asv_table_samples <- asv_table[, colnames(asv_table) %!in% cntrls]
asv_table_samples <- asv_table_samples[, colnames(asv_table_samples) %!in% duplicates]
reads_per_sample <- apply(asv_table_samples, 2, sum)
asv_table_cntrls <- asv_table[, colnames(asv_table) %in% cntrls]
reads_per_cntrls <- apply(asv_table_cntrls, 2, sum)
# save the number of reads to a file
# set a path for outputs
output_path <- "results/output/directory"
sample_ids <- paste("sample", 1:171, sep = "_")
df_seq_depth <- data.frame(Sample_ID = sample_ids, Reads = c(as.vector(reads_per_sample), as.vector(reads_per_cntrls)))
write.xlsx(df_seq_depth, file.path(output_path, "seq_depth.xlsx"))

## Create phyloseq object
# correct taxa names in the asv table
rownames(asv_table) <- gsub("\n", "", rownames(asv_table))
ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <-  tax_table(taxa_table)
samples <-  sample_data(metadata)
ps <- phyloseq(ASV, TAX, samples)

## Samples and taxa filtering
# remove all taxa with read count for all samples == 0
ps_without_zeros <- filter_taxa(ps, function(x) sum(x > 0) > 0, prune = TRUE)
# select only control samples
ps_controls <- subset_samples(ps_without_zeros, rownames(ps_without_zeros@sam_data) %in% cntrls)
# remove the control samples from the object
ps_without_controls <- subset_samples(ps_without_zeros, rownames(ps_without_zeros@sam_data) %!in% cntrls)

# remove those samples that have less than 1000 total reads
ps_without_controls <- phyloseq::subset_samples(ps_without_controls, phyloseq::sample_sums(ps_without_controls) > 1000)
ps <- phyloseq::subset_samples(ps_without_zeros, phyloseq::sample_sums(ps_without_zeros) > 1000)

# taxa filtering with decontam
# create the object with relative abundances which is the one to be used with decontam
ps_rel_abund <-  phyloseq::transform_sample_counts(ps, function(x){x / sum(x)})
# create a column in the metadata needed by decontam
sample_data(ps_rel_abund)$is.neg <- rownames(sample_data(ps_rel_abund)) %in% cntrls
contamdf.prev <- isContaminant(ps_rel_abund, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
# get the names of the taxa identified as contamination
contamination <- rownames(contamdf.prev)[contamdf.prev$contaminant==TRUE]
# determine which taxa are not contaminants
taxa_to_keep <- !(taxa_names(ps_without_controls) %in% contamination)
ps_filtered <- prune_taxa(taxa_to_keep, ps_without_controls)

# remove taxa with a read count < 30 and present in less than 2 samples
ps <- filter_taxa(ps_filtered, function(x) sum(x > 30) > 1, prune=TRUE)  # 243 taxa * 162 samples

# remove unclassified.P1 taxa, mitochondria and chloroplasts from the object
ps <- subset_taxa(ps, rownames(otu_table(ps)) != "unclassified unclassified.S160")
ps <- subset_taxa(ps, rownames(otu_table(ps)) != "unclassified unclassified.S137")
ps <- subset_taxa(ps, rownames(otu_table(ps)) != "unclassified unclassified.S1024")

# Recalculate the relative abundance
ps_rel_abund <-  phyloseq::transform_sample_counts(ps, function(x){x / sum(x)})

# save the result of the setup
save.image(file.path(output_path, "data_after_setup.RData"))