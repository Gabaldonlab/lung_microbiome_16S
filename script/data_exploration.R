# load the needed libraries
library(dada2); packageVersion("dada2")            # sequence processing
library(phyloseq); packageVersion("phyloseq")     # v. 1.44.0 # data visualization and analysis
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")        # data visualization and analysis
library(DESeq2)
library(microbiome); packageVersion("microbiome")     
library(vegan); packageVersion("vegan")                          
library(picante); packageVersion("picante")                       
#library(ALDEx2); packageVersion("ALDEx2")                        
library(metagenomeSeq); packageVersion("metagenomeSeq")          
library(HMP); packageVersion("HMP")                              
library(dendextend); packageVersion("dendextend")                
library(selbal); packageVersion("selbal")   
library(rms); packageVersion("rms")
library(breakaway); packageVersion("breakaway")
library(xlsx)
library(tidyverse)
library(naniar)
library(ggpubr)
library(ggsci)
library(dplyr)
"%!in%" <- Negate("%in%")
library(zCompositions)
library(CoDaSeq)
library(robCompositions)

## Upload data and tables
taxa_table <- readRDS("./results/taxTables.rds")
asv_table <- readRDS("./results/asvs.rds")
# metadata
metadata <- read.xlsx("./data/metadata.xlsx", sheetIndex = 1)

# taxa table
taxa_table <- taxa_table$Species 
asv_table <- asv_table$Species

# save sample name to re-add as row name after substituting the "." by NA since the function remove rownames
ID_sample <- metadata$Microomics.ID# metadata$Patient.ID
# re-assign rownames
rownames(metadata) <- ID_sample
metadata <- metadata[,-1]   # remove ID_sample column

# modify the colnames by selecting only the middle part of the column names
colnames(asv_table) <- str_replace(colnames(asv_table), pattern=".*-([^-]*)-.*", replacement = "\\1")
# prepare regular expression for digits
regexp <- "[[:digit:]]+"
# process string 
prova <- as.numeric(str_extract(colnames(asv_table)[startsWith(colnames(asv_table), "FLU")], regexp)) #%>%

length(print(paste0("FLU", sprintf("%03d", prova))))
# adjust "FLU" samples names by adding leading zeros
colnames(asv_table)[startsWith(colnames(asv_table), "FLU")] <- paste0("FLU", sprintf("%03d", prova))
rownames(metadata)[rownames(metadata) %!in% colnames(asv_table)] # "IAP135"
# add control samples to metadata
cntrls <- colnames(asv_table)[colnames(asv_table) %!in% rownames(metadata)]
dim(metadata) # 167 * 125
# create a dataframe with rownames = control samples and add 125 columns of NAs
cntrls_df <- as.data.frame(matrix(ncol = 125, nrow = 4), row.names = cntrls)
colnames(cntrls_df) <- colnames(metadata)
metadata <- rbind(metadata, cntrls_df)
dim(metadata) # 171*125
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

# correct taxa names in asv_table
rownames(asv_table) <- gsub("\n", "", rownames(asv_table))
ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <-  tax_table(taxa_table)
samples <-  sample_data(metadata)
ps <- phyloseq(ASV, TAX, samples)

## Taxonomic filtering of the taxa
# select only control samples
ps_controls <- subset_samples(ps, rownames(ps@sam_data) %in% cntrls)
# select all but control samples
ps_without_controls <- subset_samples(ps, rownames(ps@sam_data) %!in% cntrls)
# filter taxa
ps_controls <- filter_taxa(ps_controls, function(x) sum(x > 0) > 0, prune = TRUE)
ps_without_controls <- filter_taxa(ps_without_controls, function(x) sum(x > 0) > 0, prune = TRUE)
table(rownames(ps_controls@otu_table) %in% rownames(ps_without_controls@otu_table))

# Create dataframes with each group
otus_control <- as.data.frame(ps_controls@otu_table)
# Add the maximum number of reads per each taxa in the negative controls
otus_control$max <- apply(otus_control, 1, max)
otus_control
otus_samples <- as.data.frame(ps_without_controls@otu_table)

# Filtering of the taxa according to the maximum number of reads in the negative controls
for (e in 1:nrow(otus_control)) {
  tax_id_c <- rownames(otus_control)[e]
  taxa_reads_c <- otus_control[e,"max"]
  
  for (i in 1:nrow(otus_samples)) {
    tax_id_s <- rownames(otus_samples)[i]
    
    if (tax_id_c == tax_id_s) {
      
      for (s in colnames(otus_samples)) {
        taxa_reads_s <- otus_samples[tax_id_s,s]
        # if reads in samples is < than reads in cntrl, change it to 0
        if (taxa_reads_s < taxa_reads_c) {
          otus_samples[tax_id_s,s] <- 0
          #  print ("changed to 0")
        }
      }
    }
  }
}


#Remove the taxa that are 0 in all the rows:
taxa_samples_f <- otus_samples
taxa_samples_f <- taxa_samples_f[apply(taxa_samples_f[,-1], 1, function(x) !all(x==0)),]
otus_samples[rownames(otus_samples) %!in% rownames(taxa_samples_f),]


# Change the the asv and taxtables from the phyloseq object with the filtered ones
ps_without_controls@otu_table <- otu_table(taxa_samples_f, taxa_are_rows = TRUE)
ps_without_controls@tax_table <- ps_without_controls@tax_table[rownames(ps_without_controls@tax_table) %in% rownames(ps_without_controls@otu_table), ]

# remove taxa with less than 30 reads appearing in less than 2 samples
ps <- filter_taxa(ps_without_controls, function(x) sum(x > 30) > 1, prune=TRUE)

# remove those samples that have less than 1000 total reads
ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 1000)

# Remove unclassified phyla and mitochondria
ps <- subset_taxa(ps, rownames(otu_table(ps)) != "unclassified unclassified.S160")
ps <- subset_taxa(ps, rownames(otu_table(ps)) != "unclassified unclassified.S137")

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# List of phyloseq objects agglomerated at different ranks 
# Species is excluded because it raises an error, will add it manually
ps_list_glom <- sapply(tax_ranks,
                       FUN = function(rank) {
                         new_ps <- tax_glom(ps_filt, rank)
                         # Rename taxa with the according rank
                         taxa_names(new_ps) <- unlist(as.vector(
                           tax_table(new_ps)[, rank]))
                         new_ps
                       }, simplify = FALSE)

# List of phyloseq objects with 0-replaced data at different ranks
ps_list_czm <- lapply(ps_list_glom[-1],
                      FUN = function(ps) {
                        # Impute with the CZM method
                        # Transpose bc cmultRepl assumes samples on rows
                        # and transpose back to avoid errors later
                        old_otu_tab <- otu_table(ps)
                        czm_test <- try(t(cmultRepl(t(old_otu_tab),
                                                    method = "CZM", label = 0)))
                        # Errors arise if there are no 0s
                        if (inherits(czm_test, "try-error")) {
                          # Return original data if no 0s
                          new_otu_tab <- old_otu_tab
                        } else {
                          # Return imputed data
                          new_otu_tab <- czm_test
                        }
                        new_ps <- phyloseq(
                          otu_table(new_otu_tab, taxa_are_rows = TRUE),
                          tax_table(ps),
                          sample_data(ps)
                        )
                        new_ps
                      }
)

# List of phyloseq objects with CLR-transformed data at different ranks
ps_list_clr <- lapply(ps_list_czm,
                      FUN = function(ps) {
                        # Impute with the CZM method
                        # Transpose bc cmultRepl assumes samples on rows
                        # and transpose back to avoid errors later
                        old_otu_tab <- otu_table(ps)
                        # CLR-transform the imputed data
                        new_otu_tab <- codaSeq.clr(data.frame(old_otu_tab, 
                                                              check.names = FALSE), samples.by.row = FALSE)
                        new_ps <- phyloseq(
                          otu_table(new_otu_tab, taxa_are_rows = TRUE),
                          tax_table(ps),
                          sample_data(ps)
                        )
                        new_ps
                      }
)


# Plots of abundances with CLR data
my_comparisons <- list(c("CAPA","COVID-19"), c("CAPA","IAPA"),
                       c("COVID-19","Influenza"),
                       c("IAPA", "Influenza")
)
# REORDER FACTOR LEVELS (for plot)
ps_list_clr$Phylum@sam_data$Disease.group <- 
  factor(ps_list_clr$Phylum@sam_data$Disease.group, c("Influenza", "IAPA", "COVID-19", "CAPA"))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
my_colors <- gg_color_hue(length(unique(ps@sam_data$Disease.group)))
names(my_colors) <- unique(ps@sam_data$Disease.group)

# Alpha-diversity measures
#Generate a data.frame with adiv measures
a_div <- data.frame(
  "Observed" = phyloseq::estimate_richness(ps, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(ps, measures = "Simpson"),
  "Disease.group" = phyloseq::sample_data(ps)$Disease.group)

# Aitchistion distance:
f.n0 <- cmultRepl(t(ps@otu_table), method = "CZM", label = 0)
f.clr <- codaSeq.clr(f.n0)
aitch <- as.matrix(aDist(f.n0))
diag(aitch) <- NA

#Add it into the metadata:
ps@sam_data[rownames(aitch), "Aitchison.distance"] <- rowMeans(aitch, na.rm=T)

aitch_d <- as.data.frame(aitch)
aitch_d[is.na(aitch_d)] <- 0

dis_group <- cmdscale(aitch_d,eig=TRUE, k=2) # k is the number of dim
ord_plot <- plot_ordination(ps_clr, dis_group, color = "Disease.group") +
  labs(col = "Disease group")+ scale_size_manual(12) +
  stat_ellipse(aes(group = Disease.group), linetype = 2) +
  scale_color_manual(values = my_colors) 

### PERMANOVA

# Generate distance matrix
# clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
# clr_dist_matrix <- as.matrix(clr_dist_matrix)
# diag(clr_dist_matrix) <- NA

get_adonis <- function(trait, covs, dist_meas, mTab)#, strataVar= FALSE) 
{
  
  
  #Definition of the different distances to be used (names in the meta data)
  b_dists <- #c("JSD","Weighted_Unifrac","Unweighted_Unifrac","Bray.Curtis","Canberra", "Jaccard",
    "Aitchison.distance"#)
  
  #Vector of the data related to the distances caclulated before 
  bdist.codes <- structure(#c("jsd.distance_samples","weighted_Unifrac_samples","unweighted_Unifrac_samples","bray_distance_samples","canberra_distance_samples", "jaccard_distance_samples",
    "aitch"#)
    , .Names = b_dists)
  
  
  #exclude trait from covs vector if already there to avoid duplicating
  if (trait %in% covs) covs <- covs[ covs != trait ]
  
  # get only those columns that are needed
  mTab <- mTab[ , c(trait, covs)]#, strataVar) ]
  # then keep only rows without NAs, or will get an error from adonis function
  mTab <- mTab[ sapply(rownames(mTab), function(x) sum(is.na(mTab[x, ]))==0), ]
  # sample names that are kept for this analysis
  sn <- rownames(mTab) #contain for each variable, the name of the samples 
  
  # get string of fixed effects
  fixed_effs <- paste(c(trait, covs), collapse = " + ")
  
  # adonis analysis
  ado <- adonis(as.formula(sprintf('as.dist(%s[sn, sn]) ~ %s', 
                                   bdist.codes[ dist_meas ], fixed_effs)), 
                data=mTab)#, strata = mTab[ , strataVar])
  return( as.data.frame(ado$aov.tab) )
}

# the names should correspond to the colnames in the metadata
adonis_vars <- c("Disease.group", "AB.treatment.14.days.before.BAL",
                 "Days.on.antibiotics.14.days.before.BAL",
                 "Age.at.ICU.admission","Gender")#,"run")

adonis_vars <- c("Disease.group", "EORTC.risk.factor", "Age.at.ICU.admission", "Days.between.start.MV.en.study.BAL")
covs <- c("Gender", "Age.at.ICU.admission")

library(microbiome)
ps.meta <- meta(ps)

adonis_list <- list()
for (trait in adonis_vars) {
  print(trait)
  adonis_list[[ trait ]] <- list()
  for (dist_meas in  #c("Bray-Curtis", "Jaccard", "Canberra", "JSD", 
       "Aitchison.distance") {
    
    print(dist_meas)
    
    adonis_list[[ trait ]][[ dist_meas ]] <- get_adonis(trait, cov,
                                                        dist_meas, ps.meta)#,
    #strataVar= "run")
  }
}

adonis_list$Disease.group
adonis_list$Age.at.ICU.admission


#Get the pvals for each trait with each dist_meas

adonis.pvals <- sapply(names(adonis_list), function(x) 
  lapply(adonis_list[[ x ]], function(y)
    y[x, "Pr(>F)"]))


# Write only those pvals < 0.05
 
adonis.signifs <- adonis.pvals
adonis.signifs[ adonis.signifs > 0.05 ] <- ""
adonis.signifs


# For each dist_meas, get p-vals of covariates, to see what covaries with each trait
 
adonis.covs.pvals <- lapply(adonis_list, function(y) {
  cov.ps <- sapply(names(y), function(x) {
    y[[ x ]][ covs, "Pr(>F)"]
  })
  rownames(cov.ps) <- covs
  cov.ps
})

# Then write just those which are significant, ignoring the values for same trait (if one of the covs)

 
adonis.covs.signifs <- lapply(names(adonis.covs.pvals), function(x) {
  covsTab <- adonis.covs.pvals[[ x ]]
  covsTab[ covsTab > 0.1 ] <- ""
  
  if (x %in% covs) {
    covsTab[x, covsTab[x, ] != "" ] <- ""
  }
  return(covsTab)
}); names(adonis.covs.signifs) <- names(adonis.covs.pvals) # lapply does not give the names to the list items


 
adonis.covs.signifs
adonis.signifs