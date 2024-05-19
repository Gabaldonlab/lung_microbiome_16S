# Load the libraries
library(phyloseq)
library(xlsx)
library(zCompositions)
library(CoDaSeq)
library(robCompositions)
library(microbiome)

# Load the data
output_path <- "results/output/directory"

load(file.path(output_path, "data_after_betadiversity.RData"))

# Filter out the taxa present in less than 10% of the samples (~16 samples)
ps_strict_filter <- filter_taxa(ps, function(x) sum(x > 0) > 15, prune = TRUE) # from 238 taxa to 64
ps_strict_filter

# Subset ps object by disease
ps_FLU_strictfilter <- subset_samples(ps_strict_filter, startsWith(sample_names(ps_strict_filter), "FLU"))
ps_IAPA_strictfilter <- subset_samples(ps_strict_filter, startsWith(sample_names(ps_strict_filter), "IAP"))
ps_COV_strictfilter <- subset_samples(ps_strict_filter, startsWith(sample_names(ps_strict_filter), "COV"))
ps_CAPA_strictfilter <- subset_samples(ps_strict_filter, startsWith(sample_names(ps_strict_filter), "CAP"))

# Merge phyloseq objects
ps_FLU_IAPA_strictfilter <- merge_phyloseq(ps_FLU_strictfilter, ps_IAPA_strictfilter)
ps_COV_CAPA_strictfilter <- merge_phyloseq(ps_COV_strictfilter, ps_CAPA_strictfilter)
ps_FLU_COV_strictfilter <- merge_phyloseq(ps_FLU_strictfilter, ps_COV_strictfilter)
ps_IAPA_CAPA_strictfilter <- merge_phyloseq(ps_IAPA_strictfilter, ps_CAPA_strictfilter)

# Function to calculate Aitchison distance and CLR transformation for new filtered objects
aitch_clr <- function(ps) {
  f.n0 <- cmultRepl(t(otu_table(ps)), method = "CZM", label = 0)
  f.clr <- codaSeq.clr(f.n0)
  aitch <- as.matrix(aDist(f.n0))
  diag(aitch) <- NA
  sample_data(ps)[rownames(aitch), "Aitchison.distance"] <- rowMeans(aitch, na.rm = TRUE)
  phyloseq(otu_table(t(f.clr), taxa_are_rows = TRUE),
           sample_data(ps),
           tax_table(ps))
}

#### MANN-WHITNEY TEST

perform_differential_abundance <- function(ps, disease_groups) {
  # Define the ranks to aggregate at
  all_ranks <- c("Genus", "Species")
  
  # Taxonomic aggregation
  ps_list_glom <- sapply(all_ranks, function(rank) {
    new_ps <- tax_glom(ps, rank)
    taxa_names(new_ps) <- unlist(as.vector(tax_table(new_ps)[, rank]))
    new_ps
  }, simplify = FALSE)
  
  # Perform Wilcoxon test across all ranks and taxa
  wlcx_test_results <- sapply(all_ranks, function(rank) {
    sapply(taxa_names(ps_list_glom[[rank]]), function(taxon) {
      wilcox.test(
        otu_table(ps_list_glom[[rank]])[taxon, meta(ps)$Disease.group == disease_groups[1]]@.Data,
        otu_table(ps_list_glom[[rank]])[taxon, meta(ps)$Disease.group == disease_groups[2]]@.Data
      )
    }, simplify = FALSE)
  }, simplify = FALSE)
  
  # Extract and adjust p-values
  wlcx_pvals <- sapply(wlcx_test_results, function(rank_list) {
    sapply(rank_list, function(test) test$p.value)
  })
  
  list(test_results = wlcx_test_results, p_values = wlcx_pvals)
}

# Function to save Wilcoxon results
save_wilcoxon_results <- function(p_values, output_path, file_prefix) {
  lapply(names(p_values), function(rank) {
    df <- as.data.frame(p_values[[rank]][p_values[[rank]] < 0.05], stringsAsFactors = FALSE)
    if (nrow(df) > 0) {
      write.xlsx(df, file.path(output_path, paste0("MW_", file_prefix, "_", rank, ".xlsx")))
    }
  })
}

# Main analysis function
perform_analysis <- function(ps, disease_groups, file_prefix) {
  # Apply Aitchison CLR transformation
  ps_clr <- aitch_clr(ps)
  
  # Wilcoxon analysis
  results <- perform_differential_abundance(ps, disease_groups)
  
  # Save Wilcoxon results
  save_wilcoxon_results(results$p_values, output_path, file_prefix)
  
  # Adjust p-values for multiple testing
  p_adj <- lapply(results$p_values, function(p) p.adjust(p, method = "fdr"))
  sum(unlist(p_adj) < 0.05)
  
  # Return a list of results
  return(list(p_values = results$p_values, 
              p_adj = p_adj))
}

# Run the main analysis function
wilcoxon_FLU_IAPA <- perform_analysis(ps_FLU_IAPA_strictfilter, c("Influenza", "IAPA"), "FLU_IAPA")
wilcoxon_COV_CAPA <- perform_analysis(ps_COV_CAPA_strictfilter, c("COVID-19", "CAPA"), "COV_CAPA")
wilcoxon_FLU_COV <- perform_analysis(ps_FLU_COV_strictfilter, c("Influenza", "COVID-19"), "FLU_COV")
wilcoxon_IAPA_CAPA <- perform_analysis(ps_IAPA_CAPA_strictfilter, c("IAPA", "CAPA"), "IAPA_CAPA")


#### LINEAR MODELS

# Aggregate at different taxonomic levels the phyloseq object
get_gloms <- function(phy, tl, tTabs) {
  # agglomerate taxa at given level
  tg <- tax_glom(phy, taxrank = tl)@otu_table
  return(tg)
}


agglomerate_fun <- function(ps_obj){
  gloms_clr <- list()
  otu_clr.types <- list()
  for (tl in c("Genus", "Species")) {
    print(tl)
    gloms_clr [["aspergillosis"]] [[tl]] <- get_gloms(ps_obj, tl, taxTables)
  }
  return(gloms_clr)
}

ps_FLU_IAPA_sf_clr <- aitch_clr(ps_FLU_IAPA_strictfilter)
ps_COV_CAPA_sf_clr <- aitch_clr(ps_COV_CAPA_strictfilter)
ps_FLU_COV_sf_clr <- aitch_clr(ps_FLU_COV_strictfilter)
ps_IAPA_CAPA_sf_clr <- aitch_clr(ps_IAPA_CAPA_strictfilter)

gloms_FLU_IAPA_clr <- agglomerate_fun(ps_FLU_IAPA_sf_clr)
gloms_COV_CAPA_clr <- agglomerate_fun(ps_COV_CAPA_sf_clr)
gloms_FLU_COV_clr <- agglomerate_fun(ps_FLU_COV_sf_clr)
gloms_IAPA_CAPA_clr <- agglomerate_fun(ps_IAPA_CAPA_sf_clr)

# Apply multiple string replacements on row names of specific phyloseq object levels
apply_name_corrections <- function(phylo_obj, levels, corrections) {
  for (level in levels) {
    current_names <- rownames(phylo_obj[[level]])
    for (correction in corrections) {
      current_names <- gsub(correction$pattern, correction$replacement, current_names, fixed = TRUE)
    }
    rownames(phylo_obj[[level]]) <- current_names
  }
}

# Define the levels and corrections to apply
levels_to_correct <- c("Genus", "Species")
corrections <- list(
  list(pattern = " ", replacement = "."),
  list(pattern = "[", replacement = ""),
  list(pattern = "]", replacement = ""),
  list(pattern = "-", replacement = "."),
  list(pattern = "(", replacement = "."),
  list(pattern = ")", replacement = ".")
)

# Apply transformations to FLU-IAPA and COV-CAPA phyloseq objects
apply_name_corrections(gloms_FLU_IAPA_clr[["aspergillosis"]], levels_to_correct, corrections)
apply_name_corrections(gloms_COV_CAPA_clr[["aspergillosis"]], levels_to_correct, corrections)
apply_name_corrections(gloms_FLU_COV_clr[["aspergillosis"]], levels_to_correct, corrections)
apply_name_corrections(gloms_IAPA_CAPA_clr[["aspergillosis"]], levels_to_correct, corrections)

# Function to perform linear modeling and ANOVA analysis
lin_mod_fun <- function(ps_obj, fixed_effects, glom_data){
  # Combine fixed effects into a formula part
  fs_full <- paste(fixed_effects, collapse = " + ")

  # Initialize lists to store p-values and adjusted p-values
  Anova_pvals <- list()
  Anova_padj <- list()

  # Metadata table
  mtab <- meta(ps_obj)

  # Loop through taxonomy levels
  for (tl in c("Genus", "Species")) {
    # Extract and prepare the taxonomy table for analysis
    glomTab <- glom_data[["aspergillosis"]][[tl]]
    glomTab <- glomTab[, colnames(glomTab) %in% rownames(mtab)]
    rownames(glomTab) <- gsub("[()-]", ".", rownames(glomTab))  # Sanitize row names

    # Bind the glom data to metadata
    mtab <- cbind(mtab, t(glomTab))

    # Initialize matrix to store p-values
    pval_mat <- matrix(NA, nrow = nrow(glomTab), ncol = length(fixed_effects))
    rownames(pval_mat) <- rownames(glomTab)
    colnames(pval_mat) <- fixed_effects

    # Loop through each dependent variable
    for (dependVar in rownames(pval_mat)) {
      formula <- as.formula(sprintf("%s ~ %s", dependVar, fs_full))
      mod_lm <- glm(formula, data = mtab)

      # Perform ANOVA and handle potential errors
      am <- try(as.data.frame(Anova(mod_lm)))
      if (class(am) == "try-error") {
        next
      } else {
        pval_mat[dependVar, fixed_effects] <- am[fixed_effects, "Pr(>Chisq)"]
      }
    }

    # Filter out rows with only NA values
    valid_rows <- rowSums(is.na(pval_mat)) < ncol(pval_mat)
    Anova_pvals[[tl]] <- pval_mat[valid_rows, ]

    # Adjust p-values for multiple testing
    Anova_padj[[tl]] <- apply(Anova_pvals[[tl]], 2, function(p) p.adjust(p, method = "fdr"))
  }

  # Return a list containing ANOVA results and adjusted p-values
  list(Anova_results = Anova_pvals, Anova_padj = Anova_padj)
}

# Example usage of the function
fs <- c("Disease.group", "EORTC.risk.factor", "Age.at.ICU.admission", "Days.between.start.MV.en.study.BAL")

results_FLU_IAPA <- lin_mod_fun(ps_FLU_IAPA_sf_clr, fs, gloms_FLU_IAPA_clr)
results_COV_CAPA <- lin_mod_fun(ps_COV_CAPA_sf_clr, fs, gloms_COV_CAPA_clr)
results_FLU_COV <- lin_mod_fun(ps_FLU_COV_sf_clr, fs, gloms_FLU_COV_clr)
results_IAPA_CAPA <- lin_mod_fun(ps_IAPA_CAPA_sf_clr, fs, gloms_IAPA_CAPA_clr)


#### SELBAL

# Data preparation
# modify into numeric the metadata EORTC
for (el in 1:length(ps_strict_filter@sam_data$EORTC.risk.factor)){
  ps_strict_filter@sam_data$EORTC.risk.factor[el] = as.numeric(ps_strict_filter@sam_data$EORTC.risk.factor[el])
}
# modify into numeric the metadata Age
for (el in 1:length(ps_strict_filter@sam_data$Age.at.ICU.admission)){
  ps_strict_filter@sam_data$Age.at.ICU.admission[el] = as.numeric(ps_strict_filter@sam_data$Age.at.ICU.admission[el])
}
# modify into numeric the metadata MV-BAL
for (el in 1:length(ps_strict_filter@sam_data$Days.between.start.MV.en.study.BAL)){
  ps_strict_filter@sam_data$MV_BAL.numeric[el] = as.numeric(ps_strict_filter@sam_data$Days.between.start.MV.en.study.BAL[el])
}
# change the NAs into another number, otherwise selbal gives error
ps_strict_filter@sam_data$MV_BAL.numeric[is.na(ps_strict_filter@sam_data$MV_BAL.numeric)] = 40

ps_strict_filter@sam_data$Patient.ID <- gsub("VARIOMIC_", "", ps_strict_filter@sam_data$Patient.ID)
ps_strict_filter@sam_data$Patient.ID <- gsub("Variomic_", "", ps_strict_filter@sam_data$Patient.ID)
ps_strict_filter@sam_data$Patient.ID <- gsub("FLU_", "FLU", ps_strict_filter@sam_data$Patient.ID)
ps_strict_filter@sam_data$Patient.ID <- gsub("IAPA_", "IAP", ps_strict_filter@sam_data$Patient.ID)
ps_strict_filter@sam_data$Patient.ID <- gsub("COVID_", "COV", ps_strict_filter@sam_data$Patient.ID)
ps_strict_filter@sam_data$Patient.ID <- gsub("CAPA_", "CAP", ps_strict_filter@sam_data$Patient.ID)

# subset the phyloseq objects
# FLU
ps_FLU <- subset_samples(ps_strict_filter, startsWith(sample_names(ps_strict_filter), "FLU"))
# IAPA
ps_IAPA <- subset_samples(ps_strict_filter, startsWith(sample_names(ps_strict_filter), "IAP"))
# COV
ps_COV <- subset_samples(ps_strict_filter, startsWith(sample_names(ps_strict_filter), "COV"))
# CAPA
ps_CAPA <- subset_samples(ps_strict_filter, startsWith(sample_names(ps_strict_filter), "CAP"))

# merge phyloseq objects
ps_FLU_IAPA <- merge_phyloseq(ps_FLU, ps_IAPA)
ps_COV_CAPA <- merge_phyloseq(ps_COV, ps_CAPA)
ps_FLU_COV <- merge_phyloseq(ps_FLU, ps_COV)
ps_IAPA_CAPA <- merge_phyloseq(ps_IAPA, ps_CAPA)

# filter taxa appearing in less than 2 samples
taxa_filtering <- function(ps_object){
  filter_taxa(ps_object, function(x) sum(x > 0) > 1, prune=TRUE)
}

ps_FLU_IAPA_selbal <- taxa_filtering(ps_FLU_IAPA)
ps_FLU_COV_selbal <- taxa_filtering(ps_FLU_COV)
ps_COV_CAPA_selbal <- taxa_filtering(ps_COV_CAPA)
ps_IAPA_CAPA_selbal <- taxa_filtering(ps_IAPA_CAPA)

# Transform the disease.group condition into numeric (as character -> for selbal if it is dichotoumous it needs to be a factor)
for (el in 1:length(ps_FLU_IAPA_selbal@sam_data$Disease.group)){
  if (ps_FLU_IAPA_selbal@sam_data$Disease.group[el] == "Influenza"){
    ps_FLU_IAPA_selbal@sam_data$Disease.group.numeric[el] = "0"
  }
  else{
    ps_FLU_IAPA_selbal@sam_data$Disease.group.numeric[el] = "1"
  }
}

for (el in 1:length(ps_FLU_COV_selbal@sam_data$Disease.group)){
  if (ps_FLU_COV_selbal@sam_data$Disease.group[el] == "Influenza"){
    ps_FLU_COV_selbal@sam_data$Disease.group.numeric[el] = "0"
  }
  else{
    ps_FLU_COV_selbal@sam_data$Disease.group.numeric[el] = "2"
  }
}

for (el in 1:length(ps_COV_CAPA_selbal@sam_data$Disease.group)){
  if (ps_COV_CAPA_selbal@sam_data$Disease.group[el] == "COVID-19"){
    ps_COV_CAPA_selbal@sam_data$Disease.group.numeric[el] = "2"
  }
  else{
    ps_COV_CAPA_selbal@sam_data$Disease.group.numeric[el] = "3"
  }
}

for (el in 1:length(ps_IAPA_CAPA_selbal@sam_data$Disease.group)){
  if (ps_IAPA_CAPA_selbal@sam_data$Disease.group[el] == "IAPA"){
    ps_IAPA_CAPA_selbal@sam_data$Disease.group.numeric[el] = "1"
  }
  else{
    ps_IAPA_CAPA_selbal@sam_data$Disease.group.numeric[el] = "3"
  }
}

# Function to create data frame with covariates
covars_df <- function(ps){
  df <- data.frame(
    sample_ID <- ps@sam_data$Patient.ID,
    Age <- ps@sam_data$Age.at.ICU.admission,
    EORTC <- ps@sam_data$EORTC.risk.factor,
    MV_BAL <- ps@sam_data$MV_BAL.numeric
  )
  # make sure the rownames match with asv_df
  rownames(df) <- df$sample_ID
  df <- df[,-1]
  return(df)
}

cv_df_FLU_IAPA <- covars_df(ps_FLU_IAPA_strictfilter)
cv_df_COV_CAPA <- covars_df(ps_COV_CAPA_strictfilter)
cv_df_FLU_COV <- covars_df(ps_FLU_COV_strictfilter)
cv_df_IAPA_CAPA <- covars_df(ps_IAPA_CAPA_strictfilter)


# Function to run selbal
run_selbal_with_covar <- function(ps_object, covars_df){
  asv_df <- as.matrix(t(ps_object@otu_table)) # species
  disease_group <- as.vector(ps_object@sam_data$Disease.group.numeric)
  x <- asv_df
  y <- factor(disease_group)
  z <- covars_df
  #rownames(z) <- rownames(x)
  CV.bal <- selbal.cv(x = x, y = y, 
                      n.iter=10, covar = z)
}

