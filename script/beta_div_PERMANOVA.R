# Load libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(zCompositions)
library(CoDaSeq)
library(robCompositions)
library(microbiome)

# Load the data
output_path <- "results/output/directory"
load(file.path(output_path, "data_after_setup.RData"))

## Obtain the Aitchison distances
f.n0 <- cmultRepl(t(ps@otu_table), method = "CZM", label = 0)
f.clr <- codaSeq.clr(f.n0)
aitch <- as.matrix(aDist(f.n0))
diag(aitch) <- NA
# add them into the metadata:
ps@sam_data[rownames(aitch), "Aitchison.distance"] <- rowMeans(aitch, na.rm=T)

# Aitchison distances dataframe
aitch_d <- as.data.frame(aitch)
aitch_d[is.na(aitch_d)] <- 0

dis_group <- cmdscale(aitch_d,eig=TRUE, k=2) # k is the number of dim

# normalize the data with centered log-ratio transformation
ps_clr <- phyloseq(otu_table(t(f.clr), taxa_are_rows = T),
                   sample_data(ps),
                   tax_table(ps))

# reorder factor levels (for plot)
ps_clr@sam_data$Disease.group <- 
  factor(ps_clr@sam_data$Disease.group, c("Influenza", "IAPA", "COVID-19", "CAPA"))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
my_colors <- gg_color_hue(length(unique(ps@sam_data$Disease.group)))
names(my_colors) <- unique(ps@sam_data$Disease.group)

ord_plot <- plot_ordination(ps_clr, dis_group, color = "Disease.group") +
  labs(col = "Disease group")+ scale_size_manual(12) +
  stat_ellipse(aes(group = Disease.group), linetype = 2) +
  scale_color_manual(values = my_colors) 
ord_plot

ggsave(file.path(output_path, "2D_plot.png"))
write.csv(aitch_d, file.path(output_path, "beta_div_matrix.csv"))

# Pairwise PERMANOVA

# Define the variables
# categorical variables
cat_vars   <- c("Disease.group", "EORTC.risk.factor", "Gender",
                "AB.treatment.14.days.before.BAL")
names(cat_vars) <- cat_vars
# continuous variables
cont_vars  <- c("Age.at.ICU.admission", "Days.between.start.MV.en.study.BAL")
names(cont_vars) <- cont_vars
# consider all the variables 
vars <- c("Disease.group", "EORTC.risk.factor", "Gender",
          "AB.treatment.14.days.before.BAL", "Age.at.ICU.admission", "Days.between.start.MV.en.study.BAL.numeric")

# Subset ps object by disease
ps_FLU <- subset_samples(ps, startsWith(sample_names(ps), "FLU"))
ps_IAPA <- subset_samples(ps, startsWith(sample_names(ps), "IAP"))
ps_COV <- subset_samples(ps, startsWith(sample_names(ps), "COV"))
ps_CAPA <- subset_samples(ps, startsWith(sample_names(ps), "CAP"))

# merge ps objects and remove taxa with count == 0
ps_FLU_IAPA <- merge_phyloseq(ps_FLU, ps_IAPA)
ps_FLU_IAPA <- filter_taxa(ps_FLU_IAPA, function(x) sum(x > 0) > 0, prune = TRUE)
ps_COV_CAPA <- merge_phyloseq(ps_COV, ps_CAPA)
ps_COV_CAPA <- filter_taxa(ps_COV_CAPA, function(x) sum(x > 0) > 0, prune = TRUE)
ps_FLU_COV <- merge_phyloseq(ps_FLU, ps_COV)
ps_FLU_COV <- filter_taxa(ps_FLU_COV, function(x) sum(x > 0) > 0, prune = TRUE)
ps_IAPA_CAPA <- merge_phyloseq(ps_IAPA, ps_CAPA)
ps_IAPA_CAPA <- filter_taxa(ps_IAPA_CAPA, function(x) sum(x > 0) > 0, prune = TRUE)

# change the variable "sample_data(ps_filt)$Days.between.start.MV.en.study.BAL" from character into numeric
sample_data(ps_FLU_IAPA)$Days.between.start.MV.en.study.BAL.numeric <- as.numeric(sample_data(ps_FLU_IAPA)$Days.between.start.MV.en.study.BAL)
sample_data(ps_COV_CAPA)$Days.between.start.MV.en.study.BAL.numeric <- as.numeric(sample_data(ps_COV_CAPA)$Days.between.start.MV.en.study.BAL)
sample_data(ps_FLU_COV)$Days.between.start.MV.en.study.BAL.numeric <- as.numeric(sample_data(ps_FLU_COV)$Days.between.start.MV.en.study.BAL)
sample_data(ps_IAPA_CAPA)$Days.between.start.MV.en.study.BAL.numeric <- as.numeric(sample_data(ps_IAPA_CAPA)$Days.between.start.MV.en.study.BAL)

# List of ps objects agglomerated at different ranks 
# Taxonomic ranks
all_ranks <- c("Genus", "Species")

pairwise_permanova <- function(ps, all_ranks, vars) {
  # Agglomerate at different taxonomic ranks
  ps_list_glom <- sapply(all_ranks, function(rank) {
    new_ps <- tax_glom(ps, rank)
    # Rename taxa with the according rank
    taxa_names(new_ps) <- unlist(as.vector(tax_table(new_ps)[, rank]))
    return(new_ps)
  }, simplify = FALSE)
  
  # Imputation with CZM method and calculate Aitchison distances
  ps_list_czm <- lapply(ps_list_glom, function(ps_item) {
    # Impute with the CZM method
    # Transpose bc cmultRepl assumes samples on rows
    # and transpose back to avoid errors later
    otu_tab <- otu_table(ps_item)
    try_czm <- try(t(cmultRepl(t(otu_tab), method = "CZM", label = 0)), silent = TRUE)
    if (inherits(try_czm, "try-error")) {
      new_otu_tab <- otu_tab
    } else {
      new_otu_tab <- try_czm
    }
    new_ps <- phyloseq(otu_table(new_otu_tab, taxa_are_rows = TRUE), 
                       tax_table(ps_item), 
                       sample_data(ps_item))
    return(new_ps)
  })
  
  aitch_list <- sapply(ps_list_czm, function(ps_czm) {
    # Calculate Aitchison distance (DOES CLR INTERNALLY)
    # and assumes samples are rows
    as.matrix(aDist(t(otu_table(ps_czm))))
  }, simplify = FALSE)
  
  # Prepare metadata and perform PERMANOVA tests
  meta_data <- meta(ps)[, vars]
  meta_data <- meta_data[sapply(rownames(meta_data), 
                     function(x) sum(is.na(meta_data[x, ]))==0),]
  sample_names <- rownames(meta_data)
  fixed_effects <- paste(vars, collapse = " + ")
  
  perm_results <- list()
  for (rank in names(aitch_list)) {
    formula_str <- sprintf("as.dist(aitch_list$%s[sample_names, sample_names]) ~ %s", rank, fixed_effects)
    perm_results[[rank]] <- adonis2(as.formula(formula_str), data = meta_data)
  }
  # Adjust p-values for multiple comparisons
  for (rank in names(aitch_list)) {
    p_adj_results <- lapply(perm_results, function(x) {
      p.adjust(x$`Pr(>F)`, method = "fdr")
    })
  }
    
  # Return a list of results
  return(list(glommed = ps_list_glom, 
              czm = ps_list_czm, 
              aitchison = aitch_list, 
              permanova = perm_results, 
              p_adjusted = p_adj_results))
}

permanova_FLU_IAPA <- pairwise_permanova(ps_FLU_IAPA, all_ranks, vars)
permanova_COV_CAPA <- pairwise_permanova(ps_COV_CAPA, all_ranks, vars)
permanova_FLU_COV <- pairwise_permanova(ps_FLU_COV, all_ranks, vars)
permanova_IAPA_CAPA <- pairwise_permanova(ps_IAPA_CAPA, all_ranks, vars)