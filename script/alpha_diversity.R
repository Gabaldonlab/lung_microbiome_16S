# Load the libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the data
output_path <- "results/output/directory"
load(file.path(output_path, "data_after_setup.RData"))

# Generate a data frame with alpha diversity measures
a_div <- data.frame(
  "Observed" = phyloseq::estimate_richness(ps, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(ps, measures = "Simpson")
)
# Save the data to csv file
write.csv(a_div, file.path(output_path, "a_div.csv"))

# Define comparisons
my_comparisons <- list(c("CAPA","COVID-19"), 
                       c("CAPA","IAPA"),
                       c("COVID-19","Influenza"),
                       c("IAPA", "Influenza")
)

# Function to perform pairwise Wilcoxon tests
run_wilcox_tests <- function(data, variable) {
  lapply(my_comparisons, function(pair) {
    test_results <- wilcox.test(reformulate("Disease.group", response = variable), 
                                data = data, subset = Disease.group %in% pair, exact = FALSE)
    test_results$p.value <- format.pval(test_results$p.value, digits = 3)
    test_results
  })
}

# Apply function for Observed, Shannon, and Simpson and extract the p-values
p_values_observed <- sapply(run_wilcox_tests(a_div, "Observed"), `[[`, "p.value")
p_values_shannon  <- sapply(run_wilcox_tests(a_div, "Shannon"), `[[`, "p.value")
p_values_simpson  <- sapply(run_wilcox_tests(a_div, "Simpson"), `[[`, "p.value")

# Print the p-values
print(p_values_observed)
print(p_values_shannon)
print(p_values_simpson)

## Melt data for ggplot and plot
a_div_long <- a_div %>%
  pivot_longer(cols = c("Observed", "Shannon", "Simpson"), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson")))

ggplot(a_div_long, aes(x = Disease.group, y = value, color = Disease.group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  facet_wrap(~ metric, scales = "free") +
  labs(x = "Disease group", y = "") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", label = "p.signif", hide.ns = TRUE, exact = FALSE) +
  theme_bw(base_size = 9) +
  labs(colour = "Disease group")

ggsave(file.path(output_path, "alpha-div.png"), height = 7, width = 10)

save.image(file.path(output_path, "temp_data_till_adiv.RData"))