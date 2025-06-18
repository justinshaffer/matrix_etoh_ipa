####################################################################################################
# 2025
# Matrix project
# Justin Shaffer
# justinparkshaffer@gmail.com
####################################################################################################

####################################################################################################
# Set working directory
####################################################################################################
getwd()
setwd("/matrix/r")

install.packages("tidyverse")
remotes::install_github("jbisanz/qiime2R")

library(tidyverse)
library(qiime2R)


####################################################################################################
# Merge alpha-diversity metrics with sample metadata
####################################################################################################

# Read in metadata
md <- read_tsv("metadata_samples_qiita_20250326a.txt", na = c("", "NA", "not applicable"))


# Read in alpha-diversity metrics for 16S data
alpha_16S_lbm_shannon <- read_qza("matrix_16s_deblur_gg2_biom_silva_noMit_noChl_noUnassigned_noEuk_noDomain_noControls_noSpike_lbm_noSingletons_rar277_alpha_shannon.qza")
alpha_16S_lbm_richness <- read_qza("matrix_16s_deblur_gg2_biom_silva_noMit_noChl_noUnassigned_noEuk_noDomain_noControls_noSpike_lbm_noSingletons_rar277_alpha_richness.qza")
alpha_16S_hbm_shannon <- read_qza("matrix_16s_deblur_gg2_biom_silva_noMit_noChl_noUnassigned_noEuk_noDomain_noControls_noSpike_hbm_noSingletons_rar20636_alpha_shannon.qza")
alpha_16S_hbm_richness <- read_qza("matrix_16s_deblur_gg2_biom_silva_noMit_noChl_noUnassigned_noEuk_noDomain_noControls_noSpike_hbm_noSingletons_rar20636_alpha_richness.qza")


# Concatenate alpha-diversity metrics from low- and high-biomass samples for 16S data
alpha_16S_shannon <- rbind((rownames_to_column(as.data.frame(alpha_16S_hbm_shannon$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_16S_lbm_shannon$data), "sample_name")))
alpha_16S_richness <- rbind((rownames_to_column(as.data.frame(alpha_16S_hbm_richness$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_16S_lbm_richness$data), "sample_name")))


# Merge alpha-diversity metrics with sample metadata
alpha_16S_shannon %>%
  right_join(md) -> md
colnames(md)[2] <- "alpha_16S_shannon"

alpha_16S_richness %>%
  right_join(md) -> md
colnames(md)[2] <- "alpha_16S_richness"

# Export new metadata file
write_tsv(md, file = "metadata_samples_qiita_20250319h.txt")


# Read in metadata for shotgun data
md <- read_tsv("metadata_samples_qiita_20250319d.txt", na = c("", "NA", "not applicable"))


# Read in alpha-diversity metrics for shotgun data
alpha_shotgun_lbm_shannon <- read_qza("matrix_shotgun_wolr2pe_biom_lbm_noControls_noSpike_noSingletons_rar55K_alpha_shannon.qza")
alpha_shotgun_lbm_richness <- read_qza("matrix_shotgun_wolr2pe_biom_lbm_noControls_noSpike_noSingletons_rar55K_alpha_richness.qza")
alpha_shotgun_hbm_shannon <- read_qza("matrix_shotgun_wolr2pe_biom_hbm_noControls_noSpike_noSingletons_rar1515K_alpha_shannon.qza")
alpha_shotgun_hbm_richness <- read_qza("matrix_shotgun_wolr2pe_biom_hbm_noControls_noSpike_noSingletons_rar1515K_alpha_richness.qza")


# Concatenate alpha-diversity metrics from low- and high-biomass samples for shotgun data
alpha_shotgun_shannon <- rbind((rownames_to_column(as.data.frame(alpha_shotgun_hbm_shannon$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_shotgun_lbm_shannon$data), "sample_name")))
alpha_shotgun_richness <- rbind((rownames_to_column(as.data.frame(alpha_shotgun_hbm_richness$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_shotgun_lbm_richness$data), "sample_name")))


# Merge alpha-diversity metrics with sample metadata for shotgun data
alpha_shotgun_shannon %>%
  right_join(md) -> md
colnames(md)[2] <- "alpha_shotgun_shannon"

alpha_shotgun_richness %>%
  right_join(md) -> md
colnames(md)[2] <- "alpha_shotgun_richness"

# Export new metadata file
write_tsv(md, file = "metadata_samples_qiita_20250319d.txt")


# Read in metadata for metabolomics data
md <- read_tsv("metadata_samples_qiita_20250305e.txt", na = c("", "NA", "not applicable"))


# Read in alpha-diversity metrics for metab data
alpha_metab_lbm_shannon <- read_qza("matrix_lcms_lbm_biom_qiita_ids_alpha_shannon.qza")
alpha_metab_lbm_richness <- read_qza("matrix_lcms_lbm_biom_qiita_ids_alpha_richness.qza")
alpha_metab_hbm_shannon <- read_qza("matrix_lcms_hbm_biom_qiita_ids_alpha_shannon.qza")
alpha_metab_hbm_richness <- read_qza("matrix_lcms_hbm_biom_qiita_ids_alpha_richness.qza")
alpha_metab_noSingletons_shannon <- read_qza("matrix_lcms_merged_biom_qiita_ids_noSingletons_alpha_shannon.qza")
alpha_metab_noSingletons_richness <- read_qza("matrix_lcms_merged_biom_qiita_ids_noSingletons_alpha_richness.qza")


# Concatenate alpha-diversity metrics from low- and high-biomass samples for metab data
alpha_metab_shannon <- rbind((rownames_to_column(as.data.frame(alpha_metab_hbm_shannon$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_metab_lbm_shannon$data), "sample_name")))
alpha_metab_richness <- rbind((rownames_to_column(as.data.frame(alpha_metab_hbm_richness$data), "sample_name")), (rownames_to_column(as.data.frame(alpha_metab_lbm_richness$data), "sample_name")))


# Merge alpha-diversity metrics with sample metadata for metab data
alpha_metab_shannon %>%
  right_join(md) -> md
colnames(md)[2] <- "alpha_lcmsms_shannon"

alpha_metab_richness %>%
  right_join(md) -> md
colnames(md)[2] <- "alpha_lcmsms_richness"

alpha_metab_noSingletons_shannon_data <- rownames_to_column(alpha_metab_noSingletons_shannon$data, var = "sample_name")
alpha_metab_noSingletons_shannon_data %>%
  right_join(md) -> md
colnames(md)[2] <- "alpha_lcmsms_shannon"

alpha_metab_noSingletons_richness_data <- rownames_to_column(alpha_metab_noSingletons_richness$data, var = "sample_name")
alpha_metab_noSingletons_richness_data %>%
  right_join(md) -> md
colnames(md)[2] <- "alpha_lcmsms_richness"

# Export new metadata file
write_tsv(md, file = "metadata_samples_qiita_20250326b.txt")


#######################################################################################################
# Tidy up metadata for plotting
#######################################################################################################

# Import metadata
md <- read_tsv("metadata_samples_qiita_20250326b.txt")


# Subset metadata to exclude pilot and control samples
md_subset <- subset(md, md$extraction_protocol == "Matrix" & md$empo_1 != "Control" & md$round == "1_2" | md$extraction_protocol == "MagMax" & md$empo_1 != "Control" & md$round == "1_2")
md_subset <- droplevels(md_subset)


# Re-order levels for sample_type2
md_subset$sample_type2 <- factor(md_subset$sample_type2, levels = c("human_feces", 
                                                                    "mouse_feces", 
                                                                    "human_saliva_before_brushing", 
                                                                    "human_saliva_after_brushing", 
                                                                    "skin_armpit", "skin_hand",  
                                                                    "surface_swab_floor", 
                                                                    "surface_swab_keyboard"))


# Re-name levels for sample_type2
levels(md_subset$sample_type2) <- c("Human feces", 
                                    "Mouse feces", 
                                    "Human saliva, before", 
                                    "Human saliva, after", 
                                    "Skin, armpit", 
                                    "Skin, hand",  
                                    "Surface, floor", 
                                    "Surface, keyboard")


# Make vector for colors for sample type 2
sample_type_2_colors <- c("#F48022", "#A6CFE5", "#8AC752", "#FCBF6D", "#22201F", "#F9F49D", "#B15829", "#6A3F98")


# Re-name levels for biomass_plate
md_subset$biomass <- as.factor(md_subset$biomass)
levels(md_subset$biomass) <- c("High biomass", "Low biomass")


# Create factor variable for storage_solution
md_subset$storage_solution <- as.factor(md_subset$storage_solution)


# Make vector for colors for storage_solution
storage_solution_colors <- c("black", "white")


# Re-order levels for host_subject_id
md_subset$host_subject_id <- factor(md_subset$host_subject_id, levels = c("A", "B", "C", "D", 
                                                                          "mouse1", "mouse2", "mouse3", "mouse4",
                                                                          "bayd", "room1247",
                                                                          "pc_bayd", "pc_bayg"))


# Re-name levels for host_subject_id
levels(md_subset$host_subject_id) <- c("Human subject A", "Human subject B", "Human subject C", "Human subject D",
                                       "Mouse subject 1", "Mouse subject 2", "Mouse subject 3", "Mouse subject 4",
                                       "Floor, Bay D", "Floor, Room 1247",
                                       "Keyboard, Bay D", "Keyboard, Bay G")


# Make vector for colors for host_subject_id
host_subject_colors <- c("orange", "#A6CFE5", "blue", "red",
                         "orange", "#A6CFE5", "blue", "red",
                         "orange", "#A6CFE5",
                         "orange", "#A6CFE5")


#######################################################################################################
# Plot DNA concentrations as boxplots
#######################################################################################################

# Boxplot - four groups - protocol x storage solution (PAPER FIGURE 1)
#######################################################################################################
# Create new variable
md_subset$storage_solution_extraction_protocol <- paste(md_subset$storage_solution, md_subset$extraction_protocol, sep = "-")


# Re-order levels for new variable
md_subset$storage_solution_extraction_protocol <- factor(md_subset$storage_solution_extraction_protocol, levels = c("etoh-MagMax",
                                                                                                                    "isopropanol-MagMax",
                                                                                                                    "etoh-Matrix",
                                                                                                                    "isopropanol-Matrix"))

# Re-name levels for new variable
levels(md_subset$storage_solution_extraction_protocol) <- c("EtOH-Plate",
                                                            "IPA-Plate",
                                                            "EtOH-Matrix",
                                                            "IPA-Matrix")


# Obtain sample sizes per group
summary(md_subset$sample_type2)


# Plot
ggplot(md_subset, aes(x = storage_solution_extraction_protocol, y = dna_conc)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.5) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 2) +
  ylab("DNA concentration (ng/µL)") +
  xlab("storage solution-extraction protocol") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)


# Kruskal-Wallis tests
## Subset data
md_st_feces_human <- subset(md_subset, md_subset$sample_type2 == "Human feces")
md_st_feces_mouse <- subset(md_subset, md_subset$sample_type2 == "Mouse feces")
md_st_saliva_before <- subset(md_subset, md_subset$sample_type2 == "Human saliva, before")
md_st_saliva_after <- subset(md_subset, md_subset$sample_type2 == "Human saliva, after")
md_st_skin_armpit <- subset(md_subset, md_subset$sample_type2 == "Skin, armpit")
md_st_skin_hand <- subset(md_subset, md_subset$sample_type2 == "Skin, hand")
md_st_surface_floor <- subset(md_subset, md_subset$sample_type2 == "Surface, floor")
md_st_surface_keyboard <- subset(md_subset, md_subset$sample_type2 == "Surface, keyboard")

## Perform tests
kruskal.test(md_st_feces_human$dna_conc ~ md_st_feces_human$extraction_protocol)
kruskal.test(md_st_feces_mouse$dna_conc ~ md_st_feces_mouse$storage_solution_extraction_protocol)
kruskal.test(md_st_saliva_before$dna_conc ~ md_st_saliva_before$storage_solution_extraction_protocol)
kruskal.test(md_st_saliva_after$dna_conc ~ md_st_saliva_after$storage_solution_extraction_protocol)
kruskal.test(md_st_skin_armpit$dna_conc ~ md_st_skin_armpit$storage_solution_extraction_protocol)
kruskal.test(md_st_skin_hand$dna_conc ~ md_st_skin_hand$storage_solution_extraction_protocol)
kruskal.test(md_st_surface_floor$dna_conc ~ md_st_surface_floor$storage_solution_extraction_protocol)
kruskal.test(md_st_surface_keyboard$dna_conc ~ md_st_surface_keyboard$storage_solution_extraction_protocol)

# Subset skin_armpit dataset to perform Wilcoxon tests between storage solutions for each extraction protocol
md_st_skin_armpit_plate <- subset(md_st_skin_armpit, md_st_skin_armpit$extraction_protocol == "MagMax")
md_st_skin_armpit_matrix <- subset(md_st_skin_armpit, md_st_skin_armpit$extraction_protocol == "Matrix")

# Perform Wilcoxon tests
wilcox.test(md_st_skin_armpit_plate$dna_conc ~ md_st_skin_armpit_plate$storage_solution)
wilcox.test(md_st_skin_armpit_matrix$dna_conc ~ md_st_skin_armpit_matrix$storage_solution)


#######################################################################################################
# Plot DNA yield, read counts, and alpha-diversity as scatter plots, comparing extraction protocols
#######################################################################################################
# Create new variable
md_subset$storage_solution_extraction_protocol <- paste(md_subset$storage_solution, md_subset$extraction_protocol, sep = "-")


# Re-order levels for new variable
md_subset$storage_solution_extraction_protocol <- factor(md_subset$storage_solution_extraction_protocol, levels = c("etoh-MagMax",
                                                                                                                    "isopropanol-MagMax",
                                                                                                                    "etoh-Matrix",
                                                                                                                    "isopropanol-Matrix"))

# Re-name levels for new variable
levels(md_subset$storage_solution_extraction_protocol) <- c("EtOH-Plate",
                                                            "IPA-Plate",
                                                            "EtOH-Matrix",
                                                            "IPA-Matrix")


# Obtain sample sizes per group
summary(md_subset$sample_type2)


install.packages('ggpubr')
install.packages('cowplot')
install.packages('svglite')
install.packages('scales')

library(tidyverse)
library(ggpubr)
library(cowplot)
library(svglite)
library(scales)


# Tidy-up metadata for running correlations and plotting
#######################################################################################################

# Subset metadata to exclude spike-in samples for shotgun analysis
md_subset_noSpike <- subset(md_subset, md_subset$metagenomic_spike_in != TRUE)


# Subset metadata to include DNA concentrations, read counts, alpha-diversity metrics, host subject id, new variable
md_long <- md_subset_noSpike[,c(1:26, 64, 96)]


# Create wide form data for comparing extraction protocols (all samples)
md_wide <- pivot_wider(md_long,
                                       id_cols = c(sample_name_mantel_protocol,
                                                   biomass,
                                                   storage_solution, # remove this line if issue arises
                                                   sample_type,
                                                   sample_type2,
                                                   host_subject_id),
                                       names_from = extraction_protocol,
                                       values_from = c(dna_conc,
                                                       read_count_16s_deblur,
                                                       alpha_16S_richness,
                                                       alpha_16S_shannon,
                                                       alpha_16S_faithspd,
                                                       read_count_shotgun_wolr2pe,
                                                       alpha_shotgun_richness,
                                                       alpha_shotgun_shannon,
                                                       alpha_shotgun_faithspd,
                                                       alpha_lcmsms_richness,
                                                       alpha_lcmsms_shannon),
                                       values_fn = mean)


# Plot DNA concentrations
#######################################################################################################

# Create new variable combining biomass and sample type, and reorder levels
md_wide$biomass_storage_solution <- factor(paste(md_wide$biomass, md_wide$storage_solution, sep = "_"),
                                      levels = c("High biomass_EtOH",
                                                 "Low biomass_EtOH",
                                                 "High biomass_IPA",
                                                 "Low biomass_IPA"))

# Rename levels of new variable
levels(md_wide$biomass_storage_solution) <- c("High biomass, EtOH",
                                      "Low biomass, EtOH",
                                      "High biomass, IPA",
                                      "Low biomass, IPA")

# Create shape vector for new variable
biomass_storage_solution_shape <- c(21, 22, 24, 23)

# Check for negative values
summary(md_wide$dna_conc_Matrix)
summary(md_wide$dna_conc_MagMax)

# Test for normality
ggqqplot(md_wide$dna_conc_Matrix)
ggqqplot(md_wide$dna_conc_MagMax)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
plate_vs_matrix_spearman <- cor.test(x = md_wide$dna_conc_MagMax, y = md_wide$dna_conc_Matrix, method = "spearman")
plate_vs_matrix_kendall <- cor.test(x = md_wide$dna_conc_MagMax, y = md_wide$dna_conc_Matrix, method = "kendall")

## Per storage solution
plate_vs_matrix_spearman_etoh <- cor.test(x = subset(md_wide,
                                                md_wide$storage_solution == "EtOH")$dna_conc_MagMax, 
                                     y = subset(md_wide,
                                                md_wide$storage_solution == "EtOH")$dna_conc_Matrix, 
                                     method = "spearman")
plate_vs_matrix_spearman_ipa <- cor.test(x = subset(md_wide,
                                                     md_wide$storage_solution == "IPA")$dna_conc_MagMax, 
                                          y = subset(md_wide,
                                                     md_wide$storage_solution == "IPA")$dna_conc_Matrix, 
                                          method = "spearman")
plate_vs_matrix_kendall_etoh <- cor.test(x = subset(md_wide,
                                                     md_wide$storage_solution == "EtOH")$dna_conc_MagMax, 
                                          y = subset(md_wide,
                                                     md_wide$storage_solution == "EtOH")$dna_conc_Matrix, 
                                          method = "kendall")
plate_vs_matrix_kendall_ipa <- cor.test(x = subset(md_wide,
                                                    md_wide$storage_solution == "IPA")$dna_conc_MagMax, 
                                         y = subset(md_wide,
                                                    md_wide$storage_solution == "IPA")$dna_conc_Matrix, 
                                         method = "kendall")


# View correlations
plate_vs_matrix_spearman
plate_vs_matrix_kendall
plate_vs_matrix_spearman_etoh
plate_vs_matrix_kendall_etoh
plate_vs_matrix_spearman_ipa
plate_vs_matrix_kendall_ipa

# Run paired t-test (the data are paired, and we expect more DNA from the plated based method)
t.test(x = md_wide$dna_conc_MagMax, y = md_wide$dna_conc_Matrix, paired = TRUE, alternative = "greater")

## Per storage solution
t.test(x = subset(md_wide,
                  md_wide$storage_solution == "EtOH")$dna_conc_MagMax, 
       y = subset(md_wide,
                  md_wide$storage_solution == "EtOH")$dna_conc_Matrix, 
       paired = TRUE, 
       alternative = "greater")
t.test(x = subset(md_wide,
                  md_wide$storage_solution == "IPA")$dna_conc_MagMax, 
       y = subset(md_wide,
                  md_wide$storage_solution == "IPA")$dna_conc_Matrix, 
       paired = TRUE, 
       alternative = "greater")


# Run linear model to get slope and intercept
lm_protocol <- lm(md_wide$dna_conc_Matrix ~ md_wide$dna_conc_MagMax)
lm_protocol

## Per storage solution
lm_protocol_etoh <- lm(subset(md_wide,
                         md_wide$storage_solution == "EtOH")$dna_conc_Matrix ~ 
                    subset(md_wide,
                           md_wide$storage_solution == "EtOH")$dna_conc_MagMax)
lm_protocol_etoh
lm_protocol_ipa <- lm(subset(md_wide,
                         md_wide$storage_solution == "IPA")$dna_conc_Matrix ~ 
                    subset(md_wide,
                           md_wide$storage_solution == "IPA")$dna_conc_MagMax)
lm_protocol_ipa


# Plot Plate vs. Matrix - one smooth line for each storage solution
ggplot(md_wide, aes(x = dna_conc_MagMax,
                    y = dna_conc_Matrix)) +
  xlab("Plate-based method\nlog(DNA yield [ng/µL])") +
  ylab("Matrix tubes\nlog(DNA yield [ng/µL])") +
  geom_abline(linetype = "dashed") +
  geom_smooth(data = subset(md_wide, md_wide$storage_solution == "EtOH"),
              aes(x = dna_conc_MagMax,
              y = dna_conc_Matrix),
              color = "lightblue",
              fill = "lightblue",
              linewidth = 0.5,
              method = "lm") +
  geom_smooth(data = subset(md_wide, md_wide$storage_solution == "IPA"),
              aes(x = dna_conc_MagMax,
                  y = dna_conc_Matrix),
              color = "orange",
              fill = "orange",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass_storage_solution)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values =  biomass_storage_solution_shape) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,21,21,
                                                           21,21,21,21)))) +
  guides(shape = guide_legend(override.aes = list(shape = c(21,22,24,23)))) +
  labs(shape = "Sample biomass (B, C)", fill = "Sample type") +
  scale_x_log10() +
  scale_y_log10()


# Plot read counts - 16S deblur reads
#######################################################################################################
# Check for negative values
summary(md_wide$read_count_16s_deblur_Matrix)
summary(md_wide$read_count_16s_deblur_MagMax)

# Test for normality
ggqqplot(md_wide$read_count_16s_deblur_Matrix)
ggqqplot(md_wide$read_count_16s_deblur_MagMax)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
plate_vs_matrix_spearman <- cor.test(x = md_wide$read_count_16s_deblur_MagMax, y = md_wide$read_count_16s_deblur_Matrix, method = "spearman")
plate_vs_matrix_kendall <- cor.test(x = md_wide$read_count_16s_deblur_MagMax, y = md_wide$read_count_16s_deblur_Matrix, method = "kendall")

## Per storage solution
plate_vs_matrix_spearman_etoh <- cor.test(x = subset(md_wide,
                                                     md_wide$storage_solution == "EtOH")$read_count_16s_deblur_MagMax, 
                                          y = subset(md_wide,
                                                     md_wide$storage_solution == "EtOH")$read_count_16s_deblur_Matrix, 
                                          method = "spearman")
plate_vs_matrix_spearman_ipa <- cor.test(x = subset(md_wide,
                                                    md_wide$storage_solution == "IPA")$read_count_16s_deblur_MagMax, 
                                         y = subset(md_wide,
                                                    md_wide$storage_solution == "IPA")$read_count_16s_deblur_Matrix, 
                                         method = "spearman")
plate_vs_matrix_kendall_etoh <- cor.test(x = subset(md_wide,
                                                    md_wide$storage_solution == "EtOH")$read_count_16s_deblur_MagMax, 
                                         y = subset(md_wide,
                                                    md_wide$storage_solution == "EtOH")$read_count_16s_deblur_Matrix, 
                                         method = "kendall")
plate_vs_matrix_kendall_ipa <- cor.test(x = subset(md_wide,
                                                   md_wide$storage_solution == "IPA")$read_count_16s_deblur_MagMax, 
                                        y = subset(md_wide,
                                                   md_wide$storage_solution == "IPA")$read_count_16s_deblur_Matrix, 
                                        method = "kendall")

# View correlations
plate_vs_matrix_spearman
plate_vs_matrix_kendall
plate_vs_matrix_spearman_etoh
plate_vs_matrix_kendall_etoh
plate_vs_matrix_spearman_ipa
plate_vs_matrix_kendall_ipa


# Run paired t-test (the data are paired)
t.test(x = md_wide$read_count_16s_deblur_MagMax, 
       y = md_wide$read_count_16s_deblur_Matrix, 
       paired = TRUE, 
       alternative = "two.sided")


## Per storage solution
t.test(x = subset(md_wide,
                  md_wide$storage_solution == "EtOH")$read_count_16s_deblur_MagMax, 
       y = subset(md_wide,
                  md_wide$storage_solution == "EtOH")$read_count_16s_deblur_Matrix, 
       paired = TRUE, 
       alternative = "two.sided")
t.test(x = subset(md_wide,
                  md_wide$storage_solution == "IPA")$read_count_16s_deblur_MagMax, 
       y = subset(md_wide,
                  md_wide$storage_solution == "IPA")$read_count_16s_deblur_Matrix, 
       paired = TRUE, 
       alternative = "two.sided")


# Run linear model to get slope and intercept
lm_protocol <- lm(md_wide$read_count_16s_deblur_Matrix ~ md_wide$read_count_16s_deblur_MagMax)
lm_protocol

## Per storage solution
lm_protocol_etoh <- lm(subset(md_wide,
                              md_wide$storage_solution == "EtOH")$read_count_16s_deblur_Matrix ~ 
                         subset(md_wide,
                                md_wide$storage_solution == "EtOH")$read_count_16s_deblur_MagMax)
lm_protocol_etoh
lm_protocol_ipa <- lm(subset(md_wide,
                             md_wide$storage_solution == "IPA")$read_count_16s_deblur_Matrix ~ 
                        subset(md_wide,
                               md_wide$storage_solution == "IPA")$read_count_16s_deblur_MagMax)
lm_protocol_ipa


# Plot Plate vs. Matrix - one smooth line for each storage solution
ggplot(md_wide, aes(x = read_count_16s_deblur_MagMax,
                    y = read_count_16s_deblur_Matrix)) +
  xlab("Plate-based method\n16S Deblur reads") +
  ylab("Matrix tubes\n16S Deblur reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(data = subset(md_wide, md_wide$storage_solution == "EtOH"),
              aes(x = read_count_16s_deblur_MagMax,
                  y = read_count_16s_deblur_Matrix),
              color = "lightblue",
              fill = "lightblue",
              linewidth = 0.5,
              method = "lm") +
  geom_smooth(data = subset(md_wide, md_wide$storage_solution == "IPA"),
              aes(x = read_count_16s_deblur_MagMax,
                  y = read_count_16s_deblur_Matrix),
              color = "orange",
              fill = "orange",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass_storage_solution)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values =  biomass_storage_solution_shape) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,21,21,
                                                           21,21,21,21)))) +
  guides(shape = guide_legend(override.aes = list(shape = c(21,22,24,23)))) +
  labs(shape = "Sample biomass\nStorage solution (A, B)", fill = "Sample type")


# Plot read counts - shotgun WoLr2 PE reads
#######################################################################################################
# Check for negative values
summary(md_wide$read_count_shotgun_wolr2pe_Matrix)
summary(md_wide$read_count_shotgun_wolr2pe_MagMax)

# Test for normality
ggqqplot(md_wide$read_count_shotgun_wolr2pe_Matrix)
ggqqplot(md_wide$read_count_shotgun_wolr2pe_MagMax)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
plate_vs_matrix_spearman <- cor.test(x = md_wide$read_count_shotgun_wolr2pe_MagMax, y = md_wide$read_count_shotgun_wolr2pe_Matrix, method = "spearman")
plate_vs_matrix_kendall <- cor.test(x = md_wide$read_count_shotgun_wolr2pe_MagMax, y = md_wide$read_count_shotgun_wolr2pe_Matrix, method = "kendall")

## Per storage solution
plate_vs_matrix_spearman_etoh <- cor.test(x = subset(md_wide,
                                                     md_wide$storage_solution == "EtOH")$read_count_shotgun_wolr2pe_MagMax, 
                                          y = subset(md_wide,
                                                     md_wide$storage_solution == "EtOH")$read_count_shotgun_wolr2pe_Matrix, 
                                          method = "spearman")
plate_vs_matrix_spearman_ipa <- cor.test(x = subset(md_wide,
                                                    md_wide$storage_solution == "IPA")$read_count_shotgun_wolr2pe_MagMax, 
                                         y = subset(md_wide,
                                                    md_wide$storage_solution == "IPA")$read_count_shotgun_wolr2pe_Matrix, 
                                         method = "spearman")
plate_vs_matrix_kendall_etoh <- cor.test(x = subset(md_wide,
                                                    md_wide$storage_solution == "EtOH")$read_count_shotgun_wolr2pe_MagMax, 
                                         y = subset(md_wide,
                                                    md_wide$storage_solution == "EtOH")$read_count_shotgun_wolr2pe_Matrix, 
                                         method = "kendall")
plate_vs_matrix_kendall_ipa <- cor.test(x = subset(md_wide,
                                                   md_wide$storage_solution == "IPA")$read_count_shotgun_wolr2pe_MagMax, 
                                        y = subset(md_wide,
                                                   md_wide$storage_solution == "IPA")$read_count_shotgun_wolr2pe_Matrix, 
                                        method = "kendall")

# View correlations
plate_vs_matrix_spearman
plate_vs_matrix_kendall
plate_vs_matrix_spearman_etoh
plate_vs_matrix_kendall_etoh
plate_vs_matrix_spearman_ipa
plate_vs_matrix_kendall_ipa


# Run paired t-test (the data are paired)
t.test(x = md_wide$read_count_shotgun_wolr2pe_MagMax, y = md_wide$read_count_shotgun_wolr2pe_Matrix, paired = TRUE, alternative = "two.sided")

## Per storage solution
t.test(x = subset(md_wide,
                  md_wide$storage_solution == "EtOH")$read_count_shotgun_wolr2pe_MagMax, 
       y = subset(md_wide,
                  md_wide$storage_solution == "EtOH")$read_count_shotgun_wolr2pe_Matrix, 
       paired = TRUE, 
       alternative = "two.sided")
t.test(x = subset(md_wide,
                  md_wide$storage_solution == "IPA")$read_count_shotgun_wolr2pe_MagMax, 
       y = subset(md_wide,
                  md_wide$storage_solution == "IPA")$read_count_shotgun_wolr2pe_Matrix, 
       paired = TRUE, 
       alternative = "two.sided")



# Run linear model to get slope and intercept
lm_protocol <- lm(log10(md_wide$read_count_shotgun_wolr2pe_Matrix) ~ log10(md_wide$read_count_shotgun_wolr2pe_MagMax))
lm_protocol

## Per storage solution
lm_protocol_etoh <- lm(subset(md_wide,
                              md_wide$storage_solution == "EtOH")$read_count_shotgun_wolr2pe_Matrix ~ 
                         subset(md_wide,
                                md_wide$storage_solution == "EtOH")$read_count_shotgun_wolr2pe_MagMax)
lm_protocol_etoh
lm_protocol_ipa <- lm(subset(md_wide,
                             md_wide$storage_solution == "IPA")$read_count_shotgun_wolr2pe_Matrix ~ 
                        subset(md_wide,
                               md_wide$storage_solution == "IPA")$read_count_shotgun_wolr2pe_MagMax)
lm_protocol_ipa


# Plot Plate vs. Matrix - one smooth line for each storage solution
ggplot(md_wide, aes(x = read_count_shotgun_wolr2pe_MagMax,
                    y = read_count_shotgun_wolr2pe_Matrix)) +
  xlab("Plate-based method\nlog(metagenomic reads, WoL)") +
  ylab("Matrix tubes\nlog(metagenomic reads, WoL)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(data = subset(md_wide, md_wide$storage_solution == "EtOH"),
              aes(x = read_count_shotgun_wolr2pe_MagMax,
                  y = read_count_shotgun_wolr2pe_Matrix),
              color = "lightblue",
              fill = "lightblue",
              linewidth = 0.5,
              method = "lm") +
  geom_smooth(data = subset(md_wide, md_wide$storage_solution == "IPA"),
              aes(x = read_count_shotgun_wolr2pe_MagMax,
                  y = read_count_shotgun_wolr2pe_Matrix),
              color = "orange",
              fill = "orange",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass_storage_solution)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values =  biomass_storage_solution_shape) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,21,21,
                                                           21,21,21,21)))) +
  guides(shape = guide_legend(override.aes = list(shape = c(21,22,24,23)))) +
  labs(shape = "Sample biomass\nStorage solution", fill = "Sample type") + 
  scale_y_log10() +
  scale_x_log10()


# Plot alpha-diversity for 16S data
#######################################################################################################
# Check for negative values
summary(md_wide$alpha_16S_faithspd_Matrix)
summary(md_wide$alpha_16S_faithspd_MagMax)

# Test for normality
ggqqplot(md_wide$alpha_16S_faithspd_Matrix)
ggqqplot(md_wide$alpha_16S_faithspd_MagMax)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
plate_vs_matrix_alpha_16s_faithspd_spearman <- cor.test(x = md_wide$alpha_16S_faithspd_MagMax, y = md_wide$alpha_16S_faithspd_Matrix, method = "spearman")
plate_vs_matrix_alpha_16s_faithspd_kendall <- cor.test(x = md_wide$alpha_16S_faithspd_MagMax, y = md_wide$alpha_16S_faithspd_Matrix, method = "kendall")

## Per storage solution
plate_vs_matrix_spearman_etoh <- cor.test(x = subset(md_wide,
                                                     md_wide$storage_solution == "EtOH")$alpha_16S_faithspd_MagMax, 
                                          y = subset(md_wide,
                                                     md_wide$storage_solution == "EtOH")$alpha_16S_faithspd_Matrix, 
                                          method = "spearman")
plate_vs_matrix_spearman_ipa <- cor.test(x = subset(md_wide,
                                                    md_wide$storage_solution == "IPA")$alpha_16S_faithspd_MagMax, 
                                         y = subset(md_wide,
                                                    md_wide$storage_solution == "IPA")$alpha_16S_faithspd_Matrix, 
                                         method = "spearman")
plate_vs_matrix_kendall_etoh <- cor.test(x = subset(md_wide,
                                                    md_wide$storage_solution == "EtOH")$alpha_16S_faithspd_MagMax, 
                                         y = subset(md_wide,
                                                    md_wide$storage_solution == "EtOH")$alpha_16S_faithspd_Matrix, 
                                         method = "kendall")
plate_vs_matrix_kendall_ipa <- cor.test(x = subset(md_wide,
                                                   md_wide$storage_solution == "IPA")$alpha_16S_faithspd_MagMax, 
                                        y = subset(md_wide,
                                                   md_wide$storage_solution == "IPA")$alpha_16S_faithspd_Matrix, 
                                        method = "kendall")


# View correlations
plate_vs_matrix_alpha_16s_faithspd_spearman
plate_vs_matrix_alpha_16s_faithspd_kendall
plate_vs_matrix_spearman_etoh
plate_vs_matrix_kendall_etoh
plate_vs_matrix_spearman_ipa
plate_vs_matrix_kendall_ipa


# Run paired t-test (the data are paired)
t.test(x = md_wide$alpha_16S_faithspd_MagMax, y = md_wide$alpha_16S_faithspd_Matrix, paired = TRUE, alternative = "two.sided")

## Per storage solution
t.test(x = subset(md_wide,
                  md_wide$storage_solution == "EtOH")$alpha_16S_faithspd_MagMax, 
       y = subset(md_wide,
                  md_wide$storage_solution == "EtOH")$alpha_16S_faithspd_Matrix, 
       paired = TRUE, 
       alternative = "two.sided")
t.test(x = subset(md_wide,
                  md_wide$storage_solution == "IPA")$alpha_16S_faithspd_MagMax, 
       y = subset(md_wide,
                  md_wide$storage_solution == "IPA")$alpha_16S_faithspd_Matrix, 
       paired = TRUE, 
       alternative = "two.sided")


# Run linear model to get slope and intercept
lm_protocol_alpha_16s_faithspd <- lm(md_wide$alpha_16S_faithspd_Matrix ~ md_wide$alpha_16S_faithspd_MagMax)
lm_protocol_alpha_16s_faithspd

## Per storage solution
lm_protocol_etoh <- lm(subset(md_wide,
                              md_wide$storage_solution == "EtOH")$alpha_16S_faithspd_Matrix ~ 
                         subset(md_wide,
                                md_wide$storage_solution == "EtOH")$alpha_16S_faithspd_MagMax)
lm_protocol_etoh
lm_protocol_ipa <- lm(subset(md_wide,
                             md_wide$storage_solution == "IPA")$alpha_16S_faithspd_Matrix ~ 
                        subset(md_wide,
                               md_wide$storage_solution == "IPA")$alpha_16S_faithspd_MagMax)
lm_protocol_ipa


# Plot Faith's PD - Plate vs. Matrix - one smooth line for each storage solution
ggplot(md_wide, aes(x = alpha_16S_faithspd_MagMax,
                    y = alpha_16S_faithspd_Matrix)) +
  xlab("Plate-based method\n16S Faith's PD") +
  ylab("Matrix tubes\n16S Faith's PD") +
  geom_abline(linetype = "dashed") +
  geom_smooth(data = subset(md_wide, md_wide$storage_solution == "EtOH"),
              aes(x = alpha_16S_faithspd_MagMax,
                  y = alpha_16S_faithspd_Matrix),
              color = "lightblue",
              fill = "lightblue",
              linewidth = 0.5,
              method = "lm") +
  geom_smooth(data = subset(md_wide, md_wide$storage_solution == "IPA"),
              aes(x = alpha_16S_faithspd_MagMax,
                  y = alpha_16S_faithspd_Matrix),
              color = "orange",
              fill = "orange",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass_storage_solution)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values =  biomass_storage_solution_shape) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,21,21,
                                                           21,21,21,21)))) +
  guides(shape = guide_legend(override.aes = list(shape = c(21,22,24,23)))) +
  labs(shape = "Sample biomass (C-G)", fill = "Sample type")


# Plot alpha-diversity for shotgun data
#######################################################################################################
# Check for negative values
summary(md_wide$alpha_shotgun_richness_Matrix)
summary(md_wide$alpha_shotgun_shannon_Matrix)
summary(md_wide$alpha_shotgun_faithspd_Matrix)
summary(md_wide$alpha_shotgun_richness_MagMax)
summary(md_wide$alpha_shotgun_shannon_MagMax)
summary(md_wide$alpha_shotgun_faithspd_MagMax)

# Test for normality
ggqqplot(md_wide$alpha_shotgun_richness_Matrix)
ggqqplot(md_wide$alpha_shotgun_shannon_Matrix)
ggqqplot(md_wide$alpha_shotgun_faithspd_Matrix)
ggqqplot(md_wide$alpha_shotgun_richness_MagMax)
ggqqplot(md_wide$alpha_shotgun_shannon_MagMax)
ggqqplot(md_wide$alpha_shotgun_faithspd_MagMax)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
plate_vs_matrix_alpha_shotgun_richness_spearman <- cor.test(x = md_wide$alpha_shotgun_richness_MagMax, y = md_wide$alpha_shotgun_richness_Matrix, method = "spearman")
plate_vs_matrix_alpha_shotgun_richness_kendall <- cor.test(x = md_wide$alpha_shotgun_richness_MagMax, y = md_wide$alpha_shotgun_richness_Matrix, method = "kendall")
plate_vs_matrix_alpha_shotgun_shannon_spearman <- cor.test(x = md_wide$alpha_shotgun_shannon_MagMax, y = md_wide$alpha_shotgun_shannon_Matrix, method = "spearman")
plate_vs_matrix_alpha_shotgun_shannon_kendall <- cor.test(x = md_wide$alpha_shotgun_shannon_MagMax, y = md_wide$alpha_shotgun_shannon_Matrix, method = "kendall")
plate_vs_matrix_alpha_shotgun_faithspd_spearman <- cor.test(x = md_wide$alpha_shotgun_faithspd_MagMax, y = md_wide$alpha_shotgun_faithspd_Matrix, method = "spearman")
plate_vs_matrix_alpha_shotgun_faithspd_kendall <- cor.test(x = md_wide$alpha_shotgun_faithspd_MagMax, y = md_wide$alpha_shotgun_faithspd_Matrix, method = "kendall")

## Per storage solution
plate_vs_matrix_spearman_etoh <- cor.test(x = subset(md_wide,
                                                     md_wide$storage_solution == "EtOH")$alpha_shotgun_richness_MagMax, 
                                          y = subset(md_wide,
                                                     md_wide$storage_solution == "EtOH")$alpha_shotgun_richness_Matrix, 
                                          method = "spearman")
plate_vs_matrix_spearman_ipa <- cor.test(x = subset(md_wide,
                                                    md_wide$storage_solution == "IPA")$alpha_shotgun_richness_MagMax, 
                                         y = subset(md_wide,
                                                    md_wide$storage_solution == "IPA")$alpha_shotgun_richness_Matrix, 
                                         method = "spearman")
plate_vs_matrix_kendall_etoh <- cor.test(x = subset(md_wide,
                                                    md_wide$storage_solution == "EtOH")$alpha_shotgun_richness_MagMax, 
                                         y = subset(md_wide,
                                                    md_wide$storage_solution == "EtOH")$alpha_shotgun_richness_Matrix, 
                                         method = "kendall")
plate_vs_matrix_kendall_ipa <- cor.test(x = subset(md_wide,
                                                   md_wide$storage_solution == "IPA")$alpha_shotgun_richness_MagMax, 
                                        y = subset(md_wide,
                                                   md_wide$storage_solution == "IPA")$alpha_shotgun_richness_Matrix, 
                                        method = "kendall")

# View correlations
plate_vs_matrix_alpha_shotgun_richness_spearman
plate_vs_matrix_alpha_shotgun_richness_kendall
plate_vs_matrix_alpha_shotgun_shannon_spearman
plate_vs_matrix_alpha_shotgun_shannon_kendall
plate_vs_matrix_alpha_shotgun_faithspd_spearman
plate_vs_matrix_alpha_shotgun_faithspd_kendall
plate_vs_matrix_spearman_etoh
plate_vs_matrix_kendall_etoh
plate_vs_matrix_spearman_ipa
plate_vs_matrix_kendall_ipa


# Run paired t-test (the data are paired)
t.test(x = md_wide$alpha_shotgun_richness_MagMax, y = md_wide$alpha_shotgun_richness_Matrix, paired = TRUE, alternative = "two.sided")
t.test(x = md_wide$alpha_shotgun_shannon_MagMax, y = md_wide$alpha_shotgun_shannon_Matrix, paired = TRUE, alternative = "two.sided")
t.test(x = md_wide$alpha_shotgun_faithspd_MagMax, y = md_wide$alpha_shotgun_faithspd_Matrix, paired = TRUE, alternative = "two.sided")

## Per storage solution
t.test(x = subset(md_wide,
                  md_wide$storage_solution == "EtOH")$alpha_shotgun_richness_MagMax, 
       y = subset(md_wide,
                  md_wide$storage_solution == "EtOH")$alpha_shotgun_richness_Matrix, 
       paired = TRUE, 
       alternative = "two.sided")
t.test(x = subset(md_wide,
                  md_wide$storage_solution == "IPA")$alpha_shotgun_richness_MagMax, 
       y = subset(md_wide,
                  md_wide$storage_solution == "IPA")$alpha_shotgun_richness_Matrix, 
       paired = TRUE, 
       alternative = "two.sided")


# Run linear model to get slope and intercept
lm_protocol_alpha_shotgun_richness <- lm(md_wide$alpha_shotgun_richness_Matrix ~ md_wide$alpha_shotgun_richness_MagMax)
lm_protocol_alpha_shotgun_richness
lm_protocol_alpha_shotgun_shannon <- lm(md_wide$alpha_shotgun_shannon_Matrix ~ md_wide$alpha_shotgun_shannon_MagMax)
lm_protocol_alpha_shotgun_shannon
lm_protocol_alpha_shotgun_faithspd <- lm(md_wide$alpha_shotgun_faithspd_Matrix ~ md_wide$alpha_shotgun_faithspd_MagMax)
lm_protocol_alpha_shotgun_faithspd

## Per storage solution
lm_protocol_etoh <- lm(subset(md_wide,
                              md_wide$storage_solution == "EtOH")$alpha_shotgun_richness_Matrix ~ 
                         subset(md_wide,
                                md_wide$storage_solution == "EtOH")$alpha_shotgun_richness_MagMax)
lm_protocol_etoh
lm_protocol_ipa <- lm(subset(md_wide,
                             md_wide$storage_solution == "IPA")$alpha_shotgun_richness_Matrix ~ 
                        subset(md_wide,
                               md_wide$storage_solution == "IPA")$alpha_shotgun_richness_MagMax)
lm_protocol_ipa


# Plot Faith's PD - Plate vs. Matrix - one smooth line for each storage solution
ggplot(md_wide, aes(x = alpha_shotgun_faithspd_MagMax,
                    y = alpha_shotgun_faithspd_Matrix)) +
  xlab("Plate-based method\nmetagenomic Faith's PD") +
  ylab("Matrix tubes\nmetagenomic Faith's PD") +
  geom_abline(linetype = "dashed") +
  geom_smooth(data = subset(md_wide, md_wide$storage_solution == "EtOH"),
              aes(x = alpha_shotgun_faithspd_MagMax,
                  y = alpha_shotgun_faithspd_Matrix),
              color = "lightblue",
              fill = "lightblue",
              linewidth = 0.5,
              method = "lm") +
  geom_smooth(data = subset(md_wide, md_wide$storage_solution == "IPA"),
              aes(x = alpha_shotgun_faithspd_MagMax,
                  y = alpha_shotgun_faithspd_Matrix),
              color = "orange",
              fill = "orange",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass_storage_solution)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values =  biomass_storage_solution_shape) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,21,21,
                                                           21,21,21,21)))) +
  guides(shape = guide_legend(override.aes = list(shape = c(21,22,24,23)))) +
  labs(shape = "Sample biomass\nStorage solution", fill = "Sample type")


#########################################################################################################
# Scatter plots comparing DNA yield, read counts, and alpha-diversity between storage solutions for each method
#########################################################################################################
install.packages('ggpubr')
install.packages('cowplot')
install.packages('svglite')
install.packages('scales')

library(tidyverse)
library(ggpubr)
library(cowplot)
library(svglite)
library(scales)


# Tidy-up metadata for running correlations and plotting
#######################################################################################################

# Subset metadata to exclude spike-in samples for shotgun analysis
md_subset_noSpike <- subset(md_subset, md_subset$metagenomic_spike_in != TRUE)


# Subset metadata to include DNA concentrations, read counts, alpha-diversity metrics, host subject id, new variable
md_long <- md_subset_noSpike[,c(1:26, 64, 96)]


# Find sample sizes
md_long_subset_plate <- subset(md_long, md_long$extraction_protocol == "MagMax") # 132
md_long_subset_matrix <- subset(md_long, md_long$extraction_protocol == "Matrix") # 131
md_long_shotgun_subset_plate <- subset(md_long_shotgun, md_long_shotgun$extraction_protocol == "MagMax") # 156
md_long_shotgun_subset_matrix <- subset(md_long_shotgun, md_long_shotgun$extraction_protocol == "Matrix") # 155


# Create wide form data for comparing storage solutions within each extraction protocol
md_wide_storage <- pivot_wider(md_long,
                       id_cols = c(sample_name_mantel_solution,
                                   biomass,
                                   sample_type,
                                   sample_type2,
                                   host_subject_id),
                       names_from = c(extraction_protocol,
                                      storage_solution),
                       values_from = c(dna_conc,
                                       read_count_16s_deblur,
                                       alpha_16S_richness,
                                       alpha_16S_shannon,
                                       alpha_16S_faithspd,
                                       read_count_shotgun_wolr2pe,
                                       alpha_shotgun_richness,
                                       alpha_shotgun_shannon,
                                       alpha_shotgun_faithspd,
                                       alpha_lcmsms_richness,
                                       alpha_lcmsms_shannon),
                       values_fn = mean)

# Plot DNA concentrations
#######################################################################################################

# Check for negative values
summary(md_wide_storage$dna_conc_Matrix_EtOH)
summary(md_wide_storage$dna_conc_Matrix_IPA)


# Test for normality
ggqqplot(md_wide_storage$dna_conc_Matrix_EtOH)
ggqqplot(md_wide_storage$dna_conc_Matrix_IPA)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
matrix_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage$dna_conc_Matrix_EtOH, y = md_wide_storage$dna_conc_Matrix_IPA, method = "spearman")
matrix_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage$dna_conc_Matrix_EtOH, y = md_wide_storage$dna_conc_Matrix_IPA, method = "kendall")

plate_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage$dna_conc_MagMax_EtOH, y = md_wide_storage$dna_conc_MagMax_IPA, method = "spearman")
plate_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage$dna_conc_MagMax_EtOH, y = md_wide_storage$dna_conc_MagMax_IPA, method = "kendall")


# View correlations
matrix_etoh_vs_ipa_spearman
matrix_etoh_vs_ipa_kendall

plate_etoh_vs_ipa_spearman
plate_etoh_vs_ipa_kendall


# Run paired t-test (the data are paired, and we expect more DNA from the plated based method)
t.test(x = md_wide_storage$dna_conc_Matrix_EtOH, y = md_wide_storage$dna_conc_Matrix_IPA, paired = TRUE, alternative = "two.sided")
t.test(x = md_wide_storage$dna_conc_MagMax_EtOH, y = md_wide_storage$dna_conc_MagMax_IPA, paired = TRUE, alternative = "two.sided")


# Run linear model to get slope and intercept
lm_protocol_matrix <- lm(md_wide_storage$dna_conc_Matrix_IPA ~ md_wide_storage$dna_conc_Matrix_EtOH)
lm_protocol_plate <- lm(md_wide_storage$dna_conc_MagMax_IPA ~ md_wide_storage$dna_conc_MagMax_EtOH)
lm_protocol_matrix
lm_protocol_plate


# Plot IPA vs. EtOH for the Matrix method
ggplot(md_wide_storage, aes(x = dna_conc_Matrix_EtOH,
                            y = dna_conc_Matrix_IPA)) +
  #geom_point() +
  #facet_wrap(~sample_type_2,
  #           scales = "free") +
  xlab("Matrix tubes\nlog(DNA yield from EtOH [ng/µL])") +
  ylab("Matrix tubes\nlog(DNA yield from IPA [ng/µL])") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_log10() +
  scale_y_log10()


# Plot IPA vs. EtOH for the Plate-based method
ggplot(md_wide_storage, aes(x = dna_conc_MagMax_EtOH,
                            y = dna_conc_MagMax_IPA)) +
  xlab("Plate-based\nlog(DNA yield from EtOH [ng/µL])") +
  ylab("Plate-based\nlog(DNA yield from IPA [ng/µL])") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_log10() +
  scale_y_log10()


# Plot read counts for 16S data
#######################################################################################################

# Check for negative values
summary(md_wide_storage$read_count_16s_deblur_Matrix_EtOH)
summary(md_wide_storage$read_count_16s_deblur_Matrix_IPA)
summary(md_wide_storage$read_count_16s_deblur_MagMax_EtOH)
summary(md_wide_storage$read_count_16s_deblur_MagMax_IPA)


# Test for normality
ggqqplot(md_wide_storage$read_count_16s_deblur_Matrix_EtOH)
ggqqplot(md_wide_storage$read_count_16s_deblur_Matrix_IPA)
ggqqplot(md_wide_storage$read_count_16s_deblur_MagMax_EtOH)
ggqqplot(md_wide_storage$read_count_16s_deblur_MagMax_IPA)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
matrix_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage$read_count_16s_deblur_Matrix_EtOH, y = md_wide_storage$read_count_16s_deblur_Matrix_IPA, method = "spearman")
matrix_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage$read_count_16s_deblur_Matrix_EtOH, y = md_wide_storage$read_count_16s_deblur_Matrix_IPA, method = "kendall")

plate_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage$read_count_16s_deblur_MagMax_EtOH, y = md_wide_storage$read_count_16s_deblur_MagMax_IPA, method = "spearman")
plate_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage$read_count_16s_deblur_MagMax_EtOH, y = md_wide_storage$read_count_16s_deblur_MagMax_IPA, method = "kendall")


# View correlations
matrix_etoh_vs_ipa_spearman
matrix_etoh_vs_ipa_kendall

plate_etoh_vs_ipa_spearman
plate_etoh_vs_ipa_kendall


# Run paired t-test (the data are paired, and we expect more DNA from the plated based method)
t.test(x = md_wide_storage$read_count_16s_deblur_Matrix_EtOH, y = md_wide_storage$read_count_16s_deblur_Matrix_IPA, paired = TRUE, alternative = "two.sided")
t.test(x = md_wide_storage$read_count_16s_deblur_MagMax_EtOH, y = md_wide_storage$read_count_16s_deblur_MagMax_IPA, paired = TRUE, alternative = "two.sided")


# Run linear model to get slope and intercept
lm_protocol_matrix <- lm(md_wide_storage$read_count_16s_deblur_Matrix_IPA ~ md_wide_storage$read_count_16s_deblur_Matrix_EtOH)
lm_protocol_plate <- lm(md_wide_storage$read_count_16s_deblur_MagMax_IPA ~ md_wide_storage$read_count_16s_deblur_MagMax_EtOH)
lm_protocol_matrix
lm_protocol_plate


# Plot IPA vs. EtOH for the Matrix method
ggplot(md_wide_storage, aes(x = read_count_16s_deblur_Matrix_EtOH,
                            y = read_count_16s_deblur_Matrix_IPA)) +
  xlab("EtOH, matrix tubes\n16S Deblur reads") +
  ylab("IPA, matrix tubes\n16S Deblur reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Plot IPA vs. EtOH for the Plate-based method
ggplot(md_wide_storage, aes(x = read_count_16s_deblur_MagMax_EtOH,
                            y = read_count_16s_deblur_MagMax_IPA)) +
  xlab("EtOH, plate-based\n16S Deblur reads") +
  ylab("IPA, plate-based\n16S Deblur reads") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Plot read counts for shotgun data
#######################################################################################################

# Check for negative values
summary(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_EtOH)
summary(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_IPA)
summary(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_EtOH)
summary(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_IPA)


# Test for normality
ggqqplot(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_EtOH)
ggqqplot(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_IPA)
ggqqplot(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_EtOH)
ggqqplot(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_IPA)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
matrix_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_EtOH, y = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_IPA, method = "spearman")
matrix_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_EtOH, y = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_IPA, method = "kendall")

plate_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_EtOH, y = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_IPA, method = "spearman")
plate_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_EtOH, y = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_IPA, method = "kendall")


# View correlations
matrix_etoh_vs_ipa_spearman
matrix_etoh_vs_ipa_kendall

plate_etoh_vs_ipa_spearman
plate_etoh_vs_ipa_kendall


# Run paired t-test (the data are paired, and we expect more DNA from the plated based method)
t.test(x = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_EtOH, y = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_IPA, paired = TRUE, alternative = "two.sided")
t.test(x = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_EtOH, y = md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_IPA, paired = TRUE, alternative = "two.sided")


# Run linear model to get slope and intercept
lm_protocol_matrix <- lm(log10(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_IPA) ~ log10(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_Matrix_EtOH))
lm_protocol_plate <- lm(log10(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_IPA) ~ log10(md_wide_storage_shotgun$read_count_shotgun_wolr2pe_MagMax_EtOH))
lm_protocol_matrix
lm_protocol_plate


# Plot IPA vs. EtOH for the Matrix method
ggplot(md_wide_storage_shotgun, aes(x = read_count_shotgun_wolr2pe_Matrix_EtOH,
                                    y = read_count_shotgun_wolr2pe_Matrix_IPA)) +
  xlab("EtOH, matrix tubes\nlog(metagenomic reads, WoL)") +
  ylab("IPA, matrix tubes\nlog(metagenomic reads, WoL)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_log10() +
  scale_y_log10()


# Plot IPA vs. EtOH for the Plate-based method
ggplot(md_wide_storage_shotgun, aes(x = read_count_shotgun_wolr2pe_MagMax_EtOH,
                                    y = read_count_shotgun_wolr2pe_MagMax_IPA)) +
  xlab("EtOH, plate-based\nlog(metagenomic reads, WoL)") +
  ylab("IPA, plate-based\nlog(metagenomic reads, WoL)") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type") +
  scale_x_log10() +
  scale_y_log10()


# Plot alpha-diversity for 16S data
#######################################################################################################

# Check for negative values
summary(md_wide_storage$alpha_16S_faithspd_Matrix_EtOH)
summary(md_wide_storage$alpha_16S_faithspd_Matrix_IPA)
summary(md_wide_storage$alpha_16S_faithspd_MagMax_EtOH)
summary(md_wide_storage$alpha_16S_faithspd_MagMax_IPA)


# Test for normality
ggqqplot(md_wide_storage$alpha_16S_faithspd_Matrix_EtOH)
ggqqplot(md_wide_storage$alpha_16S_faithspd_Matrix_IPA)
ggqqplot(md_wide_storage$alpha_16S_faithspd_MagMax_EtOH)
ggqqplot(md_wide_storage$alpha_16S_faithspd_MagMax_IPA)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
matrix_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage$alpha_16S_faithspd_Matrix_EtOH, y = md_wide_storage$alpha_16S_faithspd_Matrix_IPA, method = "spearman")
matrix_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage$alpha_16S_faithspd_Matrix_EtOH, y = md_wide_storage$alpha_16S_faithspd_Matrix_IPA, method = "kendall")

plate_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage$alpha_16S_faithspd_MagMax_EtOH, y = md_wide_storage$alpha_16S_faithspd_MagMax_IPA, method = "spearman")
plate_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage$alpha_16S_faithspd_MagMax_EtOH, y = md_wide_storage$alpha_16S_faithspd_MagMax_IPA, method = "kendall")


# View correlations
matrix_etoh_vs_ipa_spearman
matrix_etoh_vs_ipa_kendall

plate_etoh_vs_ipa_spearman
plate_etoh_vs_ipa_kendall


# Run paired t-test (the data are paired, and we expect more DNA from the plated based method)
t.test(x = md_wide_storage$alpha_16S_faithspd_Matrix_EtOH, y = md_wide_storage$alpha_16S_faithspd_Matrix_IPA, paired = TRUE, alternative = "two.sided")
t.test(x = md_wide_storage$alpha_16S_faithspd_MagMax_EtOH, y = md_wide_storage$alpha_16S_faithspd_MagMax_IPA, paired = TRUE, alternative = "two.sided")


# Run linear model to get slope and intercept
lm_protocol_matrix <- lm(md_wide_storage$alpha_16S_faithspd_Matrix_IPA ~ md_wide_storage$alpha_16S_faithspd_Matrix_EtOH)
lm_protocol_plate <- lm(md_wide_storage$alpha_16S_faithspd_MagMax_IPA ~ md_wide_storage$alpha_16S_faithspd_MagMax_EtOH)
lm_protocol_matrix
lm_protocol_plate


# Plot IPA vs. EtOH for the Matrix method
ggplot(md_wide_storage, aes(x = alpha_16S_faithspd_Matrix_EtOH,
                            y = alpha_16S_faithspd_Matrix_IPA)) +
  xlab("EtOH, matrix tubes\n16S Faith's PD") +
  ylab("IPA, matrix tubes\n16S Faith's PD") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Plot IPA vs. EtOH for the Plate-based method
ggplot(md_wide_storage, aes(x = alpha_16S_faithspd_MagMax_EtOH,
                            y = alpha_16S_faithspd_MagMax_IPA)) +
  xlab("EtOH, plate-based\n16S Faith's PD") +
  ylab("IPA, plate-based\n16S Faith's PD") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Plot alpha-diversity for shotgun data
#######################################################################################################

# Check for negative values
summary(md_wide_storage$alpha_shotgun_faithspd_Matrix_EtOH)
summary(md_wide_storage$alpha_shotgun_faithspd_Matrix_IPA)
summary(md_wide_storage$alpha_shotgun_faithspd_MagMax_EtOH)
summary(md_wide_storage$alpha_shotgun_faithspd_MagMax_IPA)


# Test for normality
ggqqplot(md_wide_storage$alpha_shotgun_faithspd_Matrix_EtOH)
ggqqplot(md_wide_storage$alpha_shotgun_faithspd_Matrix_IPA)
ggqqplot(md_wide_storage$alpha_shotgun_faithspd_MagMax_EtOH)
ggqqplot(md_wide_storage$alpha_shotgun_faithspd_MagMax_IPA)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
matrix_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage$alpha_shotgun_faithspd_Matrix_EtOH, y = md_wide_storage$alpha_shotgun_faithspd_Matrix_IPA, method = "spearman")
matrix_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage$alpha_shotgun_faithspd_Matrix_EtOH, y = md_wide_storage$alpha_shotgun_faithspd_Matrix_IPA, method = "kendall")

plate_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage$alpha_shotgun_faithspd_MagMax_EtOH, y = md_wide_storage$alpha_shotgun_faithspd_MagMax_IPA, method = "spearman")
plate_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage$alpha_shotgun_faithspd_MagMax_EtOH, y = md_wide_storage$alpha_shotgun_faithspd_MagMax_IPA, method = "kendall")


# View correlations
matrix_etoh_vs_ipa_spearman
matrix_etoh_vs_ipa_kendall

plate_etoh_vs_ipa_spearman
plate_etoh_vs_ipa_kendall


# Run paired t-test (the data are paired, and we expect more DNA from the plated based method)
t.test(x = md_wide_storage$alpha_shotgun_faithspd_Matrix_EtOH, y = md_wide_storage$alpha_shotgun_faithspd_Matrix_IPA, paired = TRUE, alternative = "two.sided")
t.test(x = md_wide_storage$alpha_shotgun_faithspd_MagMax_EtOH, y = md_wide_storage$alpha_shotgun_faithspd_MagMax_IPA, paired = TRUE, alternative = "two.sided")


# Run linear model to get slope and intercept
lm_protocol_matrix <- lm(log10(md_wide_storage$alpha_shotgun_faithspd_Matrix_IPA) ~ log10(md_wide_storage$alpha_shotgun_faithspd_Matrix_EtOH))
lm_protocol_plate <- lm(log10(md_wide_storage$alpha_shotgun_faithspd_MagMax_IPA) ~ log10(md_wide_storage$alpha_shotgun_faithspd_MagMax_EtOH))
lm_protocol_matrix
lm_protocol_plate


# Plot IPA vs. EtOH for the Matrix method
ggplot(md_wide_storage, aes(x = alpha_shotgun_faithspd_Matrix_EtOH,
                                    y = alpha_shotgun_faithspd_Matrix_IPA)) +
  xlab("EtOH, matrix tubes\nmetagenomic Faith's PD") +
  ylab("IPA, matrix tubes\nmetagenomic Faith's PD") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Plot IPA vs. EtOH for the Plate-based method
ggplot(md_wide_storage, aes(x = alpha_shotgun_faithspd_MagMax_EtOH,
                                    y = alpha_shotgun_faithspd_MagMax_IPA)) +
  xlab("EtOH, plate-based\nmetagenomic Faith's PD") +
  ylab("IPA, plate-based\nmetagenomic Faith's PD") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


# Plot alpha-diversity for metabolomics data
#######################################################################################################

# Check for negative values
summary(md_wide_storage$alpha_lcmsms_richness_Matrix_EtOH)
summary(md_wide_storage$alpha_lcmsms_richness_Matrix_IPA)
summary(md_wide_storage$alpha_lcmsms_richness_MagMax_EtOH)
summary(md_wide_storage$alpha_lcmsms_richness_MagMax_IPA)


# Test for normality
ggqqplot(md_wide_storage$alpha_lcmsms_richness_Matrix_EtOH)
ggqqplot(md_wide_storage$alpha_lcmsms_richness_Matrix_IPA)
ggqqplot(md_wide_storage$alpha_lcmsms_richness_MagMax_EtOH)
ggqqplot(md_wide_storage$alpha_lcmsms_richness_MagMax_IPA)

## Confirmed the data are non-normal - use non-parametric correlation test


# Run correlations
matrix_etoh_vs_ipa_spearman <- cor.test(x = md_wide_storage$alpha_lcmsms_richness_Matrix_EtOH, y = md_wide_storage$alpha_lcmsms_richness_Matrix_IPA, method = "spearman")
matrix_etoh_vs_ipa_kendall <- cor.test(x = md_wide_storage$alpha_lcmsms_richness_Matrix_EtOH, y = md_wide_storage$alpha_lcmsms_richness_Matrix_IPA, method = "kendall")


# View correlations
matrix_etoh_vs_ipa_spearman
matrix_etoh_vs_ipa_kendall


# Run paired t-test (the data are paired, and we expect more DNA from the plated based method)
t.test(x = md_wide_storage$alpha_lcmsms_richness_Matrix_EtOH, y = md_wide_storage$alpha_lcmsms_richness_Matrix_IPA, paired = TRUE, alternative = "two.sided")


# Run linear model to get slope and intercept
lm_protocol_matrix <- lm(log10(md_wide_storage$alpha_lcmsms_richness_Matrix_IPA) ~ log10(md_wide_storage$alpha_lcmsms_richness_Matrix_EtOH))
lm_protocol_matrix


# Plot IPA vs. EtOH for the Matrix method - Shannon
ggplot(md_wide_storage, aes(x = alpha_lcmsms_shannon_Matrix_EtOH,
                            y = alpha_lcmsms_shannon_Matrix_IPA)) +
  xlab("EtOH, matrix tubes\nmetabolite Shannon diversity") +
  ylab("IPA, matrix tubes\nmetabolite Shannon diversity") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")

# Plot IPA vs. EtOH for the Matrix method - richness
ggplot(md_wide_storage, aes(x = alpha_lcmsms_richness_Matrix_EtOH,
                            y = alpha_lcmsms_richness_Matrix_IPA)) +
  xlab("EtOH, matrix tubes\nmetabolite richness") +
  ylab("IPA, matrix tubes\nmetabolite richness") +
  geom_abline(linetype = "dashed") +
  geom_smooth(color = "gray60",
              linewidth = 0.5,
              method = "lm") +
  geom_point(aes(fill = sample_type2,
                 shape = biomass)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.background =element_rect(fill="white")) +
  scale_fill_manual(values = sample_type_2_colors) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill=guide_legend(override.aes=list(shape=c(21,21,21,21,
                                                     22,22,22,22)))) +
  labs(shape="Sample biomass", fill="Sample type")


####################################################################################################
# Plot within sample distances (technical replicates)
####################################################################################################

# Read in data files
tech_rep_16s_hbm_jaccard <- read_tsv("matrix_tech_rep_16S_hbm_jaccard.txt", na = c("", "NA", "not applicable"))
tech_rep_16s_hbm_unifrac <- read_tsv("matrix_tech_rep_16S_hbm_unifrac.txt", na = c("", "NA", "not applicable"))
tech_rep_16s_hbm_wunifrac <- read_tsv("matrix_tech_rep_16S_hbm_weighted_unifrac.txt", na = c("", "NA", "not applicable"))
tech_rep_16s_hbm_rpca <- read_tsv("matrix_tech_rep_16S_hbm_rpca.txt", na = c("", "NA", "not applicable"))
tech_rep_16s_hbm_phylo_rpca <- read_tsv("matrix_tech_rep_16S_hbm_phylo_rpca.txt", na = c("", "NA", "not applicable"))
tech_rep_16s_lbm_jaccard <- read_tsv("matrix_tech_rep_16S_lbm_jaccard.txt", na = c("", "NA", "not applicable"))
tech_rep_16s_lbm_unifrac <- read_tsv("matrix_tech_rep_16S_lbm_unifrac.txt", na = c("", "NA", "not applicable"))
tech_rep_16s_lbm_wunifrac <- read_tsv("matrix_tech_rep_16S_lbm_weighted_unifrac.txt", na = c("", "NA", "not applicable"))
tech_rep_16s_lbm_rpca <- read_tsv("matrix_tech_rep_16S_lbm_rpca.txt", na = c("", "NA", "not applicable"))
tech_rep_16s_lbm_phylo_rpca <- read_tsv("matrix_tech_rep_16S_lbm_phylo_rpca.txt", na = c("", "NA", "not applicable"))
tech_rep_shotgun_hbm_jaccard <- read_tsv("matrix_tech_rep_shotgun_hbm_jaccard.txt", na = c("", "NA", "not applicable"))
tech_rep_shotgun_hbm_unifrac <- read_tsv("matrix_tech_rep_shotgun_hbm_unifrac.txt", na = c("", "NA", "not applicable"))
tech_rep_shotgun_hbm_wunifrac <- read_tsv("matrix_tech_rep_shotgun_hbm_weighted_unifrac.txt", na = c("", "NA", "not applicable"))
tech_rep_shotgun_hbm_rpca <- read_tsv("matrix_tech_rep_shotgun_hbm_rpca.txt", na = c("", "NA", "not applicable"))
tech_rep_shotgun_hbm_phylo_rpca <- read_tsv("matrix_tech_rep_shotgun_hbm_phylo_rpca.txt", na = c("", "NA", "not applicable"))
tech_rep_shotgun_lbm_jaccard <- read_tsv("matrix_tech_rep_shotgun_lbm_jaccard.txt", na = c("", "NA", "not applicable"))
tech_rep_shotgun_lbm_unifrac <- read_tsv("matrix_tech_rep_shotgun_lbm_unifrac.txt", na = c("", "NA", "not applicable"))
tech_rep_shotgun_lbm_wunifrac <- read_tsv("matrix_tech_rep_shotgun_lbm_weighted_unifrac.txt", na = c("", "NA", "not applicable"))
tech_rep_shotgun_lbm_rpca <- read_tsv("matrix_tech_rep_shotgun_lbm_rpca.txt", na = c("", "NA", "not applicable"))
tech_rep_shotgun_lbm_phylo_rpca <- read_tsv("matrix_tech_rep_shotgun_lbm_phylo_rpca.txt", na = c("", "NA", "not applicable"))
tech_rep_metab_jaccard <- read_tsv("matrix_tech_rep_metab_jaccard.txt", na = c("", "NA", "not applicable"))
tech_rep_metab_rpca <- read_tsv("matrix_tech_rep_metab_rpca.txt", na = c("", "NA", "not applicable"))
tech_rep_metab_cosine <- read_tsv("matrix_tech_rep_metab_cosine.txt", na = c("", "NA", "not applicable"))
tech_rep_metab_canberra_adkins <- read_tsv("matrix_tech_rep_metab_canberra_adkins.txt", na = c("", "NA", "not applicable"))


# Add dataset specific columns for merging
tech_rep_16s_hbm_jaccard$dataset = "16s"
tech_rep_16s_hbm_jaccard$distance_metric = "jaccard"
tech_rep_16s_hbm_unifrac$dataset = "16s"
tech_rep_16s_hbm_unifrac$distance_metric = "unifrac"
tech_rep_16s_hbm_wunifrac$dataset = "16s"
tech_rep_16s_hbm_wunifrac$distance_metric = "weighted_unifrac"
tech_rep_16s_hbm_rpca$dataset = "16s"
tech_rep_16s_hbm_rpca$distance_metric = "rpca"
tech_rep_16s_hbm_phylo_rpca$dataset = "16s"
tech_rep_16s_hbm_phylo_rpca$distance_metric = "phylo_rpca"
tech_rep_16s_lbm_jaccard$dataset = "16s"
tech_rep_16s_lbm_jaccard$distance_metric = "jaccard"
tech_rep_16s_lbm_unifrac$dataset = "16s"
tech_rep_16s_lbm_unifrac$distance_metric = "unifrac"
tech_rep_16s_lbm_wunifrac$dataset = "16s"
tech_rep_16s_lbm_wunifrac$distance_metric = "weighted_unifrac"
tech_rep_16s_lbm_rpca$dataset = "16s"
tech_rep_16s_lbm_rpca$distance_metric = "rpca"
tech_rep_16s_lbm_phylo_rpca$dataset = "16s"
tech_rep_16s_lbm_phylo_rpca$distance_metric = "phylo_rpca"
tech_rep_shotgun_hbm_jaccard$dataset = "shotgun"
tech_rep_shotgun_hbm_jaccard$distance_metric = "jaccard"
tech_rep_shotgun_hbm_unifrac$dataset = "shotgun"
tech_rep_shotgun_hbm_unifrac$distance_metric = "unifrac"
tech_rep_shotgun_hbm_wunifrac$dataset = "shotgun"
tech_rep_shotgun_hbm_wunifrac$distance_metric = "weighted_unifrac"
tech_rep_shotgun_hbm_rpca$dataset = "shotgun"
tech_rep_shotgun_hbm_rpca$distance_metric = "rpca"
tech_rep_shotgun_hbm_phylo_rpca$dataset = "shotgun"
tech_rep_shotgun_hbm_phylo_rpca$distance_metric = "phylo_rpca"
tech_rep_shotgun_lbm_jaccard$dataset = "shotgun"
tech_rep_shotgun_lbm_jaccard$distance_metric = "jaccard"
tech_rep_shotgun_lbm_unifrac$dataset = "shotgun"
tech_rep_shotgun_lbm_unifrac$distance_metric = "unifrac"
tech_rep_shotgun_lbm_wunifrac$dataset = "shotgun"
tech_rep_shotgun_lbm_wunifrac$distance_metric = "weighted_unifrac"
tech_rep_shotgun_lbm_rpca$dataset = "shotgun"
tech_rep_shotgun_lbm_rpca$distance_metric = "rpca"
tech_rep_shotgun_lbm_phylo_rpca$dataset = "shotgun"
tech_rep_shotgun_lbm_phylo_rpca$distance_metric = "phylo_rpca"
tech_rep_metab_jaccard$dataset = "metab"
tech_rep_metab_jaccard$distance_metric = "jaccard"
tech_rep_metab_rpca$dataset = "metab"
tech_rep_metab_rpca$distance_metric = "rpca"
tech_rep_metab_cosine$dataset = "metab"
tech_rep_metab_cosine$distance_metric = "cosine"
tech_rep_metab_canberra_adkins$dataset = "metab"
tech_rep_metab_canberra_adkins$distance_metric = "canberra_adkins"


# Concatenate files
tech_rep_master <- rbind(tech_rep_16s_hbm_jaccard,
                       tech_rep_16s_hbm_unifrac,
                       tech_rep_16s_hbm_wunifrac,
                       tech_rep_16s_hbm_rpca,
                       tech_rep_16s_hbm_phylo_rpca,
                       tech_rep_16s_lbm_jaccard,
                       tech_rep_16s_lbm_unifrac,
                       tech_rep_16s_lbm_wunifrac,
                       tech_rep_16s_lbm_rpca,
                       tech_rep_16s_lbm_phylo_rpca,
                       tech_rep_shotgun_hbm_jaccard,
                       tech_rep_shotgun_hbm_unifrac,
                       tech_rep_shotgun_hbm_wunifrac,
                       tech_rep_shotgun_hbm_rpca,
                       tech_rep_shotgun_hbm_phylo_rpca,
                       tech_rep_shotgun_lbm_jaccard,
                       tech_rep_shotgun_lbm_unifrac,
                       tech_rep_shotgun_lbm_wunifrac,
                       tech_rep_shotgun_lbm_rpca,
                       tech_rep_shotgun_lbm_phylo_rpca,
                       tech_rep_metab_jaccard,
                       tech_rep_metab_rpca,
                       tech_rep_metab_cosine,
                       tech_rep_metab_canberra_adkins)

# Re-order levels for sample_type2
tech_rep_master$sample_type2 <- factor(tech_rep_master$sample_type2, levels = c("human_feces", 
                                                                    "mouse_feces", 
                                                                    "human_saliva_before_brushing", 
                                                                    "human_saliva_after_brushing", 
                                                                    "skin_armpit", "skin_hand",  
                                                                    "surface_swab_floor", 
                                                                    "surface_swab_keyboard"))


# Re-name levels for sample_type2
levels(tech_rep_master$sample_type2) <- c("Human feces", 
                                    "Mouse feces", 
                                    "Human saliva, before", 
                                    "Human saliva, after", 
                                    "Skin, armpit", 
                                    "Skin, hand",  
                                    "Surface, floor", 
                                    "Surface, keyboard")


# Re-order levels for host_subject_id
tech_rep_master$host_subject_id <- factor(tech_rep_master$host_subject_id, levels = c("A", "B", "C", "D", 
                                                                          "mouse1", "mouse2", "mouse3", "mouse4",
                                                                          "bayd", "room1247",
                                                                          "pc_bayd", "pc_bayg"))


# Re-name levels for host_subject_id
levels(tech_rep_master$host_subject_id) <- c("Human subject A", "Human subject B", "Human subject C", "Human subject D",
                                       "Mouse subject 1", "Mouse subject 2", "Mouse subject 3", "Mouse subject 4",
                                       "Floor, Bay D", "Floor, Room 1247",
                                       "Keyboard, Bay D", "Keyboard, Bay G")


# Re-order levels for extraction_protocol_storage_solution
tech_rep_master$extraction_protocol_storage_solution <- factor(tech_rep_master$extraction_protocol_storage_solution, levels = c("EtOH-MagMax",
                                                                                                                    "IPA-MagMax",
                                                                                                                    "EtOH-Matrix",
                                                                                                                    "IPA-Matrix"))
# Re-name levels for extraction_protocol_storage_solution
levels(tech_rep_master$extraction_protocol_storage_solution) <- c("EtOH-Plate",
                                                                       "IPA-Plate",
                                                                       "EtOH-Matrix",
                                                                       "IPA-Matrix")

# Subset datasets for plotting
tech_rep_16s <- subset(tech_rep_master, tech_rep_master$dataset == "16s")
tech_rep_shotgun <- subset(tech_rep_master, tech_rep_master$dataset == "shotgun")
tech_rep_metab <- subset(tech_rep_master, tech_rep_master$dataset == "metab")



# Plot data and run stats
## 16S
### Jaccard
ggplot(subset(tech_rep_16s, 
              tech_rep_16s$distance_metric == "jaccard"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nJaccard distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.1, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_16s$sample_type2)) {
  data = subset(tech_rep_16s,
                tech_rep_16s$distance_metric == "jaccard" &
                  sample_type2 == i)
  result <- kruskal.test(value ~ extraction_protocol_storage_solution,
               data = data)
  print(i)
  print (nrow(data))
  print(result)
}


### RPCA
ggplot(subset(tech_rep_16s, 
              tech_rep_16s$distance_metric == "rpca"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nRPCA distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_16s$sample_type2)) {
  data = subset(tech_rep_16s,
                tech_rep_16s$distance_metric == "rpca" &
                  sample_type2 == i)
  result <- kruskal.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


### UniFrac
ggplot(subset(tech_rep_16s, 
              tech_rep_16s$distance_metric == "unifrac"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nUniFrac distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_16s$sample_type2)) {
  data = subset(tech_rep_16s,
                tech_rep_16s$distance_metric == "unifrac" &
                  sample_type2 == i)
  result <- kruskal.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


### Weighted UniFrac
ggplot(subset(tech_rep_16s, 
              tech_rep_16s$distance_metric == "weighted_unifrac"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nweighted UniFrac distance") +
  xlab("storage solution-extraction protocol") +
  #ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_16s$sample_type2)) {
  data = subset(tech_rep_16s,
                tech_rep_16s$distance_metric == "weighted_unifrac" &
                  sample_type2 == i)
  result <- kruskal.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


### Phylo-RPCA
ggplot(subset(tech_rep_16s, 
              tech_rep_16s$distance_metric == "phylo_rpca"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nphylo-RPCA distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_16s$sample_type2)) {
  data = subset(tech_rep_16s,
                tech_rep_16s$distance_metric == "phylo_rpca" &
                  sample_type2 == i)
  result <- kruskal.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


# Shotgun
### Jaccard
ggplot(subset(tech_rep_shotgun, 
              tech_rep_shotgun$distance_metric == "jaccard"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nJaccard distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_shotgun$sample_type2)) {
  data = subset(tech_rep_shotgun,
                tech_rep_shotgun$distance_metric == "jaccard" &
                  sample_type2 == i)
  result <- kruskal.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


### RPCA
ggplot(subset(tech_rep_shotgun, 
              tech_rep_shotgun$distance_metric == "rpca"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nRPCA distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_shotgun$sample_type2)) {
  data = subset(tech_rep_shotgun,
                tech_rep_shotgun$distance_metric == "rpca" &
                  sample_type2 == i)
  result <- kruskal.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


### UniFrac
ggplot(subset(tech_rep_shotgun, 
              tech_rep_shotgun$distance_metric == "unifrac"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nUniFrac distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_shotgun$sample_type2)) {
  data = subset(tech_rep_shotgun,
                tech_rep_shotgun$distance_metric == "unifrac" &
                  sample_type2 == i)
  result <- kruskal.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


### Weighted UniFrac
library(scales)

ggplot(subset(tech_rep_shotgun, 
              tech_rep_shotgun$distance_metric == "weighted_unifrac"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nweighted UniFrac distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 1.5)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4) + 
  scale_y_continuous(labels = label_number(accuracy = 0.01),
                     limits = c(0.0, 1.5))

for (i in levels(tech_rep_shotgun$sample_type2)) {
  data = subset(tech_rep_shotgun,
                tech_rep_shotgun$distance_metric == "weighted_unifrac" &
                  sample_type2 == i)
  result <- kruskal.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}

### Phylo-RPCA
ggplot(subset(tech_rep_shotgun, 
              tech_rep_shotgun$distance_metric == "phylo_rpca"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nphylo-RPCA distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_shotgun$sample_type2)) {
  data = subset(tech_rep_shotgun,
                tech_rep_shotgun$distance_metric == "phylo_rpca" &
                  sample_type2 == i)
  result <- kruskal.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


# Metabolomics

tech_rep_metab <- droplevels(tech_rep_metab)

### Jaccard
ggplot(subset(tech_rep_metab, 
              tech_rep_metab$distance_metric == "jaccard"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nJaccard distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_metab$sample_type2)) {
  data = subset(tech_rep_metab,
                tech_rep_metab$distance_metric == "jaccard" &
                  sample_type2 == i)
  result <- wilcox.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


### RPCA
ggplot(subset(tech_rep_metab, 
              tech_rep_metab$distance_metric == "rpca"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nRPCA distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 5.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4) + 
  scale_y_continuous(labels = label_number(accuracy = 0.01),
                     limits = c(0.0, 5.0))

for (i in levels(tech_rep_metab$sample_type2)) {
  data = subset(tech_rep_metab,
                tech_rep_metab$distance_metric == "rpca" &
                  sample_type2 == i)
  result <- wilcox.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


### Cosine
ggplot(subset(tech_rep_metab, 
              tech_rep_metab$distance_metric == "cosine"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\ncosine distance") +
  xlab("storage solution-extraction protocol") +
  #ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_metab$sample_type2)) {
  data = subset(tech_rep_metab,
                tech_rep_metab$distance_metric == "cosine" &
                  sample_type2 == i)
  result <- wilcox.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}


### Canberra-Adkins
ggplot(subset(tech_rep_metab, 
              tech_rep_metab$distance_metric == "canberra_adkins"), 
       aes(x = extraction_protocol_storage_solution, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 0.25) +
  facet_wrap(~sample_type2,
             scales = "free_y",
             ncol = 4) +
  ylab("within-sample\nCanberra-Adkins distance") +
  xlab("storage solution-extraction protocol") +
  ylim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.title=element_blank()) +
  stat_summary(fun = "mean",
               color = "red",
               size = 0.1,
               pch = 4)

for (i in levels(tech_rep_metab$sample_type2)) {
  data = subset(tech_rep_metab,
                tech_rep_metab$distance_metric == "canberra_adkins" &
                  sample_type2 == i)
  result <- wilcox.test(value ~ extraction_protocol_storage_solution,
                         data = data)
  print(i)
  print (nrow(data))
  print(result)
}
