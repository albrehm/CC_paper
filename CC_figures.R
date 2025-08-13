library(tidyverse)
library(readxl)
library(ggpubr)
library(ggh4x)
library(rstatix)
library(binom)
library(ggsignif)
library(patchwork)

#####Figure 3
######PBS dip prevalence and load 

# Read in PBS_Ooze_Dip.xlsx and convert Sample to a factor
# --- Set consistent Sample order right after reading the data ---
desired_order <- c("DGRP-517-", "DGRP-517+", "PBS Control")

ooze_dip <- ooze_dip %>%
  mutate(Sample = factor(Sample, levels = desired_order))


# Count total number of samples per group
df_group_sizes_dip <- ooze_dip %>%
  group_by(Sample) %>% 
  summarize(n = n(), .groupd = "drop")

# Count number of positives per group
df_positives_per_group_dip <- ooze_dip %>%
  filter(Positive == "Yes") %>%
  group_by(Sample) %>% 
  summarize(n_positive = n(), .groups = "drop")

# Merge total sample counts with positive counts
df_prevalence_by_group_dip <- left_join(df_group_sizes_dip, df_positives_per_group_dip)

# Replace NA positive counts with zero (for groups with no positives)
df_prevalence_by_group_dip$n_positive <- replace_na(df_prevalence_by_group_dip$n_positive, 0)

# Calculate binomial confidence intervals for prevalence
prev_conf_int_dip <- binom.confint(
  df_prevalence_by_group_dip$n_positive, 
  df_prevalence_by_group_dip$n, 
  method = "exact"
)

# Add percent positive and confidence interval bounds to dataframe
df_prevalence_by_group_dip$percent_positive <- prev_conf_int_dip$mean * 100
df_prevalence_by_group_dip$lower_conf_int   <- prev_conf_int_dip$lower * 100
df_prevalence_by_group_dip$upper_conf_int   <- prev_conf_int_dip$upper * 100

# Format percent positive values as whole numbers (for labeling)
df_prevalence_by_group_dip <- df_prevalence_by_group_dip %>% 
  mutate(percent_positive_text = sprintf("%.0f", percent_positive))

# Set y-axis limits for viral load plots
y_axis_min <- 1e-5
y_axis_max <- 1e+4

# For viral load: replace non-positives with a minimal value for plotting log scales
df1 <- ooze_dip %>% mutate(dCt = if_else(Positive != "Yes", y_axis_min, dCt))  

# Merge viral load data with prevalence data
df1 <- left_join(df1, df_prevalence_by_group_dip)

# Reset df1 to just prevalence data (likely for the prevalence plot)
df1 <- df_prevalence_by_group_dip

# Define comparison groups for statistical tests
my_comparisons <- list(
  c("DGRP-517-", "DGRP-517+"), 
  c("DGRP-517+", "PBS Control"), 
  c("DGRP-517-", "PBS Control")
)

# --- Prevalence plot ---
prevalence_dip <- ggplot(df1, aes(x = Sample, y = percent_positive, fill = Sample)) +
  geom_point(shape = 21, size = 3, stroke = 0.25) +  # Prevalence points
  geom_errorbar(aes(ymin = lower_conf_int, ymax = upper_conf_int, color = Sample), width = 0.1) +
  theme_bw(base_size = 13) +
  xlab("") + ylab("") +
  scale_fill_manual(values = c("DGRP-517-" = "lightblue", "DGRP-517+" = "orange", "PBS Control" = "red")) +
  scale_color_manual(values = c("DGRP-517-" = "lightblue", "DGRP-517+" = "orange", "PBS Control" = "red")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

prevalence_dip


# Define significance symbol mapping
my_symnum <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
  symbols = c("****", "***", "**", "*", "NS")  # Force uppercase
)

# --- Galbut virus Ct plot ---
galbut_ct <- ggplot(ooze_dip, aes(x = Sample, y = Gal, fill = Sample)) +
  geom_jitter(data = filter(ooze_dip, Positive == "Yes"),
              shape = 21, size = 3, stroke = 0.25, height = 0, width = 0.1) +
  geom_jitter(data = filter(ooze_dip, Positive != "Yes"),
              shape = 21, size = 3, fill = "grey90", alpha = 0.5,
              stroke = 0.25, height = 0, width = 0.1) +
  theme_bw(base_size = 13) +
  scale_fill_manual(values = c("DGRP-517-" = "lightblue",
                               "DGRP-517+" = "orange",
                               "PBS Control" = "red")) +
  xlab("") + ylab("Galbut virus Ct") +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text.x = element_text(),
        axis.text.y = element_text()) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     symnum.args = my_symnum,
                     paired = FALSE,
                     step.increase = .2,
                     method = "wilcox.test") +
  ylim(0, 50) +
  scale_x_discrete(drop = FALSE)

# --- RpL32 Ct plot ---
rpl_ct <- ggplot(ooze_dip, aes(x = Sample, y = RpL, fill = Sample)) +
  geom_jitter(data = filter(ooze_dip, RpL_Positive == "Yes"),
              shape = 21, size = 3, stroke = 0.25, height = 0, width = 0.1) +
  geom_jitter(data = filter(ooze_dip, RpL_Positive != "Yes"),
              shape = 21, size = 3, fill = "grey90", alpha = 0.5,
              stroke = 0.25, height = 0, width = 0.1) +
  theme_bw(base_size = 13) +
  scale_fill_manual(values = c("DGRP-517-" = "lightblue",
                               "DGRP-517+" = "orange",
                               "PBS Control" = "red")) +
  xlab("") + ylab("RpL32 Ct") +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text.x = element_text(),
        axis.text.y = element_text()) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     symnum.args = my_symnum,
                     paired = FALSE,
                     method = "wilcox.test") +
  ylim(0, 50) +
  scale_x_discrete(drop = FALSE)




# --- Add y-axis labels to leftmost plots ---
prevalence_dip <- prevalence_dip + ylab("Prevalence \n(% positive PBS wells)") + theme(axis.text.y = element_text())
galbut_ct <- galbut_ct + ylab("Galbut virus Ct") + theme(axis.text.y = element_text())
rpl_ct <- rpl_ct + ylab("RpL32 Ct") + theme(axis.text.y = element_text())

# --- Combine plots into one figure ---
full_plot_dip <- (prevalence_dip / galbut_ct / rpl_ct)+plot_annotation(tag_levels = "A")
full_plot_dip

# --- Save to PDF ---
plot_filename <- paste0("PBS_ooze_new.pdf")
ggsave(filename = plot_filename, plot = big_p1, width = 10, height = 12, units = "in")

  



### Figure 2
##### Ziploc bag experiment

# Load w1118 ooze experiment data
w1118 <- read_excel("./w1118_ooze.xlsx")

# --- Calculate sample sizes ---
df_group_sizes <- w1118 %>%
  group_by(Sample) %>%
  summarize(n = n(), .groupd = "drop")

# --- Count positive samples in each group ---
df_positives_per_group <- w1118 %>%
  filter(Positive == "Yes") %>%
  group_by(Sample) %>% 
  summarize(n_positive = n(), .groups = "drop")

# --- Merge total sample sizes and positive counts ---
df_prevalence_by_group <- left_join(df_group_sizes, df_positives_per_group)

# Replace NA positives (no positives found) with 0
df_prevalence_by_group$n_positive <- replace_na(df_prevalence_by_group$n_positive, 0)

# --- Calculate prevalence and 95% confidence intervals ---
prev_conf_int <- binom.confint(df_prevalence_by_group$n_positive, 
                               df_prevalence_by_group$n, 
                               method = "exact")

# Add prevalence and CI columns to data
df_prevalence_by_group$percent_positive <- prev_conf_int$mean * 100
df_prevalence_by_group$lower_conf_int   <- prev_conf_int$lower * 100
df_prevalence_by_group$upper_conf_int   <- prev_conf_int$upper * 100

# Store prevalence as formatted text for plotting labels
df_prevalence_by_group <- df_prevalence_by_group %>% 
  mutate(percent_positive_text = sprintf("%.0f", percent_positive))

# Define min/max for dCt axis scaling
y_axis_min <- 1e-5
y_axis_max <- 1e+4

# --- Replace negative samples' dCt with very small value for plotting ---
df <- w1118 %>% mutate(dCt = if_else(Positive != "Yes", y_axis_min, dCt))  

# Merge prevalence data into main dataframe
df <- left_join(df, df_prevalence_by_group)


### --- PREVALENCE PLOT ---
df_to_plot <- df_prevalence_by_group

# Define pairwise comparisons for significance testing
comp <- list(
  c("W1118 experimental", "W1118 uninfected control"),
  c("W1118 uninfected control", "Positive control"),
  c("W1118 experimental", "Positive control")
)

# Plot percent positive with error bars
percent_infected_bag <- ggplot(df_to_plot, aes(x = Sample, y = percent_positive, fill = Sample)) +
  geom_point(shape = 21, size = 3, stroke = 0.25) +
  geom_errorbar(aes(ymin = lower_conf_int, ymax = upper_conf_int, color = Sample), width = 0.1) +
  scale_fill_manual(values = c(
    "Positive control" = "lightblue",
    "W1118 experimental" = "orange",
    "W1118 uninfected control" = "red"
  )) +
  scale_color_manual(values = c(
    "Positive control" = "lightblue",
    "W1118 experimental" = "orange",
    "W1118 uninfected control" = "red"
  )) +
  theme_bw(base_size = 13) +
  ylim(c(0, 100)) +
  xlab("") +
  ylab("") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )

percent_infected_bag


### --- LOAD PLOT ---
ymin = y_axis_min
ymax = y_axis_max

df_to_plot <- df

# Plot viral load (dCt) with positives and negatives distinguished
viral_load_bag <- ggplot(df_to_plot, aes(x = Sample, y = dCt, fill = Sample)) +
  # Positive samples
  geom_jitter(data = filter(df_to_plot, Positive == "Yes"), 
              shape = 21, size = 3, stroke = 0.25, height = 0, width = 0.1) +
  # Negative samples (greyed out)
  geom_jitter(data = filter(df_to_plot, Positive != "Yes"), 
              aes(y = dCt), shape = 21, size = 3, fill = "grey90", alpha = 0.5, stroke = 0.25, height = 0, width = 0.1) +
  theme_bw(base_size = 13) +
  scale_fill_manual(values = c(
    "Positive control" = "lightblue",
    "W1118 experimental" = "orange",
    "W1118 uninfected control" = "red"
  )) +
  scale_y_log10() +
  stat_compare_means(comparisons = comp, label = "p.signif", method = "wilcox.test", paired = FALSE) +
  xlab("") +
  ylab("") +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    axis.text.x = element_text(),
    axis.text.y = element_blank()
  )

viral_load_bag


### --- MEDIAN CALCULATIONS ---
# Calculate median Ct values per sample type for Galbut virus and RpL32
tapply(w1118$Gal, w1118$Sample, median)
tapply(w1118$RpL, w1118$Sample, median)

# Ratio calculations (example)
(16.87895 / 12.89359)
10^1.309096
(29.62207 / 16.97860)
10^1.744671


### --- ADD LABELS TO PLOTS ---
percent_infected_bag  <- percent_infected_bag  + ylab("Galbut virus prevalence \n(% Positive Flies)") + theme(axis.text.y = element_text())
viral_load_bag <- viral_load_bag + ylab("Galbut virus RNA levels \nrelative to RpL32 mRNA") + theme(axis.text.y = element_text())

# Combine prevalence and load plots vertically
big_p <- percent_infected_bag / viral_load_bag+plot_annotation(tag_levels = "A")
big_p

# --- SAVE FIGURE ---
plot_filename1 <- paste0("Ziploc_ooze_new.pdf")
ggsave(filename = plot_filename1, plot = big_p, width = 10, height = 12, units = "in")
