library(tidyverse)
library(readxl)
library(ggpubr)
library(ggh4x)
library(rstatix)
library(binom)
library(ggsignif)
library(patchwork)
#####Figure 2
######PBS dip prevalence and load 

ooze_dip<- read_excel("./PBS_Ooze_Dip.xlsx") %>%
  mutate(Sample = factor(Sample))

df_group_sizes1<- ooze_dip%>%
  group_by(Sample) %>% summarize(n =n(), .groupd="drop")

df_positives_per_group1 <- ooze_dip %>%
  filter(Positive == "Yes") %>%
  group_by(Sample) %>% 
  summarize(n_positive = n(), .groups="drop")

df_prevalence_by_group1 <- left_join(df_group_sizes1, df_positives_per_group1)

df_prevalence_by_group1$n_positive <- replace_na(df_prevalence_by_group1$n_positive, 0)

prev_conf_int1 <- binom.confint(df_prevalence_by_group1$n_positive, 
                               df_prevalence_by_group1$n, 
                               method="exact")

df_prevalence_by_group1$percent_positive  <- prev_conf_int1$mean * 100
df_prevalence_by_group1$lower_conf_int    <- prev_conf_int1$lower * 100
df_prevalence_by_group1$upper_conf_int    <- prev_conf_int1$upper * 100

df_prevalence_by_group1 <- df_prevalence_by_group1 %>% 
  mutate(percent_positive_text = sprintf("%.0f", percent_positive))

y_axis_min <- 1e-5
y_axis_max <- 1e+4

df1 <- ooze_dip %>% mutate(dCt = if_else(Positive != "Yes", y_axis_min, dCt))  

df1 <- left_join(df1, df_prevalence_by_group1)


df1 <- df_prevalence_by_group1
my_comparisons <- list( c("DGRP-517-", "DGRP-517+"), c("DGRP-517+", "PBS Control"), c("DGRP-517-", "PBS Control"))
  
p5 <- ggplot(df1, aes(x = Sample, y = percent_positive, fill = Sample)) +
    geom_point(shape = 21, size = 3, stroke = 0.25) +
    geom_errorbar(aes(x = Sample, ymin = lower_conf_int, ymax = upper_conf_int, color = Sample), width = 0.1) +
    theme_bw(base_size = 13) +
    xlab("") +
    ylab("") +
    scale_fill_manual(values = c("DGRP-517-" = "lightblue", "DGRP-517+" = "orange", "PBS Control" = "red"))+
    scale_color_manual(values = c("DGRP-517-" = "lightblue", "DGRP-517+" = "orange", "PBS Control" = "red"))+
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank())+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, step_increase = .12)+ylim(0,140)
  
  
  
  p5




p1 <- ggplot(df1, aes(x = Sample, y = Gal, fill = Sample)) +
  geom_jitter(data = filter(df1, Positive == "Yes"), aes(x=Sample, y =Gal),
              shape = 21, size = 3, stroke = 0.25, height = 0, width = 0.1) +
  geom_jitter(data = filter(df1, Positive != "Yes"),
              aes(x = Sample, y = Gal), shape = 21, size = 3, fill = "grey90", alpha = 0.5, stroke = 0.25, height = 0, width = 0.1) +
  theme_bw(base_size = 13) +
  scale_fill_manual(values = c("DGRP-517-" = "lightblue", "DGRP-517+" = "orange", "PBS Control" = "red"))+
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text.x = element_text(),
        axis.text.y = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label="p.signif", hide.ns = TRUE)+ylim(0,45)



p1


  
  my_comparisons <- list( c("DGRP-517-", "DGRP-517+"), c("DGRP-517+", "PBS Control"), c("DGRP-517-", "PBS Control"))
  
  p2 <- ggplot(df1, aes(x = Sample, y = RpL, fill = Sample)) +
    geom_jitter(data = filter(df1, RpL_Positive == "Yes"), aes(x=Sample, y =RpL),
                shape = 21, size = 3, stroke = 0.25, height = 0, width = 0.1) +
    geom_jitter(data = filter(df1, RpL_Positive != "Yes"),
                aes(x = Sample, y = RpL), shape = 21, size = 3, fill = "grey90", alpha = 0.5, stroke = 0.25, height = 0, width = 0.1) +
    theme_bw(base_size = 13) +
    scale_fill_manual(values = c("DGRP-517-" = "lightblue", "DGRP-517+" = "orange", "PBS Control" = "red"))+
    xlab("") +
    ylab("") +
    theme(legend.position = "none",
          strip.text = element_blank(),
          axis.text.x = element_text(),
          axis.text.y = element_blank())+
    stat_compare_means(comparisons = my_comparisons, label="p.signif")+ ylim(0,45)



p2


  # add labels to left-most plots
  p5 <- p5 + ylab("Prevalence \n(% Positive Individuals)") + theme(axis.text.y = element_text())
  p1 <- p1 + ylab("Galbut virus Ct") + theme(axis.text.y = element_text())
  p2 <- p2 + ylab("RpL32 Ct") + theme(axis.text.y = element_text())
  
  # make a big plot with patchwork
  big_p1 <- (p5/ p1/p2)
  
  big_p1
  
  # save as PDF
  plot_filename <- paste0("PBS_ooze.pdf")
  ggsave(filename = plot_filename, plot = big_p1,  width=10, height=12, units="in")
  




###Figure 3
#####Ziploc bag experiment
w1118<-read_excel("./w1118_ooze.xlsx")

df_group_sizes<- w1118%>%
group_by(Sample) %>% summarize(n =n(), .groupd="drop")

df_positives_per_group <- w1118 %>%
  filter(Positive == "Yes") %>%
  group_by(Sample) %>% 
  summarize(n_positive = n(), .groups="drop")

df_prevalence_by_group <- left_join(df_group_sizes, df_positives_per_group)

df_prevalence_by_group$n_positive <- replace_na(df_prevalence_by_group$n_positive, 0)

prev_conf_int <- binom.confint(df_prevalence_by_group$n_positive, 
                               df_prevalence_by_group$n, 
                               method="exact")
df_prevalence_by_group$percent_positive  <- prev_conf_int$mean * 100
df_prevalence_by_group$lower_conf_int    <- prev_conf_int$lower * 100
df_prevalence_by_group$upper_conf_int    <- prev_conf_int$upper * 100

df_prevalence_by_group <- df_prevalence_by_group %>% 
  mutate(percent_positive_text = sprintf("%.0f", percent_positive))

y_axis_min <- 1e-5
y_axis_max <- 1e+4

df <- w1118 %>% mutate(dCt = if_else(Positive != "Yes", y_axis_min, dCt))  

df <- left_join(df, df_prevalence_by_group)




###PREVALENCE PLOT
  df_to_plot <- df_prevalence_by_group
  comp<- list(c("W1118 experimental", "W1118 uninfected control"), c("W1118 uninfected control", "Positive control"), c("W1118 experimental", "Positive control"))
  
  
  
  p <- ggplot(df_to_plot, (aes(x = Sample, y = percent_positive, fill = Sample))) +
    geom_point(shape = 21, size = 3, stroke = 0.25) +
    geom_errorbar(aes(x = Sample, ymin = lower_conf_int, ymax = upper_conf_int, color = Sample), width = 0.1) +
    scale_fill_manual(values = c("Positive control" = "lightblue", "W1118 experimental" = "orange", "W1118 uninfected control" = "red"))+
    scale_color_manual(values = c("Positive control" = "lightblue", "W1118 experimental" = "orange", "W1118 uninfected control" = "red"))+
    theme_bw(base_size = 13) +
    ylim(c(0, 100)) +
    xlab("") +
    ylab("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank())+
    geom_signif(comparisons = comp, map_signif_level = TRUE, step_increase = .12)+ ylim(0,130)
    
    
p


ymin = y_axis_min
ymax = y_axis_max
  
###LOAD PLOT
  df_to_plot <- df

  
  p4 <- ggplot(df_to_plot, aes(x = Sample, y = dCt, fill = Sample)) +
    geom_jitter(data = filter(df_to_plot, Positive == "Yes"), 
                 shape = 21, size = 3, stroke = 0.25, height = 0, width = 0.1) +
    geom_jitter(data = filter(df_to_plot, Positive != "Yes"), 
                aes(x = Sample, y = dCt), shape = 21, size = 3, fill = "grey90", alpha = 0.5, stroke = 0.25, height = 0, width = 0.1) +
    theme_bw(base_size = 13) +
    scale_fill_manual(values = c("Positive control" = "lightblue", "W1118 experimental" = "orange", "W1118 uninfected control" = "red"))+
    scale_y_log10()+
    stat_compare_means(comparisons = comp, label = "p.signif")+
    xlab("") +
    ylab("") +
    theme(legend.position = "none",
          strip.text = element_blank(),
          axis.text.x = element_text(),
          axis.text.y = element_blank())

  
 p4

#####Calculate the median of each of my Ct values per sample type and then determine amount greater or less than for each group
 tapply(w$Gal, w$Sample, median)
 
tapply(w$RpL, w$Sample, median)
 
(16.87895/12.89359)

10^1.309096

(29.62207/16.97860)

10^1.744671

   # add labels to left-most plots
   p <- p + ylab("Prevalence \n(% Positive Individuals)") + theme(axis.text.y = element_text())
   p4 <- p4 + ylab("RNA Levels \n(Viral RNA Relative to RpL32)") + theme(axis.text.y = element_text())
   
   # make a big plot with patchwork
   big_p <- p / p4
   
   big_p
   # save a
   
   plot_filename1 <- paste0("Ziploc_ooze.pdf")
   ggsave(filename = plot_filename1, plot = big_p,  width=10, height=12, units="in")
   
   big_p

 