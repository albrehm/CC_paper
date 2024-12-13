library(tidyverse)
library(readxl)
library(ggpubr)
install.packages("ggh4x")
library(ggh4x)

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


create_prev_plot1 <- function (a_Name = NA, 
                              a_Strain = NA, 
                              a_Sex = NA, 
                              a_Stage = NA, 
                              a_Transmission_mode = NA, 
                              a_Virus_target = NA,
                              facet_by = NA,
                              positive_color = "orange",
                              negative_color = "lightblue",
                              exp_color= "red") {
  
  df_to_plot1 <- df_prevalence_by_group1
  
  p1 <- ggplot(df_to_plot1) +
    geom_point(aes(x = Sample, y = percent_positive, fill = Sample), shape = 21, size = 3, stroke = 0.25) +
    geom_errorbar(aes(x = Sample, ymin = lower_conf_int, ymax = upper_conf_int, color = Sample), width = 0.1) +
    theme_bw(base_size = 13) +
    ylim(c(0, 100)) +
    xlab("") +
    ylab("") +
    scale_fill_manual(values = c("DGRP-517-" = negative_color, "DGRP-517+" = positive_color, "PBS Control" = exp_color))+
    scale_color_manual(values = c("DGRP-517-" = negative_color, "DGRP-517+" = positive_color, "PBS Control" = exp_color))+
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  if (!is.na(facet_by)) {
    p1 <- p1 + facet_wrap(vars(.data[[facet_by]]))
  }
  
  
  
  p1
}

prev1<-create_prev_plot1()

prev1

create_levels_plot1 <- function (a_Name = NA, 
                                a_Strain = NA, 
                                a_Sex = NA, 
                                a_Stage = NA, 
                                a_Transmission_mode = NA, 
                                a_Virus_target = NA,
                                facet_by = NA,
                                ymin = y_axis_min,
                                ymax = y_axis_max,
                                positive_color = "orange",
                                negative_color = "lightblue",
                                exp_color= "red") {
  
  df_to_plot1 <- df1
  
  
  p1 <- ggplot() +
    geom_jitter(data = filter(df_to_plot1, Positive == "Yes"), 
                aes(x = Sample, y = dCt, fill = Sample), shape = 21, size = 3, stroke = 0.25, height = 0, width = 0.1) +
    geom_jitter(data = filter(df_to_plot1, Positive != "Yes"), 
                aes(x = Sample, y = dCt), shape = 21, size = 3, fill = "grey90", alpha = 0.5, stroke = 0.25, height = 0, width = 0.1) +
    theme_bw(base_size = 13) +
    scale_y_log10(limits = c(ymin, ymax)) +
    scale_fill_manual(values = c("DGRP-517-" = negative_color, "DGRP-517+" = positive_color, "PBS Control" = exp_color))+
    xlab("") +
    ylab("") +
    
    theme(legend.position = "none",
          strip.text = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_blank())
  
  if (!is.na(facet_by)) {
    p1 <- p1 + facet_wrap(vars(.data[[facet_by]]))
  }
  
  p1
}


level1<-create_levels_plot1()
  
  
level1
  # add labels to left-most plots
  prev1 <- prev1 + ylab("Prevalence (% Positive Individuals)") + theme(axis.text.y = element_text())+
    ggtitle("Prevalence and viral load of flies dipped in PBS")
  level1 <- level1 + ylab("Viral RNA Levels in Infected Individuals\n(Viral RNA Relative to RpL32)") + theme(axis.text.y = element_text())
  
  # make a big plot with patchwork
  big_p1 <- (prev1 / level1)
  
  big_p1
  
  # save as PDF
  plot_filename <- paste0("PBS_ooze.pdf")
  ggsave(filename = plot_filename, plot = big_p1,  width=12, height=8, units="in")
  




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


create_prev_plot <- function (a_Name = NA, 
                              a_Strain = NA, 
                              a_Sex = NA, 
                              a_Stage = NA, 
                              a_Transmission_mode = NA, 
                              a_Virus_target = NA,
                              facet_by = NA,
                              positive_color = "orange",
                              negative_color = "lightblue",
                              exp_color= "red") {
  
  df_to_plot <- df_prevalence_by_group
  
  p <- ggplot(df_to_plot) +
    geom_point(aes(x = Sample, y = percent_positive, fill = Sample), shape = 21, size = 3, stroke = 0.25) +
    geom_errorbar(aes(x = Sample, ymin = lower_conf_int, ymax = upper_conf_int, color = Sample), width = 0.1) +
    scale_fill_manual(values = c("Positive control" = negative_color, "W1118 experimental" = positive_color, "W1118 uninfected control" = exp_color))+
    scale_color_manual(values = c("Positive ontrol" = negative_color, "W1118 experimental" = positive_color, "W1118 uninfected control" = exp_color))+
    theme_bw(base_size = 13) +
    ylim(c(0, 100)) +
    xlab("") +
    ylab("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  if (!is.na(facet_by)) {
    p <- p + facet_wrap(vars(.data[[facet_by]]))
  }
  
  
  
  p
}

prev<-create_prev_plot()



create_levels_plot <- function (a_Name = NA, 
                                a_Strain = NA, 
                                a_Sex = NA, 
                                a_Stage = NA, 
                                a_Transmission_mode = NA, 
                                a_Virus_target = NA,
                                facet_by = NA,
                                ymin = y_axis_min,
                                ymax = y_axis_max,
                                positive_color = "orange",
                                negative_color = "lightblue",
                                exp_color= "red"){ 
  
  df_to_plot <- df

  
  p <- ggplot() +
    geom_jitter(data = filter(df_to_plot, Positive == "Yes"), 
                aes(x = Sample, y = dCt, fill = Sample), shape = 21, size = 3, stroke = 0.25, height = 0, width = 0.1) +
    geom_jitter(data = filter(df_to_plot, Positive != "Yes"), 
                aes(x = Sample, y = dCt), shape = 21, size = 3, fill = "grey90", alpha = 0.5, stroke = 0.25, height = 0, width = 0.1) +
    theme_bw(base_size = 13) +
    scale_fill_manual(values = c("Positive control" = negative_color, "W1118 experimental" = positive_color, "W1118 uninfected control" = exp_color))+
    scale_y_log10(limits = c(ymin, ymax)) +
    xlab("") +
    ylab("") +
    theme(legend.position = "none",
          strip.text = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_blank())
  
  if (!is.na(facet_by)) {
    p <- p + facet_wrap(vars(.data[[facet_by]]))
  }
  
  p
}


 level<-create_levels_plot()

 

   
 
   # add labels to left-most plots
   prev <- prev + ylab("Prevalence (% Positive Individuals)") + theme(axis.text.y = element_text())+
     ggtitle("Prevalence and viral load of flies co-frozen")
   level <- level + ylab("Viral RNA Levels in Infected Individuals\n(Viral RNA Relative to RpL32)") + theme(axis.text.y = element_text())
   
   # make a big plot with patchwork
   big_p <- (prev / level)
   
   # save as PDF
   plot_filename1 <- paste0("Ziploc_ooze.pdf")
   ggsave(filename = plot_filename1, plot = big_p,  width=12, height=8, units="in")
   
   big_p

 