#Analysis for experimental data in Chow et al (2024) 

library(tidyverse)
library(ggplot2)
library(emmeans)
library(afex)
afex_options(emmeans_model = "multivariate")
# library(KSgeneral)
library(Matching)
library(dplyr)
library(here)
library(rstatix)
library(knitr)

file_name_root <- paste0("fig/")
file_name_root_data <- paste0("output/")

graph_file_type <- ".tif"
gg_width <- 10
gg_height <- 8
dpi <- 600

#functions ----
Descriptives <- function(.data, .dv, ...) {
  out <- .data %>%
    group_by(...) %>%
    dplyr::summarise(mean = mean({{ .dv }}),
                     n = n(),
                     sd = sd({{ .dv }}),
                     se = sd({{ .dv }})/sqrt(n))
}

Partial_Eta2 <- function(F, df) {
  t <- sqrt(F)
  pes <- t^2/(t^2 + df)
  return(pes) 
}

Contrast_Table <- function(out) { 
  out <- as.data.frame(out) %>%
    mutate(F.ratio = t.ratio^2, 
           # lower.CI = estimate - 1.96*SE, 
           # upper.CI = estimate + 1.96*SE,
           pes = Partial_Eta2(F.ratio, df)
    )
  print(kable(out))
}

################################
experiment <- 1
Exp1_data <- read.csv(here("original_data/Exp1_compiled_data.csv"))

#check their rating of the last 10 fixed cues-----
Exp1_fixedcheck <- Exp1_data%>%
  filter(cueDescription%in%c("fixed0","fixed100"))%>%
  filter(phase == "train" & trialNumPerCue > 11) #21 presentations each

Exp1_process <- Descriptives(Exp1_fixedcheck,response,participant,cueDescription)
Exp1_process <- as.data.frame(Exp1_process)
Exp1_process <- Exp1_process%>%
  mutate(fixed_check = case_when(
    cueDescription == "fixed0" & mean > 10 ~ "fail",
    cueDescription == "fixed100" & mean < 90 ~ "fail",
    TRUE ~ ""
  ))
extract_subj <- Exp1_process%>%
  filter(fixed_check == "fail")
extract_subj <- unique(extract_subj$participant)

Exp1_data <- Exp1_data %>%
  mutate(fixedcheck = ifelse(
    participant%in%extract_subj, "fail","pass"
  ))
n_distinct(Exp1_data$participant)

write_csv(Exp1_data, paste0(file_name_root_data, "Exp1_ProcessedData.csv"))

#update Exp1 to exclude participants who failed manipulation check
Exp1_data <- Exp1_data %>%
  filter(!participant%in%extract_subj)
n_distinct(Exp1_data$participant)

#check if cueDescription matches outcome
Exp1_description_check <- Exp1_data%>%
  filter(!cueDescription%in%c("fixed0","fixed100"))

Exp1_descrip_check <- Descriptives(Exp1_description_check,outcome,cueDescription)
kable(Exp1_descrip_check)

Exp1_data_target <- Exp1_description_check%>%
  mutate(cueName = case_when(
    cueDescription == "fixed25" ~ "Fixed Value\nMean = 25",
    cueDescription == "fixed50" ~ "Fixed Value\nMean = 50",
    cueDescription == "fixed75" ~ "Fixed Value\nMean = 75",
    cueDescription == "norm25" ~ "Normal Distribution\nMean = 25",
    cueDescription == "norm50" ~ "Normal Distribution\nMean = 50",
    cueDescription == "norm75" ~ "Normal Distribution\nMean = 75",
    cueDescription == "rect25" ~ "Uniform Distribution\nMean = 25",
    cueDescription == "rect50" ~ "Uniform Distribution\nMean = 50",
    cueDescription == "rect75" ~ "Uniform Distribution\nMean = 75"
  ))

Exp1_description_plot <- ggplot(data = Exp1_data_target, mapping = aes(x=outcome))+facet_wrap(~cueName)+
  geom_histogram(data = Exp1_data_target, mapping = aes(x=outcome,after_stat(count)), binwidth = 1, fill = "pink", colour = "red")+
  labs(x = "Outcome Value", y = "N trials")+
  theme_classic()+
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text =  element_text(size = 13, colour = "black"), axis.title = element_text(size = 15),
        strip.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10,r=0,b= 0,l=0)),plot.title = element_text(size = 20),
        axis.title.y = element_text(margin = margin(t = 0,r=10,b= 0,l=0)),
        legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,1,0,0.5), "cm"), #add this to avoid the plot margin going all the way to the right and cutting off the scale
        panel.spacing = unit(1.1, "lines"))+
  scale_x_continuous()
  # theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+

ggsave(paste0(file_name_root, "Exp1_distribution_plot",graph_file_type),Exp1_description_plot,"jpg",
       height = gg_width*1.8, width = gg_width*2, units="cm",dpi = dpi)
include_graphics(paste0(file_name_root,"Exp1_distribution_plot",graph_file_type))

#analyse
Exp1_train <- Exp1_data_target %>%
  filter(trialNumPerCue > 11)

#target 25
fixed_25 <- Exp1_train %>%
  filter(cueDescription == "fixed25")

variable_25 <- Exp1_train%>%
  filter(cueDescription%in%c("norm25","rect25"))

normal_25 <- Exp1_train %>%
  filter(cueDescription == "norm25")

rect_25 <- Exp1_train %>%
  filter(cueDescription == "rect25")

Exp1_fixedvsvariable_25 <- ks.boot(fixed_25$response, variable_25$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp1_normalvsUniform_25 <- ks.boot(normal_25$response, rect_25$response, nboots=1000, alternative = "two.sided", print.level=0)


#target 50
fixed_50 <- Exp1_train %>%
  filter(cueDescription == "fixed50")

variable_50 <- Exp1_train%>%
  filter(cueDescription%in%c("norm50","rect50"))

normal_50 <- Exp1_train %>%
  filter(cueDescription == "norm50")

rect_50 <- Exp1_train %>%
  filter(cueDescription == "rect50")

Exp1_fixedvsvariable_50 <- ks.boot(fixed_50$response, variable_50$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp1_normalvsUniform_50 <- ks.boot(normal_50$response, rect_50$response, nboots=1000, alternative = "two.sided", print.level=0)

#mean 75
fixed_75 <- Exp1_train %>%
  filter(cueDescription == "fixed75")

variable_75 <- Exp1_train%>%
  filter(cueDescription%in%c("norm75","rect75"))

normal_75 <- Exp1_train %>%
  filter(cueDescription == "norm75")

rect_75 <- Exp1_train %>%
  filter(cueDescription == "rect75")

Exp1_fixedvsvariable_75 <- ks.boot(fixed_75$response, variable_75$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp1_normalvsUniform_75 <- ks.boot(normal_75$response, rect_75$response, nboots=1000, alternative = "two.sided", print.level=0)


#compare to observed outcomes
Exp1_fixed25_observed <- ks.boot(fixed_25$response, fixed_25$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp1_fixed50_observed <- ks.boot(fixed_50$response, fixed_50$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp1_fixed75_observed <- ks.boot(fixed_75$response, fixed_75$outcome, nboots=1000, alternative = "two.sided", print.level=0)

Exp1_rect25_observed <- ks.boot(rect_25$response, rect_25$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp1_rect50_observed <- ks.boot(rect_50$response, rect_50$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp1_rect75_observed <- ks.boot(rect_75$response, rect_75$outcome, nboots=1000, alternative = "two.sided", print.level=0)

Exp1_normal25_observed <- ks.boot(normal_25$response, normal_25$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp1_normal50_observed <- ks.boot(normal_50$response, normal_50$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp1_normal75_observed <- ks.boot(normal_75$response, normal_75$outcome, nboots=1000, alternative = "two.sided", print.level=0)

################################
experiment <- 2
Exp2_data <- read.csv(here("original_data/Exp2_compiled_data.csv"))

#check their rating of the last 10 fixed cues-----
Exp2_fixedcheck <- Exp2_data%>%
  filter(cueDescription%in%c("fixed0","fixed100"))%>%
  filter(phase == "train" & trialNumPerCue > 16)
n_distinct(Exp2_data$participant)

Exp2_process <- Descriptives(Exp2_fixedcheck,response,participant,cueDescription)
Exp2_process <- as.data.frame(Exp2_process)
Exp2_process <- Exp2_process%>%
  mutate(fixed_check = case_when(
    cueDescription == "fixed0" & mean > 10 ~ "fail",
    cueDescription == "fixed100" & mean < 90 ~ "fail",
    TRUE ~ ""
  ))
extract_subj <- Exp2_process%>%
  filter(fixed_check == "fail")
extract_subj <- unique(extract_subj$participant)

Exp2_data <- Exp2_data %>%
  mutate(fixedcheck = ifelse(
    participant%in%extract_subj, "fail","pass"
  ))
n_distinct(Exp2_data$participant)

write_csv(Exp2_data, paste0(file_name_root_data, "Exp2_ProcessedData.csv"))


#update Exp2 to exclude participants who failed manipulation check
Exp2_data <- Exp2_data %>%
  filter(!participant%in%extract_subj)
n_distinct(Exp2_data$participant)

Exp2_groupN <- Exp2_data%>%
  group_by(participant)%>%
  slice(1)%>%
  group_by(group)%>%
  dplyr::summarise(count=n())

#check if cueDescription matches outcome
Exp2_description_check <- Exp2_data%>%
  filter(cueDescription%in%c("target15","target30","target50","target70","target85"))

Exp2_description_check <- Descriptives(Exp2_description_check,outcome,cueDescription)
kable(Exp2_description_check)

Exp2_data_target <- Exp2_data %>%
  filter(cueDescription%in%c("target15","target30","target50","target70","target85"))

Exp2_train <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 15)

#target cue, outcome mean = 15
target15_group30 <- Exp2_train %>%
  filter(group=="30" & cueDescription == "target15")

target15_group50 <- Exp2_train %>%
  filter(group=="50" & cueDescription == "target15")

target15_group70 <- Exp2_train %>%
  filter(group=="70" & cueDescription == "target15")

Exp2_30vs50_target15 <- ks.boot(target15_group30$response, target15_group50$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_30vs70_target15 <- ks.boot(target15_group30$response, target15_group70$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_50vs70_target15 <- ks.boot(target15_group50$response, target15_group70$response, nboots=1000, alternative = "two.sided", print.level=0)

#target cue, outcome mean = 30
target30_group30 <- Exp2_train %>%
  filter(group=="30" & cueDescription == "target30")

target30_group50 <- Exp2_train %>%
  filter(group=="50" & cueDescription == "target30")

target30_group70 <- Exp2_train %>%
  filter(group=="70" & cueDescription == "target30")

Exp2_30vs50_target30 <- ks.boot(target30_group30$response, target30_group50$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_30vs70_target30 <- ks.boot(target30_group30$response, target30_group70$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_50vs70_target30 <- ks.boot(target30_group50$response, target30_group70$response, nboots=1000, alternative = "two.sided", print.level=0)

#target cue, outcome mean = 50
target50_group30 <- Exp2_train %>%
  filter(group=="30" & cueDescription == "target50")

target50_group50 <- Exp2_train %>%
  filter(group=="50" & cueDescription == "target50")

target50_group70 <- Exp2_train %>%
  filter(group=="70" & cueDescription == "target50")

Exp2_30vs50_target50 <- ks.boot(target50_group30$response, target50_group50$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_30vs70_target50 <- ks.boot(target50_group30$response, target50_group70$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_50vs70_target50 <- ks.boot(target50_group50$response, target50_group70$response, nboots=1000, alternative = "two.sided", print.level=0)

#target cue, outcome mean = 70
target70_group30 <- Exp2_train %>%
  filter(group=="30" & cueDescription == "target70")

target70_group50 <- Exp2_train %>%
  filter(group=="50" & cueDescription == "target70")

target70_group70 <- Exp2_train %>%
  filter(group=="70" & cueDescription == "target70")

Exp2_30vs50_target70 <- ks.boot(target70_group30$response, target70_group50$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_30vs70_target70 <- ks.boot(target70_group30$response, target70_group70$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_50vs70_target70 <- ks.boot(target70_group50$response, target70_group70$response, nboots=1000, alternative = "two.sided", print.level=0)

#target cue, outcome mean = 85
target85_group30 <- Exp2_train %>%
  filter(group=="30" & cueDescription == "target85")

target85_group50 <- Exp2_train %>%
  filter(group=="50" & cueDescription == "target85")

target85_group70 <- Exp2_train %>%
  filter(group=="70" & cueDescription == "target85")

Exp2_30vs50_target85 <- ks.boot(target85_group30$response, target85_group50$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_30vs70_target85 <- ks.boot(target85_group30$response, target85_group70$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_50vs70_target85 <- ks.boot(target85_group50$response, target85_group70$response, nboots=1000, alternative = "two.sided", print.level=0)

#compare to observed outcomes
target15_group30 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target15"& group == "30")


target30_group30 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target30" & group == "30")

target50_group30 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target50" & group == "30")

target70_group30 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target70" & group == "30")

target85_group30 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target85" & group == "30")

Exp2_30target15_observed <- ks.boot(target15_group30$response, target15_group30$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_30target30_observed <- ks.boot(target30_group30$response, target30_group30$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_30target50_observed <- ks.boot(target50_group30$response, target50_group30$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_30target70_observed <- ks.boot(target70_group30$response, target70_group30$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_30target85_observed <- ks.boot(target85_group30$response, target85_group30$outcome, nboots=1000, alternative = "two.sided", print.level=0)

target15_group50 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target15"& group == "50")

target30_group50 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target30" & group == "50")

target50_group50 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target50" & group == "50")

target70_group50 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target70" & group == "50")

target85_group50 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target85" & group == "50")

Exp2_50target15_observed <- ks.boot(target15_group50$response, target15_group50$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_50target30_observed <- ks.boot(target30_group50$response, target30_group50$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_50target50_observed <- ks.boot(target50_group50$response, target50_group50$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_50target70_observed <- ks.boot(target70_group50$response, target70_group50$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_50target85_observed <- ks.boot(target85_group50$response, target85_group50$outcome, nboots=1000, alternative = "two.sided", print.level=0)

target15_group70 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target15"& group == "70")

target30_group70 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target30" & group == "70")

target50_group70 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target50" & group == "70")

target70_group70 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target70" & group == "70")

target85_group70 <- Exp2_data_target %>%
  filter(phase == "train" & trialNumPerCue > 12 & cueDescription == "target85" & group == "70")

Exp2_70target15_observed <- ks.boot(target15_group70$response, target15_group70$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_70target30_observed <- ks.boot(target30_group70$response, target30_group70$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_70target50_observed <- ks.boot(target50_group70$response, target50_group70$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_70target70_observed <- ks.boot(target70_group70$response, target70_group70$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp2_70target85_observed <- ks.boot(target85_group70$response, target85_group70$outcome, nboots=1000, alternative = "two.sided", print.level=0)

#test data
#analyse test ratings 
Exp2_mean <- Exp2_data_target %>%
  filter(phase == "meantest")%>%
  mutate(cueMean = case_when(
    cueDescription == "target15" ~ "15",
    cueDescription == "target30" ~ "30",
    cueDescription == "target50" ~ "50",
    cueDescription == "target70" ~ "70",
    cueDescription == "target85" ~ "85"
  ),
  groupName = case_when(
    group == "30" ~ "Base rate = 30",
    group == "50" ~ "Base rate = 50",
    group == "70" ~ "Base rate = 70"
  ))

#plot
# cueDescrip_labels <- c("Normal","Positive Skew","Negative Skew","Normal","Positive Skew","Negative Skew")
fig_cols <- c("blue","orange","red","darkgreen","black")

rainplot_mean_Exp2 <- ggplot(Exp2_mean, aes(x = cueMean, y = response))+facet_grid(~groupName)+
  geom_boxplot(
    aes(color=cueMean),
    width = .3,
    outlier.shape = NA
  )+
  geom_point(
    aes(color = cueMean),
    fill = "white",
    shape = 21,
    stroke = .4,
    size = 2,
    position = position_jitter(seed = 1, width = .12)
  )+ylim(0,100)+
  # geom_hline(yintercept=0, linetype = "dashed")+
  stat_summary(fun.data=mean_se, geom = "errorbar",width = 0.2, position = position_dodge(width = 0.2))+
  stat_summary(fun.data=mean_se, geom = "point",size = 3, shape = 21, fill = "white", position = position_dodge(width = 0.2))+
  scale_color_manual(values = fig_cols)+
  scale_fill_manual(values = fig_cols)+
  # scale_x_discrete(labels = cueDescrip_labels)+
  labs(x = "Outcome Mean", y = "Estimate")+
  theme_classic()+
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text =  element_text(size = 14, colour = "black"), axis.title = element_text(size = 20),
        strip.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10,r=0,b= 0,l=0)),plot.title = element_text(size = 20),
        legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave(paste0(file_name_root, "Exp2_test_mean",graph_file_type),rainplot_mean_Exp2,"jpg",
       height = gg_width*1.8, width = gg_width*2, units="cm",dpi = dpi)
include_graphics(paste0(file_name_root,"Exp2_test_mean",graph_file_type))

#analyse mean by group for each cue

Exp2.mean.aov <- anova_test(data=Exp2_mean, dv = response, wid= participant, between = group,within = cueMean,effect.size = "pes")

Exp2_meanEstimates <- Descriptives(Exp2_mean,response,cueMean)


################################
experiment <- 3
Exp3_data <- read.csv(here("original_data/Exp3_compiled_data.csv"))
n_distinct(Exp3_data$participant)

#check their rating of the last 10 fixed cues-----
Exp3_fixedcheck <- Exp3_data%>%
  filter(cueDescription%in%c("fixed0","fixed100"))%>%
  filter(phase == "train" & trialNumPerCue > 10)

Exp3_process <- Descriptives(Exp3_fixedcheck,response,participant,cueDescription)
Exp3_process <- as.data.frame(Exp3_process)
Exp3_process <- Exp3_process%>%
  mutate(fixed_check = case_when(
    cueDescription == "fixed0" & mean > 10 ~ "fail",
    cueDescription == "fixed100" & mean < 90 ~ "fail",
    TRUE ~ ""
  ))
extract_subj <- Exp3_process%>%
  filter(fixed_check == "fail")
extract_subj <- unique(extract_subj$participant)

Exp3_data <- Exp3_data %>%
  mutate(fixedcheck = ifelse(
    participant%in%extract_subj, "fail","pass"
  ))

# write_csv(Exp3_data, paste0(file_name_root_data, "Exp3_ProcessedData.csv"))

#update Exp3 to exclude participants who failed manipulation check
Exp3_data <- Exp3_data %>%
  filter(!participant%in%extract_subj)
n_distinct(Exp3_data$participant)

#check if cueDescription matches outcome
Exp3_description_check <- Exp3_data%>%
  filter(!cueDescription%in%c("fixed0","fixed100"))

Exp3_description_check <- Descriptives(Exp3_description_check,outcome,cueDescription)
kable(Exp3_description_check)

Exp3_description_plot <- Exp3_data%>%
  filter(!cueDescription%in%c("fixed0","fixed100"))

Exp3_description_plot <- ggplot(data = Exp3_description_plot, mapping = aes(x=outcome))+facet_wrap(~cueDescription)+
  # geom_histogram(data = Exp3_description_plot, mapping = aes(x=response, y = ..density..), binwidth = 1, fill = "pink", colour = "red")+
  geom_freqpoly(data = Exp3_description_plot, colour = "red",size = 1, stat ="density",linetype = "dashed")+
  geom_vline(xintercept=35)+geom_vline(xintercept=65)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(limits = c(0,100))

Exp3_data <- Exp3_data%>%
  mutate(
  distribution = case_when(
    cueDescription == "negskew35" | cueDescription == "negskew65" ~ "negative-skew",
    cueDescription == "posskew35" | cueDescription == "posskew65" ~ "positive-skew",
    cueDescription == "norm35" | cueDescription == "norm65" ~ "normal",
    cueDescription == "fixed0" | cueDescription == "fixed100" ~ "fixed"
  ),
  mean = case_when(
    cueDescription == "negskew35" | cueDescription == "posskew35"| cueDescription == "norm35" ~ "35",
    cueDescription == "negskew65" | cueDescription == "posskew65"| cueDescription == "norm65" ~ "65",
    cueDescription == "fixed0"  ~ "0",
    cueDescription == "fixed100"  ~ "100"
  ))
  
Exp3_data_target <- Exp3_data %>%
  filter(distribution!="fixed")

#ks test on last 20 training predictions
Exp3_train <- Exp3_data_target %>%
  filter(phase == "train" & trialNumPerCue > 20)

#ksTest
test_posskew35 <- Exp3_train %>%
  filter(mean == "35" & distribution == "positive-skew")

test_normal35 <- Exp3_train %>%
  filter(mean == "35" & distribution == "normal")


test_negskew35 <- Exp3_train %>%
  filter(mean == "35" & distribution == "negative-skew")

Exp3_NormalvsPosskew_35 <- ks.boot(test_posskew35$response, test_normal35$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp3_NormalvsNegskew_35 <- ks.boot(test_negskew35$response, test_normal35$response, nboots=1000, alternative = "two.sided", print.level=0)

#65
test_posskew65 <- Exp3_train %>%
  filter(mean == "65" & distribution == "positive-skew")

test_normal65 <- Exp3_train %>%
  filter(mean == "65" & distribution == "normal")

test_negskew65 <- Exp3_train %>%
  filter(mean == "65" & distribution == "negative-skew")

Exp3_NormalvsPosskew_65 <- ks.boot(test_posskew65$response, test_normal65$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp3_NormalvsNegskew_65 <- ks.boot(test_negskew65$response, test_normal65$response, nboots=1000, alternative = "two.sided", print.level=0)

#compare to observed outcomes
Exp3_norm35_observed <- ks.boot(test_normal35$response, test_normal35$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp3_negskew35_observed <- ks.boot(test_negskew35$response, test_negskew35$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp3_posskew35_observed <- ks.boot(test_posskew35$response, test_posskew35$outcome, nboots=1000, alternative = "two.sided", print.level=0)

Exp3_norm65_observed <- ks.boot(test_normal65$response, test_normal65$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp3_negskew65_observed <- ks.boot(test_negskew65$response, test_negskew65$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp3_posskew65_observed <- ks.boot(test_posskew65$response, test_posskew65$outcome, nboots=1000, alternative = "two.sided", print.level=0)


#analyse test ratings 
Exp3_mean <- Exp3_data_target %>%
  filter(phase == "meantest")

#plot
# cueDescrip_labels <- c("Normal","Positive Skew","Negative Skew","Normal","Positive Skew","Negative Skew")
fig_cols <- c("blue","orange","red","blue","orange","red")

rainplot_mean_Exp3 <- ggplot(Exp3_mean, aes(x = distribution, y = response))+facet_grid(~ mean)+
  geom_boxplot(
    aes(color=distribution),
    width = .3,
    outlier.shape = NA
  )+
  geom_point(
    aes(color = distribution),
    fill = "white",
    shape = 21,
    stroke = .4,
    size = 2,
    position = position_jitter(seed = 1, width = .12)
  )+ylim(0,100)+
  # geom_hline(yintercept=0, linetype = "dashed")+
  stat_summary(fun.data=mean_se, geom = "errorbar",width = 0.2, position = position_dodge(width = 0.2))+
  stat_summary(fun.data=mean_se, geom = "point",size = 3, shape = 21, fill = "white", position = position_dodge(width = 0.2))+
  scale_color_manual(values = fig_cols)+
  scale_fill_manual(values = fig_cols)+
  # scale_x_discrete(labels = cueDescrip_labels)+
  labs(x = "Distribution", y = "Ratings")+ggtitle("(a) Mean Outcome")+
  theme_classic()+
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text =  element_text(size = 14, colour = "black"), axis.title = element_text(size = 20),
        strip.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10,r=0,b= 0,l=0)),plot.title = element_text(size = 20),
        legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#contrasts
mean.aov <- anova_test(data=Exp3_mean, dv = response, wid= participant, within = c("distribution","mean"),effect.size = "pes")

Exp3_meanEstimates <- Descriptives(Exp3_mean,response,cueDescription)

#planned comparison normal vs pos, normal vs neg ratings to determine whether shape influence estimates
Exp3_meanModel <- aov_ez("participant","response",Exp3_mean, within = c("distribution","mean"))
knitr::kable(nice(Exp3_meanModel))

#contrasts
Meanmodel <- emmeans(Exp3_meanModel,c("distribution","mean"))
Meanmodel

#maineffects & interactions
model_contrast <- list(
  "main_NormalvsNeg" = c(1,0,-1,1,0,-1),
  "main_NormalvsPos" = c(0,-1,1,0,-1,1),
  "main_35vs65" = c(1,1,1,-1,-1,-1)
  # "NPintx" = main_NormalvsPos*main_35vs65,
  # "NNintx" = main_NormalvsNeg*main_35vs65
)

meanmodel_out <- contrast(Meanmodel,model_contrast)
Contrast_Table(meanmodel_out)

#outcome mode
Exp3_mode <- Exp3_data_target%>%
  filter(phase=="modetest")

rainplot_mode_Exp3 <- ggplot(Exp3_mode, aes(x = distribution, y = response))+facet_grid(~ mean)+
  geom_boxplot(
    aes(color=distribution),
    width = .3,
    outlier.shape = NA
  )+
  geom_point(
    aes(color = distribution),
    fill = "white",
    shape = 21,
    stroke = .4,
    size = 2,
    position = position_jitter(seed = 1, width = .12)
  )+ylim(0,100)+
  # geom_hline(yintercept=0, linetype = "dashed")+
  stat_summary(fun.data=mean_se, geom = "errorbar",width = 0.2, position = position_dodge(width = 0.2))+
  stat_summary(fun.data=mean_se, geom = "point",size = 3, shape = 21, fill = "white", position = position_dodge(width = 0.2))+
  scale_color_manual(values = fig_cols)+
  scale_fill_manual(values = fig_cols)+
  # scale_x_discrete(labels = cueDescrip_labels)+
  labs(x = "Distribution", y = "Ratings")+ggtitle("(b) Most Common Outcome")+
  theme_classic()+
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text =  element_text(size = 14, colour = "black"), axis.title = element_text(size = 20),
        strip.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10,r=0,b= 0,l=0)),plot.title = element_text(size = 20),
        legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())


Exp3_modeModel <- aov_ez("participant","response",Exp3_mode, within = c("distribution","mean"))
knitr::kable(nice(Exp3_modeModel))

#contrasts
Modemodel <- emmeans(Exp3_modeModel,c("distribution","mean"))
Modemodel

model_contrast <- list(
  "main_NormalvsPos" = c(1,0,-1,1,0,-1),
  "main_NormalvsNeg" = c(0,-1,1,0,-1,1),
  "main_35vs65" = c(1,1,1,-1,-1,-1),
  "NPintx" = main_NormalvsPos*main_35vs65,
  "NNintx" = main_NormalvsNeg*main_35vs65
)
#mode estimate contrast
modemodel_out <- contrast(Modemodel,model_contrast)
Contrast_Table(modemodel_out)

#causal rating
Exp3_causal <- Exp3_data_target%>%
  filter(phase=="causaltest")

rainplot_causal_Exp3 <- ggplot(Exp3_causal, aes(x = distribution, y = response))+facet_grid(~ mean)+
  geom_boxplot(
    aes(color=distribution),
    width = .3,
    outlier.shape = NA
  )+
  geom_point(
    aes(color = distribution),
    fill = "white",
    shape = 21,
    stroke = .4,
    size = 2,
    position = position_jitter(seed = 1, width = .12)
  )+ylim(-100,100)+
  # geom_hline(yintercept=0, linetype = "dashed")+
  stat_summary(fun.data=mean_se, geom = "errorbar",width = 0.2, position = position_dodge(width = 0.2))+
  stat_summary(fun.data=mean_se, geom = "point",size = 3, shape = 21, fill = "white", position = position_dodge(width = 0.2))+
  scale_color_manual(values = fig_cols)+
  scale_fill_manual(values = fig_cols)+geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_x_discrete(labels = cueDescrip_labels)+
  labs(x = "Distribution", y = "Ratings")+ggtitle("(c) Causal rating")+
  theme_classic()+
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text =  element_text(size = 14, colour = "black"), axis.title = element_text(size = 20),
        strip.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10,r=0,b= 0,l=0)),plot.title = element_text(size = 20),
        legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())


test_combine <- ggarrange(rainplot_mean_Exp3,rainplot_mode_Exp3,rainplot_causal_Exp3, nrow=3)
ggsave(paste0(file_name_root, "Exp3_test",graph_file_type),test_combine,"jpg",
       height = gg_width*5, width = gg_width*3, units="cm",dpi = dpi)
include_graphics(paste0(file_name_root,"Exp3_test",graph_file_type))


Exp3_CausalModel <- aov_ez("participant","response",Exp3_causal, within = c("distribution","mean"))
knitr::kable(nice(Exp3_CausalModel))

#contrasts
Causalmodel <- emmeans(Exp3_CausalModel,c("distribution","mean"))
Causalmodel

model_contrast <- list(
  "main_NormalvsNeg" = c(1,0,-1,1,0,-1),
  "main_NormalvsPos" = c(-1,1,0,-1,1,0),
  "main_35vs65" = c(1,1,1,-1,-1,-1),
  "NPintx" = main_NormalvsPos*main_35vs65,
  "NNintx" = main_NormalvsNeg*main_35vs65
)
#mode estimate contrast
causalmodel_out <- contrast(Causalmodel,model_contrast)
Contrast_Table(causalmodel_out)

Exp3_descrip_causal <- Descriptives(Exp3_causal, response,mean)
kable(Exp3_descrip_causal)

################################
experiment <- 4
Exp4_data <- read.csv(here("original_data/Exp4_compiled_data.csv"))
n_distinct(Exp4_data$participant)

#check their rating of the last 10 fixed cues-----
Exp4_fixedcheck <- Exp4_data%>%
  filter(cueDescription%in%c("fixed0","fixed100"))%>%
  filter(phase == "train" & trialNumPerCue > 10)

Exp4_process <- Descriptives(Exp4_fixedcheck,response,participant,cueDescription)
Exp4_process <- as.data.frame(Exp4_process)
Exp4_process <- Exp4_process%>%
  mutate(fixed_check = case_when(
    cueDescription == "fixed0" & mean > 10 ~ "fail",
    cueDescription == "fixed100" & mean < 90 ~ "fail",
    TRUE ~ ""
  ))
extract_subj <- Exp4_process%>%
  filter(fixed_check == "fail")
extract_subj <- unique(extract_subj$participant)

Exp4_data <- Exp4_data %>%
  mutate(fixedcheck = ifelse(
    participant%in%extract_subj, "fail","pass"
  ))

#check if cueDescription matches outcome
Exp4_description_check <- Exp4_data%>%
  filter(!cueDescription%in%c("fixed0","fixed100"))

Exp4_description_check <- Descriptives(Exp4_data,outcome,cueDescription)
kable(Exp4_description_check)

Exp4_description_plot <- Exp4_data%>%
  filter(!cueDescription%in%c("fixed0","fixed100"))

Exp4_description_plot <- ggplot(data = Exp4_description_plot, mapping = aes(x=outcome))+facet_wrap(~cueDescription)+
  # geom_histogram(data = Exp4_description_plot, mapping = aes(x=response, y = ..density..), binwidth = 1, fill = "pink", colour = "red")+
  geom_freqpoly(data = Exp4_description_plot, colour = "red",size = 1, stat ="density",linetype = "dashed")+
  geom_vline(xintercept=35)+geom_vline(xintercept=65)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(limits = c(0,100))

Exp4_data <- Exp4_data%>%
  mutate(
  distribution = case_when(
    cueDescription == "normwide65" | cueDescription == "normwide35" ~ "normal-wide",
    cueDescription == "normnarrow65" | cueDescription == "normnarrow35" ~ "normal-narrow",
    cueDescription == "bimodal35" | cueDescription == "bimodal65" ~ "bimodal",
    cueDescription == "fixed0" | cueDescription == "fixed100" ~ "fixed"
  ),
  mean = case_when(
    cueDescription == "normwide35" | cueDescription == "normnarrow35"| cueDescription == "bimodal35" ~ "35",
    cueDescription == "normwide65" | cueDescription == "normnarrow65"| cueDescription == "bimodal65" ~ "65",
    cueDescription == "fixed0"  ~ "0",
    cueDescription == "fixed100"  ~ "100"
  ))

write_csv(Exp4_data, paste0(file_name_root_data, "Exp4_ProcessedData.csv"))

#update Exp4 to exclude participants who failed manipulation check
Exp4_data <- Exp4_data %>%
  filter(!participant%in%extract_subj)
n_distinct(Exp4_data$participant)

Exp4_data_target <- Exp4_data %>%
  filter(distribution!="fixed")

#ks test on last 20 training predictions
Exp4_train <- Exp4_data_target %>%
  filter(phase == "train" & trialNumPerCue > 20)

#ksTest
test_bimodal35 <- Exp4_train %>%
  filter(mean == 35 & distribution == "bimodal")

test_normal <- Exp4_train %>%
  mutate(contrastGroup = case_when(
    distribution == "normal-narrow" | distribution == "normal-wide" ~ "normal",
    distribution == "bimodal" ~ "bimodal"
  ))%>%
  filter(mean == "35" & contrastGroup == "normal")

Exp4_NormalvsBimodal_35 <- ks.boot(test_bimodal35$response, test_normal$response, nboots=1000, alternative = "two.sided", print.level=0)

test_bimodal65 <- Exp4_train %>%
  filter(mean == "65" & distribution == "bimodal")

test_normal65 <- Exp4_train %>%
  mutate(contrastGroup = case_when(
    distribution == "normal-narrow" | distribution == "normal-wide" ~ "normal",
    distribution == "bimodal" ~ "bimodal"
  ))%>%
  filter(mean == "65" & contrastGroup == "normal")

Exp4_NormalvsBimodal_65 <- ks.boot(test_bimodal65$response, test_normal65$response, nboots=1000, alternative = "two.sided", print.level=0)

#contrast 2
test_narrow35 <- Exp4_train %>%
  filter(mean == 35 & distribution == "normal-narrow")

test_wide35 <- Exp4_train %>%
  filter(mean == 35 & distribution == "normal-wide")


test_narrow65 <- Exp4_train %>%
  filter(mean ==  65 & distribution == "normal-narrow")

test_wide65 <- Exp4_train %>%
  filter(mean == 65  & distribution == "normal-wide")

Exp4_NarrowvsWide_35 <- ks.boot(test_narrow35$response, test_wide35$response, nboots=1000, alternative = "two.sided", print.level=0)
Exp4_NarrowvsWide_65 <- ks.boot(test_narrow65$response, test_wide65$response, nboots=1000, alternative = "two.sided", print.level=0)


#predictions vs observed outcomes
Exp4_narrow35_observed <- ks.boot(test_narrow35$response, test_narrow35$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp4_wide35_observed <- ks.boot(test_wide35$response, test_wide35$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp4_bimodal35_observed <- ks.boot(test_bimodal35$response, test_bimodal35$outcome, nboots=1000, alternative = "two.sided", print.level=0)

Exp4_narrow65_observed <- ks.boot(test_narrow65$response, test_narrow65$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp4_wide65_observed <- ks.boot(test_wide65$response, test_wide65$outcome, nboots=1000, alternative = "two.sided", print.level=0)
Exp4_bimodal65_observed <- ks.boot(test_bimodal65$response, test_bimodal65$outcome, nboots=1000, alternative = "two.sided", print.level=0)

#example plot
Exp4_plot_narrow35 <- ggplot(data = test_narrow35, mapping = aes(x=outcome))+
  geom_histogram(data = test_narrow35, mapping = aes(x=response, y = ..density..), binwidth = 1, fill = "pink", colour = "red")+
  geom_freqpoly(data = test_narrow35, colour = "maroon",size = 1, stat ="density",linetype = "solid", adjust = 2)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 13, colour = "black"), 
        axis.title = element_text(size = 15, colour = "black"))+
  scale_x_continuous(limits = c(0,100))+ggtitle("(a) Normal-Narrow Distribution (35)")

Exp4_plot_wide35 <- ggplot(data = test_wide35, mapping = aes(x=outcome))+
  geom_histogram(data = test_wide35, mapping = aes(x=response, y = ..density..), binwidth = 1, fill = "lightblue", colour = "blue")+
  geom_freqpoly(data = test_wide35, colour = "darkblue",size = 1, stat ="density",linetype = "solid",adjust = 2)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 13, colour = "black"), 
        axis.title = element_text(size = 15, colour = "black"))+  
  scale_x_continuous(limits = c(0,100))+ggtitle("(b) Normal-Wide Distribution (35)")


Exp4_plot_bimodal35 <- ggplot(data = test_bimodal35, mapping = aes(x=outcome))+
  geom_histogram(data = test_bimodal35, mapping = aes(x=response, y = ..density..), binwidth = 1, fill = "#9DC183", colour = "darkgreen")+
  geom_freqpoly(data = test_bimodal35, colour = "darkgreen",size = 1, stat ="density",linetype = "solid",adjust = 1.5)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 13, colour = "black"), 
        axis.title = element_text(size = 15, colour = "black"))+
  scale_x_continuous(limits = c(0,100))+ggtitle("(c) Bimodal Distribution (35)")

Exp4_35s <- ggarrange(Exp4_plot_narrow35,Exp4_plot_wide35,Exp4_plot_bimodal35,ncol=3) 

Exp4_plot_narrow65 <- ggplot(data = test_narrow65, mapping = aes(x=outcome))+
  geom_histogram(data = test_narrow65, mapping = aes(x=response, y = ..density..), binwidth = 1, fill = "pink", colour = "red")+
  geom_freqpoly(data = test_narrow65, colour = "maroon",size = 1, stat ="density",linetype = "solid", adjust = 2)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 13, colour = "black"), 
        axis.title = element_text(size = 15, colour = "black"))+
  scale_x_continuous(limits = c(0,100))+ggtitle("(d) Normal-Narrow Distribution (65)")

Exp4_plot_wide65 <- ggplot(data = test_wide65, mapping = aes(x=outcome))+
  geom_histogram(data = test_wide65, mapping = aes(x=response, y = ..density..), binwidth = 1, fill = "lightblue", colour = "blue")+
  geom_freqpoly(data = test_wide65, colour = "darkblue",size = 1, stat ="density",linetype = "solid",adjust = 2)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 13, colour = "black"), 
        axis.title = element_text(size = 15, colour = "black"))+
  scale_x_continuous(limits = c(0,100))+ggtitle("(e) Normal-Wide Distribution (65)")


Exp4_plot_bimodal65 <- ggplot(data = test_bimodal65, mapping = aes(x=outcome))+
  geom_histogram(data = test_bimodal65, mapping = aes(x=response, y = ..density..), binwidth = 1, fill = "#9DC183", colour = "darkgreen")+
  geom_freqpoly(data = test_bimodal65, colour = "darkgreen",size = 1, stat ="density",linetype = "solid",adjust = 1.5)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 13, colour = "black"), 
        axis.title = element_text(size = 15, colour = "black"))+
  scale_x_continuous(limits = c(0,100))+ggtitle("(f) Bimodal Distribution (65)")

Exp4_65s <- ggarrange(Exp4_plot_narrow65,Exp4_plot_wide65,Exp4_plot_bimodal65,ncol=3) 

Exp4_predOutcomes <- ggarrange(Exp4_35s,Exp4_65s,nrow = 2)

ggsave(paste0(file_name_root, "Exp4_predictionsVoutcomes",graph_file_type),Exp4_predOutcomes,"jpg",
       height = gg_width*2, width = gg_width*3.5, units="cm",dpi = dpi)
include_graphics(paste0(file_name_root,"Exp4_predictionsVoutcomes",graph_file_type))


#analyse test ratings 
Exp4_mean <- Exp4_data_target %>%
  filter(phase == "meantest")

#plot
rainplot_mean_Exp4 <- ggplot(Exp4_mean, aes(x = distribution, y = response))+facet_grid(~ mean)+
  geom_boxplot(
    aes(color=distribution),
    width = .3,
    outlier.shape = NA
  )+
  geom_point(
    aes(color = distribution),
    fill = "white",
    shape = 21,
    stroke = .4,
    size = 2,
    position = position_jitter(seed = 1, width = .12)
  )+ylim(0,100)+
  # geom_hline(yintercept=0, linetype = "dashed")+
  stat_summary(fun.data=mean_se, geom = "errorbar",width = 0.2, position = position_dodge(width = 0.2))+
  stat_summary(fun.data=mean_se, geom = "point",size = 3, shape = 21, fill = "white", position = position_dodge(width = 0.2))+
  scale_color_manual(values = fig_cols)+
  scale_fill_manual(values = fig_cols)+
  # scale_x_discrete(labels = cueDescrip_labels)+
  labs(x = "Distribution", y = "Ratings")+ggtitle("(a) Mean Outcome")+
  theme_classic()+
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text =  element_text(size = 14, colour = "black"), axis.title = element_text(size = 20),
        strip.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10,r=0,b= 0,l=0)),plot.title = element_text(size = 20),
        legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#contrasts
mean.aov <- anova_test(data=Exp4_mean, dv = response, wid= participant, within = c("distribution","mean"),effect.size = "pes")

Exp4_meanEstimates <- Descriptives(Exp4_mean,response,cueDescription)

#planned comparison normal vs pos, normal vs neg ratings to determine whether shape influence estimates
Exp4_meanModel <- aov_ez("participant","response",Exp4_mean, within = c("distribution","mean"))
knitr::kable(nice(Exp4_meanModel))

#contrasts
Meanmodel <- emmeans(Exp4_meanModel,c("distribution","mean"))
Meanmodel

#maineffects & interactions
model_contrast <- list(
  "main_NormalvsBimodal" = c(1,-2,1,1,-2,1),
  "main_NarrowvsWide" = c(1,0,-1,1,0,-1),
  "main_35vs65" = c(1,1,1,-1,-1,-1),
  "NBintx" = main_NormalvsBimodal*main_35vs65,
  "NNintx" = main_NarrowvsWide*main_35vs65
)

meanmodel_out <- contrast(Meanmodel,model_contrast)
Contrast_Table(meanmodel_out)

#mode
Exp4_mode <- Exp4_data_target%>%
  filter(phase=="modetest")

rainplot_mode_Exp4 <- ggplot(Exp4_mode, aes(x = distribution, y = response))+facet_grid(~ mean)+
  geom_boxplot(
    aes(color=distribution),
    width = .3,
    outlier.shape = NA
  )+
  geom_point(
    aes(color = distribution),
    fill = "white",
    shape = 21,
    stroke = .4,
    size = 2,
    position = position_jitter(seed = 1, width = .12)
  )+ylim(0,100)+
  # geom_hline(yintercept=0, linetype = "dashed")+
  stat_summary(fun.data=mean_se, geom = "errorbar",width = 0.2, position = position_dodge(width = 0.2))+
  stat_summary(fun.data=mean_se, geom = "point",size = 3, shape = 21, fill = "white", position = position_dodge(width = 0.2))+
  scale_color_manual(values = fig_cols)+
  scale_fill_manual(values = fig_cols)+
  # scale_x_discrete(labels = cueDescrip_labels)+
  labs(x = "Distribution", y = "Ratings")+ggtitle("(b) Most Common Outcome")+
  theme_classic()+
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text =  element_text(size = 14, colour = "black"), axis.title = element_text(size = 20),
        strip.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10,r=0,b= 0,l=0)),plot.title = element_text(size = 20),
        legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())



Exp4_modeModel <- aov_ez("participant","response",Exp4_mode, within = c("distribution","mean"))
knitr::kable(nice(Exp4_modeModel))

#contrasts
Modemodel <- emmeans(Exp4_modeModel,c("distribution","mean"))
Modemodel

#mode estimate contrast
modemodel_out <- contrast(Modemodel,model_contrast)
Contrast_Table(modemodel_out)

#causal rating
Exp4_causal <- Exp4_data_target%>%
  filter(phase=="causaltest")

rainplot_causal_Exp4 <- ggplot(Exp4_causal, aes(x = distribution, y = response))+facet_grid(~ mean)+
  geom_boxplot(
    aes(color=distribution),
    width = .3,
    outlier.shape = NA
  )+
  geom_point(
    aes(color = distribution),
    fill = "white",
    shape = 21,
    stroke = .4,
    size = 2,
    position = position_jitter(seed = 1, width = .12)
  )+ylim(-100,100)+
  # geom_hline(yintercept=0, linetype = "dashed")+
  stat_summary(fun.data=mean_se, geom = "errorbar",width = 0.2, position = position_dodge(width = 0.2))+
  stat_summary(fun.data=mean_se, geom = "point",size = 3, shape = 21, fill = "white", position = position_dodge(width = 0.2))+
  scale_color_manual(values = fig_cols)+
  scale_fill_manual(values = fig_cols)+
  # scale_x_discrete(labels = cueDescrip_labels)+
  labs(x = "Distribution", y = "Ratings")+ggtitle("(c) Causal Ratings")+geom_hline(yintercept = 0, linetype = "dashed")+
  theme_classic()+
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text =  element_text(size = 14, colour = "black"), axis.title = element_text(size = 20),
        strip.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10,r=0,b= 0,l=0)),plot.title = element_text(size = 20),
        legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())

Exp4_test_combine <- ggarrange(rainplot_mean_Exp4,rainplot_mode_Exp4,rainplot_causal_Exp4,nrow =3)
ggsave(paste0(file_name_root, "Exp4_test_combine",graph_file_type),Exp4_test_combine,"tiff",
       height = gg_width*5, width = gg_width*3, units="cm",dpi = dpi)
include_graphics(paste0(file_name_root,"Exp4_test_combine",graph_file_type))

Exp4_CausalModel <- aov_ez("participant","response",Exp4_causal, within = c("distribution","mean"))
knitr::kable(nice(Exp4_CausalModel))

#contrasts
Causalmodel <- emmeans(Exp4_CausalModel,c("distribution","mean"))
Causalmodel

model_contrast <- list(
  "main_NormalvsBimodal" = c(-2,1,1,-2,1,1),
  "main_NarrowvsWide" = c(0,1,-1,0,1,-1),
  "main_35vs65" = c(1,1,1,-1,-1,-1),
  "NBintx" = main_NormalvsBimodal*main_35vs65,
  "NNintx" = main_NarrowvsWide*main_35vs65
)
#mode estimate contrast
causalmodel_out <- contrast(Causalmodel,model_contrast)
Contrast_Table(causalmodel_out)

Exp4_descrip_causal <- Descriptives(Exp4_causal, response,mean)
kable(Exp4_descrip_causal)


#Are participants using a matching strategy?----------------
Exp1_modelOutput <- read.csv(here("output/Exp1_SimulatedPredictions.csv"))

Exp1_matching_df <- Exp1_modelOutput %>%
  select(participant,trialNumPerCue,cueDescription,response,outcome)%>%
  group_by(participant,cueDescription) %>%  # Ensure calculations are done per participant
  arrange(trialNumPerCue) %>%         # Ensure data is sorted by trial
  mutate(outcome_prev = lag(outcome, n = 1),   # Shift outcomes by one trial to get previous outcome
         prediction_prev = lag(response, n = 1))%>%
  filter(!is.na(outcome_prev))

#could probably have added this to the code above
Exp1_matching_df <- Exp1_matching_df %>%
  group_by(participant, cueDescription) %>%
  mutate(cumulative_mean_outcome = cumsum(ifelse(!is.na(outcome_prev), outcome_prev, 0)) / 
           seq_along(outcome_prev)) %>% # Use ifelse to handle NA values
  ungroup() %>% 
  filter(!is.na(outcome_prev))

# Step 2: Calculate the correlation between prediction on trial n+1 and outcome on trial n
Exp1_correlation_results <- Exp1_matching_df %>%
  group_by(participant,cueDescription)%>%
  filter(!str_starts(cueDescription,"fixed"))%>% #remove fixed cues since can't compute correlation without variance so all fixed trials will produce NA
  summarise(correlation = cor(response, outcome_prev, use = "complete.obs"),
            p_value = cor.test(response, outcome_prev)$p.value,  # Perform correlation test
            .groups = "drop")
print(Exp1_correlation_results)

Exp1_correlation_results$participant <- as.factor(Exp1_correlation_results$participant)
Exp1_correlation_results <- Exp1_correlation_results %>%
  mutate(Experiment = 1)

#plot correlation matrix
Exp1_corrPlot <- ggplot(Exp1_correlation_results, aes(x = participant, y = cueDescription, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-1, 1), name = "Correlation") +
  labs(title = "Experiment 1",
       x = "Participant", y = "Cue type") +
  scale_x_discrete()+
  theme_classic()+
  theme(legend.position = "right", 
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
        panel.background = element_rect())
#one black grid because the participant produced the same prediction on every trial

#
library(dplyr)

# Create a data frame to store results
Exp1_likelihood_ratio_results <- data.frame(cueDescription = character(), 
                                            log_likelihood_ratio = numeric(), 
                                            BayesFactor = numeric(), 
                                            stringsAsFactors = FALSE)

# Loop over each cue type
for (cue in unique(Exp1_matching_df$cueDescription)) {
  # Subset data for the current cue type
  data_subset <- Exp1_matching_df %>% filter(cueDescription == cue)
  
  # Calculate cumulative mean of outcome_prev for each participant up to the current trial
  data_subset <- data_subset %>%
    group_by(participant) %>%
    arrange(trialNumPerCue)
  
  # Model 1: Matching behavior (linear regression model using outcome_prev)
  matching_model <- lm(response ~ outcome_prev, data = data_subset, na.action = na.exclude)
  
  # Model 2: Non-matching behavior (linear regression model using cumulative_mean_outcome)
  non_matching_model <- lm(response ~ cumulative_mean_outcome, data = data_subset, na.action = na.exclude)
  
  # Extract log-likelihoods of both models
  logLik_matching <- logLik(matching_model)
  logLik_non_matching <- logLik(non_matching_model)
  
  # Calculate the log-likelihood ratio
  log_likelihood_ratio <- as.numeric(logLik_matching - logLik_non_matching)
  
  # Calculate the Bayes Factor
  BayesFactor <- exp(log_likelihood_ratio)
  
  # Store the results
  Exp1_likelihood_ratio_results <- rbind(Exp1_likelihood_ratio_results, 
                                         data.frame(cueDescription = cue, 
                                                    log_likelihood_ratio = log_likelihood_ratio, 
                                                    BayesFactor = BayesFactor))
}

Exp1_likelihood_ratio_results$cueDescription <- factor(Exp1_likelihood_ratio_results$cueDescription, 
                                                       levels = c("fixed0", "fixed100", "fixed25", "fixed50", "fixed75", 
                                                                  "rect25","rect50","rect75",
                                                                  "norm25","norm50","norm75"))

# Display the results
Exp1_matching <- Exp1_likelihood_ratio_results%>%
  mutate(Experiment = 1)


#Experiment 2
Exp2_modelOutput <- read.csv(here("output/Exp2_SimulatedPredictions.csv"))

Exp2_matching_df <- Exp2_modelOutput %>%
  select(participant,trialNumPerCue,cueDescription,response,outcome)%>%
  group_by(participant,cueDescription) %>%  # Ensure calculations are done per participant
  arrange(trialNumPerCue) %>%         # Ensure data is sorted by trial
  mutate(outcome_prev = lag(outcome, n = 1),   # Shift outcomes by one trial to get previous outcome
         prediction_prev = lag(response, n = 1))%>%
  filter(!is.na(outcome_prev))%>%  # Remove rows where previous outcome is NA (i.e., first trial)
  filter(str_starts(cueDescription, "target"))

Exp2_matching_df <- Exp2_matching_df %>%
  group_by(participant, cueDescription) %>%
  mutate(cumulative_mean_outcome = cumsum(ifelse(!is.na(outcome_prev), outcome_prev, 0)) / 
           seq_along(outcome_prev)) %>% # Use ifelse to handle NA values
  ungroup() %>%  # Ungroup to avoid issues with further calculations
  filter(!is.na(outcome_prev))

# Step 2: Calculate the correlation between prediction on trial n+1 and outcome on trial n
Exp2_correlation_results <- Exp2_matching_df %>%
  group_by(participant,cueDescription)%>%
  summarise(correlation = cor(response, outcome_prev, use = "complete.obs"),
            p_value = cor.test(response, outcome_prev)$p.value,  # Perform correlation test
            .groups = "drop")
print(Exp2_correlation_results)

Exp2_correlation_results$participant <- as.factor(Exp2_correlation_results$participant)
Exp2_correlation_results <- Exp2_correlation_results %>%
  mutate(Experiment = 2)
#plot correlation matrix
Exp2_corrPlot <- ggplot(Exp2_correlation_results, aes(x = participant, y = cueDescription, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-1, 1), name = "Correlation") +
  labs(title = "Experiment 2",
       x = "Participant", y = "Cue type") +
  scale_x_discrete()+
  theme_classic()+
  theme(legend.position = "right", 
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
        panel.background = element_rect())

#
library(dplyr)

# Create a data frame to store results
Exp2_likelihood_ratio_results <- data.frame(cueDescription = character(), 
                                       log_likelihood_ratio = numeric(), 
                                       BayesFactor = numeric(), 
                                       stringsAsFactors = FALSE)

# Loop over each cue type
for (cue in unique(Exp2_matching_df$cueDescription)) {
  # Subset data for the current cue type
  data_subset <- Exp2_matching_df %>% filter(cueDescription == cue)
  
  # Calculate cumulative mean of outcome_prev for each participant up to the current trial
  data_subset <- data_subset %>%
    group_by(participant) %>%
    arrange(trialNumPerCue)
  
  # Model 1: Matching behavior (linear regression model using outcome_prev)
  matching_model <- lm(response ~ outcome_prev, data = data_subset, na.action = na.exclude)
  
  # Model 2: Non-matching behavior (linear regression model using cumulative_mean_outcome)
  non_matching_model <- lm(response ~ cumulative_mean_outcome, data = data_subset, na.action = na.exclude)
  
  # Extract log-likelihoods of both models
  logLik_matching <- logLik(matching_model)
  logLik_non_matching <- logLik(non_matching_model)
  
  # Calculate the log-likelihood ratio
  log_likelihood_ratio <- as.numeric(logLik_matching - logLik_non_matching)
  
  # Calculate the Bayes Factor
  BayesFactor <- exp(log_likelihood_ratio)
  
  # Store the results
  Exp2_likelihood_ratio_results <- rbind(Exp2_likelihood_ratio_results, 
                                    data.frame(cueDescription = cue, 
                                               log_likelihood_ratio = log_likelihood_ratio, 
                                               BayesFactor = BayesFactor))
}

# Display the results
Exp2_matching <- Exp2_likelihood_ratio_results%>%
  mutate(Experiment = 2)


#use Experiment 3 as test data first
Exp3_modelOutput <- read.csv(here("output/Exp3_SimulatedPredictions.csv"))

Exp3_matching_df <- Exp3_modelOutput %>%
  select(participant,trialNumPerCue,cueDescription,response,outcome)%>%
  group_by(participant,cueDescription) %>%  # Ensure calculations are done per participant
  arrange(trialNumPerCue) %>%         # Ensure data is sorted by trial
  mutate(outcome_prev = lag(outcome, n = 1),   # Shift outcomes by one trial to get previous outcome
         prediction_prev = lag(response, n = 1))%>%
  filter(!is.na(outcome_prev))  # Remove rows where previous outcome is NA (i.e., first trial)

Exp3_matching_df <- Exp3_matching_df %>%
  group_by(participant, cueDescription) %>%
  mutate(cumulative_mean_outcome = cumsum(ifelse(!is.na(outcome_prev), outcome_prev, 0)) / 
           seq_along(outcome_prev)) %>% # Use ifelse to handle NA values
  ungroup() %>%  # Ungroup to avoid issues with further calculations
  filter(!is.na(outcome_prev))

# Step 2: Calculate the correlation between prediction on trial n+1 and outcome on trial n
Exp3_correlation_results <- Exp3_matching_df %>%
  filter(!cueDescription%in%c("fixed0","fixed100"))%>%
  group_by(participant,cueDescription)%>%
  summarise(correlation = cor(response, outcome_prev, use = "complete.obs"),
            p_value = cor.test(response, outcome_prev)$p.value,  # Perform correlation test
            .groups = "drop")
print(Exp3_correlation_results)

Exp3_correlation_results$participant <- as.factor(Exp3_correlation_results$participant)
Exp3_correlation_results <- Exp3_correlation_results %>%
  mutate(Experiment = 3)
#plot correlation matrix
Exp3_corrPlot <- ggplot(Exp3_correlation_results, aes(x = participant, y = cueDescription, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-1, 1), name = "Correlation") +
  geom_text(aes(label = ifelse(p_value < 0.05, sprintf("p=\n%.3f", p_value), "")), # Display p-value only if < 0.05
            color = "black", size = 3) +  
  labs(title = "Participant Prediction (trial n+1) and Observed Outcome (trial n)",
       x = "Participant", y = "Cue type") +
  scale_y_discrete(labels = c("negative skew (35)", "negative skew (65)",
                              "normal (35)", "normal (65)",
                                "positive skew (35)", "positive skew (65)")) +
  theme_classic() +
  theme(legend.position = "right", 
        plot.title = element_text(size = 16),legend.title = element_text(size = 16),
        legend.text = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),size = 15),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size = 15),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(size = 14, colour = "black"), axis.title = element_text(size = 14),
        panel.background = element_rect())

# Display the plot
print(Exp3_corrPlot)

ggsave(paste0(file_name_root, "Exp3_matchingcorr",graph_file_type),Exp3_corrPlot,"tiff",
       height = gg_width*1.5, width = gg_width*3, units="cm",dpi = dpi)
include_graphics(paste0(file_name_root,"Exp3_matchingcorr",graph_file_type))

#correlate predictions from the delta and distributed model against the previous outcome
sim_df <- Exp3_modelOutput
sim_df <- sim_df%>%select(, 1:13) #not sure why this is giving me an error when it works

matching_sim_df <- sim_df %>%
  select(participant,trialNumPerCue,cueDescription,Delta_Predictions,Distributed_Predictions,outcome)%>%
  group_by(participant,cueDescription) %>%  # Ensure calculations are done per participant
  arrange(trialNumPerCue) %>%         # Ensure data is sorted by trial
  mutate(outcome_prev = lag(outcome, n = 1),   # Shift outcomes by one trial to get previous outcome
         delta_prev = lag(Delta_Predictions, n = 1),
         dist_prev = lag(Distributed_Predictions, n = 1))%>%
  filter(!is.na(outcome_prev))

#delta corr prior predictions
correlation_delta <- matching_sim_df %>%
  filter(!cueDescription%in%c("fixed0","fixed100"))%>%
  group_by(participant,cueDescription)%>%
  summarise(correlation = cor(Delta_Predictions, outcome_prev, use = "complete.obs"),
            p_value = cor.test(Delta_Predictions, outcome_prev)$p.value,  # Perform correlation test
            .groups = "drop")
print(correlation_delta)

correlation_delta$participant <- as.factor(correlation_delta$participant)
Exp3_delta_matching <- ggplot(correlation_delta, aes(x = participant, y = cueDescription, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-1, 1), name = "Correlation") +
  geom_text(aes(label = ifelse(p_value < 0.05, sprintf("p=\n%.3f", p_value), "")), # Display p-value only if < 0.05
            color = "black", size = 3) +  
  labs(title = "Simple Delta Prediction (trial n+1) and Observed Outcome (trial n) ",
       x = "Participant", y = "Cue type") +
  scale_y_discrete(labels = c("negative skew (35)", "negative skew (65)",
                              "normal (35)", "normal (65)",
                              "positive skew (35)", "positive skew (65)")) +
  theme_classic() +
  theme(legend.position = "right", 
        plot.title = element_text(size = 16),legend.title = element_text(size = 16),
        legend.text = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),size = 15),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size = 15),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(size = 14, colour = "black"), axis.title = element_text(size = 14),
        panel.background = element_rect())

correlation_dist <- matching_sim_df %>%
  filter(!cueDescription%in%c("fixed0","fixed100"))%>%
  group_by(participant,cueDescription)%>%
  summarise(correlation = cor(Distributed_Predictions, outcome_prev, use = "complete.obs"),
            p_value = cor.test(Distributed_Predictions, outcome_prev)$p.value,  # Perform correlation test
            .groups = "drop")
print(correlation_dist)

correlation_dist$participant <- as.factor(correlation_dist$participant)

Exp3_dist_matching <- ggplot(correlation_dist, aes(x = participant, y = cueDescription, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-1, 1), name = "Correlation") +
  geom_text(aes(label = ifelse(p_value < 0.05, sprintf("p=\n%.3f", p_value), "")), # Display p-value only if < 0.05
            color = "black", size = 3) +  
  labs(title = "Distributed Model Prediction (trial n+1) and Observed Outcome (trial n) ",
       x = "Participant", y = "Cue type") +
  scale_y_discrete(labels = c("negative skew (35)", "negative skew (65)",
                              "normal (35)", "normal (65)",
                              "positive skew (35)", "positive skew (65)")) +
  theme_classic() +
  theme(legend.position = "right", 
        plot.title = element_text(size = 16),legend.title = element_text(size = 16),
        legend.text = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),size = 15),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size = 15),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(size = 14, colour = "black"), axis.title = element_text(size = 14),
        panel.background = element_rect())

combine_correlations <- ggarrange(Exp3_delta_matching,Exp3_dist_matching, nrow=2)
ggsave(paste0(file_name_root, "Exp3_matchingcorr_models",graph_file_type),combine_correlations,"jpg",
       height = gg_width*3, width = gg_width*3, units="cm",dpi = dpi)

library(dplyr)

# Create a data frame to store results
likelihood_ratio_results <- data.frame(cueDescription = character(), 
                                       log_likelihood_ratio = numeric(), 
                                       BayesFactor = numeric(), 
                                       stringsAsFactors = FALSE)

# Loop over each cue type
for (cue in unique(Exp3_matching_df$cueDescription)) {
  # Subset data for the current cue type
  data_subset <- Exp3_matching_df %>% filter(cueDescription == cue)
  
  # Calculate cumulative mean of outcome_prev for each participant up to the current trial
  data_subset <- data_subset %>%
    group_by(participant) %>%
    arrange(trialNumPerCue)
  
  # Model 1: Matching behavior (linear regression model using outcome_prev)
  matching_model <- lm(response ~ outcome_prev, data = data_subset, na.action = na.exclude)
  
  # Model 2: Non-matching behavior (linear regression model using cumulative_mean_outcome)
  non_matching_model <- lm(response ~ cumulative_mean_outcome, data = data_subset, na.action = na.exclude)
  
  # Extract log-likelihoods of both models
  logLik_matching <- logLik(matching_model)
  logLik_non_matching <- logLik(non_matching_model)
  
  # Calculate the log-likelihood ratio
  log_likelihood_ratio <- as.numeric(logLik_matching - logLik_non_matching)
  
  # Calculate the Bayes Factor
  BayesFactor <- exp(log_likelihood_ratio)
  
  # Store the results
  likelihood_ratio_results <- rbind(likelihood_ratio_results, 
                                    data.frame(cueDescription = cue, 
                                               log_likelihood_ratio = log_likelihood_ratio, 
                                               BayesFactor = BayesFactor))
}

# Display the results
Exp3_matching <- likelihood_ratio_results%>%
  mutate(Experiment = 3)

#Experiment 4
Exp4_modelOutput <- read.csv(here("output/Exp4_SimulatedPredictions.csv"))

Exp4_matching_df <- Exp4_modelOutput %>%
  select(participant,trialNumPerCue,cueDescription,response,outcome)%>%
  group_by(participant,cueDescription) %>%  # Ensure calculations are done per participant
  arrange(trialNumPerCue) %>%         # Ensure data is sorted by trial
  mutate(outcome_prev = lag(outcome, n = 1),   # Shift outcomes by one trial to get previous outcome
         prediction_prev = lag(response, n = 1))%>%
  filter(!is.na(outcome_prev))%>%  # Remove rows where previous outcome is NA (i.e., first trial)
  filter(!str_starts(cueDescription, "fixed"))

Exp4_matching_df <- Exp4_matching_df %>%
  group_by(participant, cueDescription) %>%
  mutate(cumulative_mean_outcome = cumsum(ifelse(!is.na(outcome_prev), outcome_prev, 0)) / 
           seq_along(outcome_prev)) %>% # Use ifelse to handle NA values
  ungroup() %>%  # Ungroup to avoid issues with further calculations
  filter(!is.na(outcome_prev))

# Step 2: Calculate the correlation between prediction on trial n+1 and outcome on trial n
Exp4_correlation_results <- Exp4_matching_df %>%
  group_by(participant,cueDescription)%>%
  summarise(correlation = cor(response, outcome_prev, use = "complete.obs"),
            p_value = cor.test(response, outcome_prev)$p.value,  # Perform correlation test
            .groups = "drop")
print(Exp4_correlation_results)

Exp4_correlation_results$participant <- as.factor(Exp4_correlation_results$participant)
Exp4_correlation_results <- Exp4_correlation_results %>%
  mutate(Experiment = 4)
#plot correlation matrix
Exp4_corrPlot <- ggplot(Exp4_correlation_results, aes(x = participant, y = cueDescription, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-1, 1), name = "Correlation") +
  geom_text(aes(label = ifelse(p_value < 0.05, sprintf("p=\n%.3f", p_value), "")), # Display p-value only if < 0.05
            color = "black", size = 3) +  
  labs(title = "Experiment 4",
       x = "Participant", y = "Cue type") +
  scale_x_discrete()+
  scale_y_discrete(labels = c("bimodal (35)","bimodal (65)",
                              "normal narrow (35)", "normal narrow (65)",
                              "normal wide (35)", "normal wide (65)"))+
  theme_classic()+
  theme(legend.position = "right", 
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
        panel.background = element_rect())

combine_corrPlots <- ggarrange(Exp3_corrPlot,Exp4_corrPlot,nrow=2)
ggsave(paste0(file_name_root, "matchingcorrPlots",graph_file_type),combine_corrPlots,"jpg",
       height = gg_width*3, width = gg_width*3.5, units="cm",dpi = dpi)
include_graphics(paste0(file_name_root,"matchingcorrPlots",graph_file_type))
#
library(dplyr)

# Create a data frame to store results
Exp4_likelihood_ratio_results <- data.frame(cueDescription = character(), 
                                            log_likelihood_ratio = numeric(), 
                                            BayesFactor = numeric(), 
                                            stringsAsFactors = FALSE)

# Loop over each cue type
for (cue in unique(Exp4_matching_df$cueDescription)) {
  # Subset data for the current cue type
  data_subset <- Exp4_matching_df %>% filter(cueDescription == cue)
  
  # Calculate cumulative mean of outcome_prev for each participant up to the current trial
  data_subset <- data_subset %>%
    group_by(participant) %>%
    arrange(trialNumPerCue)
  
  # Model 1: Matching behavior (linear regression model using outcome_prev)
  matching_model <- lm(response ~ outcome_prev, data = data_subset, na.action = na.exclude)
  
  # Model 2: Non-matching behavior (linear regression model using cumulative_mean_outcome)
  non_matching_model <- lm(response ~ cumulative_mean_outcome, data = data_subset, na.action = na.exclude)
  
  # Extract log-likelihoods of both models
  logLik_matching <- logLik(matching_model)
  logLik_non_matching <- logLik(non_matching_model)
  
  # Calculate the log-likelihood ratio
  log_likelihood_ratio <- as.numeric(logLik_matching - logLik_non_matching)
  
  # Calculate the Bayes Factor
  BayesFactor <- exp(log_likelihood_ratio)
  
  # Store the results
  Exp4_likelihood_ratio_results <- rbind(Exp4_likelihood_ratio_results, 
                                         data.frame(cueDescription = cue, 
                                                    log_likelihood_ratio = log_likelihood_ratio, 
                                                    BayesFactor = BayesFactor))
}

# Display the results
Exp4_matching <- Exp4_likelihood_ratio_results%>%
  mutate(Experiment = 4)

#combine all matching output
matching_output <- rbind(Exp1_matching,Exp2_matching,Exp3_matching,Exp4_matching)%>%
  select(Experiment,cueDescription,log_likelihood_ratio,BayesFactor)%>%
  filter(!cueDescription%in%c("fixed0","fixed100"))

write_csv(matching_output, paste0(file_name_root_data, "MatchingStrategy_ModelOutput.csv"))

matching_correlation <- rbind(Exp1_correlation_results,Exp2_correlation_results,correlation_results,Exp4_correlation_results)%>%
  select(Experiment,participant,cueDescription,correlation,p_value)

write_csv(matching_correlation, paste0(file_name_root_data, "MatchingStrategy_Correlations.csv"))
