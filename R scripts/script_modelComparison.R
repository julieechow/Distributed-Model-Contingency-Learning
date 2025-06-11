#compile model outputs and do model comparison of parameter fits

#load individual experiment data output from ModelFits
library(here)
library(tidyverse)

dat1 <- read.csv(here("output/Experiment_1_FitData_single_MLE_All.csv"))
dat2 <- read.csv(here("output/Experiment_2_FitData_single_MLE_All.csv"))
dat3 <- read.csv(here("output/Experiment_3_FitData_single_MLE_All.csv"))
dat4 <- read.csv(here("output/Experiment_4_FitData_single_MLE_All.csv"))

dat1 <- dat1 %>%
  mutate(Experiment = 1,group_array = NA)
dat2 <- dat2 %>%
  mutate(Experiment = 2)
dat3 <- dat3 %>%
  mutate(Experiment = 3,group_array = NA)
dat4 <- dat4 %>%
  mutate(Experiment = 4,group_array = NA)

model_output_df <- rbind(dat1,dat2,dat3,dat4)
model_output_df <- model_output_df%>%
  mutate(BICdiff = deltaBIC - distBIC)

#calculate bayes factor to determine evidence for better model fit with the Distributed Model relative to the Delta Model
#function to calculate means---
Descriptives <- function(.data, .dv, ...) {
  out <- .data %>%
    group_by(...) %>%
    dplyr::summarise(mean = mean({{ .dv }}),
                     n = n(),
                     sd = sd({{ .dv }}),
                     se = sd({{ .dv }})/sqrt(n))
}

#average and compute BF---
descrip_BIC <- Descriptives(model_output_df,BICdiff,Experiment)
descrip_BIC <- as.data.frame(descrip_BIC)

descrip_BIC <- descrip_BIC%>%
  mutate(BF = exp(mean/2))

write.csv(model_output_df, "output/ALL_ModelOutputs.csv")

