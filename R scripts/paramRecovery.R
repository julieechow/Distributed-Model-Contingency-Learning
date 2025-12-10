#Parameter Recovery analysis
#This script requires model recovery fit parameter values from modelRecovery.R (also found in output/)
#Note that this code allows you to run parameter recovery for simulations generated with both data-informed and uninformed process; only data-informed sim parameter recovery are reported in Supp Materials
#This script will save the dataframe generated with iterative sampling, and the average correlation coefficients between true and recovered parameter values for each iteration

#load packages
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(here)
library(dplyr)
library(tidyr)
library(purrr)

graph_file_type<- ".jpeg"
experiment <- 4 #change this to call the appropriate experiment

generatedDataType <- "informed"


#Load the param values from simulated data fit
if (generatedDataType=="informed"){
  recovery_df <- read.csv(paste0("output/Experiment_", experiment, "_Informed_ModelRecovery_AllFits.csv"))
  
  # load the RData file
  loaded_name <- load(paste0("simulated_data/Exp", experiment, "_Informed_allSims.RData"))
  

  assign("allSimData", get(loaded_name))
  rm(list = loaded_name)     # so you don't have to add "informed" to subsequent code

  
}else{
  recovery_df <- read.csv(paste0("output/Experiment_", experiment, "_ModelRecovery_AllFits.csv"))
  load(paste0("simulated_data/Exp", experiment, "_allSims.RData"))
}

if(experiment ==2){
  allSimData <- allSimData%>%
    filter(simulation <=1000)
}

#this gets the true param values used to generate the simulated data
true_params <- allSimData %>%
  select(simulation, sampled_from_participant, model,
         param_T,param_A,param_ContextA, param_NoiseSD, param_OutcomeSigma)%>%
  distinct()%>%
  mutate(
    true_param4 = ifelse(model == "Delta",
                         param_NoiseSD, param_OutcomeSigma)
  )%>%
  select(simulation, sampled_from_participant,model,
         true_T = param_T,
         true_A = param_A,
         true_ContextA = param_ContextA,
         true_param4)

#this gets the recovered parameter values
recovery_fits <- recovery_df %>%
  group_by(simulation, generated_model, fit_model) %>%
  slice_min(bic, n = 1) %>% #this takes the lowest BIC for each simulation so you dont have 50 values to 1 true_param value
  ungroup()%>%
  rename(rec_T = param1,
         rec_A = param2,
         rec_ContextA = param3,
         rec_param4 = param4) #noise if Delta, OutcomeVariance if Distributed


param_recovery <- recovery_fits%>%
  left_join(true_params, by = c("simulation", "generated_model" = "model"))

#only look at cases where the generating and fit model are the same
matched_models <- param_recovery %>%
  filter(generated_model == fit_model)

#Parameter Recovery Analysis----

#First need to select one sim data per participant for each iteration 
#Then get the correlation between true and recovered parameter values

# Function to compute correlations for one random sample
compute_sample_correlations <- function(data, iteration_number) {
  # Sample one simulation per unique participant
  sampled_data <- data %>%
    group_by(sampled_from_participant, generated_model) %>%
    slice_sample(n = 1) %>% #randomly picks one row for each participant and model combination
    ungroup() %>%
    mutate(iteration = iteration_number)  # Add iteration column
  
  # Compute correlations for each parameter
  correlations <- sampled_data %>%
    group_by(generated_model) %>%
    summarise(
      cor_T = cor(true_T, rec_T, use = "complete.obs", method = "pearson"),
      cor_A = cor(true_A, rec_A, use = "complete.obs", method = "pearson"),
      cor_ContextA = cor(true_ContextA, rec_ContextA, use = "complete.obs", method = "pearson"),
      cor_param4 = cor(true_param4, rec_param4, use = "complete.obs", method = "pearson"),
      .groups = "drop"
    ) %>%
    mutate(iteration = iteration_number)  
  
  # Return both as a list
  return(list(
    sampled_data = sampled_data,
    correlations = correlations
  ))
}

# Run the sampling process n_iteration times
n_iterations <- 100
set.seed(123)  # For reproducibility

# Collect results
all_results <- map(1:n_iterations, ~compute_sample_correlations(matched_models, .x))

# Extract correlation results
correlation_results <- map_dfr(all_results, ~.x$correlations)

# Extract all sampled data from each iteration
all_sampled_iterations <- map_dfr(all_results, ~.x$sampled_data)
#save data from each iteration
save(all_sampled_iterations, file = paste0("output/Exp", experiment, "_", generatedDataType, "_IterativeSampling.RData"))
output_filename <- paste0("output/Experiment_", experiment, "_", generatedDataType, "_IterativeSampling.csv")
write.csv(all_sampled_iterations, file = output_filename, row.names = FALSE)
#save correlations
output_filename <- paste0("output/Experiment_", experiment, "_", generatedDataType, "_Correlations.csv")
write.csv(correlation_results, file = output_filename, row.names = FALSE)

# Compute summary statistics: mean and 95% confidence intervals by model
summary_stats <- correlation_results %>%
  pivot_longer(cols = starts_with("cor_"), 
               names_to = "parameter", 
               values_to = "correlation") %>%
  group_by(generated_model, parameter) %>%
  summarise(
    median_r = median(correlation, na.rm = TRUE),
    mean_r = mean(correlation, na.rm = TRUE),
    lower_95 = quantile(correlation, 0.025, na.rm = TRUE),
    upper_95 = quantile(correlation, 0.975, na.rm = TRUE),
    sd_r = sd(correlation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    parameter = gsub("cor_", "", parameter)
  )

#this gives us some indication of the % of recovered param values that falls within the range of true values
rec_in_true_range <- matched_models %>%
  group_by(generated_model) %>%
  summarise(
    # For each parameter, calculate % within range
    pct_T_in_range = mean(rec_T >= min(true_T) & rec_T <= max(true_T)) * 100,
    pct_A_in_range = mean(rec_A >= min(true_A) & rec_A <= max(true_A)) * 100,
    pct_ContextA_in_range = mean(rec_ContextA >= min(true_ContextA) & rec_ContextA <= max(true_ContextA)) * 100,
    pct_param4_in_range = mean(rec_param4 >= min(true_param4) & rec_param4 <= max(true_param4)) * 100
  )

assign(paste0("Exp", experiment, "_rec_in_true_range"), rec_in_true_range)

#what proportion of recovered params are in range?
experiment <- 4
load(paste0("output/Exp", experiment, "_", generatedDataType, "_", "IterativeSampling.RData"))

true_param_range <- all_sampled_iterations %>%
  group_by(generated_model) %>%
  summarise(
    min_param1 = min(true_T),
    max_param1 = max(true_T),
    min_param2 = min(true_A),
    max_param2 = max(true_A),
    min_param3 = min(true_ContextA),
    max_param3 = max(true_ContextA),
    min_param4 = min(true_param4),
    max_param4 = max(true_param4)
  )

# Add a flag for whether parameters are in range
recovery_range <- all_sampled_iterations %>%
  left_join(true_param_range, by = "generated_model") %>%
  mutate(
    T_in_range = rec_T >= min_param1 & rec_T <= max_param1,
    A_in_range = rec_A >= min_param2 & rec_A <= max_param2,
    ContextA_in_range = rec_ContextA >= min_param3 & rec_ContextA <= max_param3,
    param4_in_range = rec_param4 >= min_param4 & rec_param4 <= max_param4,
    all_in_range = T_in_range & A_in_range & ContextA_in_range & param4_in_range
  ) %>%
  select(-starts_with("min_param"), -starts_with("max_param"))

get_range <- recovery_range %>%
  group_by(generated_model) %>%
    summarise(
      pct_all_in_range = mean(all_in_range) * 100,
      pct_T_in_range = mean(T_in_range) * 100,
      pct_A_in_range = mean(A_in_range) * 100,
      pct_ContextA_in_range = mean(ContextA_in_range) * 100,
      pct_param4_in_range = mean(param4_in_range) * 100
    )%>%
  ungroup()%>%
  mutate(Experiment = experiment)

assign(paste0("Exp", experiment, "_rec_in_true_range"), get_range)

rec_in_true_range <- rbind(Exp1_rec_in_true_range,Exp2_rec_in_true_range,Exp3_rec_in_true_range,Exp4_rec_in_true_range)

rec_in_true_range%>%
  group_by(generated_model)%>%
  summarise(      
            pct_all_in_range = mean(pct_all_in_range), 
            pct_T_in_range = mean(pct_T_in_range),
            pct_A_in_range = mean(pct_A_in_range),
            pct_ContextA_in_range = mean(pct_ContextA_in_range),
            pct_param4_in_range = mean(pct_param4_in_range)
  )


# Plotting-----
# Visualize the distribution
correlation_results %>%
  pivot_longer(cols = starts_with("cor_"), 
               names_to = "parameter", 
               values_to = "correlation") %>%
  filter(generated_model == "Delta", parameter %in% c("cor_A", "cor_ContextA","cor_T","cor_param4")) %>%
  ggplot(aes(x = correlation)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept = median(correlation)), color = "red", linetype = "dashed") +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of correlations across 100 iterations")

#plot scatterplot from data in a single random iteration
plot_corr <- all_sampled_iterations%>%
  filter(iteration == 37)%>% #this could be any number
  select(sampled_from_participant, generated_model,
         starts_with("true_"), starts_with("rec_")) %>%
  pivot_longer(
    cols = starts_with(c("true_", "rec_")),
    names_to = c("type", "parameter"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(
    id_cols = c(sampled_from_participant, generated_model, parameter),
    names_from = type,
    values_from = value) %>%
  mutate(
    parameter_label = case_when(
      parameter == "A" ~ "Alpha",
      parameter == "ContextA" ~ "Context Alpha",
      parameter == "T" ~ "Softmax Temperature k",
      parameter == "param4" & generated_model == "Delta" ~ "Noise SD",
      parameter == "param4" & generated_model == "Distributed" ~ "Outcome Sigma",
      TRUE ~ parameter
    ),
    # Create combined facet label
    facet_label = paste(generated_model, "-", parameter_label)
  ) %>%
  ggplot(aes(x = true, y = rec, color = generated_model, fill = generated_model)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +  # Adds correlation line with 95% CI
  facet_wrap(~facet_label, scales = "free", ncol = 4) +  # Changed to facet_wrap
  labs(title = "True vs Recovered Values",
       y = "Recovered Parameter Values",
       x = "True Parameter Values") +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        plot.subtitle = element_text(size = 15))

#saved as pdf and strip labels fixed in Adobe Photoshop
graph_file_type <- ".pdf"
ggsave(paste0(file_name_root, "Exp", experiment, "_paramRecovery", graph_file_type), plot_corr, "pdf",
       height = gg_height*1.5, width = gg_width*2.5, units = "cm", dpi = dpi)

