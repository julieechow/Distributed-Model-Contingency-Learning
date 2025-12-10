#Model recovery analysis
#This script requires simulated data generated from the generate_predictions and generate_predictions_informed (can be found in output/simulated_data/)
#This script requires parameter values from empirical model fit from ModelFitsAll_Single_MLE_Exp[1-4] (also found in output/ALL_ModelOutputs.csv)
#This script will save the parameter values used in each run of model recovery on the same simulated data

library(tidyverse)
library(here)

# --- Load input ---
# redefine if not following directly from generate predictions
experiment <- 1
generatedDataType <- "informed"

if (generatedDataType=="informed"){
  load(paste0("simulated_data/Exp", experiment, "_Informed_allSims.RData"))  
  allSimData <- allSimData_Informed #rename to allSimData so the following code is identical for both informed and uninformed sims
}else{
  load(paste0("simulated_data/Exp", experiment, "_allSims.RData"))
}

dat <- allSimData %>%
  filter(simulation <= 1000) # Select the first 1000 simulations per model if simulated more than 1000 datasets

#get relevant parameter values 
best_params <- read.csv(here("output/ALL_ModelOutputs.csv")) %>%
  filter(Experiment == experiment)

# --- Parameter bounds (from real fitted data) ---
param_bounds <- list(
  DeltaTVal = range(best_params$DeltaTVal),
  DeltaAVal = range(best_params$DeltaAVal),
  DeltaContextAVal = range(best_params$DeltaContextAVal),
  DeltaNoiseSDVal = range(best_params$DeltaNoiseSDVal),
  
  DistTVal = range(best_params$DistributedTVal),
  DistAVal = range(best_params$DistributedAVal),
  DistContextAVal = range(best_params$DistributedContextAVal),
  DistNoiseSDVal = range(best_params$DistributedOutcomeSigmaVal)
)

# --- Prepare to fit ---
npars <- 4
nTrials <- max(dat$trial)

# Get unique simulation IDs for each model
delta_sims <- dat %>% filter(model == "Delta") %>% distinct(simulation) %>% pull(simulation)
dist_sims <- dat %>% filter(model == "Distributed") %>% distinct(simulation) %>% pull(simulation)

fit_results <- list()

# --- Loop over simulations ---
for (sim in delta_sims) {
  
  #this only prints every 50 sims to avoid showing me every single run
  if (sim %% 50 == 0) cat("Fitting Delta simulation", sim, "/", length(delta_sims), "\n")
  
  generated_model <- "Delta"
  
  sim_data <- dat %>%
    filter(simulation == sim, model == generated_model)
    
  simulation <- unique(sim_data$simulation)
  Response <- sim_data$predicted_response
  Outcome  <- sim_data$outcome
  Phase    <- sim_data$phase
  cueNum   <- sim_data$cueNum
  cueLabel <- sim_data$cueLabel
  trialNumPerCue <- sim_data$trialNumPerCue
  
  # Store outcome of every fit
  fit_run_results <- data.frame()
  
  # Fit with both models
  for (fit_model in c("Delta", "Distributed")) {
    
    bestfit <- 1e6
    numstarts <- 50
    
    for (i in 1:numstarts) {
      if (fit_model == "Delta") {
        startPars <- c(
          runif(1, param_bounds$DeltaTVal[1], param_bounds$DeltaTVal[2]),
          runif(1, param_bounds$DeltaAVal[1], param_bounds$DeltaAVal[2]),
          runif(1, param_bounds$DeltaContextAVal[1], param_bounds$DeltaContextAVal[2]),
          runif(1, param_bounds$DeltaNoiseSDVal[1], param_bounds$DeltaNoiseSDVal[2])
        )
      } else {
        startPars <- c(
          runif(1, param_bounds$DistTVal[1], param_bounds$DistTVal[2]),
          runif(1, param_bounds$DistAVal[1], param_bounds$DistAVal[2]),
          runif(1, param_bounds$DistContextAVal[1], param_bounds$DistContextAVal[2]),
          runif(1, param_bounds$DistNoiseSDVal[1], param_bounds$DistNoiseSDVal[2])
        )
      }
      
      # Run optimizer
      res <- optim(
        startPars, LL,
        nTrials = nTrials, nTrain = nTrials,
        model = fit_model,
        experiment = experiment,
        output = "fit",
        lower = c(0.0001, 0.0001, 0.0001, 0.0001),
        upper = c(5, 1, 1, 1),
        method = 'L-BFGS-B',
        control = list(trace = 0)
      )

      
      # Compute BIC for this run
      bic <- res$value * 2 + npars * log(nTrials - 1)
      
      # Store every run result
      fit_run_results <- rbind(fit_run_results, data.frame(
        simulation = sim,
        generated_model = generated_model,
        fit_model = fit_model,
        start = i,
        negLL = res$value,
        bic = bic,
        param1 = res$par[1], #saves the best fitting param values for each run; #Temp 
        param2 = res$par[2], #Cue Alpha
        param3 = res$par[3], #Context Alpha
        param4 = res$par[4] #Prediction Noise
      ))
    } # end numstarts
  } # end fit_model loop
  
  # Store to global list
  fit_results[[length(fit_results) + 1]] <- fit_run_results
  
} # end Delta simulations

# Fit Distributed-generated data
for (sim in dist_sims) {
  if (sim %% 50 == 0) cat("Fitting Distributed simulation", sim, "/", length(dist_sims), "\n")
  
  generated_model <- "Distributed"
  
  sim_data <- dat %>%
    filter(simulation == sim, model == generated_model)
  
  simulation <- unique(sim_data$simulation)
  Response <- sim_data$predicted_response
  Outcome  <- sim_data$outcome
  Phase    <- sim_data$phase
  cueNum   <- sim_data$cueNum
  cueLabel <- sim_data$cueLabel
  trialNumPerCue <- sim_data$trialNumPerCue
  
  # Store outcome of every fit
  fit_run_results <- data.frame()
  
  # Fit with both models
  for (fit_model in c("Delta", "Distributed")) {
    
    bestfit <- 1e6
    numstarts <- 50
    
    for (i in 1:numstarts) {
      if (fit_model == "Delta") {
        startPars <- c(
          runif(1, param_bounds$DeltaTVal[1], param_bounds$DeltaTVal[2]),
          runif(1, param_bounds$DeltaAVal[1], param_bounds$DeltaAVal[2]),
          runif(1, param_bounds$DeltaContextAVal[1], param_bounds$DeltaContextAVal[2]),
          runif(1, param_bounds$DeltaNoiseSDVal[1], param_bounds$DeltaNoiseSDVal[2])
        )
      } else {
        startPars <- c(
          runif(1, param_bounds$DistTVal[1], param_bounds$DistTVal[2]),
          runif(1, param_bounds$DistAVal[1], param_bounds$DistAVal[2]),
          runif(1, param_bounds$DistContextAVal[1], param_bounds$DistContextAVal[2]),
          runif(1, param_bounds$DistNoiseSDVal[1], param_bounds$DistNoiseSDVal[2])
        )
      }
      
      # Run optimizer
      res <- optim(
        startPars, LL,
        nTrials = nTrials, nTrain = nTrials,
        model = fit_model,
        experiment = experiment,
        output = "fit",
        lower = c(0.0001, 0.0001, 0.0001, 0.0001),
        upper = c(5, 1, 1, 1),
        method = 'L-BFGS-B',
        control = list(trace = 0)
      )
      
      # Compute BIC for this run
      bic <- res$value * 2 + npars * log(nTrials - 1)
      
      # Store every run result
      fit_run_results <- rbind(fit_run_results, data.frame(
        simulation = sim,
        generated_model = generated_model,
        fit_model = fit_model,
        start = i,
        negLL = res$value,
        bic = bic,
        param1 = res$par[1], #Temperature
        param2 = res$par[2], #Cue Alpha
        param3 = res$par[3], #Contex Alpha
        param4 = res$par[4] #Node Activation Variance
      ))
    } # end numstarts
  } # end fit_model loop
  
  # Store to global list
  fit_results[[length(fit_results) + 1]] <- fit_run_results
  
} # end Distributed simulations

cat("Model recovery fitting complete!\n")

# --- Combine results ---
recovery_df <- bind_rows(fit_results)

# Save all fit results
if (generatedDataType=="informed"){
  write.csv(recovery_df, 
            file = paste0("output/Experiment_", experiment, "_Informed_ModelRecovery_AllFits.csv"),
            row.names = FALSE) 
}else{
  write.csv(recovery_df, 
            file = paste0("output/Experiment_", experiment, "_ModelRecovery_AllFits.csv"),
            row.names = FALSE)
}


# Summary: find the best fit of all runs for each simulated dataset
recovery_summary <- recovery_df %>%
  group_by(simulation, generated_model, fit_model) %>%
  summarise(
    minNegLL = min(negLL),
    minBIC = min(bic),
    best_param1 = param1[which.min(bic)],
    best_param2 = param2[which.min(bic)],
    best_param3 = param3[which.min(bic)],
    best_param4 = param4[which.min(bic)],
    .groups = "drop"
  ) %>%
  group_by(simulation, generated_model) %>%
  mutate(
    chosen_model = if_else(
      minBIC[fit_model == "Delta"] < minBIC[fit_model == "Distributed"],
      "Delta",
      "Distributed"
    ),
    deltaBIC_minus_distBIC = minBIC[fit_model == "Delta"] - minBIC[fit_model == "Distributed"]
  ) %>%
  ungroup()


summary_overall <- recovery_summary %>%
  group_by(generated_model) %>%
  summarise(
    n = n(),
    percent_Delta = mean(chosen_model == "Delta") * 100,
    percent_Distributed = mean(chosen_model == "Distributed") * 100,
    .groups = "drop"
  )


# --- Save outputs ---
if (generatedDataType=="informed"){
  write.csv(recovery_df, paste0("output/Experiment", experiment, "_Informed_ModelRecovery_AllSims.csv"), row.names = FALSE)
}else{
  write.csv(recovery_df, paste0("output/Experiment", experiment, "_ModelRecovery_AllSims.csv"), row.names = FALSE)
}

#Code to generate individual experiment tile plots
Exp1_generateMatrix <- summary_overall %>%
  pivot_longer(
    cols = c("percent_Delta","percent_Distributed"),
    names_to = "chosen_model",
    names_prefix = "percent_",
    values_to = "percent"
  )

# Heatmap plot
Exp_modelR_matrix <- ggplot(Exp1_generateMatrix, aes(x = generated_model, y = chosen_model, fill = percent)) +
  geom_tile(color="black") +
  geom_text(aes(label = round(percent,1)), color="white", size=5) +
  scale_fill_continuous(limits= c(0,100),low="lightblue", high="darkblue")+
  # scale_fill_gradient(low="lightblue", high="darkblue") +
  theme_minimal(base_size = 14) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        panel.grid = element_blank())+
  labs(title="Experiment 1", x="Generating Model", y="Recovered Model")

file_name_root <- 'fig/'
graph_file_type <- '.jpeg'
dpi <- 96
gg_width <- 10
gg_height <- 10
if (generatedDataType=="informed"){
  ggsave(paste0(file_name_root, "Exp", experiment, "_Informed_modelR_matrix", graph_file_type), Exp_modelR_matrix, "jpeg",
         height = gg_height*1.5, width = gg_width*1.5, units = "cm", dpi = dpi)
}else{
  ggsave(paste0(file_name_root, "Exp", experiment, "_modelR_matrix", graph_file_type), Exp_modelR_matrix, "jpeg",
         height = gg_height*1.5, width = gg_width*1.5, units = "cm", dpi = dpi)
}
