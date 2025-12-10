#script to generate simulated data using best fitting params from model fit (data-informed procedure)
###################################################
# Generate a large number of simulated datasets for each model
# (1) randomly select a participant and adopt their specific trial order,
# (2) get the ppt's best fitting param values from model fit
###################################################
library(here)
library(tidyverse)
source(here("ModelFunctionsLL.R"))

# Load file with best-fitting parameter results
experiment <- 1 #change this for separate experiment files
read_params <- read.csv(file= paste("output/Experiment",experiment,"FitData_single_MLE_All.csv",sep='_'))

dat <- read.csv(file= paste("output/Exp",experiment,"_ProcessedData.csv",sep=''))
dat <- dat %>%
  filter(fixedcheck=="pass")

nCues = max(dat$cueNum)
dat$cueNum <- as.factor(dat$cueNum)
dat <- subset(dat, dat$phase == "train") # only fit training and test phase
split_data <- split(x = dat, dat$participant)  #splits data by participant to fit on participant level

# Number of synthetic datasets to generate per MODEL
nSimulations <- 1000  

# Number of trials per participant (same as in original data)
nTrials <- max(dat$trialNumber)
nTrain <- nTrials
nParticipants <- length(split_data)

# Extract parameter matrices for Delta and Distributed models
delta_best_pars <- read_params %>%
  select(DeltaTVal, DeltaAVal, DeltaContextAVal, DeltaNoiseSDVal) %>%
  as.matrix()

dist_best_pars <- read_params %>%
  select(DistributedTVal, DistributedAVal, DistributedContextAVal, DistributedOutcomeSigmaVal) %>%
  as.matrix()

sub_array <- read_params$sub_array  


# Storage for simulation results
deltaSimData <- list()
distSimData  <- list()

# --- GENERATE SIMULATIONS ---
cat("Generating", nSimulations, "simulations per model...\n")

for (sim in 1:nSimulations) {
  if (sim %% 100 == 0) cat("Simulation:", sim, "/", nSimulations, "\n")
  
  set.seed(sim)  # for reproducibility
  
  # --- Randomly sample a participant's trial structure (with replacement) ---
  sampled_participant_idx <- sample(1:nParticipants, 1)
  
  # Extract trial structure from randomly sampled participant
  Outcome <- split_data[[sampled_participant_idx]][1:nTrials, "outcome"]
  Phase <- split_data[[sampled_participant_idx]][1:nTrials, "phase"]
  cueNum <- split_data[[sampled_participant_idx]][1:nTrials, "cueNum"]
  cueLabel <- split_data[[sampled_participant_idx]][1:nTrials, "cueLabel"]
  trialNumPerCue <- split_data[[sampled_participant_idx]][1:nTrials, "trialNumPerCue"]
  cueNums <- as.numeric(cueNum)
  #if experiment = 2
  # Group <- split_data[[sampled_participant_idx]][1:nTrials, "group"]

  #Informed procedure: use participants' best fitting parameter values from model fit
  sim_output_delta <- LL([sampled_participant_idx,], 
                         nTrials = nTrials, nTrain = nTrain,
                         model = "Delta", experiment = experiment, output = "full")
  
  deltaSimData[[sim]] <- data.frame(
    simulation = sim,
    sampled_from_participant = sub_array[sampled_participant_idx],
     # group = Group, #only relevant for experiment 2
    trial = 1:nTrials,
    model = "Delta",
    phase = Phase,
    cueNum = cueNum,
    cueLabel = cueLabel,
    trialNumPerCue = trialNumPerCue,
    outcome = Outcome,
    predicted_response = sim_output_delta$discretePredictions,
    param_T = delta_best_pars[sampled_participant_idx,1],
    param_A = delta_best_pars[sampled_participant_idx,2],
    param_ContextA = delta_best_pars[sampled_participant_idx,3],
    param_NoiseSD = delta_best_pars[sampled_participant_idx,4]
  )
  
  # --- Distributed model simulation ---

  sim_output_dist <- LL(dist_best_pars[sampled_participant_idx, ], 
                        nTrials = nTrials, nTrain = nTrain,
                        model = "Distributed", experiment = experiment, output = "full")
  
  distSimData[[sim]] <- data.frame(
    simulation = sim,
    sampled_from_participant = sub_array[sampled_participant_idx],
    # group = Group, #only relevant for experiment 2
    trial = 1:nTrials,
    model = "Distributed",
    phase = Phase,
    cueNum = cueNum,
    cueLabel = cueLabel,
    trialNumPerCue = trialNumPerCue,
    outcome = Outcome,
    predicted_response = sim_output_dist$discretePredictions,
    param_T = dist_best_pars[sampled_participant_idx,1],
    param_A = dist_best_pars[sampled_participant_idx,2],
    param_ContextA = dist_best_pars[sampled_participant_idx,3],
    param_OutcomeSigma = dist_best_pars[sampled_participant_idx,4]
  )
}

# Combine all simulations into single dataframes
deltaSimData <- bind_rows(deltaSimData)
distSimData <- bind_rows(distSimData)

cat("Simulation complete!\n")
cat("Delta model simulations:", nrow(deltaSimData), "rows\n")
cat("Distributed model simulations:", nrow(distSimData), "rows\n")

# --- COMBINE BOTH MODELS INTO A SINGLE DATAFRAME ---
allSimData_Informed <- bind_rows(deltaSimData, distSimData)

# --- SAVE TO CSV ---
output_filename <- paste0("simulated_data/Experiment_", experiment, "_Informed_SimulatedData_", nSimulations,"sims.csv")
write.csv(allSimData_Informed, file = output_filename, row.names = FALSE)

cat("Data saved to:", output_filename, "\n")

#save as R data object for model recovery
save(allSimData_Informed, file = paste0("simulated_data/Exp", experiment, "_Informed_allSims.RData"))
