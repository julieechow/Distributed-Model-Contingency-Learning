#script to generate simulated data for model recovery
###################################################
# Generate a large number of simulated datasets for each model 
# Decouples empirical fits from simulated data by (1) randomly selecting participant trial order, and 
# (2) randomly selecting starting parameter values from range of best fitting parameter values from empirical model fit
###################################################

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

# --- DETERMINE PARAMETER RANGES ACROSS ALL PARTICIPANTS ---
# For Delta model
delta_param_ranges <- list(
  T = range(delta_best_pars[, 1]),
  A = range(delta_best_pars[, 2]),
  ContextA = range(delta_best_pars[, 3]),
  NoiseSD = range(delta_best_pars[, 4])
)

# For Distributed model
dist_param_ranges <- list(
  T = range(dist_best_pars[, 1]),
  A = range(dist_best_pars[, 2]),
  ContextA = range(dist_best_pars[, 3]),
  OutcomeSigma = range(dist_best_pars[, 4])
)

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
  # participant_label <- paste0("Sim_", sim)  # create unique simulation ID
  
  # Extract trial structure from randomly sampled participant
  Outcome <- split_data[[sampled_participant_idx]][1:nTrials, "outcome"]
  Phase <- split_data[[sampled_participant_idx]][1:nTrials, "phase"]
  cueNum <- split_data[[sampled_participant_idx]][1:nTrials, "cueNum"]
  cueLabel <- split_data[[sampled_participant_idx]][1:nTrials, "cueLabel"]
  trialNumPerCue <- split_data[[sampled_participant_idx]][1:nTrials, "trialNumPerCue"]
  cueNums <- as.numeric(cueNum)
  
  
  #if experiment = 2
  # Group <- split_data[[sampled_participant_idx]][1:nTrials, "group_array"]
  
  # --- Delta model simulation ---
  # Randomly sample parameters from the observed range
  delta_params <- c(
    runif(1, delta_param_ranges$T[1], delta_param_ranges$T[2]),
    runif(1, delta_param_ranges$A[1], delta_param_ranges$A[2]),
    runif(1, delta_param_ranges$ContextA[1], delta_param_ranges$ContextA[2]),
    runif(1, delta_param_ranges$NoiseSD[1], delta_param_ranges$NoiseSD[2])
  )
  
  sim_output_delta <- LL(delta_params, 
                         nTrials = nTrials, nTrain = nTrain,
                         model = "Delta", experiment = experiment, output = "full")
  
  deltaSimData[[sim]] <- data.frame(
    #group = Group, #only relevant for experiment 2
    simulation = sim,
    sampled_from_participant = sub_array[sampled_participant_idx],
    trial = 1:nTrials,
    model = "Delta",
    phase = Phase,
    cueNum = cueNum,
    cueLabel = cueLabel,
    trialNumPerCue = trialNumPerCue,
    outcome = Outcome,
    predicted_response = sim_output_delta$discretePredictions,
    param_T = delta_params[1],
    param_A = delta_params[2],
    param_ContextA = delta_params[3],
    param_NoiseSD = delta_params[4]
  )
  
  # --- Distributed model simulation ---
  # Randomly sample parameters from the observed range
  dist_params <- c(
    runif(1, dist_param_ranges$T[1], dist_param_ranges$T[2]),
    runif(1, dist_param_ranges$A[1], dist_param_ranges$A[2]),
    runif(1, dist_param_ranges$ContextA[1], dist_param_ranges$ContextA[2]),
    runif(1, dist_param_ranges$OutcomeSigma[1], dist_param_ranges$OutcomeSigma[2])
  )
  
  sim_output_dist <- LL(dist_params, 
                        nTrials = nTrials, nTrain = nTrain,
                        model = "Distributed", experiment = experiment, output = "full")
  
  distSimData[[sim]] <- data.frame(
    simulation = sim,
    #group = Group, #only relevant for experiment 2
    sampled_from_participant = sub_array[sampled_participant_idx],
    trial = 1:nTrials,
    model = "Distributed",
    phase = Phase,
    cueNum = cueNum,
    cueLabel = cueLabel,
    trialNumPerCue = trialNumPerCue,
    outcome = Outcome,
    predicted_response = sim_output_dist$discretePredictions,
    param_T = dist_params[1],
    param_A = dist_params[2],
    param_ContextA = dist_params[3],
    param_OutcomeSigma = dist_params[4]
  )
}

# Combine all simulations into single dataframes
deltaSimData <- bind_rows(deltaSimData)
distSimData <- bind_rows(distSimData)

cat("Simulation complete!\n")
cat("Delta model simulations:", nrow(deltaSimData), "rows\n")
cat("Distributed model simulations:", nrow(distSimData), "rows\n")

# --- COMBINE BOTH MODELS INTO A SINGLE DATAFRAME ---
allSimData <- bind_rows(deltaSimData, distSimData)

# --- SAVE TO CSV ---
output_filename <- paste0("simulated_data/Experiment_", experiment, "_SimulatedData_", nSimulations,"sims.csv")
write.csv(allSimData, file = output_filename, row.names = FALSE)

cat("Data saved to:", output_filename, "\n")

#save as R data object for model recovery
save(allSimData, file = paste0("simulated_data/Exp", experiment, "_allSims.RData"))
