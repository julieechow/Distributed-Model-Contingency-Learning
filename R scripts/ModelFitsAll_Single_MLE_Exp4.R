
source("ModelFunctionsLL.R")

library(tidyverse)

#variables for saving output
experiment <- 4
fitTo <- "All"

dat <- read.csv("output/Exp4_ProcessedData.csv", header = TRUE) #need long format trial by trial data
dat <- dat %>%
  filter(fixedcheck=="pass")

nCues = max(dat$cueNum)

dat$cueNum <- as.factor(dat$cueNum)
dat <- subset(dat, dat$phase == "train") # only fit training and test phase
split_data <- split(x = dat,dat$participant)

#initialise variables
nTrials <- max(dat$trialNumber) 
nTrain <- nTrials
nParticipants <- length(split_data)

delta_best_pars <- array(rep(NA), dim = c(nParticipants,4))
delta_bestValue <- array()

dist_best_pars <- array(rep(NA), dim = c(nParticipants,4))
dist_bestValue <- array()

sub_array <- array()
group_array <- array()

deltaFullOutput <- list()
distFullOutput <- list()

#Fit data
###############################################################################################
for(p in 1:nParticipants){
  #this is getting arrays of all the things we need from each trial
  Response <- split_data[[p]][1:nTrials,"response"]
  Outcome <- split_data[[p]][1:nTrials,"outcome"]
  Phase <- split_data[[p]][1:nTrials,"phase"]
  cueNum <- split_data[[p]][1:nTrials,"cueNum"]
  cueLabel <- split_data[[p]][1:nTrials,"cueLabel"]
  trialNumPerCue <- split_data[[p]][1:nTrials,"trialNumPerCue"]
  cueNums <- as.numeric(cueNum)#split_data[[p]][1:nTrials,"AllTrial"]
  #group = split_data[[p]][1,"group"]
  participant = split_data[[p]][1,"participant"]
  
  #this is saving subject number and group to attach to fit data later
  sub_array[p] <- participant
  #group_array[p] <- group
  
  numstarts <- 100 
  bestfit <- 10000
  for (i in 1:numstarts) {
    #sets the starting point for each parameter as a random point within parameter range
    startt <- runif(1,0.0001,5) #temperature parameter for softmax
    starta <- runif(1,0.0001,1) #cue learning rate
    startxa <- runif(1,0.0001,1) #context learning rate
    startsd <- runif(1,0.0001,1) #noisy prediction around EV
    
    startPars <- c(startt,starta,startxa,startsd)  # puts them in an array to feed into optim
    
    # optim function: need to define any parameters you need in model function. Lower and upper are lower and upper ranges of params.
    results <- optim(startPars, LL, nTrials = nTrials, nTrain = nTrain, model = "Delta", experiment = experiment, output = "fit",
                     lower = c(0.0001,0.0001,0.0001,0.0001), upper = c(5,1,1,1), method = c('L-BFGS-B'), control = "trace") 
    #saves all the output you need from the model
    if (results$value < bestfit) {
      bestfit <- results$value
      delta_bestValue[p] <- results$value
      delta_best_pars[p,1] <- results$par[1]
      delta_best_pars[p,2] <- results$par[2]
      delta_best_pars[p,3] <- results$par[3]
      delta_best_pars[p,4] <- results$par[4]
    }
    print(p)
    print(bestfit)
    
  }
  print("Delta done")
  #print(delta_best_pars)
  
  #gets choice probabilities on each trial
  deltaFullOutput[[p]] <- LL(delta_best_pars[p,], nTrials = nTrials, nTrain = nTrain, model = "Delta", experiment = experiment, output = "full")
  # delta_output <- cbind(p,cueNum,probability)
    
  ### Fit Distributed model
  
  bestfit <- 10000
  for (i in 1:numstarts) {    
    startt <- runif(1,0.0001,5) #temperature parameter for softmax
    starta <- runif(1,0.0001,1) #cue learning rate
    startxa <- runif(1,0.0001,1) #context learning rate
    startsd <- runif(1,0.0001,1) #outcome node sigma
    
    startPars <- c(startt,starta,startxa,startsd)  # puts them in an array to feed into optim
    
    results <- optim(startPars, LL, nTrials = nTrials, nTrain = nTrain, model = "Distributed", experiment = experiment, output = "fit",
                     lower = c(0.0001,0.0001,0.0001,0.0001), upper = c(5,1,1,1), method = c('L-BFGS-B'))
    
    if (results$value < bestfit) {
      bestfit <- results$value    
      dist_bestValue[p] <- results$value
      dist_best_pars[p,1] <- results$par[1]
      dist_best_pars[p,2] <- results$par[2]
      dist_best_pars[p,3] <- results$par[3]
      dist_best_pars[p,4] <- results$par[4]
    }
    print(p)
    print(bestfit)
    
  }
  print("Distributed done")
  
  #gets choice probabilities on each trial
  distFullOutput[[p]] <- LL(dist_best_pars[p,], nTrials = nTrials, nTrain = nTrain, model = "Distributed",  experiment = experiment, output = "full")
  # dist_output <- cbind( p,cueNum,probability)
  
} #end the entire subject loop


#save them for posthoc predictions
delta_best_pars <- as.data.frame(delta_best_pars)
names(delta_best_pars) <- c("DeltaTVal","DeltaAVal","DeltaContextAVal","DeltaNoiseSDVal")
delta_best_pars

dist_best_pars <- as.data.frame(dist_best_pars)
names(dist_best_pars) <- c("DistributedTVal","DistributedAVal","DistributedContextAVal","DistributedOutcomeSigmaVal")
dist_best_pars

npars <- 4
deltaBIC <- delta_bestValue*2 + npars*log(nTrials-1) 
distBIC <- dist_bestValue*2 + npars*log(nTrials-1)

#combine the parameters and lnL's into columns and save
AllFits <- cbind(sub_array, delta_best_pars,delta_bestValue,deltaBIC, dist_best_pars,dist_bestValue,distBIC)


#save to csv
write.csv(AllFits, file= paste(file_name_root_data,"Experiment",experiment,"FitData_single_MLE", fitTo,Sys.Date(),".csv",sep='_'))

#extract delta and distribution full outputs and save----
Exp4_data <- read.csv(here("original_data/Exp4_compiled_data.csv"))

#function to do this for every participant----
# Function to add "discretePredictions" to Exp4_data for each participant
add_discrete_predictions <- function(data, modeloutput) {
  
  # Get unique participants from the data
  participants <- unique(data$participant)
  
  # Initialize an empty list to store the combined data for all participants
  combined_data <- list()
  
  # Loop through each participant by index
  for (i in seq_along(participants)) {
    ppt <- participants[i]  # Current participant ID
    
    # Filter data for the current participant and training phase
    data_ppt <- data %>% filter(phase == "train", participant == ppt)
    
    # Extract the corresponding discrete predictions from modeloutput using the index 'i'
    discretePredictions <- modeloutput[[i]][["discretePredictions"]]
    
    # Add the discretePredictions as a new column
    data_ppt <- cbind(data_ppt, discretePredictions)
    
    # Append the combined data to the list
    combined_data[[length(combined_data) + 1]] <- data_ppt
  }
  
  # Combine all participant data into one data frame
  combined_data_df <- do.call(rbind, combined_data)
  
  return(combined_data_df)
}

# Extract trial by trial discrete model predictions----
Exp4_modelOutput <- add_discrete_predictions(Exp4_data, deltaFullOutput)
names(Exp4_modelOutput)[names(Exp4_modelOutput) == 'discretePredictions'] <- 'Delta_Predictions'

Exp4_modelOutput <- add_discrete_predictions(Exp4_modelOutput, distFullOutput)
names(Exp4_modelOutput)[names(Exp4_modelOutput) == 'discretePredictions'] <- 'Distributed_Predictions'

Exp4_description_check <- Descriptives(Exp4_modelOutput,outcome,cueDescription)
kable(Exp4_description_check)

Exp4_description_plot <- ggplot(data = Exp4_modelOutput, mapping = aes(x=outcome))+facet_wrap(~cueDescription)+
  geom_histogram(data = Exp4_modelOutput, mapping = aes(x=outcome,after_stat(count)), binwidth = 1, fill = "pink", colour = "red")+
  # geom_freqpoly(data = Exp3_data, colour = "red",size = 1, stat ="density",linetype = "dashed")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

write_csv(Exp4_modelOutput, paste0(file_name_root_data, "Exp4_SimulatedPredictions.csv"))


#modified function to add 280 x 101 matrix
add_response_probabilities <- function(data, modeloutput) {
  
  # Get unique participants from the data
  participants <- unique(data$participant)
  
  # Initialize an empty list to store the combined data for all participants
  combined_data <- list()
  
  # Loop through each participant by index
  for (i in seq_along(participants)) {
    ppt <- participants[i]  # Current participant ID
    
    # Takes the dataframe from above so already filtered for train phase
    data_ppt <- data %>% filter(participant == ppt)
    
    # Extract the corresponding response probabilities from modeloutput using the index 'i'
    responseProbabilities <- modeloutput[[i]][["responseProbabilities"]]

    # Convert the response probabilities matrix to a data frame
    responseProbabilities_df <- as.data.frame(responseProbabilities)
    
    # Generate column names for the response probabilities
    colnames(responseProbabilities_df) <- paste0("responseProb_", 0:(ncol(responseProbabilities_df)-1))
    
    # Combine the data_ppt with the response probabilities
    data_ppt <- cbind(data_ppt, responseProbabilities_df)
    
    # Append the combined data to the list
    combined_data[[length(combined_data) + 1]] <- data_ppt
  }
  
  # Combine all participant data into one data frame
  combined_data_df <- do.call(rbind, combined_data)
  
  return(combined_data_df)
}

# Example usage
Exp4_modelOutput <- add_response_probabilities(Exp4_modelOutput, deltaFullOutput)
colnames(Exp4_modelOutput) <- sub("^responseProb_", "delta_responseProb_", colnames(Exp4_modelOutput))

Exp4_modelOutput <- add_response_probabilities(Exp4_modelOutput, distFullOutput)
colnames(Exp4_modelOutput) <- sub("^responseProb_", "dist_responseProb_", colnames(Exp4_modelOutput))

#extract terminal Vs using a modified function 
#cannot merge delta and dist model terminal Vs because one is 1 value the other is 101 values

#delta model
extract_terminalV <- function(data, modeloutput) {
  
  # Get unique participants from the data
  participants <- unique(data$participant)
  
  # Initialize an empty list to store the combined data for all participants
  combined_data <- list()
  
  # Loop through each participant by index
  for (i in seq_along(participants)) {
    ppt <- participants[i]  # Current participant ID
    
    # Extract the corresponding discrete predictions from modeloutput using the index 'i'
    terminalV <- modeloutput[[i]][["terminalVs"]]
    
    participant <- rep(ppt,each = length(terminalV))
    cueNum <- rep(1:9) #context is 9
    
    
    # Add the discretePredictions as a new column
    data_ppt <- cbind(participant, cueNum, terminalV)
    
    # Append the combined data to the list
    combined_data[[length(combined_data) + 1]] <- data_ppt
  }
  
  # Combine all participant data into one data frame
  combined_data_df <- do.call(rbind, combined_data)
  
  return(combined_data_df)
}


Exp4_delta_terminalV <- extract_terminalV(Exp4_data, deltaFullOutput)
Exp4_delta_terminalV <- as.data.frame(Exp4_delta_terminalV)

test <- Exp4_data %>%
  select(cueNum, cueDescription) %>%
  distinct() %>%  # Use distinct to get unique cueNum and cueDescription pairs
  mutate(cueNum = as.integer(cueNum))  # Ensure cueNum is of integer type

# Add the additional row for cueNum 9
test <- test %>% add_row(cueNum = 9, cueDescription = "context")

# Merge cueDescription into Exp4_delta_terminalV based on cueNum
Exp4_delta_terminalV <- Exp4_delta_terminalV %>%
  left_join(test, by = "cueNum")

Exp4_delta_terminalV <- Exp4_delta_terminalV %>%
  mutate(terminalV_plusContext = ifelse(cueNum == 9, 
                                        terminalV, 
                                        terminalV + terminalV[cueNum == 9][1]))

Exp4_delta_terminalV <- Exp4_delta_terminalV %>%
  select(participant,cueNum,cueDescription,terminalV,terminalV_plusContext)

#sanity checks
Exp4_description_check <- Descriptives(Exp4_delta_terminalV,terminalV_plusContext,cueDescription)
kable(Exp4_description_check)

Exp4_description_plot <- Exp4_data%>%
  filter(!cueDescription%in%c("fixed0","fixed100"))

Exp4_description_plot <- ggplot(data = Exp4_delta_terminalV, mapping = aes(x=terminalV_plusContext))+facet_wrap(~cueDescription)+
  # geom_histogram(data = Exp4_description_plot, mapping = aes(x=response, y = ..density..), binwidth = 1, fill = "pink", colour = "red")+
  geom_freqpoly(data = Exp4_delta_terminalV, colour = "red",size = 1, stat ="density",linetype = "dashed")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

write_csv(Exp4_delta_terminalV, paste0(file_name_root_data, "Exp4_delta_terminalV.csv"))

#distributed model
extract_distribution_terminalV <- function(data, modeloutput) {
  
  # Get unique participants from the data
  participants <- unique(data$participant)
  
  # Initialize an empty list to store the combined data for all participants
  combined_data <- list()
  
  # Loop through each participant by index
  for (i in seq_along(participants)) {
    ppt <- participants[i]  # Current participant ID

    # Extract the corresponding response probabilities from modeloutput using the index 'i'
    terminalVs <- modeloutput[[i]][["terminalVs"]]
    
    # Convert the response probabilities matrix to a data frame
    terminalVs_df <- as.data.frame(terminalVs)

    # Generate column names for the response probabilities
    colnames(terminalVs_df) <- paste0("terminalV_", 0:(ncol(terminalVs_df)-1))
    
    participant <- rep(ppt,each = nrow(terminalVs_df))
    cueNum <- rep(1:9) #context is 9
    
    # Add the discretePredictions as a new column
    data_ppt <- cbind(participant, cueNum, terminalVs_df)
    
    # Append the combined data to the list
    combined_data[[length(combined_data) + 1]] <- data_ppt
  }
  
  # Combine all participant data into one data frame
  combined_data_df <- do.call(rbind, combined_data)
  
  return(combined_data_df)
}

Exp4_dist_terminalV <- extract_distribution_terminalV(Exp4_data, distFullOutput)

test <- Exp4_data %>%
  select(cueNum, cueDescription) %>%
  distinct() %>%  # Use distinct to get unique cueNum and cueDescription pairs
  mutate(cueNum = as.integer(cueNum))  # Ensure cueNum is of integer type

# Add the additional row for cueNum 9
test <- test %>% add_row(cueNum = 9, cueDescription = "context")

# Merge cueDescription into Exp4_delta_terminalV based on cueNum
Exp4_dist_terminalV <- Exp4_dist_terminalV %>%
  left_join(test, by = "cueNum")

write_csv(Exp4_dist_terminalV, paste0(file_name_root_data, "Exp4_dist_terminalV.csv"))
