#Likelihood function
###################################################
LL <- function(x,nTrials,nTrain,model,experiment,output){
  
  #gets all the starting parameters that were fed into function
  t <- x[1]
  a <- x[2]
  xa <- x[3]
  ds <- x[4]*100
  
  #initialise expected values for each option (QVs are scaled EVs)
  outcome_values <- c(0:100)
  nONodes <- 101
  nCNodes <- nCues+1 #+1 for context as the last number
  
  #DeltaModel
  if(model == "Delta"){
    QV <- replicate(101, 0)
    Vs <- replicate(nCNodes, 0)  
    Vs[nCNodes] <- 50 #associative strength starts at value of 50 assuming ppts biased towards middle of scale
  }
  
  #DistributedModel
  if(model == "Distributed"){
    QV <- replicate(101, 0)
    Vs <- array(runif(n = nONodes*nCNodes, min = -.001, max = .001), dim=c(nCNodes,nONodes)) # starting activation of each node anywhere between -.001 and .001
  }
  
  LL <- vector(length = nTrials)
  TrialType <- vector(length = nTrials)
  
  #these are important for full output return
  pred_response <- vector(length = nTrials)
  r_probabilities <- array(0, dim=c(nTrials,nONodes))
  
  for(trial in 1:nTrials){ 
    
    #if(trial>1){#fit every trial including trial 1
    
    #determine which trial chosen  
    obs_response <- Response[trial]
    obs_choiceN <- obs_response+1 
    phas <- Phase[trial]
    
    #for training trials
    if(trial <= nTrain){
      
      cue_present <- cueNums[trial]
      
      #update expected value for chosen option only
      #DeltaModel
      if(model == "Delta"){ # update by a simple delta rule
        # establish predictions by creating a distribution around the EV
        EV <- Vs[cue_present] + Vs[nCNodes] #V for context = 50
        QV <- dnorm(outcome_values, mean=EV, sd=ds)
        
        # likelihood of the prediction actually made - softmax
        c <- 3^t -1
        #allTemps <- min(700, c * QV)
        allTemps <- c * QV
        allTemps[allTemps>700] <- 700
        temp <- allTemps[obs_choiceN]
        
        #calculate probability of choosing each response
        sMax <- exp(temp) / sum(exp(allTemps))
        if(sMax==0){
          LL[trial] <- log(.Machine$double.xmin)
        }else{
          LL[trial] <- log(sMax)
        }
        
        #if doing a full output then generate a discrete prediction stochastically
        if(output == "full") {
          sAll <- exp(allTemps) / sum(exp(allTemps))
          r_probabilities[trial,] <- sAll
          #choose a response
          x = runif(n = 1, min = 0, max = 1)
          for(jj in 1:nONodes){ 
            if(x <= sum(sAll[1:jj])){
              if(jj == 1){
                pred_response[trial] <- outcome_values[jj]
              }else if(x > sum(sAll[1:jj-1])){
                pred_response[trial] <- outcome_values[jj]  
              }
            }
          }
        }

        #update
        Vs[cue_present] <- Vs[cue_present] + a * (Outcome[trial] - EV)  #update cue
        Vs[nCNodes] <- Vs[nCNodes] + xa * (Outcome[trial] - EV)  #update context
      }#end delta rule update
      
      
      #Distributed Model
      if(model == "Distributed"){ # update a distributed network via delta rule
        # establish predictions by creating a distribution around the EV
        cue_present <- cueNums[trial]
        
        EV <- Vs[cue_present,] + Vs[nCNodes,]
        QV <- EV
        
        # likelihood of the prediction actually made
        c <- 3^t -1
        #allTemps <- min(700, c * QV)
        allTemps <- c * QV
        allTemps[allTemps>700] <- 700
        temp <- allTemps[obs_choiceN]
        
        sMax <- exp(temp) / sum(exp(allTemps))
        if(sMax==0){
          LL[trial] <- log(.Machine$double.xmin)
        }else{
          LL[trial] <- log(sMax)
        }
        
        #if doing a full output then generate a discrete prediction stochastically
        if(output == "full") {
          sAll <- exp(allTemps) / sum(exp(allTemps))
          r_probabilities[trial,] <- sAll
          #choose a response
          x = runif(n = 1, min = 0, max = 1)
          for(jj in 1:nONodes){ 
            if(x <= sum(sAll[1:jj])){
              if(jj == 1){
                pred_response[trial] <- outcome_values[jj]
              }else if(x > sum(sAll[1:jj-1])){
                pred_response[trial] <- outcome_values[jj]  
              }
            }
          }
        }
        
        #update
        #this will work for the simple case where all nodes have the same sd and have mu 
        #spaced 1 integer apart - will need to do in loop if it gets more complex than this
        OnodeAct <- dnorm(outcome_values, mean=Outcome[trial], sd=ds)
        OnodeAct <- OnodeAct/sum(OnodeAct)
        Vs[cue_present,] <- Vs[cue_present,] + a * (OnodeAct - EV)  #update cue
        Vs[nCNodes,] <- Vs[nCNodes,] + xa * (OnodeAct - EV)  #update context        

      }#end distributed update      
    }
    
  } #end trial loop

  LLsum <- -sum(LL) #minus summed log likelihood - takes sum of the log likelihoods of the individual trial LLs
  
  if(output == "fit") {
    return(LLsum)
  } else if(output == "full") {
    return(list(LLfit = LLsum, discretePredictions = pred_response, responseProbabilities = r_probabilities, terminalVs = Vs ))
  }
} #end function loop

