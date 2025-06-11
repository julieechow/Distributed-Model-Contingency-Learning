#plot figures from manuscript
library(ggpubr)
graph_file_type <- ".tif"

#note this code assumes you've already done the work compiling the model_output_df dataframe from script_modelComparison

# Function to generate Gaussian distributions and plot them
plot_gaussian_distributions <- function(means, sd, black_index) {
  # Create a data frame to store the distributions
  data <- data.frame()
  
  # Generate data for each distribution
  for (i in seq_along(means)) {
    x <- seq(0, 100, length.out = 1000)  # Define the x-axis range (Outcome Value)
    y <- dnorm(x, mean = means[i], sd = sd)  # Compute the Gaussian values
    color <- ifelse(i == black_index, "black", "lightgrey")  # Assign color
    
    # Add to the data frame
    data <- rbind(data, data.frame(x = x, y = y, group = i, color = color))
  }
  
  # Separate the black distribution from the others for explicit layering
  black_data <- data[data$color == "black", ]
  grey_data <- data[data$color == "lightgrey", ]
  
  # Plot the distributions
  ggplot() +
    geom_line(data = grey_data, aes(x = x, y = y, group = group, color = color), size = 1) +
    geom_line(data = black_data, aes(x = x, y = y, group = group, color = color), size = 1) +
    scale_color_identity() +  # Use the specified colors directly
    labs(x = "Outcome Value", y = "Node Activation") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 15,colour = "black"),
      axis.title =  element_text(size = 18,colour = "black"),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_blank()
      #element_rect(color = "black", fill = NA)
    )
}

# Example usage:
# Specify the means, standard deviation, and the index of the black distribution
means <- c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95)  # Specify the means of the distributions
sd <- 10  # Specify the standard deviation
black_index <- 4  # Index of the distribution to be colored black

# Plot the figure
Distarchitecture <- plot_gaussian_distributions(means, sd, black_index)
ggsave(paste0(file_name_root, "Figure1", graph_file_type), Distarchitecture, "tiff",
       height = gg_height*1.5, width = gg_width*1.5, units = "cm", dpi = dpi)


#model outputs----
model_output_df <- read.csv(here("output/ALL_ModelOutputs.csv"))

#plots the BIC difference score for each participant, separated by Experiment
fig_BICdiff<- model_output_df %>%
  ggplot(aes(x=factor(Experiment),y=BICdiff))+
  geom_violin(aes(y = BICdiff))+geom_boxplot(width=0.2, color="red", alpha=1.5, outlier.color = NA) +
  geom_jitter(aes(y = BICdiff),height = 0, width = 0.1,alpha = 0.7)+
  labs(y = "BIC difference",x = "Experiment")+geom_hline(yintercept=0, linetype = "dashed")+
  theme_bw()+
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        panel.background = element_blank(),panel.grid = element_blank(),panel.border = element_blank(),
        axis.line.x = element_line(),axis.line.y = element_line())

ggsave(paste0(file_name_root, "Figure2", graph_file_type), fig_BICdiff, "tiff",
       height = gg_height*1.5, width = gg_width*1.5, units = "cm", dpi = dpi)


#plots terminal V
library(RColorBrewer)
library(viridis)
library(ggplot2)

nb.cols <- 54
mycolours <- colorRampPalette(brewer.pal(10,"Paired"))(nb.cols) #extends palette to 54 unique colours

#pivot longer
Exp3_dist_terminalV <- read.csv(here("output/Exp3_dist_terminalV.csv"))
Exp3_delta_terminalV <- read.csv(here("output/Exp3_delta_terminalV.csv"))

#distributed model
Exp3_distV_long <- Exp3_dist_terminalV %>%
  pivot_longer(
    col= `terminalV_0`:`terminalV_100`,
    names_prefix = "terminalV_",
    values_to = "V",
    names_to = "Nodes")

Exp3_distV_long$Nodes <- as.integer(Exp3_distV_long$Nodes)

dist_plot_point <- list()
#create distributed model plot with each participant's terminal V on each node + delta model for each participant
CueType <- unique((Exp3_distV_long$cueNum))

for (i in 1:length(CueType)){
  data = subset(Exp3_distV_long, cueNum == CueType[i])
  data2 = subset(Exp3_delta_terminalV, cueNum == CueType[i])
  Participant = as.factor(data$participant)
  Participant2 = as.factor(data2$participant)
  
  distfig <- ggplot(data, aes(x = Nodes, y = V, group = Participant, colour = Participant))+
    geom_point(size = 0.8)+scale_fill_manual(mycolours)+
    labs(title = "", x = "Nodes", y = "Terminal V")+ggtitle("", subtitle = "Distributed Model")+
    theme_classic()+
    theme(legend.position = "none", axis.line = element_line(colour = "black"),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
          plot.subtitle = element_text(size = 15),
          panel.background = element_blank())
  
  deltafig <- ggplot(data2, aes(x=Participant2, y = terminalV_plusContext, group = Participant2, colour = Participant2))+
    geom_point(size=2, alpha = 1)+scale_y_continuous(limits=c(-1,101))+
    ylab("Terminal\nValue")+xlab("Participant ID")+geom_hline(yintercept = 65, linetype = "dashed")+
    ggtitle(data$cueDescription[i],subtitle = "Simple Delta Model")+
    theme(legend.position = "none", axis.line = element_line(colour = "black"),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.text.x = element_text(colour = "black", size = 10,angle = 90),
          axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
          plot.subtitle = element_text(size = 15),
          panel.background = element_blank())
  
  fig <- ggarrange(deltafig, NULL,distfig,nrow=1, widths = c(1,0.01,1))

  dist_plot_point[[i]] <- fig
}

norm65_plot <- dist_plot_point[[1]]
posskew65_plot <- dist_plot_point[[3]]
negskew65_plot <- dist_plot_point[[5]]

individual_65_plot <- ggarrange(norm65_plot,negskew65_plot,posskew65_plot,nrow=3)

ggsave(paste0(file_name_root, "TerminalV_combine_points",graph_file_type),individual_65_plot,"tif",
       height = gg_width*2.5, width = gg_width*3, units="cm",dpi = dpi)
include_graphics(paste0(file_name_root,"TerminalV_combine_points",graph_file_type))

##average over all participants
descrip_deltaV <- Descriptives(Exp3_delta_terminalV, terminalV_plusContext,cueDescription)
kable(descrip_deltaV)
descrip_deltaV <- as.data.frame(descrip_deltaV)


descrip_distV <- Descriptives(Exp3_distV_long, V,cueDescription,Nodes)
kable(descrip_distV)
descrip_distV <- as.data.frame(descrip_distV)


test <- descrip_distV %>%
  dplyr::filter(cueDescription == "posskew65")

distfig_posskew65 <- ggplot(test, aes(x = Nodes, y = mean))+
  geom_point()+geom_line()+
  labs(title = "", x = "Nodes", y = "Average Terminal V")+ggtitle("",subtitle = "Distributed Model (Average)")+
  theme_classic()+
  theme(legend.position = "none", axis.line = element_line(colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 15),
        panel.background = element_blank())

test2 <- descrip_distV %>%
  dplyr::filter(cueDescription == "negskew65")

distfig_negskew65 <- ggplot(test2, aes(x = Nodes, y = mean))+
  geom_point()+geom_line()+
  labs(title = "", x = "Nodes", y = "Average Terminal V")+ggtitle("",subtitle = "Distributed Model (Average)")+
  theme_classic()+
  theme(legend.position = "none", axis.line = element_line(colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 15),
        panel.background = element_blank())

test3 <- descrip_distV %>%
  dplyr::filter(cueDescription == "norm65")

distfig_norm65 <- ggplot(test3, aes(x = Nodes, y = mean))+
  geom_point()+geom_line()+
  labs(title = "", x = "Nodes", y = "Average Terminal V")+ggtitle("",subtitle = "Distributed Model (Average)")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 15),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
        #rect(colour = "black", size = 1,
        #linetype = "solid", fill = NA))

combine_posskew65_plot <- ggarrange(posskew65_plot,distfig_posskew65,ncol=2,widths = c(2,1.2))
combine_negskew65_plot <- ggarrange(negskew65_plot,distfig_negskew65,ncol=2,widths = c(2,1.2))
combine_norm65_plot <- ggarrange(norm65_plot,distfig_norm65,ncol=2,widths = c(2,1.2))

combine_terminalV_plot <- ggarrange(combine_norm65_plot,combine_negskew65_plot,combine_posskew65_plot,nrow=3)

#plot is later beautified in illustrator
ggsave(paste0(file_name_root, "Figure3",graph_file_type),combine_terminalV_plot,"tif",
       height = gg_width*3, width = gg_width*3.3, units="cm",dpi = dpi)

#simulated predictions (individual example against responses)-----

#overlaying predictions
Exp3_modelPredictions <- Exp3_modelOutput %>%
  filter(trialNumPerCue > 19, !cueDescription%in%c("fixed0","fixed100"))%>%
  select(participant,cueDescription,trialNumPerCue,response,Delta_Predictions,Distributed_Predictions)

Exp3_modelPredictions <- Exp3_modelPredictions %>%
  pivot_longer(
    cols = response: Distributed_Predictions,
    names_to = "Predictions",
    values_to = "Values"
  )%>%
  mutate(Predictions = case_when(
    Predictions == "response" ~ "Participant",
    Predictions == "Delta_Predictions" ~ "Simple Delta",
    Predictions == "Distributed_Predictions" ~ "Distributed"
  ))%>%
  mutate(cueDescription = case_when(
    cueDescription=="norm65" ~ "Normal, M = 65",
    cueDescription=="norm35" ~ "Normal, M = 35",
    cueDescription=="negskew65" ~ "Negative Skew, M = 65",
    cueDescription=="negskew35" ~ "Negative Skew, M = 35",
    cueDescription=="posskew65" ~ "Positive Skew, M = 65",
    cueDescription=="posskew35" ~ "Positive Skew, M = 35"
  ))

Exp3_modelPredictions$cueDescription <- factor(Exp3_modelPredictions$cueDescription, 
                                               levels = c("Normal, M = 35","Normal, M = 65",
                                                          "Negative Skew, M = 35","Negative Skew, M = 65",
                                                          "Positive Skew, M = 35","Positive Skew, M = 65"))


Exp3_modelPredictions$Predictions <- factor(Exp3_modelPredictions$Predictions, 
                                            levels = c("Participant","Distributed","Simple Delta"))

colour_spec <- c("#00BA38","#00A5FF","red")

CueType <- levels(Exp3_modelPredictions$cueDescription)  # Ensures correct order
fitplot <- list()
for (i in 1:length(CueType)){
  
  data = subset(Exp3_modelPredictions,cueDescription == CueType[i])
  
  fig <- ggplot(data=data, aes(x=Values, group=Predictions, fill=Predictions)) +
    geom_density(adjust=1.5, alpha=.5, bw = 3) +ylim(c(0,0.045))+
    scale_fill_manual(values = colour_spec)+
    theme_minimal()+
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 13, colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14))+
    labs(title = data$cueDescription[i])
  
  fitplot[[i]] <- fig
}

Exp3_compile_fitDensity <- ggarrange(plotlist = fitplot,ncol=2,nrow = 3)

ggsave(paste0(file_name_root, "Figure4", graph_file_type), Exp3_compile_fitDensity, "tiff",
       height = gg_height*2, width = gg_width*2, units = "cm", dpi = dpi)

#pick a participant
#discrete predictions
discretePred_df <- Exp3_modelOutput%>%
  select(participant,trialNumber,cueDescription,trialNumPerCue,outcome,response,Delta_Predictions,Distributed_Predictions)%>%
  filter(participant==10 & cueDescription == "negskew65")

#pick posskew35
ppt_predictions <- discretePred_df%>%
  select(trialNumPerCue,response)

observed_outcomes <- discretePred_df%>%
  select(trialNumPerCue,outcome)

delta_predictions <- discretePred_df%>%
  select(trialNumPerCue,Delta_Predictions)

dist_predictions <- discretePred_df%>%
  select(trialNumPerCue,Distributed_Predictions)

Exp3_negskew_ppt <- ggplot(data = ppt_predictions, mapping = aes(x=response))+
  geom_histogram(data = ppt_predictions, mapping = aes(x=response, y = ..density..), binwidth = 1, fill = "#00BA38", colour = "forestgreen")+
  geom_freqpoly(data = ppt_predictions, colour = "black",size = 1, stat ="density",linetype = "dashed")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14))+  
  scale_x_continuous(limits = c(0,100))+ggtitle("Participant Response")

Exp3_negskew_outcome <- ggplot(data = observed_outcomes, mapping = aes(x=outcome))+
  geom_histogram(data = observed_outcomes, mapping = aes(x=outcome, y = ..density..), binwidth = 1, fill = "grey", colour = "black")+
  geom_freqpoly(data = observed_outcomes, colour = "black",size = 1, stat ="density",linetype = "dashed")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14))+  
  scale_x_continuous(limits = c(0,100))+ggtitle("Observed Outcomes")

Exp3_negskew_delta <- ggplot(data = delta_predictions, mapping = aes(x=Delta_Predictions))+
  geom_histogram(data = delta_predictions, mapping = aes(x=Delta_Predictions, y = ..density..), binwidth = 1, fill = "pink", colour = "red")+
  geom_freqpoly(data = delta_predictions, colour = "black",size = 1, stat ="density",linetype = "dashed")+
  xlab("Outcome Values")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14))+  
  scale_x_continuous(limits = c(0,100))+ggtitle("Simple Delta Model Predictions")

Exp3_negskew_dist <- ggplot(data = dist_predictions, mapping = aes(x=Distributed_Predictions))+
  geom_histogram(data = dist_predictions, mapping = aes(x=Distributed_Predictions, y = ..density..), binwidth = 1, fill = "#00A5FF", colour = "blue")+
  geom_freqpoly(data = dist_predictions, colour = "black",size = 1, stat ="density",linetype = "dashed")+
  xlab("Outcome Values")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(size = 12, colour = "black"), axis.title = element_text(size = 14))+  
  scale_x_continuous(limits = c(0,100))+ggtitle("Distributed Model Predictions")

combine_predictions <- ggarrange(Exp3_negskew_ppt,Exp3_negskew_outcome,Exp3_negskew_delta,Exp3_negskew_dist,ncol=2,nrow=2)

ggsave(paste0(file_name_root, "Figure5",graph_file_type),combine_predictions,"tiff",
       height = gg_width*2, width = gg_width*2, units="cm",dpi = dpi)

#figure 6 (matching correlation matrix) can be found in script_empiricalAnalysis

