# Distributed Model Contingency Learning
Code and Data to accompany Chow et al: "Predicting continuous outcomes: Some new tests of associative approaches to contingency learning."

In this paper, we introduce one approach to incorporating a distributed architecture into a prediction error model that tracks the contingency between cues and dimensional outcomes. Our Distributed Model allows associative links to form between the cue and outcome nodes that provide distributed representation depending on the magnitude of the outcome, thus enabling learning that extends beyond approximating the mean (Simple Delta Model). We provide a formalisation of both models, and fit empirical data to each model. Best fitting parameters are then used for formal model comparison using Bayesian Information Criterion (BIC).

This repository contains:<br>
1) Empirical data from 4 Experiments used for model fit (folders: original data)<br>
2) R scripts: Model functions, model fits (Experiments 1 -4), formal model comparison, analysis of empirical data, and code to reproduce figures in the paper.<br>
3) Relevant data output as csv files (folders: output)

Software Specifications:<br>
R scripts in this repository were written and run on RStudio Version: 2023.03.1+446 (R Version: 4.0.0 or higher recommended)<br>
Operating System: Windows, macOS (should also work with Linux)

Installation Instructions<br>
1) Install R: Download and install R from [CRAN](https://cran.r-project.org/)<br>
2) Install [RStudio](https://www.rstudio.com/products/rstudio/download/)<br>
3) Clone or download this repository:
Option A Using Git (command line)
```bash
   git clone https://github.com/julieechow/Distributed-Model-Contingency-Learning.git
   cd Distributed-Model-Contingency-Learning
```
Option B Download ZIP (no Git required)
    Click the green <>Code button and select "Download ZIP"
    Extract the ZIP file to your desired location
    Open extracted folder in RStudio and navigate to it in your file browser

Required R Packages<br>
The required R packages to run each script is reported at the start of each R script.<br>
To install the relevant packages, use install.packages("name_of_package")<br>
To call a package after installation, use library(name_of_package)


Project Structure<br>
.
├── ModelFunctionsLL.R                    # Core model fitting function called in multiple scripts
├── ModelFitsAll_Single_MLE_Exp1.R        # Model fitting for Experiment 1
├── ModelFitsAll_Single_MLE_Exp2.R        # Model fitting for Experiment 2
├── ModelFitsAll_Single_MLE_Exp3.R        # Model fitting for Experiment 3
├── ModelFitsAll_Single_MLE_Exp4.R        # Model fitting for Experiment 4
├── script_modelComparison.R              # Model comparison and BIC analysis - takes output from preceding model fits (saved in output/ folder)
├── plot_model_output.R                   # Generate manuscript figures
├── generate_predictions.R                # Generate simulated datasets (uninformed) for model recovery analysis
├── generate_predictions_informed.R       # Generate simulated datasets (informed) for model and parameter recovery analysis
├── modelRecovery.R                       # Model recovery analysis
├── paramRecovery.R                       # Parameter recovery analysis
├── script_empiricalAnalysis.R            # Empirical data analysis
├── original_data/                        # Empirical data used for model fit and empirical analysis are located here
├── simulated_data/                       # simulated datasets required for model and parameter recovery (output files from generate_predictions[_informed.R] scripts)
└── output/                               # output files from all scripts except generate_predictions.R and generate_predictions_informed.R


Workflow Summary<br>
Before running scripts, please ensure you unzip folders.zip and extract the three data files (original_data, simulated_data and output)<br>
Place all three folders in the same environment as R scripts.<br>
Note: All scripts save their output to the output/ folder (simulated data from generate_predictions and generate_predictions_informed are saved to simulated_data/). This allows you to run scripts out of order or in isolation, as each script reads from and writes to this shared output directory.


For complete reproduction of manuscript results, run scripts in this order:

1) Model Fitting: ModelFitsAll_Single_MLE_Exp[1-4].R <br>
2) Model Comparison: script_modelComparison.R <br>
3) Visualization: plot_model_output.R <br>

Supplementary Materials:
4) Empirical Analysis: script_empiricalAnalysis.R <br>
5) Simulated Data Generation: generate_predictions.R and generate_predictions_informed.R <br>
6) Model Recovery Analysis: modelRecovery.R <br>
7) Parameter Recovery Analysis: paramRecovery.R 



