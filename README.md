# Distributed Model Contingency Learning
Code and Data to accompany Chow et al: "Predicting continuous outcomes: Some new tests of associative approaches to contingency learning."

In this paper, we introduce one approach to incorporating a distributed architecture into a prediction error model that tracks the contingency between cues and dimensional outcomes. Our Distributed Model allows associative links to form between the cue and outcome nodes that provide distributed representation depending on the magnitude of the outcome, thus enabling learning that extends beyond approximating the mean (Simple Delta Model). We provide a formalisation of both models, and fit empirical data to each model. Best fitting parameters are then used for formal model comparison using Bayesian Information Criterion (BIC).

This repository contains:<br>
1) Empirical data from 4 Experiments used for model fit (folders: original data)<br>
2) R scripts: Model functions, model fits (Experiments 1 -4), formal model comparison, analysis of empirical data, and code to reproduce figures in the paper.<br>
3) Relevant data output as csv files (folders: output)