---
output: html_document
bibliography: Library.bib  
---


# **Graphical_Conditional_Extremes**

The repository contains R code used to model data according to the methods in the preprint “Conditional Extremes With Graphical Models” (Farrell, Eastoe, and Lee (2024)). Additionally, output (figures and tables) has been saved in the repository.

## **Dependencies**

The code has several dependencies throughout; however, these are loaded into each script.

## **Repository Overview**

The repository is rather large. Below, we detail the contents of each folder.

### `/threshold_selection_paper`

This folder contains files related to threshold selection [@Murphy_2024]. These files can also be found [here](https://github.com/conor-murphy4/automated_threshold_selection/tree/main/src).

### `/Miscellaneous_Functions`

- `General_Functions.R` includes commonly used functions: calculating the conditional covariance matrix for a multivariate Gaussian distribution; adding NAs into vectors and matrices; and calculating the mean absolute error and root mean squared error of samples.
- `Plotting_Functions.R` includes functions to create boxplots of parameter estimates.
- `MVAGG_Functions.R` includes the probability density function, cumulative distribution function, quantile functions, and random generation for the asymmetric generalised Gaussian (AGG) distribution [@Nacereddine_2019]. Functions to fit the AGG distribution to data are also provided. Additionally, functions for the density and random number generation of the multivariate asymmetric generalised Gaussian (MVAGG) distribution are included. Details can be found in Section 2.2 of the paper.
- `Transformations.R` includes functions to marginally transform data onto standard Laplace margins using the semi-parametric approach described in Section 3.1 of the paper.

### `/Model_Fitting`

This folder contains functions for fitting the conditional multivariate extreme value model (CMEVM) [@Heffernan_2004] and the structured CMEVM (SCMEVM) described in Section 3.1 of the paper. 

- `Cond_Extremes_MVN_Residuals.R` extends the CMEVM to allow the residual distribution to be a multivariate Gaussian distribution with independent, graphical, or saturated covariance structures.
- `Cond_Extremes_MVAGG_Residuals_One_Step.R`, `Cond_Extremes_MVAGG_Residuals_Two_Step.R`, and `Cond_Extremes_MVAGG_Residuals_Three_Step.R` contain code to fit the SCMEVM with independent, graphical, or saturated covariance structures for the one-, two-, and three-step methods, respectively.

### `Graphcial_Selection`

- `Graphical_Selection.R` includes functions to infer the graphical structure using the method described in Section 3.2 of the paper.
- `MVP_Example.R` provides the code for the simulation study in Section 4.2 of the paper.


### `/Prediction`

- `Sim_Surfaces.R` contains code to simulate data from the fitted model(s) using the method described in Section 3.3 of the paper.
- `Conditional_Probability_Calculations.R` provides functions to calculate the conditional distribution or survivor function for a given dataset.

### `/Simualtion_Study`

This folder contains all the scripts used to perform the simulation studies in Section 4 of the paper and in Sections 1 and 2 of the supplementary material. Unless otherwise stated, the folders correspond to different simulation studies with weak positive, strong positive, and weak negative associations between the components of the random vector.

- `/True_Distribution`: Contains simulation studies when data is simulated from the true distribution. Details can be found in Section 4.1 of the paper and Section 2 of the supplementary material.
- `/Timing`: Contains scripts for timing comparisons of the stepwise SCMEVMs in Section 4.1 of the paper.
- `Mixed_Data.R`: Contains a simulation study comparing the CMEVM, SCMEVMs, and the Engelke and Hitz model [@Engelke_2020] for data generated from a mixture distribution with various extremal behaviours. Details are in Section 4.3 of the paper.
- `CMEVM_Prediction.R`: Contains a simulation study comparing the predictive performance of the CMEVM and the three-step SCMEVM with graphical and saturated covariance structures. Details are in Section 1 of the supplementary material.
- `/MVN`: Contains simulation studies for data simulated from a multivariate Gaussian distribution. Details can be found in Section 3.1 of the supplementary material.
- `/MVL`: Contains simulation studies for data simulated from a symmetric multivariate Laplace distribution. Details can be found in Section 3.2 of the supplementary material.
- `/MVT`: Contains simulation studies for data simulated from a multivariate t-distribution. Details can be found in Section 3.3 of the supplementary material.
- `MVP.R`: Contains a simulation study for data simulated from a multivariate Pareto distribution. Details are in Section 3.4 of the supplementary material.

### `/Data`

This folder contains the data, on their original and standard Laplace margins, used in the simulation studies mentioned above.

### `/Danube_River`

This folder contains files related to the analysis of extreme river discharges in the upper Danube River basin [@Asadi_2015].

- `Danube_EDA.R`: Performs exploratory data analysis.
- `Danube_Model_Comparison_Bootstrapped.R`: Includes all other analyses on the data, such as bootstrapping the data, inferring the graph using the bootstrapped samples, fitting various models to the data, and creating comparative plots of the predictive performance of the models.

### `/Images`

The Images folder contains several subfolders related to the simulation studies and the analysis of the upper Danube River discharge dataset.

Below is a summary table that details the file used to obtain the output and a link to the either the file or output. 

| Figure/Table | Section | Code to generate Figure/Table | Link to Output | 
| :-----: | :-----: | :-----: | :-----: |
| Figure 1 <br>(left panel) | 1 | [`Danube_River_Analysis/Danube_Model_Comparison_Bootstrapped.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Danube_River_Analysis/Danube_Model_Comparison_Bootstrapped.R) | [Figure](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Images/Danube/Bootstrapped_Ouput/Danube_River.pdf) | 
| Figure 1 <br>(right panel) | 1 | [`Danube_River_Analysis/Danube_EDA.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Danube_River_Analysis/Danube_EDA.R) | [Figure](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Images/Danube/Inference_On_Data/EDA.pdf) |
| Figure 2 | 4.1 | [`Simulation_Study/True_Distribution/True_Dist_Low_Dependence.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Simulation_Study/True_Distribution/True_Dist_Low_Dependence.R) | [Figure](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Images/Simulation_Study/True_Distribution/Low_Dependence/Alpha.pdf) | 
| Table 1 | 4.1 | [`Simulation_Study/True_Distribution/True_Dist_Low_Dependence.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Simulation_Study/True_Distribution/True_Dist_Low_Dependence.R) | N/A | 
| Table 2 | 4.1 | [`Simulation_Study/Timing/Model_Fitting/Timing_Large_Sparse_Graphs.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Simulation_Study/Timing/Model_Fitting/Timing_Large_Sparse_Graphs.R) | N/A | 
| Figure 3 | 4.2 | [`Graphical_Selection/MVP_Example.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Graphical_Selection/MVP_Example.R) | [Folder](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/tree/main/Images/Simulation_Study/Graphical_Selection) | 
| Figure 4 | 4.3 | [`Simulation_Study/Mixed_Data.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Simulation_Study/Mixed_Data.R) | [Figure](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Images/Simulation_Study/Mixed_Data/Gamma.pdf) | 
| Figure 5 <br>(left panel) | 4.3 | [`Simulation_Study/Mixed_Data.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Simulation_Study/Mixed_Data.R) | [Figure](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Images/Simulation_Study/Mixed_Data/Probabilities/Site_3/Prob_5.pdf) |
| Figure 5 <br>(right panel) | 4.3 | [`Simulation_Study/Mixed_Data.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Simulation_Study/Mixed_Data.R) | [Figure](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Images/Simulation_Study/Mixed_Data/Probabilities/Site_3/Prob_6.pdf) |
| Figure 6 <br>(left panel) | 4.1 | [`Simulation_Study/Timing/Time_plot.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Simulation_Study/Timing/Time_plot.R) | [Figure](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Images/Simulation_Study/Time_Graph.pdf) | 
| Figure 6 <br>(right panel) | 5 | [`Danube_River_Analysis/Danube_Model_Comparison_Bootstrapped.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Danube_River_Analysis/Danube_Model_Comparison_Bootstrapped.R) | [Figure](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Images/Danube/Bootstrapped_Ouput/Danube_Inferred_Glasso_Residuals.pdf) | 
| Figure 7 | 5 | [`Danube_River_Analysis/Danube_Model_Comparison_Bootstrapped.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Danube_River_Analysis/Danube_Model_Comparison_Bootstrapped.R) | [Figure](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Images/Danube/Bootstrapped_Ouput/Chi_Comp.pdf) | 
| Figure 8 | 5 | [`Danube_River_Analysis/Danube_Model_Comparison_Bootstrapped.R`](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Danube_River_Analysis/Danube_Model_Comparison_Bootstrapped.R) | [Figure](https://github.com/AidenFarrell/Graphical_Conditional_Extremes/blob/main/Images/Danube/Bootstrapped_Ouput/Bias_ETA_80_CMEVM_Graph.pdf) | 

# **References**

