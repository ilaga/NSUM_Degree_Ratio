# Instructions

## Initial Setup
Set the working directory to the top-level of the folder where all R files are contained.

## Organization
- Supplementary Material:
	- Data:
		- (a) (Need to add, see Data Availability and Instructions below) twomwth2.csv: ARD from the McCarty dataset
		- (b) (Need to add, see Data Availability and Instructions below) rwanda-varnames.csv: Rwanda meal study ARD subpopulation question IDs, known population sizes, and descriptions
		- (c) (Need to add, see Data Availability and Instructions below) RWIQ6AFL.SAS7BDAT: Survey responses from Rwanda meal study
	- Figures:
		Location to save all figures presented in the manuscript
	- Latent_surface_code:
		Code provided and adapted from Breza et al. (2020), "Using aggregated relational data to feasibly identify network structure without network data"
		- (a) main.R: Main function to take input and return results from Latent Surface Models
		- (b) latent_surface_code_MH.R: Updated Metropolis-Hastings adaptation of the original Breza et al. (2020) paper that performs the MCMC for the Latent Surface Model
	- Results:
		Location to save all MCMC results

## Software
All analysis was performed in R version 4.2.3. Packages used for analysis include (and versions in parentheses):
- ggplot2 (3.4.2)
- haven (2.5.2)
- igraph (1.4.2)
- Metrics (0.1.4)
- reshape2 (1.4.4)
- ggpmisc (0.5.4-1)
- dplyr (1.1.3)

## Analysis Procedure
### Step 1. Data Preparation
Details: Process the Rwanda Meal survey to be usable by network scale-up models
- (a) Prepare_q_vs_q_data.R
	Reads in the Rwanda survey response .SAS7BDAT file, and subsets it to two files:
	1. qvsq_meal.RDS: Aggregated relational data responses for the meal question form
	2. qvsq_standard.RDS: Aggregated relational data responses for the standard question form (not used in this paper)
	
### Step 2. Fit Latent Surface Models
Details: Implement Algorithm 1 and present results for the Rwanda, McCarty, and simulated binomial and stochastic block model data.
- (a) Slope_Binomial.R
	Implements algorithm 1 with the simulated binomial model data and outputs results and figures
- (b) Slope_SBM.R
	Implements algorithm 1 with the simulated stochastic block model data and outputs results and figures
- (c) Fit_McCarty_latent_surface.R
	Implements algorithm 1 with the McCarty data and outputs results and figures
- (d) Fit_Rwanda_meal_latent_surface.R
  	Implements algorithm 1 with the Rwanda Meal data and outputs results and figures
- (e) Slope_Binomial_different_p.R
  	Implements algorithm 1 with the simulated binomial model with different p data and outputs results and figures


## Data Availability and Instructions:
- We do not have permission to distribute the McCarty ARD dataset and the Rwanda ARD dataset.
- The McCarty ARD dataset may be obtained directly by contacting Christopher McCarty (ufchris@ufl.edu)
- The Rwanda ARD dataset can be downloaded from https://dhsprogram.com/data/dataset/Rwanda_Special_2011.cfm?flag=1)[https://www.measuredhs.com/data/dataset/Rwanda_Special_2011.cfm&flag=1)[https://www.measuredhs.com/data/dataset/Rwanda_Special_2011.cfm
	- The ARD dataset used in this analysis are found in the RWIQ6ASD.ZIP SAS dataset
	- The rwanda-varnames.csv file can be downloaded from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/CCC6HF
