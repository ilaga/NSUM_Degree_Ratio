# Instructions

Initial Setup: Set the working directory to the top-level of the folder where all R files are contained.

Organization:
	Supplementary Material:
		Data:
			(a) (Need to add, see Data Availability and instructions below) twomwth2.csv: ARD from the McCarty dataset
			(b) (Need to add, see Data Availability and instructions below) rwanda-varnames.csv: Rwanda meal study ARD subpopulation question IDs, known population sizes, and descriptions
			(c) (Need to add, see Data Availability and instructions below) RWIQ6AFL.SAS7BDAT: Survey responses from Rwanda meal study
		Figures:
			Location to save all figures presented in the manuscript
		Latent_surface_code:
			Code provided and adapted from Breza et al. (2020), "Using aggregated relational data to feasibly identify network structure without network data"
			(a) main.R: Main function to take input and return results from Latent Surface Models
			(b) latent_surface_code_MH.R: Updated Metropolis-Hastings adaptation of the original Breza et al. (2020) paper that performs the MCMC for the Latent Surface Model
		Results:
			Location to save all MCMC results

Software:
	All analysis was performed in R version 4.2.3. Packages used for analysis include (and versions in parentheses):
		(a) ggplot2 (3.4.2)
		(b) haven (2.5.2)
		(c) igraph (1.4.2)
		(d) Metrics (0.1.4)
		(e) movMF (0.2-7)
		(f) reshape2 (1.4.4)

Analysis Procedure:
	Step 1. Data Preparation
		Details: Process the Rwanda Meal survey to be usable by network scale-up models
		(a) Prepare_q_vs_q_data.R
			Reads in the Rwanda survey response .SAS7BDAT file, and subsets it to two files:
			1. qvsq_meal.RDS: Aggregated relational data responses for the meal question form
			2. qvsq_standard.RDS: Aggregated relational data responses for the standard question form (not used in this paper)
			
	Step 2. Fit Latent Surface Models
		Details: Fit the latent space models to the Rwanda, McCarty, and simulated stochastic block model data.
		(a) Fit_SBM_latent_surface.R
			Fits the latent surface model to the simulated stochastic block model data and saves the results to SBM_LS.RData
		(b) Fit_McCarty_latent_surface.R
			Fits the latent surface model to the simulated stochastic block model data and saves the results to McCarty_LS.RData
		(c) Fit_Rwanda_meal_latent_surface.R
			Fits the latent surface model to the simulated stochastic block model data and saves the results to Rwanda_meal_latent_fit.RData
			
	Step 3. Analyze Results
		Details: Analyze the output from the fitted models and generate all result plots and tables in the manuscript
		(a) Results_SBM_scaling.R
			Fits the latent surface model to the simulated stochastic block model data and saves the results to SBM_LS.RData
		(b) Results_McCarty_scaling.R
			Fits the latent surface model to the McCarty data and saves the results to McCarty_LS.RData
		(c) Results_Rwanda_meal_scaling.R
			Fits the latent surface model to the Rwanda meal data and saves the results to Rwanda_meal_latent_fit.RData

Data Availability and Instructions:
	We do not have permission to distribute the McCarty ARD dataset and the Rwanda ARD dataset.
	The McCarty ARD dataset may be obtained directly by contacting Christopher McCarty (ufchris@ufl.edu)
	The Rwanda ARD dataset can be downloaded from https://dhsprogram.com/data/dataset/Rwanda_Special_2011.cfm?flag=1)[https://www.measuredhs.com/data/dataset/Rwanda_Special_2011.cfm&flag=1)[https://www.measuredhs.com/data/dataset/Rwanda_Special_2011.cfm
		The ARD dataset used in this analysis are found in the RWIQ6ASD.ZIP SAS dataset
		The rwanda-varnames.csv file can be downloaded from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/CCC6HF
