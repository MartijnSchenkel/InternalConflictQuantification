# ================================================================================================ #

Schenkel, M.A., Patten, M.M. & Ågren, J.A. (2025). Internal conflicts and the measurement of 
 evolutionary individuality. Journal TBD.

# ================================================================================================ #

CONTACT: Martijn A. Schenkel
EMAIL: maschenkel@gmail.com; martijn.schenkel@wur.nl

# ================================================================================================ # 

/figures/
  Contains PDF and PNG versions of Figure 3 and Supplementary Figures 1 and 2, which are produced
  using the analysis_script.R file. 

/source_code/
  Source code for the model. Is used to solve differential equations for fitness using a bisection
  method (see Supplementary Material for full information). 

/analysis_script.R
  Script for data analysis. Reads in the data.txt file in data.zip (need to unpack this first!), and 
  produces Figure 3 and Supplementary Figures 1 and 2. Can be run in R/RStudio, but might require
  manual installation of several packages (as per usual).

/data.zip
  Zip archive containing the (compressed) data set used for the analysis of the genomic imprinting
  model. Variables are as follows:
  
  iteration	= Iteration at which function converged on 0. Values of 10002 (higher than i_max) 
		  represent that the function never fully converged on 0 due to monotonicity.
  i_max		= How many iterations of the bisection method can be performed at most. 
  monotonic	= Fitness function increases monotonically (yes = 1, no = 0). If yes, optimal trait
		  value is equal to all-in investment (i.e., invest R resources).
  R		= Amount of resources available for reproduction
  r		= Relatedness factor used for calculating inclusive fitness
  c		= Scaling factor for survival as function of investment. Higher values result in 
		  larger increases in survival for each additional unit of resources invested.
  k		= Maximum survival rate, unused in model (always set to 1).
  z_0		= Minimal investment rate below which survival is 0. Values above z_0 result in positive
		  survival rates.
  z_m		= Solution for dw/dz = 0, based on remaining parameters. Equivalent to optimal investment
		  of resources.
  w_pat		= Relative fitness of paternal genome at z_m.
  w_mat		= Relative fitness of maternal genome at z_m.
  w_org		= Relative fitness of organism at z_m.