# ================================================================================================ #

Schenkel, M.A., Patten, M.M. & Ågren, J.A. (2025). Internal conflicts and the measurement of 
 evolutionary individuality. Journal TBD.

# ================================================================================================ #

CONTACT: Martijn A. Schenkel
EMAIL: maschenkel@gmail.com; martijn.schenkel@wur.nl

# ================================================================================================ # 

/figures/
  Contains PDF and PNG versions of Figure 4 and Supplementary Figures 3 and 4, which are produced
  using the analysis_script.R file. 

/source_code/
  Source code for the model, implemented in C++. Is used to solve differential equations for fitness 
  using a bisection method (see Supplementary Material for full information).

/analysis_script.R
  Script for data analysis. Reads in the data.txt file in data.zip (need to unpack this first!), and 
  produces Figure 3 and Supplementary Figures 3 and 4. Can be run in R/RStudio, but might require
  manual installation of several packages (as per usual). Note that there is a lengthy loop in this
  script to calculate values at a large number of values for p (fraction of genome that is sex-
  chromosomal); it might be advantageous to save the resulting data-frames (d4z_... and d4w_... 
  variables in the script).

/data.zip
  Zip archive containing the (compressed) data set used for the analysis of the genomic imprinting
  model. Variables are as follows:

  n		= Number of foundresses per deme
  zI		= Optimal sex ratio for individual (z_I^*)
  zG_X_XX	= Optimal sex ratio for X chromosome in XX individuals
  zG_X_XY	= Optimal sex ratio for X chromosome in XY individuals
  zG_Y_XY	= Optimal sex ratio for Y chromosome in XY individuals
  zG_Z_ZZ	= Optimal sex ratio for Z chromosome in ZZ individuals
  zG_Z_ZW	= Optimal sex ratio for Z chromosome in ZW individuals
  zG_W_ZW	= Optimal sex ratio for W chromosome in ZW individuals
  wI_star	= Fitness of individual using strategy z_I (equivalent to w_I^*)
  wG_X_XX	= Fitness of individual using strategy z_G_X_XX 
  wG_X_XY	= Fitness of individual using strategy z_G_X_XY
  wG_Y_XY	= Fitness of individual using strategy z_G_Y_XY
  wG_Z_ZZ	= Fitness of individual using strategy z_G_Z_ZZ
  wG_Z_ZW	= Fitness of individual using strategy z_G_Z_ZW
  wG_W_ZW	= Fitness of individual using strategy z_G_W_ZW