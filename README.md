# ordinal_summary_measures
Tools for calculating interpretable summary measures with ordinal outcomes

The R programs in this directory are used for implementing the proposed interpretable treatment effect summary measures for ordinal outcomes. 

covid_out_functions.R contain the functions used in the applied adjusted analysis of the covid out dataset. 

simulation_functions.R contain the functions used in the simulation study for an unadjusted analysis. 

covid_out.R uses the covid_out_functions.R to run the applied analysis of the publicly available COVID-OUT dataset. To request access to the COVID-OUT dataset, see https://covidout.umn.edu/accessing-data

simulation.R can be used to replicate the simulation completed in the manuscript. 

simulation_analysis.Rmd can be used to analyze the results from the simulation study. 
