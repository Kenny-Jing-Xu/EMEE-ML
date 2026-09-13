# EMEE-ML

The R code for paper "Causal Excursion Effects Estimation with Ordinal Outcomes under Micro-Randomized Trials" by Jing Xu, Tianchen Qian, Tanujit Chakraborty, Saeed AL Dahmani, Adrian Aguilera, Bibhas Chakraborty.

## R code files

- Simulation_EMEE-ML.R: generate the simulation results in Table 4 when K=6 and q=101
- Data_MRT_MC_Ord.R: generate the Monte Carlo MRT data that includes the multiple-category interventions and the ordinal proximal outcomes
- EMEE_mc_ord.R: estimate the causal effect using the EMEE method with multiple-category interventions and the ordinal proximal outcomes
- EMEE_ML.R: estimate the causal effect using the EMEE-ML method with multiple-category interventions and the ordinal proximal outcomes

## How to use the code to generate the results 

- For result of K=6 and q=101 in Table 4 in the "Simulation Study" section, run the R scripts in the Simulation_EMEE-ML.R file.




