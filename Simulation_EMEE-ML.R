# Title: "Simulation Study for EMEE-ML Method"
# Author: "Jing Xu"
# Data: "07/09/2026"

## Package installation
install.packages("foreach")
install.packages("rootSolve")
install.packages("doParallel")
install.packages("parallel")
install.packages("Formula")
install.packages("earth")
install.packages("BayesTree") # The original BART implementation
install.packages("dbarts") # The newer, faster implementation (recommended)
install.packages("bartMachine") # Alternative: BART via the bartMachine package
install.packages("truncnorm")
install.packages("brnn")
install.packages("brant")
install.packages("plotmo")
install.packages("plotrix")
install.packages("glmnet")
install.packages("glm2")
install.packages("readxl")


## Managing big dataset
options( max.print = 1000000 ) 
memory.size(1000000)

rm(list = ls()) 

## Load R packages
library(foreach)
library(rootSolve)
library(doParallel)
library(parallel)
library(Formula)
library(earth)
library(BayesTree)
library(dbarts)
library(bartMachine)
library(truncnorm)
library(brnn)
library(brant)
library(plotmo)
library(plotrix)
library(glmnet)
library(glm2)
library(readxl)


## Directory
setwd("~/Dropbox/Ordinal MRT/Git")

## R function
source("EMEE_mc_ord.R")
source("EMEE_ML.R")
source("Data_MRT_MC_Ord.R")




## Data generation
NoTrial <- 500 # The total number of simulation repetitions 
N <- 50 # Sample size
T.tot <- 60 # The total number of decision time points
M <- 2 # The number of intervention categories
K <- 6 # The number of ordinal outcome levels without the reference level
p <- 2  # The number of causal effect parameters include the intercept for each intervention category
q <- 101 # The number of nuisance parameters include the intercenpt for each ordinal outcome level

# Randomization probabilities (uniform)
prob_true <- matrix(rep(1/(M+1), (M+1) * T.tot), nrow = T.tot, ncol = M + 1) # Randomization probability of intervention message matrix
tau_t <- 1  # Availability rate at time t equals one

# Alpha_true: (q x K) matrix for H_it
alpha_true <- matrix(0, nrow = q, ncol = K) # True value of alpha: nuisance parameter
alpha_true[1, 1] <- log(6/7)   
alpha_true[1, 2] <- log(5/7)
alpha_true[1, 3] <- log(4/7)   
alpha_true[1, 4] <- log(3/7)   
alpha_true[1, 5] <- log(2/7)   
alpha_true[1, 6] <- log(1/7) 
alpha_true[2:q, ] <- 0.01  # small effect for covariates


# Beta_true: (p x M) matrix for S_it
beta_true <- matrix(0, nrow = p, ncol = M) # True value of beta: causal effect size parameter
beta_true[1, 1] <- 0.1   # intercept for treatment 1
beta_true[1, 2] <- -0.1  # intercept for treatment 2
beta_true[2:p, ] <- 0.05 # small effects for covariates

set.seed(123)
sim_H_A_Y <- foreach( s = 1:NoTrial ) %do% simulate_mrt_data(
    s, # Simulation repetition s
    N, # Sample size
    T.tot, # Number of decision time points
    alpha_true, # True value of alpha: nuisance parameter
    beta_true, # True value of beta: causal effect size parameter
    prob_true, # Randomization probability of intervention message matrix
    tau_t # Availability rate at time t
)
dim(sim_H_A_Y[[1]])
names(sim_H_A_Y[[1]])
colMeans(sim_H_A_Y[[1]]["Yit_1"])
colMeans(sim_H_A_Y[[1]]["Yit_2"])
colMeans(sim_H_A_Y[[1]]["Yit_3"])
colMeans(sim_H_A_Y[[1]]["Yit_4"])
colMeans(sim_H_A_Y[[1]]["Yit_5"])
colMeans(sim_H_A_Y[[1]]["Yit_6"])

## Parameter estimation
#debug( fit_EMEE_ML )
#undebug( fit_EMEE_ML )
est_beta <- array( NA, dim = c( dim(beta_true)[1], M, NoTrial ) )
se_beta <- array( NA, dim = c( dim(beta_true)[1], M, NoTrial ) )
for( s in 1:NoTrial)
{
  EMEE_sim <- fit_EMEE_ML( 
    dta = sim_H_A_Y[[s]], # Data used for analysis
    id_varname = "ID", # Variable name for ID
    decision_time_varname = "DecPt", # Variable name for decision time point
    treatment_varname = c("Ait_1", "Ait_2"), # Variable name for treatment
    outcome_varname = c( paste("Yit_", 1:6, sep = "") ), # Variable name for outcome
    nuisance_method = "earth", # Estimation method for the nuisance parameters
    control_varname = c( paste("x", 1:100, sep = "") ), # Variable names for the nuisance parameters Ht
    control_varname2 = c( "DecPt" ), # Variable names for the nuisance part when the 'gam' method is used
    moderator_varname = c( "x1" ), # Variable names for the moderator St
    rand_prob_varname = paste("prob_it", 1:M, sep = ""), # Variable names for the random probabilities P(Ht)
    avail_varname = NULL, # Variable name for the availability indicator It
    rand_prob_tilde_varname = NULL, # Variable names for the random probabilities \tilde{P}(S_t) 
    rand_prob_tilde = (table(sim_H_A_Y[[s]]$A)/sum( table(sim_H_A_Y[[s]]$A) ))[-1], # Random probabilities \tilde{P}(S_t) 
    estimator_initial_value = NULL # Initial value for the estimator
  )
  est_beta[,,s] <- EMEE_sim$beta_hat
  se_beta[,,s] <- EMEE_sim$beta_se
  print(s)
}
UCI95 <- est_beta + 1.96*se_beta
LCI95 <- est_beta - 1.96*se_beta

beta_est_95CP <- matrix( NA, p, M )
for(m in 1:M){
  for(j in 1:p){
    beta_est_95CP[ j, m ] <- mean( ifelse( beta_true[j,m]>=LCI95[j,m,] & beta_true[j,m]<=UCI95[j,m,], 1, 0 ) )
  }
}

apply(est_beta, c(1,2), mean)
sqrt( apply(est_beta, c(1,2), var) )
apply(se_beta, c(1,2), mean)
apply(est_beta, c(1,2), mean)-beta_true

Table_sim <- apply(
    cbind(
      c(
        1, rep("",p-1),
        2, rep("",p-1)
      ),
      round(
        cbind(
          matrix(beta_true, nrow = M*p, ncol = 1),
          matrix(apply(est_beta, c(1,2), mean), nrow = M*p, ncol = 1),
          matrix(apply(se_beta, c(1,2), mean), nrow = M*p, ncol = 1),
          matrix(sqrt( apply(est_beta, c(1,2), var) ), nrow = M*p, ncol = 1),
          matrix(beta_est_95CP, nrow = M*p, ncol = 1)
        ),
        digits = 3
      )
    ), 
    2, as.numeric
  )
colnames(Table_sim) <- c( "m", "True beta_m", "Est", "SE", "SD", "CP" )
Table_sim

## Normality test for the parameter estimates
hist( est_beta[ 1, 1, ], xlab = expression(hat(beta)[11]), main = "" )
qqnorm( est_beta[ 1, 1, ], main = "", )
qqline( est_beta[ 1, 1, ] )
shapiro.test(est_beta[ 1, 1, ] )$p.value




# delete the plots
dev.off()

