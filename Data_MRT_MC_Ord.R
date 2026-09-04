simulate_mrt_data <- function(
    s, #Simulation repetition s
    N, #sample size
    T.tot, #number of decision time points
    alpha_true, #true value of alpha: nuisance parameter
    beta_true, #true value of beta: causal effect size parameter
    prob_true, #randomization probability of intervention message matrix
    tau_t #availabity rate at time t
){
  #set.seed(j)
  
  #row number
  Total <- N*( T.tot ) #total number of time points of each sample
  #participant id
  ID <- rep( 1:N, each = T.tot )
  #decision time point
  DecPt <- rep( 1:T.tot, N )
  DecPt_1 <- DecPt - 1
  DecPt_1_13max <- pmin( DecPt_1 , 13 )
  
  #number of intervention categories and ordinal outcome level except for the reference
  if(is.vector(beta_true)) beta_true <- matrix(beta_true, ncol = 1)
  if(is.vector(alpha_true)) alpha_true <- matrix(alpha_true, ncol = 1)
  M <- dim(beta_true)[2]
  K <- dim(alpha_true)[2]
  
  #randomization probability 
  prob_it <- matrix( 0,  Total, dim( prob_true )[2] )
  prob_it[ , 1 ] <- rep( prob_true[ , 1 ], N )
  for( m in 1:( dim( prob_true )[2] - 1 ) )
  {
    prob_it[ , ( m + 1) ] <- rep( prob_true[ , ( m + 1) ], N )
  }
  colnames(prob_it) <- paste("prob_it", 0:M, sep = "")
  
  #intervention categories
  #A_idt <- matrix( 0, Total, dim( prob )[2] - 1 )
  #for ( i in 1:Total ){ A_idt[ i, ] <- t( rmultinom( 1, size = 1, prob = prob_idt[ i, ] )[ -1,  ] ) }
  A <- sample( 0:M, Total, replace = T )
  Ait <- matrix( 0, Total, M )
  for( m in 1:M ){ Ait[ A %in% c( m ), m ] <- 1 }
  colnames(Ait) <- paste("Ait_", 1:M, sep = "")
  
  #availability
  u1 = runif( Total )
  It = rep( 0, Total );
  It[ u1 <= tau_t ] = 1
  
  X <- matrix(NA, nrow = N*T.tot, ncol = 100 )
  for(j in 1:100){
    X[ , j ] <- rep( rnorm( n = N, mean = 0, sd = 1 ), each = T.tot ) 
  }
  Hit <-cbind( X ) 
  colnames(Hit) <- c( paste("x", 1:100, sep = "") )
  g_Hit <- cbind( rep(1, Total), Hit )
  Sit <- cbind( X[,c( 1 )] )
  f_Sit <- cbind( rep(1, Total), Sit )
  
  # EXPECTED VALUES for each threshold
  EY <- matrix(NA, nrow=Total, ncol=K)
  for (k in 1:K) {
    log_EYit_Ait <- g_Hit%*%alpha_true[ , k ] + rowSums( cbind( Ait*( ( f_Sit )%*%beta_true ) ) )
    EYit_Ait <- exp( log_EYit_Ait )
    EY[,k] <- EYit_Ait
  }
  
  # MONOTONICITY ENFORCEMENT
  epsilon <- 1e-6
  EY[,1] <- pmin(EY[,1], 0.99)
  for (k in 2:K) { EY[,k] <- pmin(EY[,k], EY[,k-1] - epsilon) }
  
  # CATEGORY PROBABILITIES for each ordinal level
  p_mat <- matrix( NA, nrow = Total, ncol = K+1 )
  p_mat[ , 1 ] <- 1 - EY[ , 1 ]
  for ( k in 2:K ) { p_mat[ , k ] <- EY[ , k-1 ] - EY[ , k ] }
  p_mat[ , K+1 ] <- EY[ ,K ]
  p_mat <- p_mat / rowSums( p_mat )
  
  # SIMULATE ORDINAL OUTCOME Y_t
  Yit <- apply( p_mat, 1, function( p ) sample( 0:K, 1, prob=p ) )
  Yit_lvl <- matrix( 0, nrow = Total, ncol = K )
  for( k in 1:K ){
    Yit_lvl[ Yit>=k, k ] <- 1 
  }
  colnames(Yit_lvl) <- paste("Yit_", 1:K, sep = "")
  
  # FINAL DATASET
  mrt_data <- data.frame( ID, DecPt, DecPt_1, Hit, A, Ait, prob_it, Yit, Yit_lvl )
  
  print(s)
  return(mrt_data)
}










