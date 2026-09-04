library(rootSolve) # for solver function multiroot()
library(geepack) # for fitting GEE using package

get_alpha_beta_from_multiroot_result <- function(root, p, q, M, K)
{
    beta_root <- matrix( NA, nrow = p, ncol = M )
    for ( j in 1:M ) {
      beta_root[ , j ] <- root$root[ ( q*K + 1 + ( j - 1 )*p ) : ( q*K + j*p ) ]
    }
    
    alpha_root <- matrix( root$root[ 1:( q*K ) ], nrow = q, ncol = K )

  return( list( alpha = alpha_root, beta = beta_root ) )
}

find_change_location <- function( v ){
  n <- length( v )
  if ( n <= 1 ) {
    stop( "The vector need to have length > 1." )
  }
  return( c( 1, 1 + which( v[ 1:( n - 1 ) ] != v[ 2:n ] ) ) )
}

fit_EMEE <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
    rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
    rand_prob_tilde = NULL,    
    estimator_initial_value = NULL
)
{
  sample_size <- length( unique( dta[ , id_varname ] ) )
  total_person_decisionpoint <- nrow( dta )
  
  if ( is.null( avail_varname ) ) { avail <- rep( 1, total_person_decisionpoint ) } 
  else { avail <- dta[ , avail_varname ] }
  
  A <- dta[ , treatment_varname ]
  A <- as.matrix( A )
  M <- dim( A )[ 2 ] - 0
  # checking for NA in treatment indicator    
  if ( any( is.na( A[ avail == 1, ] ) ) ) { stop( "Treatment indicator is NA where availability = 1." ) }
  A[ avail == 0 ] <- 0
 
  p_t <- dta[ , rand_prob_varname ]
  p_t <- as.matrix( p_t )
  cA <- A - p_t # centered A
  
  #ordinal outcome
  Y <- as.matrix(dta[ , outcome_varname ])
  K <- dim( Y )[ 2 ]
  
  Xdm <- as.matrix( cbind( rep( 1, nrow( dta ) ), dta[ , moderator_varname ] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep( 1, nrow( dta ) ), dta[ , control_varname ] ) ) # Z (control) design matrix, intercept added
  
  if (length(rand_prob_tilde) == M) { 
    p_t_tilde <- matrix( rep( rand_prob_tilde, each=total_person_decisionpoint ), nrow = total_person_decisionpoint, ncol = M ) 
  } else { 
      p_t_tilde <- dta[, rand_prob_tilde_varname] 
      }
  p_t_tilde <- as.matrix( p_t_tilde )
  cA_tilde <- A - p_t_tilde
  
  WCLS_weight <- rowSums( A*p_t_tilde/p_t ) + ( 1 - rowSums( A ) )*( 1 - rowSums( p_t_tilde ) )/( 1 - rowSums( p_t ) )
  WCLS_weight <- as.vector( WCLS_weight )
  
  p <- dim( Xdm )[2] #length(moderator_varname) + 1 # dimension of beta
  q <- dim( Zdm )[2] #length(control_varname) + 1 # dimension of alpha
  
  Xnames <- c( "Intercept", moderator_varname)
  Znames <- c( "Intercept", control_varname )
  
  ### estimation ###
  estimating_equation <- function( theta ) {
    alpha <- matrix( theta[ 1:( q*K ) ], nrow = q, ncol = K )
    beta <- matrix( theta[ ( q*K + 1 ):( q*K + p*M ) ], nrow = p, ncol = M )
    
    exp_Zdm_alpha <- exp( Zdm %*% alpha )
    exp_AXdm_beta <- exp( rowSums( A*( Xdm%*%beta ) ) )
    exp_AXdm_beta_K <- matrix( exp_AXdm_beta, nrow = total_person_decisionpoint, ncol = K )
    residual <- Y - exp_Zdm_alpha*exp_AXdm_beta_K
    residual <- as.matrix( residual )
    weight <- exp_AXdm_beta^( -1 ) #inverse_factor
    
    ef <- rep( NA, K*q+M*p ) # value of estimating function
    for( j in 1:K ) { 
      for( i in 1:q ){
        ef[ 0 + ( j - 1 )*q + i ] <- sum( weight * residual[ , j ] * avail * WCLS_weight * Zdm[ , i ] )  
      }
    }
    for( j in 1:M ){
      for( i in 1:p ){
        ef[ q*K + ( j- 1 )*p + i ] <- sum( weight * rowSums( residual ) * avail * WCLS_weight * cA_tilde[ , j ] * Xdm[ , i ] )
      }
    }
    ef <- ef / sample_size
    return( ef )
  }
  
  if ( is.null( estimator_initial_value ) ) { estimator_initial_value <- rep( 0, length = M*p + K*q ) }
  solution <- tryCatch( { multiroot(estimating_equation, estimator_initial_value) },
    error = function(cond) { message("\nCatched error in multiroot inside weighted_centered_least_square():")
      message(cond)
      return( list( root = rep( NaN, M*p + K*q ), msg = cond, f.root = rep( NaN, M*p + K*q ) ) ) } )
  
  estimator <- get_alpha_beta_from_multiroot_result( solution, p, q, M, K )
  alpha_hat <- estimator$alpha
  beta_hat <- estimator$beta
  
  ### asymptotic variance ###
  Mn_summand <- array( NA, dim = c( total_person_decisionpoint*K, M*p + K*q, M*p + K*q ) )
  r_term_collected <- matrix( NA, nrow = M*p + K*q, ncol = total_person_decisionpoint*K ) #residual*avail
  D_term_collected <- matrix( NA, nrow = M*p + K*q, ncol = total_person_decisionpoint*K ) #weight*centred*inverse_factor
  
  for( k in 1:K ){
    for ( it in 1:total_person_decisionpoint ) {
      Xbeta <- as.numeric( t( as.matrix( Xdm[ it, ] ) ) %*% beta_hat )
      Zalpha <- as.numeric( t( as.matrix( Zdm[ it, ] ) ) %*% alpha_hat[ , k ] )
      pre_multiplier <- exp( -sum( A[ it, ]*( Xbeta ) ) ) * WCLS_weight[ it ] #inverse_factor*weight
      gHt <- rep( 0, K*q )
      gHt[ ( ( k - 1 )*q+1 ):( k*q ) ] <- Zdm[ it, ] ## The nuisance covariates for the k-th level ordinal outcome only 
      D_term <- pre_multiplier*c( gHt, rep( cA_tilde[ it , ], each = p ) * rep( Xdm[ it , ], M ) )
      D_term_collected[ , ( k - 1 )*total_person_decisionpoint+it ] <- D_term
      r_term <- ( Y[ it, k ] - exp( Zalpha + sum( A[ it, ]*Xbeta ) ) ) * avail[ it ]
      sum_r_term <- sum( r_term ) 
      r_term_collected[ , ( k - 1 )*total_person_decisionpoint+it ] <- c( rep( r_term, K*q ), rep( sum_r_term, M*p ) )
      
      partialD_partialtheta <- matrix( NA, nrow =  M*p + K*q, ncol = M*p + K*q )
      partialD_partialtheta[ ( 1 ):( K*q + M*p ), 1:( K*q ) ] <- 0
      for( j in 1:M ){
        partialD_partialtheta[ 1:( K*q ), ( ( K*q ) + ( j -1 )*p + 1 ):( ( K*q ) + j*p ) ] <- -pre_multiplier * A[ it, j ] * ( gHt %o% Xdm[ it, ] ) 
      }
      for( i in 1:M ){
        for( j in 1:M ){
          partialD_partialtheta[ ( ( K*q ) + ( i - 1 )*p + 1 ):( ( K*q ) + i*p ), ( ( K*q ) + ( j - 1 )*p + 1 ):( ( K*q ) + j*p ) ] <- 
            -pre_multiplier * A[ it, j ] * cA_tilde[ it , i ] * ( Xdm[ it, ] %o% Xdm[ it, ] )
        }
      }
      exp_Zalpha_sumAXbeta <- exp( Zalpha + sum( A[ it, ]*Xbeta ) )
      partialr_partialtheta <- -exp_Zalpha_sumAXbeta * c( gHt, rep( A[ it, ], each = p ) * rep( Xdm[ it, ], M ) ) * avail[ it ]
      
      Mn_summand[( k - 1 )*total_person_decisionpoint+it, , ] <- 
        partialD_partialtheta*t( matrix( c(  rep( r_term, q*K ), rep( sum_r_term, M*p ) ), nrow = K*q+M*p, ncol = K*q+M*p ) ) + D_term%o%partialr_partialtheta
    }
  }
  Mn <- apply( Mn_summand, c( 2, 3 ), sum ) / sample_size
  Mn_inv <- solve( Mn )
  
  Sigman_summand <- array( NA, dim = c( sample_size, M*p + K*q, M*p + K*q ) )
  ID_T_K <- rep( dta[ , id_varname ], K )
  unqiue_ID <- names( table(ID_T_K)[ !table(ID_T_K)==0 ] )
  length(unqiue_ID)
  for ( i in 1:sample_size ) {
    D_term_i <- D_term_collected[ , ID_T_K %in% unqiue_ID[ i ] ]
    r_term_i <- matrix( r_term_collected[ , ID_T_K %in% unqiue_ID[ i ] ], nrow = K*q+M*p )
    Sigman_summand[ i, , ] <- rowSums( D_term_i*r_term_i ) %*% t( rowSums( r_term_i*D_term_i ) ) 
  } 
  Sigman <- apply( Sigman_summand, c( 2, 3 ), sum ) / sample_size
  
  varcov <- Mn_inv %*% Sigman %*% t( Mn_inv ) / sample_size
  alpha_se <- sqrt( diag( varcov )[ 1:(K*q) ] )
  beta_se <- sqrt( diag( varcov )[ ( K*q + 1 ):( K*q + M * p ) ] )
  names( alpha_hat ) <- names( alpha_se )  <- rep( Znames, K )
  names( beta_hat ) <- names( beta_se ) <- rep( Xnames, M )
  
  return(
    list(beta_hat = beta_hat, alpha_hat = alpha_hat,
              beta_se = beta_se, alpha_se = alpha_se,
              varcov = varcov,
              dims = list(p = p, q = q),
              f.root = solution$f.root
         )
    )
}



