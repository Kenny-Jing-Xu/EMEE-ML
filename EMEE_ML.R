library(rootSolve) # for solver function multiroot()
library(geepack) # for fitting GEE using package

fit_EMEE_ML <- function(
    dta, # Data used for analysis
    id_varname, # Variable name for ID
    decision_time_varname, # Variable name for decision time point
    treatment_varname, # Variable name for treatment
    outcome_varname, # Variable name for outcome
    nuisance_method, # Estimation method for the nuisance parameters
    control_varname, # Variable names for the nuisance parameters Ht
    control_varname2, # Variable names for the nuisance part when the 'gam' method is used
    moderator_varname, # Variable names for the moderator St
    rand_prob_varname, # Variable names for the random probabilities P(Ht)
    avail_varname = NULL, # Variable name for the availability indicator It
    rand_prob_tilde_varname = NULL, # Variable names for the random probabilities \tilde{P}(S_t) 
    rand_prob_tilde = NULL, # Random probabilities \tilde{P}(S_t)  
    estimator_initial_value = NULL # Initial value for the estimator
)
{
  sample_size <- length( unique( dta[ , id_varname ] ) )
  total_person_decisionpoint <- nrow( dta )
  time_length <- length( unique( dta[ , decision_time_varname ] ) )
  
  if ( is.null( avail_varname ) ) { avail <- rep( 1, total_person_decisionpoint ) } else { avail <- dta[ , avail_varname ] }
  
  A <- as.matrix( dta[ , treatment_varname ] )
  # checking for NA in treatment indicator    
  if ( any( is.na( A[ avail == 1, ] ) ) ) { stop("Treatment indicator is NA where availability = 1.") }
  M <- dim( A )[ 2 ]
  A[avail == 0] <- 0
  
  p_t <- dta[ , rand_prob_varname ]
  p_t <- as.matrix (p_t)
  cA <- A - p_t # centered A 
  
  #ordinal outcome
  Y <- dta[ , outcome_varname ]
  K <- dim( Y )[ 2 ]
  
  Xdm <- as.matrix( cbind( rep( 1, nrow( dta ) ), dta[ , moderator_varname ] ) ) # X (moderator) design matrix, intercept added
  Zdm <- as.matrix( cbind( rep( 1, nrow( dta ) ), dta[ , c( control_varname, control_varname2 ) ] ) ) # Z (control) design matrix, intercept added
  
  if ( length( rand_prob_tilde ) == M ) {
    p_t_tilde <- matrix( rep( rand_prob_tilde, total_person_decisionpoint ), nrow = total_person_decisionpoint, ncol = M )
  } else { p_t_tilde <- as.matrix( dta[ , rand_prob_tilde_varname ] ) }
  cA_tilde <- A - p_t_tilde
  
  WCLS_weight <- rowSums( A * p_t_tilde / p_t ) + ( 1 - rowSums( A ) )*( 1 - rowSums( p_t_tilde ) )/( 1 - rowSums( p_t ) )
  WCLS_weight <- as.vector( WCLS_weight )
  
  p <- length( moderator_varname ) + 1 # dimension of beta
  q <- length( control_varname ) + 1 # dimension of alpha
  
  Xnames <- c( "Intercept", moderator_varname )
  Znames <- c( "Intercept", control_varname, control_varname2 )
  
  ### estimation: nuisance parameter [] ###
  
  EY <- array( data = NA, dim = c( total_person_decisionpoint, M+1, K ) )
  for(k in 1:K){
    df4 <- data.frame( Y[ , k ], Zdm )
    colnames( df4 ) <- c( outcome_varname[ k ], Znames )
    if( nuisance_method == "glm" ){
      df40 <- df4[ rowSums( A ) %in% c( 0 ), ] 
      formula1 <- as.formula( paste0( outcome_varname[ k ], "~", paste0( c( 1, control_varname ), collapse = "+" ) ) )
      fit40 <- glm2( formula1, data = df40, family = binomial( link = 'logit' ) ) #glm for less covariate, gam for unkown form of the continuous variables, and ML for high dimension  
      EY0 <- predict( fit40, newdata = df4, type = "response" )
      EY[, 1, k] <- EY0
      for(j in 1:M){
        df41 <- df4[ A[ , j ] %in% c( 1 ), ]
        fit41 <- glm( formula1, data = df41, family = binomial( link = 'logit' ) , start = fit40$coefficients )
        EY1 <- predict( fit41, newdata = df4, type = "response" )
        EY[ , ( 1 + j ), k ] <- EY1
      } 
    }
    else if ( nuisance_method == "gam" ){
      df40 <- df4[ rowSums( A ) %in% c( 0 ), ] 
      formula1 <- as.formula( paste0( outcome_varname[ k ], "~", paste0( c( 1, control_varname ),  collapse = "+" ), "+", "s(", paste0( control_varname2, collapse = "+" ), ",k=3)" ) )
      fit40 <- gam( formula1, data = df40, family = binomial( link = 'logit' ) ) #glm for less covariate, gam for unkown form of the continuous variables, and ML for high dimension  
      EY0 <- predict( fit40, newdata = df4, type = "response" )
      EY[, 1, k] <- EY0
      for(j in 1:M){
        df41 <- df4[ A[ , j ] %in% c( 1 ), ]
        fit41 <- gam( formula1, data = df41, family = binomial( link = 'logit' ) , start = fit40$coefficients )
        EY1 <- predict( fit41, newdata = df4, type = "response" )
        EY[ , ( 1 + j ), k ] <- EY1
      } 
    }
    else if ( nuisance_method == "bart" ){
      df40 <- df4[ rowSums( A ) %in% c( 0 ), ] 
      fit40 <- bart( x.train = as.matrix( df40[ control_varname ] ), 
                     y.train = as.matrix( df40[ outcome_varname[ k ] ] ), 
                     keeptrees = TRUE 
                     )
      EY0 <- colMeans( predict( fit40, newdata = as.matrix( df4[ control_varname ] ), type = "response" ) )
      EY[, 1, k] <- EY0
      for(j in 1:M){
        df41 <- df4[ A[ , j ] %in% c( 1 ), ]
        fit41 <- bart( x.train = as.matrix( df41[ control_varname ] ), y.train = as.matrix( df41[ outcome_varname[ k ] ] ), keeptrees = TRUE )
        EY1 <- colMeans( predict( fit41, newdata = as.matrix( df4[ control_varname ] ), type = "response" ) )
        EY[ , ( 1 + j ), k ] <- EY1
      } 
    }
    else if ( nuisance_method == "glmnet" ){
      df40 <- df4[ rowSums( A ) %in% c( 0 ), ] 
      fit40 <- cv.glmnet( x = as.matrix( df40[ control_varname ] ), y = as.matrix( df40[ outcome_varname[ k ] ] ), family = "binomial" , alpha = 0.5)
      EY0 <- predict( fit40, newx = as.matrix( df4[ control_varname ] ), type = "response", s = fit40$lambda.min )
      EY[, 1, k] <- EY0
      for(j in 1:M){
        df41 <- df4[ A[ , j ] %in% c( 1 ), ]
        fit41 <- cv.glmnet( x = as.matrix( df41[ control_varname ] ), y = as.matrix( df41[ outcome_varname[ k ] ] ), family = "binomial", alpha = 0.5 )
        EY1 <- predict( fit41, newx = as.matrix( df4[ control_varname ] ), type = "response", s = fit40$lambda.min )
        EY[ , ( 1 + j ), k ] <- EY1
      } 
    }
    
    else if ( nuisance_method == "earth" ){
      df40 <- df4[ rowSums( A ) %in% c( 0 ), ]
      formula1 <- as.formula( paste0( outcome_varname[ k ], "~", paste0( c( 1, control_varname ), collapse = "+" ) ) )
      fit40 <- earth( formula1, data = df40, glm = list( family = binomial ) )  
      EY0 <- predict( fit40, newdata = as.matrix( df4[ control_varname ] ), type = "response" )
      EY[, 1, k] <- EY0
      for(j in 1:M){
        df41 <- df4[ A[ , j ] %in% c( 1 ), ]
        fit41 <- earth( formula1, data = df41, glm = list( family = binomial ) )  
        EY1 <- predict( fit41, newdata = as.matrix( df4[ control_varname ] ), type = "response" )
        EY[ , ( 1 + j ), k ] <- EY1
      }
    }
    
  }
  
  ### estimation: parameter of interest [beta] ###
  estimating_equation <- function( beta ) {
    beta <- matrix( beta, nrow = p, ncol = M )
    exp_AXdm_beta <- exp( rowSums( A*( Xdm%*%beta ) ) )
    exp_AXdm_beta_K <- matrix( exp_AXdm_beta, nrow = total_person_decisionpoint, ncol = K )
    exp_Xdm_beta <- exp( -Xdm%*%beta )
    exp_Zdm_alpha <- #as.vector( EY[  , 1 ] ) * ( 1 - rowSums( p_t_tilde ) ) + rowSums( EY[  , -1 ] * p_t_tilde * exp_Xdm_beta )
      as.matrix( EY[, 1, ] )*( 1 - rowSums(p_t_tilde) ) + apply( EY[, -1,]*array( p_t_tilde*exp_Xdm_beta, dim = c(total_person_decisionpoint, M, K) ), MARGIN = c(1, 3), FUN = sum )
    residual <- Y - exp_Zdm_alpha*exp_AXdm_beta_K
    residual <- as.matrix( residual )
    weight <- exp_AXdm_beta^( -1 )
    ef <- rep( NA, M*p ) # value of estimating function
    for ( j in 1:M ) {
      for( i in 1:p ){
        ef[ ( j - 1 )*p + i ] <- sum( weight*rowSums(residual)*avail*WCLS_weight*cA_tilde[, j]*Xdm[, i] )
      }
    }
    ef <- ef / sample_size
    return( ef )
  }
  if (is.null( estimator_initial_value ) ) { estimator_initial_value <- rep( 0, length = M * p ) }
  solution <- tryCatch(
    { multiroot( estimating_equation, estimator_initial_value ) },
    error = function( cond ) {
      message( "\nCatched error in multiroot inside estimator 2():" )
      message( cond )
      return( list( root = rep(NaN, M*p), msg = cond, f.root = rep(NaN, M*p) ) )
    }
  )
  beta_hat <- matrix( solution$root, nrow = p )
  
  ### 3. asymptotic variance ###
  Mn_summand <- array( NA, dim = c(total_person_decisionpoint*K, M*p, M*p) )
  r_term_collected <- matrix( NA, nrow = M*p, ncol = total_person_decisionpoint*K ) #residual*avail
  D_term_collected <- matrix( NA, nrow = M*p, ncol = total_person_decisionpoint*K ) #weight*centred*inverse_factor
  for(k in 1:K){
    for(it in 1:total_person_decisionpoint){
      Xbeta <- as.numeric( t( as.matrix( Xdm[it, ] ) )%*%beta_hat )
      pre_multiplier <- exp( -sum( A[it, ]*( Xbeta ) ) )*WCLS_weight[ it ] #inverse_factor*weight
      D_term <- pre_multiplier*rep( cA_tilde[it, ], each = p )*rep( Xdm[it, ], M )
      D_term_collected[ , (k - 1)*total_person_decisionpoint+it ] <- D_term
      partialD_partialbeta <- matrix(NA, nrow =  M*p, ncol = M*p)
      for( i in 1:M ){
        for( j in 1:M ){
          partialD_partialbeta[( (i - 1)*p + 1 ):( i*p ), ( (j - 1)*p + 1 ):( j*p )] <- -pre_multiplier*A[it, j]*cA_tilde[it, i]*( Xdm[it, ]%o%Xdm[it, ] )
        }
      }
      EH <- EY[it, 1, k]*( 1 - sum( p_t_tilde[it, ] ) ) + sum( EY[it, -1, k]*p_t_tilde[it,]*exp( -Xbeta ) )
      r_term <- ( Y[it, k] - EH*exp( sum( A[it, ]*Xbeta ) ) )*avail[it]
      r_term_collected[ , (k - 1)*total_person_decisionpoint+it ] <- rep( r_term, M*p )
      partialr_partialbeta <- rep( EY[it, -1, k]*p_t_tilde[it, ]*exp( -Xbeta )*exp( sum( A[it, ]*Xbeta ) ), each = p )*rep( Xdm[it, ], M )*avail[it] +
        -EH*exp( sum( A[it, ]*Xbeta ) )*rep( A[it, ], each=p )*rep( Xdm[it, ], M )*avail[ it ]
      Mn_summand[(k - 1)*total_person_decisionpoint+it, , ] <- partialD_partialbeta*r_term + D_term%o%partialr_partialbeta
    }
  }
  Mn <- apply( Mn_summand, c(2, 3), sum )/sample_size/time_length
  Mn_inv <- solve( Mn )
  
  Sigman_summand <- array( NA, dim = c(total_person_decisionpoint, M*p, M*p) )
  total_person_decisionpoint_K <- rep( 1:total_person_decisionpoint, K )
  unqiue_total_person_decisionpoint <- 1:total_person_decisionpoint
  for (i in 1:total_person_decisionpoint){
    D_term_i <- D_term_collected[ , total_person_decisionpoint_K %in% unqiue_total_person_decisionpoint[ i ] ]
    r_term_i <- matrix( r_term_collected[ , total_person_decisionpoint_K %in% unqiue_total_person_decisionpoint[ i ] ], nrow = M*p )
    Sigman_summand[ i, , ] <- rowSums( D_term_i*r_term_i ) %*% t( rowSums( r_term_i*D_term_i ) ) 
  }
  Sigman <- apply( Sigman_summand, c(2, 3), sum )/sample_size/time_length
  
  varcov <- Mn_inv%*%Sigman%*%t( Mn_inv )/sample_size/time_length
  beta_se <- sqrt( diag( varcov ) )
  beta_hat <- as.vector( beta_hat )
  names(beta_hat) <- names( beta_se )  <- rep( Xnames, M )
  
  return( list( beta_hat = beta_hat, beta_se = beta_se, varcov = varcov, dims = list( p = p ), f.root = solution$f.root ) )
}


