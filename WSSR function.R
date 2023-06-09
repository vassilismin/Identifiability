# This function estimates the Weighted Sum of Squared Residuals

#===================================#
# Weighted Sum of Squared Residuals #
#===================================#

WSSR <- function(observed, predicted, weights, comp.names = NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted) || !is.list(weights)){
    stop(" The observations, predictions and weights must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted) || length(observed) != length(weights)){
    stop(" The observations, predictions and weights must have the same compartments")
  }
  
  # Define the number of observed outputs
  N_outputs <- length(predicted)
  # Define the number of observations per output
  N_obs <- rep(NA, N_outputs)
  
  # A vector to store the values of the weighted squared sum residuals of each compartment
  outputs_res <- c()
  for (i in 1:N_outputs) { # loop over the observed outputs
    N_obs[i] <- length(observed[[i]])
    
    # Check that all observed, predicted and weights vectors have the same length
    if(N_obs[i] != length(predicted[[i]]) || N_obs[i] != length(weights[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    # The number of observations for output i
    N <- N_obs[i]
    
    # Initiate a variable to estimate the sum of squared residuals for output j
    sq_weighted_res_j <- 0
    for (j in 1:N) { #loop over the experimental points i compartment i
      sq_weighted_res_j <- sq_weighted_res_j + ((observed[[i]][j] - predicted[[i]][j]) / weights[[i]][j])^2   
    }
    outputs_res[i] <- sq_weighted_res_j
  }
  
  WSSR_results <- sum(outputs_res)
  
  # Name the list of compartment discrepancy indices
  if ( !is.null(comp.names)){
    names(WSSR_results) <- comp.names
  }else if (!is.null(names(observed))){
    names(WSSR_results) <- names(observed)
  } else if (!is.null(names(predicted)) && is.null(comp.names) ){
    names(WSSR_results) <- names(predicted)
  } else if (!is.null(names(weights)) && is.null(comp.names) ){
    names(WSSR_results) <- names(weights)
  }
  
  return(WSSR_results)
}

###############################
#           Test

y_true_1 <- seq(0,1,0.1)
# add noise to the predictions
y_obs_1 <- y_true_1 + runif(length(y_true_1), 0,0.01)

y_true_2 <- seq(10,30,1)
# add noise to the predictions
y_obs_2 <- y_true_2 + runif(length(y_true_2))

exp_errors_1 <- runif(length(y_true_1))
exp_errors_2 <- runif(length(y_true_2))

observed <- list(y_true_1, y_true_2)
predicted <- list(y_obs_1, y_obs_2)
errors <- list(exp_errors_1,exp_errors_2)

WSSR(observed, predicted, errors)




