library(deSolve)
library(ggplot2)
library(parallel)
library(gridExtra)
setwd('C:/Users/vassi/Documents/GitHub/Identifiability/Daphnia_magna_example')

#=================#
# Model functions #
#=================#

# Function for estimating length of D. magna based on age (from Betini et al. (2019))
# Input: age [days], temperature [oC], food["low"/"high"]/ Output: length [mm]
# Considers female D. magna
Size_estimation <<- function(age, temperature = 22, food="high"){
  
  # T = 15 o C
  a_low_15 <-0.354
  b_low_15 <- 0.527
  a_high_15 <- 0.105
  b_high_15 <- 0.953
  
  # T = 25 o C
  a_low_25 <- 0.811
  b_low_25 <- 0.355
  a_high_25 <- 0.698
  b_high_25 <- 0.83
  
  if(food == "low"){
    if(temperature <= 15){
      a <- a_low_15
      b <- b_low_15
    }else if(temperature >= 25){  
      a <- a_low_25
      b <- b_low_25
    }else{ 
      a <- approx(c(15,25), c(a_low_15, a_low_25), temperature)$y
      b <- approx(c(15,25), c(b_low_15, b_low_25), temperature)$y
    }
  }else if (food == "high"){
    if(temperature <= 15){
      a <- a_high_15
      b <- b_high_15
    }else if(temperature >= 25){  
      a <- a_high_25
      b <- b_high_25
    }else{ 
      a <- approx(c(15,25), c(a_high_15, a_high_25), temperature)$y
      b <- approx(c(15,25), c(b_high_15, b_high_25), temperature)$y
    }
  }else{
    stop('food must be either "low" or "high" ')
  }
  return(a + b * log(age))
}

# Dumont et al. (1975)
# Input: length [mm]/ Output: dry weight[mg]
dry_weight_estimation <<- function(L){
  
  w1 = (1.89e-06*(L*1000)^2.25)/1000 #Donka Lake
  w2 = (4.88e-05*(L*1000)^1.80)/1000 #River Sambre
  # Selected w1 after validation with Martin-Creuzburg et al. (2018
  return(w1)
}

#==============
# ODEs System
#==============

ode_func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    # Units explanation:
    # ke: 1/day
    # ku: Lwater/day
    # Cw: ng PFAS/L water
    # C_daphnia_unbound/C_daphnia_bound: mol PFAS/L D.magna
    # C_daphnia_unbound_unmol/C_daphnia_bound_unmol: ng PFAS/g D.magna
    # C_prot_un: mol prot/L
    
    age <- init_age + time
    #size in mm
    size <- Size_estimation(age, temperature = 23, food="high")
    
    
    # dry weight mg
    DW <- dry_weight_estimation(size)
    # Convert DW to WW
    WW <- 15.5 * DW  # Conversion rate 11-20 from DW to WW (Garner et al., 2018)
    # Another post discussing DW to WW can be accessed through:
    #https://www.madsci.org/posts/archives/2005-09/1127049424.Zo.r.html
    
    # Water concentration in ng/L
    dCw <- 0
    # D.magna concentration in lumenconcentration in ng/g
    ku <- 10^ku
    ke <- 10^ke
    kon <- 10^kon
    Ka <- 10^Ka
    koff <- kon/Ka
    # Reported values for Kon and koff range between 1e2-1e04 and 1e-3-1e-1 in L/mol/s 
    # and s-^-1 respectively. Our concentrations are in ng/g. Assuming density = 1000g/L
    # then the concentration in ng/g is multyplied by 1000 to make it ng/L and then by
    # multiply by 1e-9 to make it grams and divide by MW. We do this directly to kon 
    # and koff to make the concentration mol/L so that we can match the literature 
    # values of kon and koff with ours
    C_daphnia_unbound_unmol <- C_daphnia_unbound*MW/(1000*1e-09)
    C_daphnia_bound_unmol <- C_daphnia_bound*MW/(1000*1e-09)
    dC_daphnia_unbound <-  ku*(Cw*1e-09/MW)/WW  - kon*C_prot_un*C_daphnia_unbound +   koff*C_daphnia_bound - ke*C_daphnia_unbound
    dC_daphnia_bound <- kon*C_prot_un*C_daphnia_unbound - koff*C_daphnia_bound
    dC_prot_un <-   koff*C_daphnia_bound -  kon*C_prot_un*C_daphnia_unbound
    C_tot <- C_daphnia_unbound_unmol + C_daphnia_bound_unmol
    return(list(c("dCw" = dCw,   "dC_daphnia_unbound" = dC_daphnia_unbound,
                  "dC_daphnia_bound" = dC_daphnia_bound, "dC_prot_un" = dC_prot_un), 
                "WW" = WW, "C_tot" = C_tot))
  })
}

#=====================================#
#  Weighted Sum of Squared Residuals  #
#=====================================#

WSSR <- function(observed, predicted, weights, comp.names =NULL){
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

# The updated objective function is a slightly changed version of the general objective function
# in order to be used inside the profile likelihood function and be able to change 
# the constant parameter easily
obj_f <- function(x, constant_theta, constant_theta_name, params_names, constant_params=NULL,
                             PFAS_data, Cwater, age,
                             temperatures, MW, weights_values){
  
  # Indexes of body burden and exposure time in data frame
  BB_index <- c(2,4,6)
  ExpTime_index <- c(1,3,5)
  Errors_index <- c(2,4,6)
  # Age of D.magna at beginning of exposure
  init_age <- age
  score <- rep(NA, length(temperatures))
  # Iterate over PFAS names 
  # Load PFAS data
  df <- PFAS_data
  # Assign the j-th value to the theta_i (the constant parameter)
  if(length(constant_theta_name) != length(constant_theta)){
    stop("The constant_theta_name vector must be of equal length with the constant_theta vector")
  }
  if(!is.null(constant_theta)){
    for (j in 1:length(constant_theta)){
      assign(constant_theta_name[j], constant_theta[[j]])
    }  
  }
  # Assign the values of the x vector to the corresponding parameters
  if(length(x) != length(params_names)){
    stop("The params_names must be of equal length with the x vector")
  }
  for (k in 1:length(x)) {
    assign(params_names[k], x[k])
  }
  
  if(!is.null(constant_params)){
    for (k in 1:length(constant_params)) {
      assign(names(constant_params)[k], constant_params[[k]])
    }
  }
  
  # Iterate over number of distinct temperature used in the experiment
  for (temp_iter in 1:length(temperatures)){
    # Initial water concentration of PFAS at selected temperature
    C_water <-  Cwater[temp_iter]
    # Temperature of experiment
    Temp <- temperatures[temp_iter]
    # Time of measurement of selected PFAS at selected temperature
    exp_time <- round(df[!is.na(df[,ExpTime_index[temp_iter]]),ExpTime_index[temp_iter]],1)
    # Body burden of selected PFAS at selected temperature
    BodyBurden <- df[!is.na(df[,BB_index[temp_iter]]),BB_index[temp_iter]]
    # Errors values of selected PFAS at selected temperature
    error_values <- weights_values[!is.na(weights_values[,Errors_index[temp_iter]]),Errors_index[temp_iter]]
    # Time used by numerical solver that integrates the system of ODE
    sol_times <- seq(0,15, 0.1)
    
    # Initial conditions
    inits <- c( "Cw" = C_water[[1]],  "C_daphnia_unbound" = 0,
                "C_daphnia_bound" = 0, "C_prot_un" = 10^C_prot_init[[1]])
    
    params <- c("init_age"=age, "Temp" = Temp, "ku"= ku[[1]], 
                "kon" = kon[[1]], "Ka" = Ka[[1]], "ke"= ke[[1]], "MW" = MW)
    
    solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                        y = inits,
                                        parms = params,
                                        method="lsodes",
                                        rtol = 1e-5, atol = 1e-5))
    
    if(sum(round(solution$time,2) %in% exp_time) == length(exp_time)){
      results <- solution[which(round(solution$time,2) %in% exp_time), 'C_tot']
    }else{
      stop(print("Length of predictions is not equal to the length of data"))
    }
    score[temp_iter] <- WSSR(list(BodyBurden), list(results), list(error_values))
  } 
  
  # Take the average score of all PFAS and temperatures 
  final_score <- mean(score)  
  return(final_score)
}

#=======================#
#   Profile Likelihood  #
#=======================#
profile_likelihood <- function(X
     # obj_f, 
     # i,
     # thetas,
     # thetas_names, 
     # constant_params = NULL,
     # data_df,
     # errors_df,
     # lb, ub, N_samples,
     # alpha, df, q, global_optimum, 
     # min_step_coef, max_step_coef,
     # break_at_bounds = FALSE,
     # # nlopt settings for the main optimization problem
     # opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
     #             "xtol_rel" = 1e-06, 
     #             "ftol_rel" = 1e-06,
     #             "ftol_abs" = 0.0,
     #             "xtol_abs" = 0.0 ,
     #             "maxeval" = 300,
     #             "print_level" = 1),
     # # nlopt settings for the estimation of theta_step
     # opts_theta_step = list("algorithm_step" = 'NLOPT_LN_SBPLX',
     #                        "xtol_rel_step" = 1e-05, 
     #                        "ftol_rel_step" = 1e-05,
     #                        "ftol_abs_step" = 0.0,
     #                        "xtol_abs_step" = 0.0 ,
     #                        "maxeval_step" = 50,
     #                        "print_level_step" = 0),
     # create_txt = TRUE
                               ){
  
  i=X$index # The index of the parameter whose PL is to be estimated
  thetas=X$thetas # The optimal values of the parameters
  thetas_names=X$thetas_names # the names of the parameters
  constant_params = X$constant_params # the constant parameters of the model
  data_df = X$data_df
  errors_df = X$errors_df
  lb=X$lb # The lower bounds of the parameters
  ub=X$ub # The upper bounds of the parameters
  N_samples=X$N_samples # The number of samples for each parameter
  alpha = X$alpha # probability of chi-squared to estimate the quantile
  df = X$df # the degrees of freedom of chi-squared to estimate the quantile
  global_optimum = X$global_optimum # the minimum objective function value (equivalent to maximized likelihood value)
  q = X$q # Parameter to estimate the theta_step
  max_step_coef <- X$max_step_coef
  min_step_coef <- X$min_step_coef
  break_at_bounds <- X$break_at_bounds
  opts <- X$opts
  opts_theta_step <- X$opts_theta_step
  create_txt <- X$create_txt
  
  
  # Estimate the Delta_alpha parameter
  Delta_alpha <- qchisq(alpha, df)
  
  theta_i_name <- names(thetas)[i] # the name of the theta_i parameter
  
  # Function to estimate the theta_step by solving the equation
  # chi^2(theta) - chi^2(theta_hat) - q*Delta_alpha = 0 
  theta_step_estimation <- function(theta_step, theta_last, obj_f, constant_params, index, global_optimum, q, Delta_alpha){
    i <- index
    x <- theta_last
    x[i] <- x[i] + theta_step
    chi2_last <- global_optimum
    
    chi2 <- obj_f(x = x, constant_theta = NULL, constant_theta_name = NULL, params_names = names(x),
                  constant_params = constant_params,
                  PFAS_data = data_df, Cwater = Cwater, age = age, temperatures = temperatures,
                  MW = MW, weights_values = errors_df)
    
    return(abs(chi2 - chi2_last - q*Delta_alpha))
  }
  
  # Set the threshold. The threshold is estimated as chi^2(theta_hat) + Delta_alpha 
  threshold =  global_optimum + Delta_alpha
  
  if(create_txt) sink(paste0("check_progress_",thetas_names[i], ".txt")) #keep a txt to check the progress while running
  
  # Forward search
  cat("Forward search begins...", "\n")
  # Initiate a counter for the number of iterations
  iter_counter <- 0
  
  # Inititate a parameter to check if the threshold is exceeded in each iteration
  current_score <- global_optimum
  
  # Create a 2-columns dataframe to store the tested theta_i values
  # and the corresponding values of the objective function
  forward_results_df <- data.frame(matrix(NA, nrow=N_samples, ncol = 2))
  colnames(forward_results_df) <- c(theta_i_name, "Likelihood")
  
  # take as constant_theta the parameter that will be fixed to a constant value
  constant_theta <- thetas[i]
  
  # Initiate theta_last and set it equal with the theta_hat
  theta_last <- thetas
  
  while (iter_counter < N_samples & current_score < threshold) {
    iter_counter <- iter_counter + 1
    # Estimation of theta_step 
    theta_step <-  nloptr::nloptr(x0 = 0.1,
                                  eval_f = theta_step_estimation,
                                  lb	= 1e-06,
                                  ub = 1,
                                  opts = opts_theta_step,
                                  theta_last = theta_last, index = i,
                                  global_optimum = global_optimum, q = q, Delta_alpha=Delta_alpha,
                                  obj_f=obj_f, 
                                  constant_params=constant_params)$solution
    
    # If theta_step exceeds the limits, take the value of the limit
    if(theta_step > max_step_coef*abs(thetas[[i]])){
      theta_step <- max_step_coef*abs(thetas[[i]])
    }else if(theta_step < min_step_coef*abs(thetas[[i]])){
      theta_step <- min_step_coef*abs(thetas[[i]])
    }
    
    # The parameter whose profile likelihood is examined
    constant_theta <- constant_theta + theta_step
    
    #Check if the constant_theta exceeded the corresponding upper boundary
    if(constant_theta > ub[i]){
      f_exit = 3
      break
    } 
    # The name of the constant theta parameter
    constant_theta_name <- theta_i_name
    # Take as x0 all the parameters except the one which will be fixed
    x0 <- theta_last[-i]
    # The names of the rest parameters
    params_names <- names(x0)
    
    # define the lower and upper bounds of the parameters
    set.seed(12312)
    optimization<- nloptr::nloptr(x0 = x0,
                                  eval_f = obj_f,
                                  lb	= lb[-i],
                                  ub = ub[-i],
                                  opts = opts,
                                  constant_theta = constant_theta,
                                  constant_theta_name = constant_theta_name,
                                  params_names = params_names,
                                  constant_params=constant_params,
                                  PFAS_data = data_df,
                                  Cwater = Cwater,
                                  age = age ,
                                  temperatures = temperatures,
                                  MW = MW,
                                  weights_values = errors_df)
    
    cat("Calculating PL for parameter ", theta_i_name, " and iter = ", iter_counter,". ", "step =", theta_step , constant_theta_name ," = ", constant_theta , "=> Likelihood = ", optimization$objective, "\n")  
    
    forward_results_df[iter_counter,] <- c(constant_theta, optimization$objective)
    
    # update current score
    current_score <- optimization$objective
    
    #update theta_last vector with the new optimized values
    for (k in 1:length(theta_last)) {
      if(k==i){
        theta_last[k] <- constant_theta
      }else if(k>i){
        theta_last[k] <- optimization$solution[k-1]   
      }else{
        theta_last[k] <- optimization$solution[k]} 
    }
    
  }
  if(iter_counter >= N_samples){
    f_exit <- 1
  }else if(current_score > threshold){
    f_exit <- 2
  }
  cat("Forward search ended after ", iter_counter, "iterations.", "\n")
  
  
  
  # Backward search
  cat("Backward search begins...", "\n")
  # Initiate a counter for the number of iterations
  iter_counter <- 0
  # Inititate a parameter to check if the threshold is exceeded in each iteration
  current_score <- global_optimum
  # Create a 2-columns dataframe to store the tested theta_i values
  # and the corresponding values of the objective function
  backward_results_df <- data.frame(matrix(NA, nrow=N_samples, ncol = 2))
  colnames(backward_results_df) <- c(theta_i_name, "Likelihood")
  # take as constant_theta the parameter that will be fixed to a constant value
  constant_theta <- thetas[i]
  # Initiate theta_last and set it equal with the theta_hat
  theta_last <- thetas
  
  while (iter_counter < N_samples & current_score < threshold) {
    iter_counter <- iter_counter + 1
    # Estimation of theta_step 
    theta_step <-  nloptr::nloptr(x0 = 0.1,
                                  eval_f = theta_step_estimation,
                                  lb	= 1e-06,
                                  ub = 1,
                                  opts = opts_theta_step,
                                  theta_last = theta_last, index = i,
                                  global_optimum = global_optimum, q = q, Delta_alpha=Delta_alpha,
                                  obj_f=obj_f, 
                                  constant_params=constant_params)$solution
    
    # If theta_step exceeds the limits, take the value of the limit
    if(theta_step > max_step_coef*abs(thetas[[i]])){
      theta_step <- max_step_coef*abs(thetas[[i]])
    }else if(theta_step < min_step_coef*abs(thetas[[i]])){
      theta_step <- min_step_coef*abs(thetas[[i]])
    }
    
    # The parameter whose profile likelihood is examined
    constant_theta <- constant_theta - theta_step
    #Check if the constant_theta exceeded the corresponding lower boundary
    if(constant_theta < lb[i]){
      b_exit <- 3
      break
    } 
    # The name of the constant theta parameter
    constant_theta_name <- theta_i_name
    # Take as x0 all the parameters except the one which will be fixed
    x0 <- theta_last[-i]
    # The names of the rest parameters
    params_names <- names(x0)
    
    set.seed(12312)
    optimization<- nloptr::nloptr(x0 = x0,
                                  eval_f = obj_f,
                                  lb	= lb[-i],
                                  ub = ub[-i],
                                  opts = opts,
                                  constant_theta = constant_theta,
                                  constant_theta_name = constant_theta_name,
                                  params_names = params_names,
                                  constant_params=constant_params,
                                  PFAS_data = data_df,
                                  Cwater = Cwater,
                                  age = age ,
                                  temperatures = temperatures,
                                  MW = MW,
                                  weights_values = errors_df)
    
    cat("Calculating PL for parameter ", theta_i_name, " and iter = ", iter_counter,". ", "step =", theta_step , constant_theta_name ," = ", constant_theta , "=> Likelihood = ", optimization$objective, "\n")  
    
    backward_results_df[iter_counter,] <- c(constant_theta, optimization$objective)
    
    # update current score
    current_score <- optimization$objective
    
    #update theta_last vector with the new optimized values
    for (k in 1:length(theta_last)) {
      if(k==i){
        theta_last[k] <- constant_theta
      }else if(k>i){
        theta_last[k] <- optimization$solution[k-1]   
      }else{
        theta_last[k] <- optimization$solution[k]} 
    }
  }
  if(iter_counter >= N_samples){
    b_exit <- 1
  }else if(current_score > threshold){
    b_exit <- 2
  }
  cat("Backward search ended after ", iter_counter, "iterations.", "\n")
  
  results_df <- rbind(backward_results_df, forward_results_df, c(thetas[i], global_optimum))
  results_df <- results_df[order(results_df[,1]),]
  results_df <- results_df[complete.cases(results_df),]
  if(create_txt) sink()
  return(list("plik"=results_df, "b_exit"=b_exit, "f_exit"=f_exit))
}

Identifiability_analysis <- function(obj_f, thetas, thetas_names, data_df, errors_df,
                                     lb, ub,
                                     N_samples = 50,
                                     alpha = 0.95, df = 1,
                                     q = 0.5,
                                     global_optimum, 
                                     min_step_coef = 1e-03, max_step_coef = 0.2,
                                     N_cores,
                                     constant_params = NULL,
                                     exported_to_cluster = NULL,
                                     break_at_bounds = FALSE,
                                     # nlopt settings for the main optimization problem
                                     opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
                                                 "xtol_rel" = 1e-06, 
                                                 "ftol_rel" = 1e-06,
                                                 "ftol_abs" = 0.0,
                                                 "xtol_abs" = 0.0 ,
                                                 "maxeval" = 300,
                                                 "print_level" = 0),
                                     # nlopt settings for the estimation of theta_step
                                     opts_theta_step = list("algorithm" = 'NLOPT_LN_SBPLX',
                                                            "xtol_rel_step" = 1e-05, 
                                                            "ftol_rel_step" = 1e-05,
                                                            "ftol_abs_step" = 0.0,
                                                            "xtol_abs_step" = 0.0 ,
                                                            "maxeval_step" = 50,
                                                            "print_level_step" = 0),
                                     create_txt = TRUE){

  
  # Number of parameters tested in the identifiability analysis
  N_parameters <- length(thetas)
  # prepare the input for parallel processing
  X <- vector("list", N_parameters)
  for(i in 1:N_parameters){
    X[[i]] <- list(index=i, thetas=thetas, thetas_names=thetas_names,
                   constant_params = constant_params,
                   data_df = data_df,
                   errors_df = errors_df,
                   lb=lb, ub=ub, N_samples=N_samples, 
                   alpha=alpha, df=df,
                   global_optimum = global_optimum,
                   q = q,
                   max_step_coef = max_step_coef,
                   min_step_coef = min_step_coef,
                   break_at_bounds = FALSE,
                   opts = opts,
                   opts_theta_step=opts_theta_step,
                   create_txt = create_txt)
  }
  
  start.time <- Sys.time()
  cluster <- makeCluster(N_cores)
  # Export to the cluster any function or parameter that the obj_f needs to work.
  clusterExport(cl=cluster, c(names(exported_to_cluster),"obj_f"))
  output <- parLapply(cluster, X, profile_likelihood)
  stopCluster(cluster)
  total.duration <- Sys.time() - start.time
  
  # Estimation of the confidence intervals
  ci_estimation <- function(pl_results, alpha=0.95, df=1, global_optimum){
    # pl_results should be a dataframe with 2 columns containing the theta values 
    # and the corresponding chi^2 values
    
    # As confidence interevals of each theta are considered the borders of the
    # following refion: 
    # {theta | chi^2(theta) - chi^2(theta_hat) < Delta_alpha} where
    Delta_alpha = qchisq(alpha, df)
    
    # Estimate the lower confidence interval
    # Check if the threshold was exceeded at the backwards search. If true, then
    # the lower bound is the last theta under the threshold, else it is considered 
    # -Inf.
    if(pl_results[1,2] > global_optimum + Delta_alpha){
      lower_bound <- pl_results[2,1]
    }else{
      lower_bound <- -Inf
    }
    
    # Estimate the upper confidence interval
    # Check if the threshold was exceeded at the forward search. If true, then
    # the lower bound is the last theta under the threshold, else it is considered 
    # +Inf.
    if(pl_results[dim(pl_results)[1],2] > global_optimum + Delta_alpha){
      upper_bound <- pl_results[(dim(pl_results)[1]-1),1]
    }else{
      upper_bound <- +Inf
    }
    
    return(list("Lower_bound" = lower_bound,
                "Upper_bound" = upper_bound)) 
  }
  
  results_df <- data.frame(matrix(NA, ncol = 6, nrow = length(thetas)))
  rownames(results_df) <- thetas_names
  colnames(results_df) <- c("Non-Identifiability" , "Optimal", "Lower_Bound", "Upper_Bound", "Exit_code_B", "Exit_code_F")
  for (i in 1:length(thetas)) {
    pl_results <- output[[i]]$'plik'
    confidence_intervals <- ci_estimation(pl_results,  alpha = 0.95, df=1, global_optimum)
    results_df$Lower_Bound[i] <- confidence_intervals$Lower_bound
    results_df$Upper_Bound[i] <- confidence_intervals$Upper_bound
    results_df$Optimal[i] <- thetas[[i]]
    results_df$Exit_code_B[i] <- output[[i]]$b_exit
    results_df$Exit_code_F[i] <- output[[i]]$f_exit
    if(confidence_intervals$Lower_bound == -Inf & confidence_intervals$Upper_bound == Inf){
      results_df$`Non-Identifiability`[i] <- "Structural"
    }else if(confidence_intervals$Lower_bound == -Inf & confidence_intervals$Upper_bound != Inf |
             confidence_intervals$Lower_bound != -Inf & confidence_intervals$Upper_bound == Inf){
      results_df$`Non-Identifiability`[i] <- "Practical"
    }else{
      results_df$`Non-Identifiability`[i] <- "Identifiable"
    }
  }
  
  return(list("Likelihood_profiles" = output, "Identiafiability_Analysis" = results_df, "Total_duration" = total.duration))
}

profile_likelihood_plots <- function(analysis_results, thetas, global_optimum, alpha = 0.95,
                                     df = 1){
  output <- analysis_results$Likelihood_profiles
  plot_list <- list()
  for (i in 1:length(output)) {
    data_to_plot <- output[[i]]$plik
    current_param <- names(data_to_plot)[1]
    names(data_to_plot)[1] <- "Parameter"
    optimal_value <- data.frame(thetas[i], global_optimum)
    names(optimal_value) <- c("Parameter", "Likelihood")
    
    plot <- ggplot()+
      geom_hline(yintercept=global_optimum + qchisq(alpha,df), linetype="dashed", color = "red", size=1)+
      #geom_hline(yintercept=global_optimization$objective + qchisq(0.95,1), linetype="dashed", color = "green", size=1)+
      geom_hline(yintercept=global_optimum , linetype="dashed", color = "blue", size=1)+
      geom_line(data = data_to_plot,  aes(x=Parameter, y=Likelihood), color = 'black', size=2)+
      #geom_smooth(data = data_to_plot,  aes(x=Parameter, y=Likelihood), method = "loess", span = 0.5, se =0, color = 'black', size=2)+
      geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=11)+
      geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, colour="pink", size=10)+
      geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=5)+
      
      #scale_y_log10()+
      #ylim(c(6,NA))+
      
      labs(title = paste0("Profile Likelihood of ", current_param),
           y = expression(paste(chi^2, "(", theta, ")")) , x = current_param)+
      theme(plot.title = element_text(hjust = 0.5,size=30), 
            axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
            axis.text.y=element_text(size=22),
            axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
            axis.text.x=element_text(size=22),
            legend.title=element_text(hjust = 0.5,size=25), 
            legend.text=element_text(size=22),
            panel.border = element_rect(colour = "black", fill=NA, size=1.0))
    
    #print(plot)
    plot_list[[i]] <- plot
  }
  # Arrange and print the plots in a grid
  grid.arrange(grobs = plot_list, nrow = 1)
}




################################################################################
################################################################################
################################################################################
#                           END OF FUNCTIONS                                   #
################################################################################
################################################################################
################################################################################



# Extra input needed to be parsed to profile_likelihood (problem-specific parameters and data)
# Load the experimental data
data_ls <- openxlsx::read.xlsx ('C:/Users/vassi/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023/Wang_data_reduced2.xlsx', sheet = 'PFDoA')
errors_df <- read.csv('C:/Users/vassi/Documents/GitHub/Identifiability/Daphnia_magna_example/PFDoA_errors.csv')

# Convert water concentration in ng/L
Cwater = c(1.44, 2.05, 2.31)*1000
names(Cwater) <- c("16oC", "20oC", "24oC")
age = 7+7 # age of D.magna at the beginning of exposure in days
temperatures <- c(16, 20, 24) #experiment temperature in oC 

MW <- 614 # molecular weight of PFDoA

# Test on ku, kon and Ke 

# lower bounds of parameters
#lb <- c(4, 4, 5, 6, 1e-06) 
lb <- c(3,3,3)
# upper bounds of parameters
#ub <- c(8, 9, 8.5, 8.5, 1e-04) 
ub <- c(10,10,10)

constant_params <- c("Ka"=7.681146, "C_prot_init" = -4.189729)

exported_to_cluster = list("Size_estimation"=Size_estimation,
                           "dry_weight_estimation"=dry_weight_estimation,
                           "ode_func"=ode_func,
                           "WSSR"=WSSR,
                           "age"=age,
                           "Cwater"=Cwater,
                           "temperatures"=temperatures,
                           "MW"=MW)

thetas <- c(7.264904, 6.94284, 7.728766)
thetas_names <- c("ku", "kon", "ke")
names(thetas) <- thetas_names
global_optimum <- 6.54365835275129 

test <- Identifiability_analysis(obj_f = obj_f,
                                 thetas=thetas,
                                 thetas_names=thetas_names,
                                 data_df=data_ls ,
                                 errors_df=errors_df,
                                 lb=lb ,
                                 ub=ub,
                                 N_samples = 60,
                                 alpha = 0.95 ,
                                 df = 1,
                                 q = 0.5,
                                 global_optimum = global_optimum ,
                                 min_step_coef = 1e-03 ,
                                 max_step_coef = 0.2,
                                 N_cores = 3,
                                 constant_params = constant_params,
                                 exported_to_cluster = exported_to_cluster,
                                 break_at_bounds = FALSE,
                                 # nlopt settings for the main optimization problem
                                 opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
                                             "xtol_rel" = 1e-06, 
                                             "ftol_rel" = 1e-06,
                                             "ftol_abs" = 0.0,
                                             "xtol_abs" = 0.0 ,
                                             "maxeval" = 300,
                                             "print_level" = 0),
                                 # nlopt settings for the estimation of theta_step
                                 opts_theta_step = list("algorithm" = 'NLOPT_LN_SBPLX',
                                                        "xtol_rel" = 1e-05 ,
                                                        "ftol_rel" = 1e-05,
                                                        "ftol_abs" = 0.0,
                                                        "xtol_abs" = 0.0 ,
                                                        "maxeval" = 50,
                                                        "print_level" = 0),
                                 create_txt = TRUE)


# PRINT PLOTS

profile_likelihood_plots(analysis_results = test, thetas, global_optimum, alpha = 0.95,
                         df = 1)



