library(deSolve)
library(ggplot2)
library(parallel)
library(gridExtra)
setwd('C:/Users/vassi/Documents/GitHub/Identifiability')
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

# This general objective function is used only to optimize all the parameters 
# in order to find the optimal solution
general_obj_func <- function(x, PFAS_data, Cwater, age,
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
  # Iterate over number of distinct temperature used in the experiment
  ku <- x[1]
  kon <- x[2]
  Ka <- x[3]
  ke <- x[4]
  C_prot_init <- x[5]
  
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
                                        rtol = 1e-6, atol = 1e-6))
    
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

# The updated objective function is a slightly changed version of the general objective function
# in order to be used inside the profile likelihood function and be able to change 
# the constant parameter easily
updated_obj_func <- function(x, constant_theta, constant_theta_name, params_names, 
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
  assign(constant_theta_name, constant_theta)
  # Assign the values of the x vector to the corresponding parameters
  for (k in 1:length(x)) {
    assign(params_names[k], x[k])
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
profile_likelihood <- function(X){
  i=X$index # The index of the parameter whose PL is to be estimated
  thetas=X$thetas # The optimal values of the parameters
  thetas_names=X$thetas_names # the names of the parameters
  lb=X$lb # The lower bounds of the parameters
  ub=X$ub # The upper bounds of the parameters
  N_samples=X$N_samples # The number of samples for each parameter
  N_iter=X$N_iter # The maximum number of iterations of the optimizer
  alpha = X$alpha # probability of chi-squared to estimate the quantile
  df = X$df # the degrees of freedom of chi-squared to estimate the quantile
  global_optimum = X$global_optimum # the minimum objective function value (equivalent to maximized likelihood value)
  q = X$q # Parameter to estimate the theta_step
  max_step_coef <- X$max_step_coef
  min_step_coef <- X$min_step_coef
  
  # Estimate the Delta_alpha parameter
  Delta_alpha <- qchisq(alpha, df)
  
  # Load the experimental data
  data_ls <- openxlsx::read.xlsx ('C:/Users/vassi/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023/Wang_data_reduced2.xlsx', sheet = 'PFDoA')
  data_plot <- openxlsx::read.xlsx ('C:/Users/vassi/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023/Wang_data.xlsx', sheet = 'PFDoA')
  errors <- read.csv('C:/Users/vassi/Documents/GitHub/Identifiability/PFDoA_errors.csv')
  
  
  
  # Convert water concentration in ng/L
  Cwater = c(1.44, 2.05, 2.31)*1000
  names(Cwater) <- c("16oC", "20oC", "24oC")
  age = 7+7 # age of D.magna at the beginning of exposure in days
  temperatures <- c(16, 20, 24) #experiment temperature in oC 
  
  MW <- 614 # molecular weight of PFDoA
  
  theta_i_name <- names(thetas)[i] # the name of the theta_i parameter
  
  # Function to estimate the theta_step by solving the equation
  # chi^2(theta) - chi^2(theta_hat) - q*Delta_alpha = 0 
  theta_step_estimation <- function(theta_step, theta_last, index, global_optimum, q, Delta_alpha){
    i <- index
    x <- theta_last
    x[i] <- x[i] + theta_step
    chi2_last <- global_optimum
    
    chi2 <- general_obj_func(x=x, data_ls, Cwater, age, temperatures, 
                             MW, errors)
    
    return(abs(chi2 - chi2_last - q*Delta_alpha))
  }
  
  # Set the threshold. The threshold is estimated as chi^2(theta_hat) + Delta_alpha 
  threshold =  global_optimum + Delta_alpha
  
  sink(paste0("check_progress_",thetas_names[i], ".txt")) #keep a txt to check the progress while running
  
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
  
  # the selected settings for the optimizer
  opts <- list( "algorithm" = "NLOPT_LN_NELDERMEAD", #"NLOPT_LN_SBPLX", #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
                "xtol_rel" = 1e-06, 
                "ftol_rel" = 1e-06,
                "ftol_abs" = 0.0,
                "xtol_abs" = 0.0 ,
                "maxeval" = N_iter,
                "print_level" = 0)
  
  opts_theta_step <- list( "algorithm" = 'NLOPT_LN_BOBYQA', #"NLOPT_LN_NELDERMEAD", #"NLOPT_LN_SBPLX", #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
                           "xtol_rel" = 1e-05, 
                           "ftol_rel" = 1e-05,
                           "ftol_abs" = 0.0,
                           "xtol_abs" = 0.0 ,
                           "maxeval" = 50,
                           "print_level" = 0)
  
  while (iter_counter < N_samples & current_score < threshold) {
    iter_counter <- iter_counter + 1
    # Estimation of theta_step 
    theta_step <-  nloptr::nloptr(x0 = 0.1,
                                  eval_f = theta_step_estimation,
                                  lb	= 1e-06,
                                  ub = 1,
                                  opts = opts_theta_step,
                                  theta_last = theta_last, index = i,
                                  global_optimum = global_optimum, q = q, Delta_alpha=Delta_alpha)$solution
    # If theta_step exceeds the limits, take the value of the limit
    if(theta_step > max_step_coef*abs(constant_theta)){
      theta_step <- max_step_coef*abs(constant_theta)
    }else if(theta_step < min_step_coef*abs(constant_theta)){
      theta_step <- min_step_coef*abs(constant_theta)
    }
    
    # The parameter whose profile likelihood is examined
    constant_theta <- constant_theta + theta_step
    #Check if the constant_theta exceeded the corresponding upper boundary
    #if(constant_theta > ub[i]) break
    # The name of the constant theta parameter
    constant_theta_name <- theta_i_name
    # Take as x0 all the parameters except the one which will be fixed
    x0 <- theta_last[-i]
    # The names of the rest parameters
    params_names <- names(x0)
    
    
    
    # define the lower and upper bounds of the parameters
    set.seed(12312)
    optimization<- nloptr::nloptr(x0 = x0,
                                  eval_f = updated_obj_func,
                                  lb	= lb[-i],
                                  ub = ub[-i],
                                  opts = opts,
                                  constant_theta = constant_theta,
                                  constant_theta_name = constant_theta_name,
                                  params_names = params_names,
                                  PFAS_data = data_ls,
                                  Cwater = Cwater,
                                  age = age ,
                                  temperatures = temperatures,
                                  MW = MW,
                                  weights_values = errors)
    
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
                                  global_optimum = global_optimum, q = q, Delta_alpha=Delta_alpha)$solution
    # If theta_step exceeds the limits, take the value of the limit
    if(theta_step > max_step_coef*abs(constant_theta)){
      theta_step <- max_step_coef*abs(constant_theta)
    }else if(theta_step < min_step_coef*abs(constant_theta)){
      theta_step <- min_step_coef*abs(constant_theta)
    }
    # The parameter whose profile likelihood is examined
    constant_theta <- constant_theta - theta_step
    #Check if the constant_theta exceeded the corresponding lower boundary
    #if(constant_theta < lb[i]) break
    # The name of the constant theta parameter
    constant_theta_name <- theta_i_name
    # Take as x0 all the parameters except the one which will be fixed
    x0 <- theta_last[-i]
    # The names of the rest parameters
    params_names <- names(x0)
    
    set.seed(12312)
    optimization<- nloptr::nloptr(x0 = x0,
                                  eval_f = updated_obj_func,
                                  lb	= lb[-i],
                                  ub = ub[-i],
                                  opts = opts,
                                  constant_theta = constant_theta,
                                  constant_theta_name = constant_theta_name,
                                  params_names = params_names,
                                  PFAS_data = data_ls,
                                  Cwater = Cwater,
                                  age = age ,
                                  temperatures = temperatures,
                                  MW = MW,
                                  weights_values = errors)
    
    cat("Calculating PL for parameter ", theta_i_name, " and iter = ", iter_counter,". ", "step =", theta_step, ", ", constant_theta_name ," = ", constant_theta , "=> Likelihood = ", optimization$objective, "\n")  
    
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
  cat("Backward search ended after ", iter_counter, "iterations.", "\n")
  
  results_df <- rbind(backward_results_df, forward_results_df, c(global_optimization$solution[i], global_optimization$objective))
  results_df <- results_df[order(results_df[,1]),]
  results_df <- results_df[complete.cases(results_df),]
  sink()
  return(results_df)
}

#############################################
# Here we optimize the value of the 5 parameters in order to find best solution
# Load the experimental data
data_ls <- openxlsx::read.xlsx ('C:/Users/vassi/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023/Wang_data_reduced2.xlsx', sheet = 'PFDoA')
data_plot <- openxlsx::read.xlsx ('C:/Users/vassi/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023/Wang_data.xlsx', sheet = 'PFDoA')
errors <- read.csv('C:/Users/vassi/Documents/GitHub/Identifiability/PFDoA_errors.csv')

# the selected settings for the optimizer
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NELDERMEAD" ,#"NLOPT_LN_SBPLX", #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 1e-08, 
              "ftol_rel" = 1e-08,
              "ftol_abs" = 0,
              "xtol_abs" = 0 ,
              "maxeval" = 3000,
              "print_level" = 1)

# Convert water concentration in ng/L
Cwater = c(1.44, 2.05, 2.31)*1000
names(Cwater) <- c("16oC", "20oC", "24oC")
age = 7+7 # age of D.magna at the beginning of exposure in days
temperatures <- c(16, 20, 24) #experiment temperature in oC 

MW <- 614 # molecular weight of PFDoA

x0 <- c(-1, 2, 5, -3, -6)
set.seed(13)
x0 <- runif(5, min = c(-2, 1, 3, -2, -9), max = c( 3, 5, 6, 0, -2))
#x0 <- c(2.591115, 4, 4.340183, 0.9271614, -3)
#x0 <- c(2.353689, 17.09018, 4.793657, 0.1949285, -3.993427)
global_optimization<- nloptr::nloptr(x0 = x0,
                                     eval_f = general_obj_func,
                                     #       ku, kon, Ka, ke, C_init_prot
                                     lb	=  c(-2, 1, 3, -2, -9),
                                     ub =  c( 3, 5, 6, 0, -2),
                                     opts = opts,
                                     PFAS_data = data_ls,
                                     Cwater = Cwater,
                                     age = age ,
                                     temperatures = temperatures,
                                     MW = MW,
                                     weights_values = errors)

###############################

optimized_params <- global_optimization$solution
names(optimized_params) <- c('ku', 'kon', 'Ka', 'ke', 'C_prot_init')

thetas <- optimized_params
thetas_names <- names(optimized_params)

# lower bounds of parameters
#lb <- c(4, 4, 5, 6, 1e-06) 
lb <- c(3,3,3,3, -8)
# upper bounds of parameters
#ub <- c(8, 9, 8.5, 8.5, 1e-04) 
ub <- c(10,10, 10, 10,  2)

# prepare the input for parallel processing
X <- vector("list", 5)
for(i in 1:5){
  X[[i]] <- list(index=i, thetas=thetas, thetas_names=thetas_names,
                 lb=lb, ub=ub, N_samples=100, N_iter=250,
                 alpha=0.95, df=1,
                 global_optimum = global_optimization$objective,
                 q = 0.5,
                 max_step_coef = 0.5,
                 min_step_coef = 0.001)
}

start.time <- Sys.time()
cluster <- makeCluster(5)
clusterExport(cl=cluster, c("Size_estimation","dry_weight_estimation", "ode_func",
                            "WSSR", "updated_obj_func","optimized_params", "global_optimization",
                            "general_obj_func"))
output <- parLapply(cluster, X, profile_likelihood)
stopCluster(cluster)
total.duration <- Sys.time() - start.time
print(total.duration)





plot_list <- list()
alpha <- 0.95
df <- 1
for (i in 1:length(output)) {
  data_to_plot <- output[[i]]
  current_param <- names(data_to_plot)[1]
  names(data_to_plot)[1] <- "Parameter"
  optimal_value <- data.frame(global_optimization$solution[i], global_optimization$objective)
  names(optimal_value) <- c("Parameter", "Likelihood")
  
  plot <- ggplot()+
    geom_hline(yintercept=global_optimization$objective + qchisq(alpha,df), linetype="dashed", color = "red", size=1)+
    #geom_hline(yintercept=global_optimization$objective + qchisq(0.95,1), linetype="dashed", color = "green", size=1)+
    geom_hline(yintercept=global_optimization$objective , linetype="dashed", color = "blue", size=1)+
    geom_line(data = data_to_plot,  aes(x=Parameter, y=Likelihood), color = 'black', size=2)+
    #geom_smooth(data = data_to_plot,  aes(x=Parameter, y=Likelihood), method = "loess", span = 0.5, se =0, color = 'black', size=2)+
    geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=11)+
    geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, colour="pink", size=10)+
    geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=5)+
    
    #scale_y_log10()+
    ylim(c(5,NA))+
    
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
grid.arrange(grobs = plot_list, nrow = 3)
