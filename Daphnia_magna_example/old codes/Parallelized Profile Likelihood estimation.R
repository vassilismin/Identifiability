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
  # Create a counter to mark the position of fitted parameters of x
  # that corresponds to a specific combination of PFAS and temperature
  counter <- 1
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
                "C_daphnia_bound" = 0, "C_prot_un" = C_prot_init[[1]])
    
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
  # Create a counter to mark the position of fitted parameters of x
  # that corresponds to a specific combination of PFAS and temperature
  counter <- 1
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
                "C_daphnia_bound" = 0, "C_prot_un" = C_prot_init[[1]])
    
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

profile_likelihood <- function(X){
  i=X$index 
  thetas=X$thetas
  thetas_names=X$thetas_names
  lb=X$lb
  ub=X$ub
  N_thetas=X$N_thetas
  N_iter=X$N_iter
  
  # Load the experimental data
  data_ls <- openxlsx::read.xlsx ('C:/Users/vassi/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023/Wang_data_reduced2.xlsx', sheet = 'PFDoA')
  data_plot <- openxlsx::read.xlsx ('C:/Users/vassi/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023/Wang_data.xlsx', sheet = 'PFDoA')
  errors <- read.csv('C:/Users/vassi/Documents/GitHub/Identifiability/PFDoA_errors.csv')
  
  # the selected settings for the optimizer
  opts <- list( "algorithm" = "NLOPT_LN_NELDERMEAD", #"NLOPT_LN_SBPLX", #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
                "xtol_rel" = 1e-06, 
                "ftol_rel" = 1e-06,
                "ftol_abs" = 0.0,
                "xtol_abs" = 0.0 ,
                "maxeval" = N_iter,
                "print_level" = 0)
  
  # Convert water concentration in ng/L
  Cwater = c(1.44, 2.05, 2.31)*1000
  names(Cwater) <- c("16oC", "20oC", "24oC")
  age = 7+7 # age of D.magna at the beginning of exposure in days
  temperatures <- c(16, 20, 24) #experiment temperature in oC 
  
  MW <- 614 # molecular weight of PFDoA
  
  
  
  sink(paste0("check_progress_",thetas_names[i], ".txt")) #keep a txt to check the progress while running
  
  theta_i_name <- names(thetas)[i] # the name of the theta_i parameter
  theta_i_optimal <- thetas[i] # the optimal value of the theta_i parameter
  # The values around the neighborhood of the theta_i optimal 
  # that will be used for the profile likelihood.
  #theta_i_step <- seq(theta_i_optimal*0.75, theta_i_optimal*1.25, length.out=N_thetas) 
  theta_i_step <- seq(lb[i], ub[i], length.out=N_thetas) 
  
  # For the i-th key create a 2-columns dataframe to store the tested theta_i values
  # and the corresponding values of the objective function
  results_df <- data.frame(matrix(NA, nrow=length(theta_i_step), ncol = 2))
  colnames(results_df) <- c(theta_i_name, "Likelihood")
  
  for (j in 1:length(theta_i_step)) { # loop over every instance of theta_i_step
    # The parameter whose profile likelihood is examined
    constant_theta <- theta_i_step[j]
    # The name of the constant theta parameter
    constant_theta_name <- theta_i_name
    # The initial values of the rest parameters
    x0 <- thetas[-i]
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
    
    cat("Calculating PL for parameter ", theta_i_name, " and j = ", j, "= ", constant_theta , ". Likelihood = ", optimization$objective, "\n")  
    
    results_df[j,] <- c(constant_theta, optimization$objective)
  }
  results_df <- rbind(results_df, c(global_optimization$solution[i], global_optimization$objective))
  results_df <- results_df[order(results_df[,1]),]
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
opts <- list( "algorithm" = "NLOPT_LN_NELDERMEAD" ,"NLOPT_LN_SBPLX", #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 0, 
              "ftol_rel" = 0,
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

#x0 <- c(7.434352e+00, 6.109240e+00, 6.825681e+00, 7.551428e+00, 1.935566e-05)
#x0 <- c(7.020703 , 7.367442 , 7.457980, 7.178138, 2.545041e-05)
#x0 <- c(7.102776, 7.102683, 7.850294, 7.459273, 3.485151e-05)
#x0 <- c(7.180861, 7.034985, 7.777999, 7.647371, 5.231248e-05)
#x0 <- c(7.228657e+00, 6.980294e+00, 7.718556e+00, 7.692513e+00, 5.928412e-05)
x0 <- c(7.26327, 6.956105, 7.697054, 7.729635, 6.275335e-05)
global_optimization<- nloptr::nloptr(x0 = x0,
                                     eval_f = general_obj_func,
                                     # lb	=  c(3,-1,1,3, 1e-07),
                                     # ub =   c(14,14, 14, 14,  1e-03),
                                     lb	=  c(7,6.5,7,7, 1e-05),
                                     ub =   c(8,8, 8, 8,  20e-05),
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
#lb <- c(4, 4, 5, 6, 1e-06) #0.7*thetas
lb <- c(7,6.5,7,7, 1e-05)
# upper bounds of parameters
#ub <- c(8, 9, 8.5, 8.5, 1e-04) #1.3*thetas
ub <- c(8,8, 8, 8,  20e-05)

# prepare the input for parallel processing
X <- vector("list", 5)
for(i in 1:5){
  X[[i]] <- list(index=i, thetas=thetas, thetas_names=thetas_names,
                 lb=lb, ub=ub, N_thetas=50, N_iter=500)
}

start.time <- Sys.time()
cluster <- makeCluster(5)
clusterExport(cl=cluster, c("Size_estimation","dry_weight_estimation", "ode_func",
                            "WSSR", "updated_obj_func","optimized_params", "global_optimization"))
output <- parLapply(cluster, X, profile_likelihood)
stopCluster(cluster)
total.duration <- Sys.time() - start.time
print(total.duration)

# A function to estimate confidence intervals

confidence_intervals_estimation <- function(profile_likelihood_results, general_optimal_parameters,
                                            general_optimal_objective_value, thetas_names,
                                            alpha=0.95, df=1){
  
  # a dataframe with the values of the theta_i and the corresponding objective function values
  # profile_likelihood_results = PL$ku
  # general_optimal_parameters = global_optimization$solution
  # general_optimal_objective_value = global_optimization$objective
  
  Delta_aplha <- qchisq(p=alpha, df)
  
  # Keep only those theta values which satisfy the condition
  # x2(theta) - x2(theta_optimal) < Delta_aplha
  keep_thetas <- data.frame(matrix(NA, ncol = 2))
  colnames(keep_thetas) <- colnames(profile_likelihood_results)
  
  for (i in 1:nrow(profile_likelihood_results)) {
    if(profile_likelihood_results$Likelihood[i] < Delta_aplha + general_optimal_objective_value){
      keep_thetas[i,] <- profile_likelihood_results[i,]
    }
  }
  keep_thetas <- keep_thetas[order(keep_thetas[,1]), ]
  keep_thetas <- na.omit(keep_thetas)
  
  conf_intervals <- c(keep_thetas[which.min(keep_thetas[,1]),1],
                      keep_thetas[which.max(keep_thetas[,1]),1])
  return(conf_intervals)
}

conf_intervals_list <- list(NA)
for (i in 1:length(output)) {
  conf_intervals_list[[i]] <- confidence_intervals_estimation(profile_likelihood_results = output[[i]],
                                                              general_optimal_parameters = global_optimization$solution,
                                                              general_optimal_objective_value = global_optimization$objective,
                                                              thetas_names=thetas_names,
                                                              alpha=0.95, df=1)
  names(conf_intervals_list)[i] <- names(output[[i]])
}



plot_list <- list()
alpha <- 0.9
df <- 1  
for (i in 1:length(output)) {
  # remove_rows <- -which(output[[i]][,2]> 50)
  # if(length(remove_rows)==0){
      data_to_plot <- output[[i]]
  # }else{
  #  data_to_plot <- output[[i]][remove_rows,]
  # }
  current_param <- names(data_to_plot)[1]
  names(data_to_plot)[1] <- "Parameter"
  optimal_value <- data.frame(global_optimization$solution[i], global_optimization$objective)
  names(optimal_value) <- c("Parameter", "Likelihood")
  plot <- ggplot()+
    geom_hline(yintercept=global_optimization$objective + qchisq(alpha,df), linetype="dashed", color = "red", size=1)+
    geom_hline(yintercept=6.57271744517408 , linetype="dashed", color = "blue", size=1)+
    # geom_smooth(data = data_to_plot,  aes(x=Parameter, y=Likelihood), method = 'loess', se=0,
    #             span=3, color = 'black', size=2)+
    geom_line(data = data_to_plot,  aes(x=Parameter, y=Likelihood), color = 'black', size=2)+
    geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=11)+
    geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, colour="pink", size=10)+
    geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=5)+
    
    #scale_y_log10()+
    ylim(c(5,15))+
    
    labs(title = paste0("Profile Likelihood of ", current_param),
         y = 'Likelihood' , x = current_param)+
    theme(plot.title = element_text(hjust = 0.5,size=30), 
          axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
          axis.text.y=element_text(size=22),
          axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
          axis.text.x=element_text(size=22),
          legend.title=element_text(hjust = 0.5,size=25), 
          legend.text=element_text(size=22),
          panel.border = element_rect(colour = "black", fill=NA, size=1.0)) #+ 
    
    # scale_color_manual("Tissues", values=color_codes)+
    # theme(legend.key.size = unit(1.5, 'cm'),  
    #       legend.title = element_text(size=14),
    #       legend.text = element_text(size=14),
    #       axis.text = element_text(size = 14))
  
  #print(plot)
  plot_list[[i]] <- plot
}
# Arrange and print the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)