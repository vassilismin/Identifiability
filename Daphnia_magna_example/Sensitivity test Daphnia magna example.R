library(deSolve)
library(ggplot2)
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

#---------------------
# Custom AUC function
#---------------------
AUC <- function(x, y){
  individual_auc <- c()
  for (i in 1:(length(x)-1)){
    individual_auc[i] <- (y[i]+y[i+1])*(x[i+1]-x[i])/2
  }
  return(sum(individual_auc))
}

#---------------------------------------------------
# SENSITIVITY ANALYSIS FUNCTION for PBPK Models
#---------------------------------------------------
# The function has the following input variables:
# - model: a function of ODEs ready to be used with the deSolve library
# - thetas: a vector with the values of the parameters, which will be analysed
# - thetas_names: a vector with the names (as characters) of the thetas
# - constant_params: a list with the constant parameters used in the model
# - ranges: the perturbation applied to the parameters
# - targets: a vector with the names of the outputs of interest
# - ode_settings: a list with the settings to be used in the deSolve solver
# - heatmap: binary variable that defines if the heatmap plots will be returned

PBPK_sensitivity <- function(model, thetas, thetas_names, constant_params, ranges, targets,
                             ode_settings, heatmap = FALSE){
  
  if(is.null(model)){
    stop("The ODEs of the PBPK model must be provided as a function compatible
         with the deSolve library")
  }
  if(is.null(thetas)){
    stop("A vector with the values of the parameters for the 
         analysis must be provided")
  }
  
  # The variation dp of each parameter will be equal to "ranges" value 
  if(!is.numeric(ranges) | length(ranges) != 1 | ranges<=0 | ranges>1){
    # "ranges" must be a single numeric value in (0,1]
    stop("For local sensitivity analysis \"ranges\"  should be 
         a single numeric value in (0,1]")
  }else{dp <- ranges}
  
  # The total number of parameters to be analysed. 
  N_thetas <- length(thetas) 
  # The number of target outputs of interest
  N_targets <- length(targets)
  
  # Take all objects from the ode_settings list 
  inits <- ode_settings$inits # Initial conditions of the ODEs
  sample_time <- ode_settings$sample_time # Time points of solution
  # Define the solver to use later
  solver <- ifelse(ode_settings$solver == "default", "bdf", ode_settings$solver)
  
  # Assign the corresponding names to the parameters
  names(thetas) <- thetas_names
  # Merge into a vector the constant parameters and the parameters of the sensitivity test
  params <- c(params_list,thetas)
  
  # Solve the ODEs for the initial values of parameters
  solution_0 <- deSolve::ode(times = sample_time,  func = model, y = inits, parms = params, 
                             method=solver, rtol = 1e-5, atol = 1e-5)
  
  if(sum(!(targets %in% colnames(solution_0))) != 0){
    stop("As \"targets\" should be provided the name(s) of one or more
    of the outputs (compertments) of the PBPK model")
  }
  
  # Calculate AUC of each target for the initial parameters
  AUC_0 <- c()
  # Keep the maximum concentration of each output variable
  Cmax_0 <- c()
  for (i in 1:N_targets) {
    AUC_0[i] <- AUC(solution_0[,"time"],solution_0[,targets[i]])
    Cmax_0[i] <- max(solution_0[,targets[i]])
  }
  
  # Assign the optimal values of thetas to thetas_0 vector
  thetas_0 <- thetas
  # Initialize a vector to store the sensitivity indexes
  SI_auc <- matrix(NA, nrow = N_thetas, ncol = N_targets)
  SI_cmax <- matrix(NA, nrow = N_thetas, ncol = N_targets)
  for (i in 1:N_thetas) {
    # Update the thetas vector with the initial values of parameters
    thetas <- thetas_0
    # Apply the desired perturbation to the i-th parameter
    thetas[i] <- thetas[i]*(1 + dp)
    # Assign names to the parameters
    names(thetas) <- thetas_names
    # Merge into a vector the constant parameters and the parameters of the sensitivity test
    params <- c(params_list,thetas)
      
    # Solve the ode system for given initial values and time
    solution <- deSolve::ode(times = sample_time,  func = model, y = inits, parms = params, 
                             method=solver, rtol = 1e-5, atol = 1e-5)
    
    for (j in 1:N_targets) {
      # Calculate AUC for the target compartment j
      AUC_j <- AUC(solution[,"time"],solution[,targets[j]])
      # Keep the max value for the target compartment j
      Cmax_j <- max(solution[,targets[j]])
      # Calculate sensitivity index of parameter i 
      # Relative Sensitivity of AUC = (dAUC/AUC)/(dp/p)
      SI_auc[i,j] <- ((AUC_j-AUC_0[j])/AUC_0[j])/((thetas[i]- thetas_0[i])/thetas_0[i])
      # Relative Sensitivity of Cmax = (dCmax/Cmax)/(dp/p)
      SI_cmax[i,j] <- ((Cmax_j-Cmax_0[j])/Cmax_0[j])/((thetas[i]- thetas_0[i])/thetas_0[i])
    }
  }
  
  rownames(SI_auc) <- thetas_names
  colnames(SI_auc) <-  paste0("AUC_", targets)
  
  rownames(SI_cmax) <- thetas_names
  colnames(SI_cmax) <- paste0("Cmax_", targets)
  
  if(heatmap){
    
    if(length(targets)==1){
      stop("Provide more than 1 targets in order to create a heatmap
           or turn \"heatmap\" to \"FALSE\".")
    }
    
    
    #########################################################################
    
    heatmap1 <- pheatmap::pheatmap(as.matrix(abs(SI_auc)),
                                   cellwidth = 60, cellheight = 60,
                                   border_color = NA,
                                   fontsize = 20,
                                   fontsize_row = 15, 
                                   fontsize_col = 15,
                                   col = colorRampPalette(RColorBrewer::brewer.pal(8, "YlOrRd"))(25), 
                                   main="AUC Senesitivity Indexes")
    
    AUC_heatmap <- recordPlot(heatmap1)
    
    heatmap2 <- pheatmap::pheatmap(as.matrix(abs(SI_cmax)),
                                   cellwidth = 60, cellheight = 60,
                                   border_color = NA,
                                   fontsize = 20,
                                   fontsize_row = 15, 
                                   fontsize_col = 15,
                                   col = colorRampPalette(RColorBrewer::brewer.pal(8, "YlOrRd"))(25),
                                   main="Cmax Senesitivity Indexes")
    Cmax_heatmap <- recordPlot(heatmap2)
    
    plot_list <- list("AUC_heatmap"=AUC_heatmap, "Cmax_heatmap"=Cmax_heatmap)
  }
  
  # Create a list to return the results 
  
  data_list <- list("Method"="Local Sensitivity", "Targets"=targets,
                    "Change of parameters" = dp,
                    "Parameters" = thetas_names, 
                    "Normalized AUC Sensitivity Coefficients"=data.frame(SI_auc), 
                    "Normalized C_max Sensitivity Coefficients"=data.frame(SI_cmax),
                    "Plot_list"=plot_list)
  
  return(data_list)
}

################################################################################
index <- 1
# Convert water concentration in ng/L
Cwater = c(1.44, 2.05, 2.31)*1000
names(Cwater) <- c("16oC", "20oC", "24oC")
age = 7+7 # age of D.magna at the beginning of exposure in days
temperatures <- c(16, 20, 24) #experiment temperature in oC 
MW <- 614 # molecular weight of PFDoA

params_list <- list('Cwater'=Cwater[[index]],
                    'init_age'=age,
                    'temperatures'=temperatures,
                    'MW'=MW)

# Initial conditions
inits <- c( "Cw" = Cwater[[index]],  "C_daphnia_unbound" = 0,
            "C_daphnia_bound" = 0, "C_prot_un" = 10^(-4.189729))

sol_times <- seq(0, 15, 0.1)

thetas <- c(7.264905, 6.942806, 7.681146, 7.728766)#, -4.189729)
thetas_names <- c('ku', 'kon', 'Ka', 'ke') #, 'C_init_prot')

targets <- c("C_daphnia_unbound", "C_daphnia_bound", "C_prot_un", "C_tot")

ode_settings <- list(inits=inits,
                     sample_time=sol_times,
                     solver= "lsodes")

sensitivity_test <- PBPK_sensitivity(model=ode_func, thetas=thetas, 
                                     thetas_names=thetas_names,
                                     constant_params=params_list,
                                     ranges=0.001, targets = targets,
                                     ode_settings=ode_settings, heatmap = TRUE) 
