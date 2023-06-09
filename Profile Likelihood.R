library(deSolve)


#=================#
# Model functions #
#=================#

#=========================
#1. Parameters of the model
#=========================

create.params  <- function(user_input){
  with( as.list(user_input),{
    #Make sure the user provided the proper input
    if(init_age<1 ){
      stop(" The model works for D.magna's with age greater than one")
    }
    if(temperature <15 | temperature >25){
      stop(" The model supports only temperatures between 15 and 25 degrees Celsius")
    }
    
    # Convert water concentration from mg/L to ng/L
    Cwater <- Cwater*1000
    
    if(food == "H"){
      food <- "high"
    }else if(food == "L"){
      food <- "low"
    }else{
      stop("Please insert either H or L for high or low food, respectively")
    }
    
    # Based on the desired pfas, select the appropriate parameters
    if(PFAS == "PFBA"){
      ku <- 8.569318
      ke <- 8.730434
      kon <-  6.147169
      Ka <- 6.316426
      C_prot_init <- 1.103837e-05 
      MW <- 214
    }else if(PFAS == "F-53B"){
      ku <- 9.632826
      ke <- 8.016869
      kon <- 5.699413
      Ka <- 6.435146
      C_prot_init <- 1.142351e-05
      MW <- 570
    }else if(PFAS == "GenX"){
      ku <- 9.022480
      ke <- 7.701288
      kon <- 5.476549
      Ka <- 6.064907
      C_prot_init <- 1.403248e-05
      MW <- 330
    }else if(PFAS == "PFBS"){
      ku <- 7.074107
      ke <- 8.191594
      kon <- 6.422378 
      Ka <- 6.192971
      C_prot_init <- 1.996858e-05
      MW <- 300
    }else if(PFAS == "PFDA"){
      ku <- 8.6780925944 
      ke <- 8.191594
      kon <- 5.4507137475 
      Ka <- 6.1024588862 
      C_prot_init <- 1.996858e-05 
      MW <- 514
    }else if(PFAS == "PFDoA"){
      ku <- 7.434352 
      ke <- 7.551428
      kon <- 6.109240
      Ka <- 6.825681
      C_prot_init <- 1.935566e-05
      MW <- 614
    }else if(PFAS == "PFHpA"){
      ku <- 7.0065455179
      ke <- 8.8746389167 
      kon <- 6.6835913063 
      Ka <- 6.8162461454 
      C_prot_init <- 0.0001052127 
      MW <- 364
    }else if(PFAS == "PFHxA"){
      ku <- 6.574575
      ke <- 7.128324
      kon <- 6.017399
      Ka <- 6.016640
      C_prot_init <- 1.047898e-05 
      MW <- 314
    }else if(PFAS == "PFNA"){
      ku <- 8.353166
      ke <- 7.414018
      kon <- 5.295278
      Ka <- 6.049464
      C_prot_init <- 1.360233e-05 
      MW <- 464
    }else if(PFAS == "PFOA"){
      ku <- 8.326005759
      ke <- 7.466661752 
      kon <- 5.188504644 
      Ka <- 5.724201956 
      C_prot_init <- 0.000015299 
      MW <- 414
    }else if(PFAS == "PFOS"){
      ku <- 8.82163
      ke <- 7.241093
      kon <- 5.055001
      Ka <- 5.819318
      C_prot_init <- 9.684750e-06
      MW <- 500
    }else if(PFAS == "PFPeA"){
      ku <- 7.065735
      ke <-  8.054163
      kon <- 6.516348
      Ka <- 6.372870
      C_prot_init <- 8.270121e-06
      MW <- 364
    }else if(PFAS == "PFUnA"){
      ku <- 8.304492
      ke <- 8.303186
      kon <- 6.069595
      Ka <- 6.785245
      C_prot_init <- 1.931527e-05
      MW <- 564
    }else{
      stop("Please select one of the 13 PFAS that the model can describe")
    }
    
    return(list("init_age" = init_age, "food" = food, "temperature" = temperature,
                "ku" = ku, "ke" = ke, "kon" = kon, "Ka" = Ka,
                "C_prot_init" = C_prot_init, "Cwater" = Cwater, "Cwater_time" = Cwater_time, "MW" = MW))
    
  })
}

#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function(parameters){
  with( as.list(parameters),{
    "Cw" = 0; "C_daphnia_unbound" = 0;
    "C_daphnia_bound" = 0; "C_prot_un" = C_prot_init
    
    return(c( "Cw" = Cw, "C_daphnia_unbound" = 0,
              "C_daphnia_bound" = 0, "C_prot_un" = C_prot_init ))
  })
}

#===================
#3. Events function
#===================

create.events <- function(parameters){
  with(as.list(parameters), {
    
    # Calculate measured concentrations and corresponding measurement time
    lcon <- length(Cwater)
    ltimes <- length(Cwater_time)
    # If not equal, then stop 
    if (ltimes != lcon){
      stop("The times of measurement should be equal in number to the concentrations")
    }else{
      events <- list(data = rbind(data.frame(var = c("Cw"),  time = Cwater_time, 
                                             value = c(Cwater), method = c("rep")) ))
    }
    return(events)
  })
}


#==================
#4. Custom function 
#==================
custom.func <- function(age, temperature, food){
  
  # Estimation of length of D. magna based on age (from Betini et al. (2019))
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
  }
  size = a + b * log(age) #mm
  w = (1.89e-06*(size*1000)^2.25)/1000 #mg DW
  
  return(list(size = size, DW = w))
}

#==============
#5. ODEs System
#==============

ode.func <- function(time, inits, params, custom.func){
  with(as.list(c(inits,params)),{
    # Units explanation:
    # ke: 1/day
    # ku: Lwater/day
    # Cw: ng PFAS/L water
    # C_daphnia_unbound/C_daphnia_bound: mol PFAS/L D.magna
    # C_daphnia_unbound_unmol/C_daphnia_bound_unmol: ng PFAS/g D.magna
    # C_prot_un: mol prot/L
    
    age <- init_age + time #in days
    from_func <- custom.func(age, temperature, food)
    
    #size in mm
    size <- from_func$size
    
    # dry weight in mg
    DW <- from_func$DW
    # Convert DW to WW
    WW <- 15.5 * DW  # Conversion rate 11-20 from DW to WW (Garner et al., 2018)
    # Another post discussing DW to WW can be accessed through:
    #https://www.madsci.org/posts/archives/2005-09/1127049424.Zo.r.html
    
    ku <- 10^ku
    ke <- 10^ke
    kon <- 10^kon
    Ka <- 10^Ka
    koff <- kon/Ka
    
    C_daphnia_unbound_unmol <- C_daphnia_unbound*MW/(1000*1e-09)
    C_daphnia_bound_unmol <- C_daphnia_bound*MW/(1000*1e-09)
    # The concentration of the water is defined through the user input
    dCw <- 0
    dC_daphnia_unbound <-  ku*(Cw*1e-09/MW)/WW  - kon*C_prot_un*C_daphnia_unbound +   koff*C_daphnia_bound - ke*C_daphnia_unbound
    dC_daphnia_bound <- kon*C_prot_un*C_daphnia_unbound - koff*C_daphnia_bound
    dC_prot_un <-   koff*C_daphnia_bound -  kon*C_prot_un*C_daphnia_unbound
    C_tot <- C_daphnia_unbound_unmol + C_daphnia_bound_unmol
    return(list(c("dCw" = dCw,   "dC_daphnia_unbound" = dC_daphnia_unbound,
                  "dC_daphnia_bound" = dC_daphnia_bound, "dC_prot_un" = dC_prot_un), 
                "WW" = WW, "C_tot" = C_tot, "C_bound" = C_daphnia_bound_unmol,  
                "C_unbound" = C_daphnia_unbound_unmol))
  })
}

#======================
#6. Objective function
#======================

obj_func <- function(x, PFAS_data, PFAS_name, Cwater, age, temperatures, MW, metric){
  
  # Indexes of body burden and exposure time in data frame
  BB_index <- c(2,4,6)
  ExpTime_index <- c(1,3,5)
  # Age of D.magna at beginning of exposure
  init_age <- age
  # Create a counter to mark the position of fitted parameters of x
  # that corresponds to a specific combination of PFAS and temperature
  counter <- 1
  score <- rep(NA, length(temperatures))
  # Iterate over PFAS names 
  # Load PFAS data
  df <- PFAS_data[[PFAS_name]]
  # Iterate over number of distinct temperature used in the experiment
  ku <- x[1]
  kon <- x[2]
  Ka <- x[3]
  ke <- x[4]
  C_prot_init <- x[5]
  for (temp_iter in 1:length(temperatures)){
    # Initial water concentration of PFAS at selected temperature
    C_water <-  Cwater[PFAS_name,temp_iter]
    # Temperature of experiment
    Temp <- temperatures[temp_iter]
    # Time of measurement of selected PFAS at selected temperature
    exp_time <- round(df[!is.na(df[,ExpTime_index[temp_iter]]),ExpTime_index[temp_iter]],1)
    # Body burden of selected PFAS at selected temperature
    BodyBurden <- df[!is.na(df[,BB_index[temp_iter]]),BB_index[temp_iter]]
    # Time used by numerical solver that integrates the system of ODE
    sol_times <- seq(0,15, 0.1 )
    # Fitted parameters
    
    
    inits <- c( "Cw" = C_water,  "C_daphnia_unbound" = 0,
                "C_daphnia_bound" = 0, "C_prot_un" = C_prot_init)
    
    params <- c("init_age"=age, "Temp" = Temp, "ku"= ku, 
                "kon" = kon, "Ka" = Ka, "ke"= ke, "MW" = MW)
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
    
    # Find the position of the current PFAS in the PFAS_names vector
    if(metric == "AAFE"){
      score[temp_iter] <- AAFE(BodyBurden, results) 
    }else if (metric =="rmse"){
      score[temp_iter] <- rmse(BodyBurden, results)
    }else if(metric == "PBKOF"){
      score[temp_iter] <- PBKOF(list(BodyBurden), list(results))
    }       
  } 
  
  # Take the average score of all PFAS and temperatures 
  final_score <- mean(score)  
  return(final_score)
}

#=====================================#
#7. Weighted Sum of Squared Residuals #
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