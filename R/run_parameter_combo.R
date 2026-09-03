#' @import dplyr
#' @import tidyr
#' @import doParallel
#' @import purrr
#' @export
run_parameter_combo <- function(params,
                                parname1, parname2,
                                parvals1, parvals2,
                                niter, Tmax,
                                continuous_par1=F, continuous_par2=F, continuous_dist=NULL) {
  # Function to simulate multiple trajectories of OPA model under different parameter combinations
  # INPUT:
  # params = list of parameter values required for simulate_OPA() function
  # parname1/2 = name of parameter(s) to vary
  # parvals1/2 = vector of parameter(s) values to consider
  # niter = number of simulations to run for each parameter combination
  # Tmax = duration of the simulation (in days)
  # continuous_par1 = Is the first parameter a continuous numerical parameter (T/F)
  # continuous_par2 = Is the second parameter a continuous numerical parameter (T/F)
  # continuous_dist = What distribution should the continuous parameter be drawn from (options are 'rnorm' or 'runif')
  #
  # OUTPUT:
  # temporal_infections = data frame containing trajectory of the outbreak for each simulation within each parameter combination
  
  # At the moment I haven't set this up to accept two continuous parameters
  if(continuous_par1 == T){
    stop("At present this function can only accept discrete values for parameter1")
  }
  
  # Set up somewhere to store the results
  temporal_infections <- NULL
  parlist <<- NULL
  
  col_names <- c(parname1, parname2)
  
  # Create a cluster with available cores minus one for the main process
  n_cores <- detectCores(logical = FALSE)
  
  # Simulate an OPA trajectory
  cl <- makeCluster(n_cores - 1)
  
  # Export necessary functions and variables to each worker node
  clusterExport(cl, varlist = c("parlist", "Tmax"))
  clusterEvalQ(cl, {
    library(TransmitTestandCull)
    library(dplyr)
    library(tidyr)
  })
  
  # Decide how to sample based on whether the parameters to vary are continuous or discrete
  if((continuous_par1 == F) & (continuous_par2 == F)){
    # Loop through the different parameter combinations
    for(v1 in 1:length(parvals1)){
      for(v2 in 1:length(parvals2)){
        # Set the parameter values to the current combination
        params[parname1] <- parvals1[v1]
        params[parname2] <- parvals2[v2]
        
        parlist <<- rep(list(params), each=niter)

        infections_list <- parLapply(cl, parlist, simulate_OPA, Tmax = Tmax)

        infections_list <- parLapply(cl, infections_list, na.omit)
        
        infections_list <- parLapply(cl, infections_list, process_results, N0 = params[["N0"]], time_between_restocking = params[["time_between_restocking"]])

        # Extracting and adding the parameter as a new column to the data frames
        infections_list <- map2(infections_list, parlist, ~ reduce(col_names, .f = function(df, col_name) add_dynamic_column(df, .y, col_name), .init = .x))
        
        infections <- bind_rows(infections_list, .id="Iteration")
        
        # Combine with previous runs
        temporal_infections <- rbind(temporal_infections,infections)
      }
    }
  } else if ((continuous_par1 == F) & (continuous_par2 == T)){
    # Loop through the different parameter combinations
    for(v1 in 1:length(parvals1)){
      # Set the parameter values to the current combination
      params[parname1] <- parvals1[v1]
        
      parlist <- rep(list(params), each=niter)
      parlist <<- parlist
      
      if (continuous_dist == "rnorm") {
        parlist <- map(parlist, function(p) {
          new_val <- abs(rnorm(1, parvals2[1], parvals2[2]))
          p[[parname2]] <- new_val
          # If parname2 is 'sensitivity_infectious', set sensitivity_latent to half
          if (parname2 == "sensitivity_infectious") {
            p[["sensitivity_latent"]] <- new_val / 2
          }
          p
        })
      } else if (continuous_dist == "runif") {
        parlist <- map(parlist, function(p) {
          new_val <- runif(1, parvals2[1], parvals2[2])
          p[[parname2]] <- new_val
          # If parname2 is 'sensitivity_infectious', set sensitivity_latent to half
          if (parname2 == "sensitivity_infectious") {
            p[["sensitivity_latent"]] <- new_val / 2
          }
          p
        })
      } else {
        stop("Specify a valid distribution ('runif' or 'rnorm') for the continuous parameter")
      }

      # Simulate an OPA trajectory
      infections_list <- parLapply(cl, parlist, simulate_OPA, Tmax = Tmax)

      infections_list <- parLapply(cl, infections_list, na.omit)
      
      # Process the results to the desired format
      infections_list <- parLapply(cl, infections_list, process_results, N0 = params[["N0"]], time_between_restocking = params[["time_between_restocking"]])

      # Extracting and adding the parameter as a new column to the data frames
      infections_list <- map2(infections_list, parlist, ~ reduce(col_names, .f = function(df, col_name) add_dynamic_column(df, .y, col_name), .init = .x))
      
      infections <- bind_rows(infections_list, .id="Iteration")
      
      # Combine with previous runs
      temporal_infections <- rbind(temporal_infections,infections)
    }
  }
  
  stopCluster(cl)

  temporal_infections_agg <- aggregate(cbind(Time,Incorrect.Culls,Removals-Incorrect.Culls) ~ ., data=temporal_infections[,-c(3,4,7)], max)
  names(temporal_infections_agg)[length(temporal_infections_agg)] <- c("Correct.Culls")
  temporal_infections_agg <- reshape(temporal_infections_agg,
                                     direction = "long",
                                     varying = c("Incorrect.Culls","Correct.Culls"),
                                     v.names = "Culls",
                                     idvar = c("Iteration",parname1,parname2),
                                     timevar = "Diagnosis",
                                     times = c("Incorrect","Correct"))
  
  temporal_infections_agg$Time[temporal_infections_agg$Time > Tmax] <- Tmax
  
  return(list(temporal_infections = temporal_infections, aggregated_infections = temporal_infections_agg))
}
