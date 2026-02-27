process_results <- function(infections, N0, time_between_restocking) {
  # This function takes the data frame of exposure, infection and removal times from simulate_OPA() and reformats it to give a single column for time with counts of exposures, infections and removals
  # INPUTS:
  # infections = the data frame exported from simulate_OPA()
  # N0 = the initial population size
  # time_between_restocking = the time between restocking flock back up to N0
  # 
  # OUTPUTS:
  # result_temporal = the new data frame in a longer format
  
  # Combine all event times into a single vector of unique times
  all_times <- sort(unique(c(infections$Exposure, infections$Infection, infections$Removal)))
  
  # Create a data frame with all_times
  result_temporal <- data.frame(
    Time = all_times,
    Exposures = 0,
    Infections = 0,
    Removals = 0,
    Incorrect.Culls = 0,
    Pop.Size = N0
  )
  
  # Count exposures, infections, and removals up to each time point
  result_temporal <- result_temporal %>%
    rowwise() %>%
    mutate(
      Exposures = sum(infections$Exposure <= Time) - sum(infections$Infection <= Time),
      Infections = sum(infections$Infection <= Time) - sum(infections$Removal <= Time),
      Removals = sum(infections$Removal <= Time),
      Incorrect.Culls = sum(infections$Incorrect[infections$Removal <= Time])
    )
  
  # Identify reset points and calculate removals since the last reset
  reset_points <- which(result_temporal$Time %% time_between_restocking == 0)
  
  # Calculate population size
  for (i in seq_along(result_temporal$Time)) {
    if (result_temporal$Time[i] == 0) {
      result_temporal$Pop.Size[i] <- N0
    } else if (result_temporal$Time[i] %% time_between_restocking == 0) {
      result_temporal$Pop.Size[i] <- N0
    } else {
      last_reset_index <- max(reset_points[reset_points < i], 1)
      result_temporal$Pop.Size[i] <- N0 - (result_temporal$Removals[i] - result_temporal$Removals[last_reset_index])
    }
  }
  
  return(result_temporal)
}
