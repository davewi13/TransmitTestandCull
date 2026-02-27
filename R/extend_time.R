extend_data <- function(df, time_interval = 182.5, max_time = 3650) {
  # Generate a complete set of time points up to the maximum time
  if(length(time_interval) == 1){
    all_times <- seq(time_interval, max_time, by = time_interval)
  } else {
    all_times <- time_interval
  }
  
  # Create a data frame with all combinations of Iteration and time points
  Iterations <- unique(df$Iteration)
  all_combinations <- expand.grid(Iteration = Iterations, Time = all_times)
  
  # Merge the original data with the complete time series
  extended_df <- all_combinations %>%
    left_join(df, by = c("Iteration", "Time"))
  
  # Get the last row for each Iteration to fill in new entries
  last_rows <- df %>%
    group_by(Iteration) %>%
    filter(Time == max(Time, na.rm = TRUE)) %>%
    slice(1) %>%
    dplyr::select(-Time, -Exposures, -Infections)
  
  # Fill in missing values for other columns
  extended_df <- extended_df %>%
    group_by(Iteration) %>%
    arrange(Iteration, Time) %>%
    fill(-Exposures, -Infections, .direction = "downup") %>%
    ungroup()
  
  # For the new time points, set Exposures and Infections to 0
  extended_df <- extended_df %>%
    mutate(Exposures = ifelse(is.na(Exposures), 0, Exposures),
           Infections = ifelse(is.na(Infections), 0, Infections))
  
  # Ensure that columns other than Exposures and Infections have the last row values
  extended_df <- extended_df %>%
    left_join(last_rows, by = "Iteration", suffix = c("", "_last")) %>%
    mutate(across(ends_with("_last"), ~ coalesce(., get(sub("_last", "", cur_column()))))) %>%
    dplyr::select(-ends_with("_last"))
  
  # Remove any rows that are exact duplicates of the original dataframe
  extended_df <- distinct(bind_rows(df, anti_join(extended_df, df, by = c("Iteration", "Time"))))

  return(extended_df)
}
