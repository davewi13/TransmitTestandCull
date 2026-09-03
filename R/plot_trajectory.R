#' @export
plot_trajectory <- function(temporal_infections,
                            Tmax,
                            time_between_scans){
  # This function takes the output from process_results() (which in turn come from simulate_OPA()) and plots the trajectory of the prevalences through time
  # INPUTS:
  # temporal_infections = the data frame outputted from process_results()
  # Tmax = The maximum amount of time the simulations were run for
  # time_between_scans = the duration between successive diagnostic scans (in days)
  
  colors <- c("Infectious" = "red", "Total" = "black")
  
  if("Iteration" %in% colnames(temporal_infections)){
    temporal_infections$Iteration <- factor(temporal_infections$Iteration)
      
    # Apply the function to extend the time series
    temporal_infections <- extend_data(temporal_infections, time_interval = time_between_scans, max_time = Tmax)
    temporal_infections <- temporal_infections[order(temporal_infections$Iteration),]
    if(length(time_between_scans) == 1){
      testtimes <- seq(0, Tmax, by = time_between_scans)
    } else {
      testtimes <- time_between_scans
    }
    medians <- numeric(length(testtimes))
    ucl <- numeric(length(testtimes))
    lcl <- numeric(length(testtimes))
    for(i in 1:length(testtimes)){
      thistime <- temporal_infections[which(temporal_infections$Time == testtimes[i])-1,]
      prevs <- (thistime$Exposures+thistime$Infections)/thistime$Pop.Size
      prevs[is.na(prevs)] <- 0
      medians[i] <- median(prevs)
      ucl[i] <- quantile(prevs, 0.975)
      lcl[i] <- quantile(prevs, 0.025)
    }
    summary_of_vals <- data.frame(Time=testtimes, medians=medians, lcl=lcl, ucl=ucl)
      
    p.temp <- ggplot(temporal_infections, aes(x=Time/365)) +
      geom_line(aes(y=(Exposures+Infections)/Pop.Size, col="Total", group=Iteration), alpha=0.05) +
      ylab("Prevalence") +
      xlab("Time (Years)") +
      xlim(-100/365,min(Tmax/365,max(temporal_infections$Time/365))) +
      ylim(0,max((temporal_infections[temporal_infections$Time>=-100,]$Exposures+temporal_infections[temporal_infections$Time>=-100,]$Infections)/temporal_infections[temporal_infections$Time>=-100,]$Pop.Size)) +
      scale_color_manual(values = colors) + 
      guides(col=guide_legend(title="State")) +  
      theme_bw(base_size=20) +
      theme(legend.position = "none")
      
    p.temp <- p.temp +
      geom_line(data=summary_of_vals, aes(y=medians), linewidth=2) +
      geom_line(data=summary_of_vals, aes(y=lcl), linetype=2, linewidth=2) +
      geom_line(data=summary_of_vals, aes(y=ucl), linetype=2, linewidth=2)
  
  } else {
    # Apply the function to extend the time series
    temporal_infections <- extend_data(temporal_infections, time_interval = time_between_scans, max_time = Tmax)
    if(length(time_between_scans) == 1){
      testtimes <- seq(0, Tmax, by = time_between_scans)
    } else {
      testtimes <- time_between_scans
    }
    for(i in 1:length(testtimes)){
      thistime <- temporal_infections[which(temporal_infections$Time == testtimes[i])-1,]
      prevs <- (thistime$Exposures+thistime$Infections)/thistime$Pop.Size
      prevs[is.na(prevs)] <- 0
    }
    
    p.temp <- ggplot(temporal_infections, aes(x=Time/365)) +
      geom_line(aes(y=(Exposures+Infections)/Pop.Size, col="Total", group=Iteration), alpha=0.05) +
      ylab("Prevalence") +
      xlab("Time (Years)") +
      xlim(-100/365,min(Tmax/365,max(temporal_infections$Time/365))) +
      ylim(0,max((temporal_infections[temporal_infections$Time>=-100,]$Exposures+temporal_infections[temporal_infections$Time>=-100,]$Infections)/temporal_infections[temporal_infections$Time>=-100,]$Pop.Size)) +
      scale_color_manual(values = colors) + 
      guides(col=guide_legend(title="State")) +  
      theme_bw(base_size=20) +
      theme(legend.position = "none")
  }  
  p.temp
}
