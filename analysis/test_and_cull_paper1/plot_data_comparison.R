plot_data_comparison = function (data, UID, params, abc_results, niter=100, Tmax=10*365) {
  
  print(UID)
  current_data <- data[data$UID == UID,]

  
  effective_Se <- params$prop_infectious * params$sensitivity_infectious + (1-params$prop_infectious) * params$sensitivity_latent
  
  params$N0 <- round(median(current_data$scans))
  params$initial_prevalence = params$initial_prevalence <- pmax(1/params$N0, (current_data$App.Prev[1] - (1-params$specificity)) / (effective_Se - (1-params$specificity)))
  params$scanning_timepoints = unique(current_data$Year * 365)
  
  current_abc <- abc_results[abc_results$UID != UID,]
  
  abc_summary <- current_abc %>%
    group_by(UID) %>%
    summarise(SE = sd(R0)/sqrt(length(R0)),
              R0 = mean(R0)              )
  
  # abc_summary must have one row per UID with:
  #   yi = posterior mean R0 for UID
  #   SE = posterior standard error for UID
  # (rename or adjust columns if necessary)
  # e.g. abc_summary <- data.frame(UID=..., yi = mean_R0, SE = sd_R0)
  
  # Fit random-effects meta-analysis (REML)
  res <- rma(yi = R0, sei = SE, method = "REML", data = abc_summary)
  
  # Extract estimated population mean (mu) and between-UID tau
  mu_hat  <- as.numeric(coef(res))     # population mean
  tau_hat <- sqrt(res$tau2)            # between-UID SD
  
  # Simulate 100 new farms' latent true R0 (theta_i)
  n_new <- niter
  theta_new <- rnorm(n_new, mean = mu_hat, sd = tau_hat)   # latent true R0 for each new farm
  
  # If you want "observed" R0 draws that include measurement noise for a new observed estimate:
  # choose an observation-level SD for a new estimate. Reasonable choices:
  #  - median SE from your summarized UIDs: median_SE <- median(abc_summary$SE)
  #  - or set to 0 if you want the latent theta only.
  obs_sd_new <- median(abc_summary$SE, na.rm = TRUE)
  
  R0_obs_new <- rnorm(n_new, mean = theta_new, sd = obs_sd_new)
  R0_obs_new[R0_obs_new < 0] <- 0.1 # This very occasionally gives negative R0s, so replace them with a small positive value
  
  
  temporal_infections <- vector("list", niter)
  for(i in 1:niter){
    params$R0 <- R0_obs_new[i]
    
    infections <- simulate_OPA(params, Tmax = max(params$scanning_timepoints) + 1)
    infections <- na.omit(infections)
    if(nrow(infections) == 0){
      infections <- data.frame(Time=0, Exposures=0, Infections=0, Removals=0,
                               Incorrect.Culls=0, Pop.Size=params$N0)
    }
    
    infections <- process_results(infections = infections,
                                  N0 = params[["N0"]],
                                  time_between_restocking = params[["time_between_restocking"]])
    infections$Iteration <- i

    temporal_infections[[i]] <- infections
  }
  
  temporal_infections <- bind_rows(temporal_infections)
  
  temporal_infections$Iteration <- factor(temporal_infections$Iteration)
  
  # Apply the function to extend the time series
  temporal_infections <- extend_data(temporal_infections, time_interval = params[["scanning_timepoints"]], max_time = params[["scanning_timepoints"]][length(params[["scanning_timepoints"]])])
  
  Positives <- rbinom(nrow(temporal_infections), size = temporal_infections$Exposures, prob = params$sensitivity_latent) +
               rbinom(nrow(temporal_infections), size = temporal_infections$Infections, prob = params$sensitivity_infectious) +
               rbinom(nrow(temporal_infections), size = temporal_infections$Pop.Size - temporal_infections$Exposures - temporal_infections$Infections, prob = 1 - params$specificity)
  
  temporal_infections$AppPrev <- Positives / temporal_infections$Pop.Size
  
  # If the first scan (at t=0) removed all the OPA cases then we end up with an issue where the population size is becoming NA and creates other NAs later. Given this only happens when exposures+infections=0 I'll
  # just make it =1000 and the prevalence will always come out as 0.
  
  #print("next")
  #if(sum(is.na(temporal_infections$Pop.Size)) > 0){
  #View(temporal_infections)}
  #print(sum(is.na(temporal_infections$Pop.Size)))
  temporal_infections[is.na(temporal_infections$Pop.Size),]$AppPrev <- 0
  temporal_infections[is.na(temporal_infections$Pop.Size),]$Pop.Size <- 1000
  
  temporal_infections <- temporal_infections[order(temporal_infections$Iteration),]
  
  testtimes <- params[["scanning_timepoints"]]
  medians <- numeric(length(testtimes))
  ucl <- numeric(length(testtimes))
  lcl <- numeric(length(testtimes))
  for(i in 1:length(testtimes)){
    thistime <- temporal_infections[which(temporal_infections$Time == testtimes[i])-1,]
    prevs <- thistime$AppPrev
    prevs[is.na(prevs)] <- 0
    medians[i] <- median(prevs)
    ucl[i] <- quantile(prevs, 0.9)
    lcl[i] <- quantile(prevs, 0.1)
  }
  
  summary_of_vals <- data.frame(Time=testtimes, medians=medians, lcl=lcl, ucl=ucl)
  
  colors <- c("Infectious" = "red", "Total" = "black")
  
  p.temp <- ggplot(temporal_infections, aes(x=Time/365)) +
    #geom_line(aes(y=AppPrev*100, col="Total", group=Iteration), alpha=0.05) +
    ylab("") +
    xlab("") +
    xlim(-10/365,max(params[["scanning_timepoints"]])/365) +
    #ylim(0,max((temporal_infections[(temporal_infections$Time>=-10) & (temporal_infections$Time<=max(params[["scanning_timepoints"]])),]$Exposures+temporal_infections[(temporal_infections$Time>=-10) & (temporal_infections$Time<=max(params[["scanning_timepoints"]])),]$Infections) / temporal_infections[(temporal_infections$Time>=-10) & (temporal_infections$Time<=max(params[["scanning_timepoints"]])),]$Pop.Size)) +
    scale_color_manual(values = colors) + 
    guides(col=guide_legend(title="State")) + 
    theme_minimal(base_size=12) +
    theme(legend.position = "top",
          panel.background = element_rect(fill = 'white', colour = 'white'),
          panel.grid.major = element_line(colour = "lightgrey"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=10),
          plot.margin = margin(0, 0, 0, 0, "cm")) +
    geom_point(data=current_data, aes(x=Year, y=App.Prev*100), col="red", size=1.4) +
    ggtitle(UID)
  
  # desired limit
  potential_lim <- c(0, max(summary_of_vals$ucl*100))
  
  # calculate natural data range
  existing_lim <- c(0, max(current_data$upr*100))
  
  # final limits: pick the wider of the two
  final_lim <- c(min(potential_lim[1], existing_lim[1]),
                 max(potential_lim[2], existing_lim[2]))
  
  p.temp <- p.temp +
    geom_line(data=summary_of_vals, aes(y=medians*100), linewidth=1) +
    geom_line(data=summary_of_vals, aes(y=lcl*100), linetype=2, linewidth=1) +
    geom_line(data=summary_of_vals, aes(y=ucl*100), linetype=2, linewidth=1) +
    ylim(final_lim)
  
  p.temp
}
