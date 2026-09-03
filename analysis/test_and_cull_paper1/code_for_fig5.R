library(TransmitTestandCull)
library(dplyr)
library(ggplot2)

source("fit_qribbons.R")
source("plot_data_comparison.R")

# Create parameter list
params = list(N0 = 1500,
              initial_prevalence = 0.03,
              sensitivity_infectious = 0.76,
              sensitivity_latent = 0.76/2,
              specificity = 0.998,
              prop_infectious = 0.33,
              scanning_timepoints = 182.5,
              time_between_restocking = 365,
              R0 = 1.5,
              mu_L = 180,
              sh_L = 30,
              mu_I = 90,
              sh_I = 10)

# Set number of iterations of the model to run and max length of simulation (in days)
niter=1000
Tmax=10*365

# Define two parameters to vary. The run_parameter_combo function is fairly hard coded to our needs
# here and care should be taken in running it. parname1 is set to be discrete and so only the two
# values passed will be used. parname2 is set to be continuous and so you can specify one of a couple
# of distributions. Here we pass a uniform distribution so the two numbers passed are the min and max
# Note this takes a couple of hours to run
parname1 = "scanning_timepoints"
parname2 = "R0"
parvals1 = c(182.5,365)
parvals2 = c(0.3,2.2)
set.seed(134905)
model_results <- run_parameter_combo(params=params,
                                     parname1=parname1,
                                     parname2=parname2,
                                     parvals1=parvals1,
                                     parvals2=parvals2,
                                     niter=niter,
                                     Tmax=Tmax,
                                     continuous_par2 = T,
                                     continuous_dist = 'runif')

# Reformat some of the output and create a variable stating whether R0 < 1 or not
aggregated_infections <- model_results$aggregated_infections
aggregated_infections$scanning_timepoints <- factor(aggregated_infections$scanning_timepoints, levels=parvals1, labels=c("6-monthly","Yearly"))
aggregated_infections$R0.gr1 <- ifelse(aggregated_infections$R0 >= 1, "R0 > 1", "R0 < 1")

# This performs quantile regression smoothing to add prediction bands to the plots
qribbons <- aggregated_infections %>%
  group_by(scanning_timepoints) %>%
  group_map(~ {
    res <- fit_qribbons(.x, xvar = "R0", yvar = "Time", lambda = 1)
    res$scanning_timepoints <- .y$scanning_timepoints
    res
  }) %>%
  bind_rows()

# Makes the top plot of Fig 5
p1 <- ggplot() +
  geom_point(data = aggregated_infections,
             aes(x = R0,
                 y = Time/365,
                 col = scanning_timepoints),
             alpha=0.4) +
  geom_ribbon(data = qribbons,
              aes(x = x,
                  ymin = lwr/365,
                  ymax = upr/365,
                  fill = scanning_timepoints),
              alpha = 0.2) +
  geom_line(data = qribbons,
            aes(x = x,
                y = fit/365,
                col = scanning_timepoints),
            size = 1) +
  theme_bw(base_size = 14) +
  ylab("Time to eradication (years)") +
  xlab("R0") +
  guides(col = guide_legend(title = "Scanning Freq."),
         fill = "none") +
  theme(legend.position = "bottom") +
  xlim(0.3, 2.2) +
  geom_vline(xintercept=1.39, linetype=2)


# This saves the plot
#png("../Figures/ScanFreq_R0.png", width=6, height=6, units="in", res=600)
p1
#dev.off()