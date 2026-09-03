library(TransmitTestandCull)
library(dplyr)
library(ggplot2)

source("fit_qribbons.R")
source("plot_data_comparison.R")

# Create parameter list
params = list(N0 = 1500,
              initial_prevalence = 0.03,
              sensitivity_infectious = 0.76,
              sensitivity_latent = 0.38,
              specificity = 0.998,
              prop_infectious = 0.33,
              scanning_timepoints = 182.5,
              time_between_restocking = 365,
              R0 = 1.39,
              mu_L = 180,
              sh_L = 30,
              mu_I = 90,
              sh_I = 10)

niter=1000
Tmax=10*365

parname1 = "scanning_timepoints"
parname2 = "sensitivity_infectious"
parvals1 = c(182.5,365)
parvals2 = c(0.5,0.95)

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


aggregated_infections <- model_results$aggregated_infections
aggregated_infections$scanning_timepoints <- factor(aggregated_infections$scanning_timepoints, levels=parvals1, labels=c("6-monthly","Yearly"))
saveRDS(aggregated_infections, file = "../Results/scan_sensitivity.rds")

qribbons <- aggregated_infections %>%
  group_by(scanning_timepoints) %>%
  group_map(~ {
    res <- fit_qribbons(.x, xvar = "sensitivity_infectious", yvar = "Time", lambda = 1)
    res$scanning_timepoints <- .y$scanning_timepoints
    res
  }) %>%
  bind_rows()

p1 <- ggplot() +
  geom_point(data = aggregated_infections,
             aes(x = sensitivity_infectious,
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
  xlab("Scan sensitivity during infectious period") +
  guides(col = guide_legend(title = "Scanning Freq."),
         fill = "none") +
  coord_cartesian(ylim = c(0, 10), xlim = c(0.5,0.95)) +
  geom_vline(xintercept=0.76, linetype=2)

#png("../Figures/ScanFreq_infsensitivity.png", width=9, height=6, units="in", res=600)
#p1
#dev.off()

# Create default parameter list
params = list(N0 = 1500,
              initial_prevalence = 0.03,
              sensitivity_infectious = 0.76,
              sensitivity_latent = 0.38,
              specificity = 0.998,
              prop_infectious = 0.33,
              scanning_timepoints = 182.5,
              time_between_restocking = 365,
              R0 = 1.39,
              mu_L = 180,
              sh_L = 30,
              mu_I = 90,
              sh_I = 10)

# Set number of iterations
niter=1000
Tmax=10*365

# Define two parameters to vary
parname1 = "scanning_timepoints"
parname2 = "specificity"
parvals1 = c(182.5,365)
parvals2 = c(0.95,1)

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


aggregated_infections <- model_results$aggregated_infections
aggregated_infections$scanning_timepoints <- factor(aggregated_infections$scanning_timepoints, levels=parvals1, labels=c("6-monthly","Yearly"))
#saveRDS(aggregated_infections, file = "../Results/scan_specificity.rds")

qribbons <- aggregated_infections[aggregated_infections$Diagnosis == "Incorrect",] %>%
  group_by(scanning_timepoints) %>%
  group_map(~ {
    res <- fit_qribbons(.x, xvar = "specificity", yvar = "Culls", lambda = 1)
    res$scanning_timepoints <- .y$scanning_timepoints
    res
  }) %>%
  bind_rows()

p2 <- ggplot() +
  geom_point(data = aggregated_infections[aggregated_infections$Diagnosis == "Incorrect",],
             aes(x = specificity,
                 y = Culls,
                 col = scanning_timepoints),
             alpha=0.4) +
  geom_ribbon(data = qribbons,
              aes(x = x,
                  ymin = lwr,
                  ymax = upr,
                  fill = scanning_timepoints),
              alpha = 0.2) +
  geom_line(data = qribbons,
            aes(x = x,
                y = fit,
                col = scanning_timepoints),
            size = 1) +
  theme_bw(base_size = 16) +
  ylab("Culls of OPA-negative animals (flock size = 1500)") +
  xlab("Scan specificity") +
  guides(col = guide_legend(title = "Scanning Freq."),
         fill = "none") +
  xlim(0.95, 1) +
  ylim(0,1300) +
  geom_vline(xintercept=0.998, linetype=2)

#png("../Figures/sens_and_spec.png", width=12, height=6, units="in", res=600)
ggarrange(p1,p2,ncol=2,common.legend=T,legend="bottom", labels = "AUTO")
#dev.off()