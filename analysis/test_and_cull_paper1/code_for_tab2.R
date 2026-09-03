library(TransmitTestandCull)
library(dplyr)

# This looks specifically at the "average" farm in terms of R0
# Create parameter list
params = list(N0 = 1500,
              initial_prevalence = 0.03,
              sensitivity_infectious = 0.76,
              sensitivity_latent = 0.76/2,
              specificity = 0.998,
              prop_infectious = 0.33,
              scanning_timepoints = seq(0,365*10,182.5),
              time_between_restocking = 365,
              R0 = 1.39,
              mu_L = 180,
              sh_L = 30,
              mu_I = 90,
              sh_I = 10)

# Set number of iterations of the model to run and max length of simulation (in days)
niter=1000
Tmax=10*365

infections_list <- vector("list", niter)
set.seed(4398)
for(i in 1:niter){
  infections_list[[i]] <- na.omit(simulate_OPA(params, Tmax=Tmax))
  infections_list[[i]]$Iter <- i
}
infections <- bind_rows(infections_list)
average_farm_summary <- infections %>%
  group_by(Iter) %>%
  summarise(Max.time = max(Removal),
            OPA_pos_culls = sum((Incorrect==0)*(Removal < 365*10)),
            OPA_neg_culls = sum((Incorrect==1)*(Removal < 365*10)))

median(average_farm_summary$Max.time)/365
quantile(average_farm_summary$Max.time,0.025)/365
quantile(average_farm_summary$Max.time,0.975)/365

median(average_farm_summary$OPA_pos_culls)
quantile(average_farm_summary$OPA_pos_culls,0.025)
quantile(average_farm_summary$OPA_pos_culls,0.975)

median(average_farm_summary$OPA_neg_culls)
quantile(average_farm_summary$OPA_neg_culls,0.025)
quantile(average_farm_summary$OPA_neg_culls,0.975)

params = list(N0 = 1500,
              initial_prevalence = 0.03,
              sensitivity_infectious = 0.76,
              sensitivity_latent = 0.76/2,
              specificity = 0.998,
              prop_infectious = 0.33,
              scanning_timepoints = seq(0,365*10,365),
              time_between_restocking = 365,
              R0 = 1.39,
              mu_L = 180,
              sh_L = 30,
              mu_I = 90,
              sh_I = 10)

# Set number of iterations of the model to run and max length of simulation (in days)
niter=1000
Tmax=10*365

infections_list <- vector("list", niter)
set.seed(4398)
for(i in 1:niter){
  infections_list[[i]] <- na.omit(simulate_OPA(params, Tmax=Tmax))
  infections_list[[i]]$Iter <- i
}
infections <- bind_rows(infections_list)
average_farm_summary <- infections %>%
  group_by(Iter) %>%
  summarise(Max.time = max(Removal),
            OPA_pos_culls = sum((Incorrect==0)*(Removal < 365*10)),
            OPA_neg_culls = sum((Incorrect==1)*(Removal < 365*10)))

median(average_farm_summary$Max.time)/365
quantile(average_farm_summary$Max.time,0.025)/365
quantile(average_farm_summary$Max.time,0.975)/365

median(average_farm_summary$OPA_pos_culls)
quantile(average_farm_summary$OPA_pos_culls,0.025)
quantile(average_farm_summary$OPA_pos_culls,0.975)

median(average_farm_summary$OPA_neg_culls)
quantile(average_farm_summary$OPA_neg_culls,0.025)
quantile(average_farm_summary$OPA_neg_culls,0.975)