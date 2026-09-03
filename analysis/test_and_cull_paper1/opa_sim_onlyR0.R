#!/usr/bin/env Rscript
library(TransmitTestandCull)
library(doParallel)
library(purrr)
library(dplyr)
library(LaplacesDemon)

# -------------------------------
# Command line arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
chunk <- as.integer(args[1])

set.seed(1000 + chunk)  # reproducible randomness

# -------------------------------
# Load and prepare data
# -------------------------------
load("../Data/opa_dat.Rdata")

opa_dat <- opa_dat[opa_dat$N.scans > 2,]


# -------------------------------
# Parameters
# -------------------------------
params <- list(
  N0 = NULL,
  initial_prevalence = NULL,
  sensitivity_infectious = 0.76,
  sensitivity_latent = 0.76/4,
  specificity = 0.998,
  prop_infectious = 0.33,
  scanning_timepoints = NULL,
  time_between_restocking = 365,
  R0 = NULL,
  mu_L = 180,
  sh_L = 30,
  mu_I = 90,
  sh_I = 10
)

# -------------------------------
# Run 100 simulations for this fold/chunk
# -------------------------------
nsim <- 1000
out <- vector("list", nsim)


effective_Se = params$prop_infectious * params$sensitivity_infectious + (1 - params$prop_infectious) * params$sensitivity_latent


for (j in 1:nsim) {
  params$R0 <- runif(1, 0, 4)
  rmse <- vector("list", 3)
  rmse.logit <- vector("list", 3)
  mae <- vector("list", 3)
  mae.logit <- vector("list", 3)
  
  farm_inputs <- opa_dat %>%
  group_by(UID) %>%
  summarise(
      N0 = round(median(scans)),
      app0 = first(App.Prev),
      scanning_timepoints = list(unique(Year * 365)),
      app_prev = list(App.Prev),
      .groups = "drop"
  ) %>%
  mutate(
      initial_prevalence = pmax(
          1 / N0,
          (app0 - (1 - params$specificity)) /
          (effective_Se - (1 - params$specificity))
      )
  )

  
  for(i in 1:3){
    infections <- Map(
      function(n0, initial_prev, timepoints) {
        temp_params <- params
        temp_params[["N0"]] <- n0
        temp_params[["initial_prevalence"]] <- initial_prev
        temp_params[["scanning_timepoints"]] <- timepoints
        return(na.omit(simulate_OPA(temp_params, Tmax = max(timepoints)+0.1)))
      },
      farm_inputs$N0,
      farm_inputs$initial_prevalence,
      farm_inputs$scanning_timepoints
    )
  
    temporal_infections <- Map(
      function(inf, n0) {
        if(nrow(inf) > 0){
          res <- process_results(
            infections = inf,
            N0 = n0,
            time_between_restocking = params[["time_between_restocking"]]
          )
        } else {
          res <- data.frame(Time=0, Exposures=0, Infections=0, Removals=0,
                            Incorrect.Culls=0, Pop.Size=n0)
        }
        return(res)
      },
      infections,
      farm_inputs$N0
    )

  
    # -------------------------------
    # RMSE
    # -------------------------------
    rmse[[i]] <- mapply(function(sim_res, obs_prev, timepoints) {

      sim_res <- sim_res[which(sim_res$Time %in% timepoints) - 1, ]
      sim_prev <- (sim_res$Exposures * params[["sensitivity_latent"]] + sim_res$Infections * params[["sensitivity_infectious"]] + (sim_res$Pop.Size - sim_res$Exposures - sim_res$Infections) * (1-params[["specificity"]])) / sim_res$Pop.Size
      
      # pad if lengths differ
      if (length(sim_prev) != length(obs_prev)) {
        sim_prev <- c(sim_prev, rep(0, length(obs_prev) - length(sim_prev)))
      }
      
      sim_prev[sim_prev >= 1] <- 1
      obs_prev[obs_prev >= 1] <- 1
      val <- sum((sim_prev - obs_prev)^2)
      
      return(val)
    },
    temporal_infections,
    farm_inputs$app_prev,
    farm_inputs$scanning_timepoints
    )

    
    rmse.logit[[i]] <- mapply(function(sim_res, obs_prev, timepoints) {

      sim_res <- sim_res[which(sim_res$Time %in% timepoints) - 1, ]
      sim_prev <- (sim_res$Exposures * params[["sensitivity_latent"]] + sim_res$Infections * params[["sensitivity_infectious"]] + (sim_res$Pop.Size - sim_res$Exposures - sim_res$Infections) * (1-params[["specificity"]])) / sim_res$Pop.Size
      
      # pad if lengths differ
      if (length(sim_prev) != length(obs_prev)) {
        sim_prev <- c(sim_prev, rep(0, length(obs_prev) - length(sim_prev)))
      }
      
      # avoid logit(0)
      sim_prev[sim_prev <= 0] <- 0.0001
      obs_prev[obs_prev <= 0] <- 0.0001
      sim_prev[sim_prev >= 1] <- 0.9999
      obs_prev[obs_prev >= 1] <- 0.9999
        
      val <- sum((logit(sim_prev) - logit(obs_prev))^2)
      
      return(val)
    },
    temporal_infections,
    farm_inputs$app_prev,
    farm_inputs$scanning_timepoints
    )
    
    
    # -------------------------------
    # MAE
    # -------------------------------
    mae[[i]] <- mapply(function(sim_res, obs_prev, timepoints) {

      sim_res <- sim_res[which(sim_res$Time %in% timepoints) - 1, ]
      sim_prev <- (sim_res$Exposures * params[["sensitivity_latent"]] + sim_res$Infections * params[["sensitivity_infectious"]] + (sim_res$Pop.Size - sim_res$Exposures - sim_res$Infections) * (1-params[["specificity"]])) / sim_res$Pop.Size
      
      # pad if lengths differ
      if (length(sim_prev) != length(obs_prev)) {
        sim_prev <- c(sim_prev, rep(0, length(obs_prev) - length(sim_prev)))
      }
      
      val <- sum(abs(sim_prev - obs_prev))
      
      return(val)
    },
    temporal_infections,
    farm_inputs$app_prev,
    farm_inputs$scanning_timepoints
    )
    
    mae.logit[[i]] <- mapply(function(sim_res, obs_prev, timepoints) {

      sim_res <- sim_res[which(sim_res$Time %in% timepoints) - 1, ]
      sim_prev <- (sim_res$Exposures * params[["sensitivity_latent"]] + sim_res$Infections * params[["sensitivity_infectious"]] + (sim_res$Pop.Size - sim_res$Exposures - sim_res$Infections) * (1-params[["specificity"]])) / sim_res$Pop.Size
      
      # pad if lengths differ
      if (length(sim_prev) != length(obs_prev)) {
        sim_prev <- c(sim_prev, rep(0, length(obs_prev) - length(sim_prev)))
      }
      
      # avoid logit(0)
      sim_prev[sim_prev <= 0] <- 0.0001
      obs_prev[obs_prev <= 0] <- 0.0001
      sim_prev[sim_prev >= 1] <- 0.9999
      obs_prev[obs_prev >= 1] <- 0.9999
      
      val <- sum(abs(logit(sim_prev) - logit(obs_prev)))
      
      return(val)
    },
    temporal_infections,
    farm_inputs$app_prev,
    farm_inputs$scanning_timepoints
    )
  }


    
  out[[j]] <- data.frame(
    Sim = (chunk - 1) * nsim + j,
    R0 = params$R0,
    sensitivity_latent = params$sensitivity_latent,
    UID = farm_inputs$UID,
    RMSE = Reduce('+', rmse),
    RMSE.logit = Reduce('+', rmse.logit),
    MAE = Reduce('+', mae),
    MAE.logit = Reduce('+', mae.logit)
  )
}

# -------------------------------
# Save results for this chunk
# -------------------------------
out_df <- do.call(rbind, out)
dir.create("../Results_onlyR0_quarterlatent", showWarnings = FALSE)
saveRDS(out_df, file = sprintf("../Results_onlyR0_quarterlatent/chunk%d.rds", chunk))
