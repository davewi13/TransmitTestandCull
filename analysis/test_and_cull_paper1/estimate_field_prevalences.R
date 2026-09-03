# This code takes the OPA scanning data from the .csv file and adds estimates of the true prevalence with uncertainty.

# Load the data
opa_dat <- read.csv("C:/Users/de43382/opa-stochastic-model/Data/processed_OPA_data.csv", header=T)

# Note that everything from HERE to END is redundant because we are using apparent prevalence. It could probably be deleted but I have left it in case its ever useful in future

# Sample nsim values for the sensitivity and specificity from the assumed distributions for these parameters
set.seed(89345897)
nsim=1000
sensitivities <- rbeta(nsim, shape1=477, shape2=150)
specificities <- rbeta(nsim, shape1=32500, shape2=50)
prev.mat <- matrix(nrow=nrow(opa_dat), ncol=nsim)

# Calculate the prevalence based on the sampled sensitivities and specificities
for(i in 1:nsim){
  prev.mat[,i] <- ((opa_dat$pos/(opa_dat$scans + round(runif(length(opa_dat$scans), -25.499, 25.499)))) + specificities[i] - 1)/(sensitivities[i] + specificities[i] - 1)
}
prev.mat[prev.mat < 0] <- 0

# Calculate an estimate and 95% credible interval for the true prevalence on each farm at each scan
opa_dat$True.Prev <- apply(prev.mat, 1, median)
opa_dat$lwr <- apply(prev.mat, 1, function(x){quantile(x,0.025, na.rm=T)})  
opa_dat$upr <- apply(prev.mat, 1, function(x){quantile(x,0.975, na.rm=T)})

# END - everything here HERE to this point is redundant

opa_dat$App.Prev <- opa_dat$pos/opa_dat$scans

opa_dat <- opa_dat %>%
  group_by(UID) %>%
  mutate(N.scans = sum(!is.na(pos)))

opa_dat <- opa_dat[complete.cases(opa_dat), ]

# Save the output
save(opa_dat, file="../Data/opa_dat.Rdata")