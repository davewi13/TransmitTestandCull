#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame simulate_OPA(List params, double Tmax) {

  // Extract parameters from the list
  int N0 = params["N0"];
  double initial_prevalence = params["initial_prevalence"];
  double sh_L = params["sh_L"];
  double mu_L = params["mu_L"];
  double sh_I = params["sh_I"];
  double mu_I = params["mu_I"];
  double R0 = params["R0"];
  double prop_infectious = params["prop_infectious"];
  double sensitivity_latent = params["sensitivity_latent"];
  double sensitivity_infectious = params["sensitivity_infectious"];
  double specificity = params["specificity"];
  NumericVector scanning_timepoints = params["scanning_timepoints"];
  double time_between_restocking = params["time_between_restocking"];

  // Handle scanning times
  // The first part of the if statement defines equally spaced scanning times if a single number is given as a scanning interval
  // If a vector is passed then the second part of the if statement sets the times from the vector to be the scanning times
  NumericVector scan_times_vec;
  if (scanning_timepoints.size() == 1) {
    int num_scans = floor(Tmax / scanning_timepoints[0]) + 1;
    scan_times_vec = NumericVector(num_scans);
    for (int i = 0; i < num_scans; ++i) {
      scan_times_vec[i] = i * scanning_timepoints[0];
    }
  } else {
    scan_times_vec = scanning_timepoints;
  }

  // Maximum initial size for vectors
  // This can be adjusted based on the simulation at hand. Too big and too much memory will be reserved but too small and repeated resizing slows things down.
  int max_size = 5000;
  
  // Initialize vectors to store infection data
  NumericVector Exposure(max_size, NA_REAL);
  NumericVector Infection(max_size, NA_REAL);
  NumericVector Removal(max_size, NA_REAL);
  NumericVector Incorrect(max_size, 0.0);
  
  // Estimate the number of individuals that are OPA positive at the start of the simulation
  int initial_infected = 0;
  while (initial_infected == 0) {
    initial_infected = R::rbinom(N0, initial_prevalence);
  }
  int initial_infectious = R::rbinom(initial_infected, prop_infectious);
  int initial_exposed = initial_infected - initial_infectious;
  
  // Initialize the number of individuals in each compartment
  int S = N0 - initial_infected;
  int E = initial_exposed;
  int I = initial_infectious;
  int R = 0;
  int N = S + E + I + R;

  // Calculate the rate parameters for the gamma distributed durations
  double rate_L = sh_L / mu_L;
  double rate_I = sh_I / mu_I;
  
  // Calculate the transmission rate
  double beta = R0 / mu_I;
  
  // Set the start time of the simulations
  double t0 = -0.5;
  
  // Loop through all of the exposed individuals and generate exposure, infection and removal times
  if (E != 0) {
    for (int i = 0; i < E; ++i) {
      double exp_dur = R::rgamma(sh_L, 1.0 / rate_L);
      double prop = R::runif(0, 1);
      Exposure[i] = t0 - prop * exp_dur;
      Infection[i] = t0 + (1 - prop) * exp_dur;
      Removal[i] = Infection[i] + R::rgamma(sh_I, 1.0 / rate_I);
    }
  }
  
  // Loop through all of the infectious individuals and generate exposure, infection and removal times
  if (I != 0) {
    for (int i = E; i < E + I; ++i) {
      double inf_dur = R::rgamma(sh_I, 1.0 / rate_I);
      double prop = R::runif(0, 1);
      Removal[i] = t0 + (1 - prop) * inf_dur;
      Infection[i] = t0 - prop * inf_dur;
      Exposure[i] = Infection[i] - R::rgamma(sh_L, 1.0 / rate_L);
    }
  }

  // Set the next free index to be the first individual that wasn't exposed or infectious
  int next_free_index = E + I; 
  
  // Set the current time to be t0
  double tnow = t0;
  std::string next_event;
  
  // This is the main loop that simulates the model. Stop if we have no infected individuals, everything dies or we reach the maximum simulation time.
  while ((E + I != 0) && (tnow < Tmax) && (N > 0)) {
    // Set the next event to be very far in the future to start
    double tnext = 999999;
    
    // Loop through all infections and find next soonest event after current time
    for (int i = 0; i < Infection.size(); ++i) {
      if (!NumericVector::is_na(Infection[i]) && Infection[i] > tnow && Infection[i] < tnext) {
        tnext = Infection[i];
        next_event = "Infection";
      }
      if (!NumericVector::is_na(Removal[i]) && Removal[i] > tnow && Removal[i] < tnext) {
        tnext = Removal[i];
        next_event = "Removal";
      }
    }
    
    // Calculate the rate of new infections
    double rate = beta * S * I / N;
    if (NumericVector::is_na(rate) || std::isnan(rate)) {
      stop("Transmission rate undefined. Likely that N=0");
    }
    
    // Define the time of the next new infection. If rate > 0 then we sample a new time and otherwise we just place it after the next progression for an infected individual
    double tnew = (rate > 0) ? tnow + R::rexp(1.0 / rate) : tnext + 1;

    // Find the next scan time after tnow
    double next_scan_time = 9999;  // Default to ages away
    for (int i = (scan_times_vec.size()-1); i >= 0; --i) {
        if (scan_times_vec[i] > tnow) {
            next_scan_time = scan_times_vec[i];
            if (scan_times_vec[i] < tnow) break;
        }
    }
    
    // If a scan time is ahead of the next event time, set tnow to the scan time
    if ((next_scan_time < tnext) && (next_scan_time < tnew)) {
      tnow = next_scan_time;
      // Set counters for what to remove to be zero
      int nrem_latent = 0;
      int nrem_infectious = 0;
      
      if (E != 0) {
        // Decide which individuals in exposed stage will be removed by scanner
        nrem_latent = R::rbinom(E, sensitivity_latent);
        IntegerVector latent_indices;
        for (int i = 0; i < Exposure.size(); ++i) {
          if (!NumericVector::is_na(Exposure[i]) && !NumericVector::is_na(Infection[i]) && Exposure[i] < tnow && Infection[i] > tnow) {
            latent_indices.push_back(i);
          }
        }

        // This does error catching to make sure that the size of the latent class vector and its counter stay in sync
        if ((int)latent_indices.size() != E) {
          // Print a clear diagnostic to stderr
          Rcpp::Rcerr << "=== COUNTER MISMATCH ERROR ===\n";
          Rcpp::Rcerr << "tnow = " << tnow << "\n";
          Rcpp::Rcerr << "Counters: E = " << E << ", I = " << I << ", S = " << S << ", R = " << R << ", N = " << N << "\n";
          Rcpp::Rcerr << "Vector-derived: latent_indices.size() = " << latent_indices.size()
               << ", Exposure.size() = " << Exposure.size() << "\n";

          for (int i = 0; i < Exposure.size(); i++) {
            if (!NumericVector::is_na(Exposure[i]) && 
                !NumericVector::is_na(Infection[i]) && 
                Exposure[i] <= tnow && Infection[i] > tnow) {
              // counted in latent_indices
              Rcpp::Rcout << "Counted individual: i=" << i 
                          << " Exposure=" << Exposure[i] 
                          << " Infection=" << Infection[i] 
                          << " Removal=" << Removal[i] << "\n";
            } else {
              Rcpp::Rcout << "Skipped individual: i=" << i 
                          << " Exposure=" << Exposure[i] 
                          << " Infection=" << Infection[i] 
                          << " Removal=" << Removal[i] << "\n";
            }
          }

            for (int i = 0; i < latent_indices.size(); i++){
            // Debug print for each latent individual
                  Rcpp::Rcout << "latent_indices: " << latent_indices[i] << " \n";
            }

          // count NA slots as additional hint
          int naExp = 0, naInf = 0, naRem = 0;
          for (int k = 0; k < Exposure.size(); ++k) {
            if (NumericVector::is_na(Exposure[k])) naExp++;
            if (NumericVector::is_na(Infection[k])) naInf++;
            if (NumericVector::is_na(Removal[k])) naRem++;
          }
          Rcpp::Rcerr << "NAs: Exposure=" << naExp << ", Infection=" << naInf << ", Removal=" << naRem << "\n";

          // Finally stop with an informative message
          Rcpp::stop("Invariant violation: counter E disagrees with vector contents. See stderr diagnostics above.");
        }

        // If any latent infections were detected then remove them. This means setting their infection and removal times to tnow.
        if (latent_indices.size() > 1) {
          IntegerVector sampled_latent_indices = sample(latent_indices, nrem_latent);
          for (int idx : sampled_latent_indices) {
            Infection[idx] = tnow;
            Removal[idx] = tnow;
          }
        } else if (latent_indices.size() == 1 && nrem_latent == 1) {
          int idx = latent_indices[0];
          Infection[idx] = tnow;
          Removal[idx] = tnow;
        }
      }
      
      // This repeats the above section but for individuals in the infectious rather than latent/exposed stage.
      if (I != 0) {
        nrem_infectious = R::rbinom(I, sensitivity_infectious);
        IntegerVector infectious_indices;
        for (int i = 0; i < Infection.size(); ++i) {
          if (!NumericVector::is_na(Infection[i]) && !NumericVector::is_na(Removal[i]) && Infection[i] < tnow && Removal[i] > tnow) {
            infectious_indices.push_back(i);
          }
        }

        if ((int)infectious_indices.size() != I) {
          // Print a clear diagnostic to stderr
          Rcpp::Rcerr << "=== COUNTER MISMATCH ERROR ===\n";
          Rcpp::Rcerr << "tnow = " << tnow << "\n";
          Rcpp::Rcerr << "Counters: E = " << E << ", I = " << I << ", S = " << S << ", R = " << R << ", N = " << N << "\n";
          Rcpp::Rcerr <<  "Vector-derived: infectious_indices.size() = " << infectious_indices.size()
               << ", Exposure.size() = " << Exposure.size() << "\n";

          // show first few infectious indices & their times
          int showN = std::min(10, (int)infectious_indices.size());
          for (int ii = 0; ii < showN; ++ii) {
            int idx = infectious_indices[ii];
            Rcpp::Rcerr << "infect idx[" << ii << "] = " << idx
                 << "  Exposure=" << Exposure[idx]
                 << "  Infection=" << Infection[idx]
                 << "  Removal=" << Removal[idx] << "\n";
          }

          // count NA slots as additional hint
          int naExp = 0, naInf = 0, naRem = 0;
          for (int k = 0; k < Infection.size(); ++k) {
            if (NumericVector::is_na(Exposure[k])) naExp++;
            if (NumericVector::is_na(Infection[k])) naInf++;
            if (NumericVector::is_na(Removal[k])) naRem++;
          }
          Rcpp::Rcerr << "NAs: Exposure=" << naExp << ", Infection=" << naInf << ", Removal=" << naRem << "\n";

          // Finally stop with an informative message
          Rcpp::stop("Invariant violation: counter I disagrees with vector contents. See stderr diagnostics above.");
        }


        if (infectious_indices.size() > 1) {
          IntegerVector sampled_infectious_indices = sample(infectious_indices, nrem_infectious);
          for (int idx : sampled_infectious_indices) {
            Removal[idx] = tnow;
          }
        } else if (infectious_indices.size() == 1 && nrem_infectious == 1) {
          int idx = infectious_indices[0];
          Removal[idx] = tnow;
        }
      }
      
      // We also (may) remove susceptibles, which is dependent on the specificity.
      int nrem_susceptible = R::rbinom(S, 1 - specificity);

      // Update all the counters
      E -= nrem_latent;
      I -= nrem_infectious;
      R += (nrem_latent + nrem_infectious);
      N -= (nrem_latent + nrem_infectious);
      
      // This deals with removing susceptibles and flags (via Incorrect) that these were animals without the disease
      if (nrem_susceptible > 0) {
        int start_idx = next_free_index;
        for (int i = 0; i < nrem_susceptible; ++i) {
          int idx = start_idx + i;
          Exposure[idx] = tnow;
          Infection[idx] = tnow;
          Removal[idx] = tnow;
          Incorrect[idx] = 1.0;
        }

        // advance next_free_index
        next_free_index += nrem_susceptible;

        S -= nrem_susceptible;
        R += nrem_susceptible;
        N -= nrem_susceptible;
      }
      
      // Restock back up to the original flock size at specified intervals
      if (std::fmod(tnow, time_between_restocking) == 0) {
        S = N0 - E - I;
        N = N0;
      }
      
    // This section deals with the case where the next event if a new exposure  
    } else if (tnew < tnext) {
      // Update time and counters
      tnow = tnew;
      S -= 1;
      E += 1;
      // Choose an individual to infect and generate a set of event times for it
      int current_idx = next_free_index++;
      Exposure[current_idx] = tnow;
      Infection[current_idx] = tnow + R::rgamma(sh_L, 1.0 / rate_L);
      Removal[current_idx]  = Infection[current_idx] + R::rgamma(sh_I, 1.0 / rate_I);
      if (S < 0) stop("Can't have negative number of susceptible individuals");
    } else {
      // This updates counters when next event is a transition from exposed to infectious class
      if (next_event == "Infection") {
        E -= 1;
        I += 1;
        if (E < 0) stop("Can't have negative number of exposed individuals");
      // This updates counters when next event is a transition from infectious to removed class
      } else if (next_event == "Removal") {
        I -= 1;
        R += 1;
        N -= 1;
        if (I < 0) stop("Can't have negative number of infectious individuals");
      } else {
      // This should never happen but prevents crashes if something completely unexpected has happened
        stop("Something went wrong");
      }
      // Update tnew to the current time following the event update
      tnow = tnext;
    }
    
    // Make the storage vectors larger if we are running out of space
    if (sum(is_na(Exposure)) < 10) {
      int old_size = Exposure.size();
      int new_size = old_size + N0;

      NumericVector new_Exposure(new_size, NA_REAL);
      NumericVector new_Infection(new_size, NA_REAL);
      NumericVector new_Removal(new_size, NA_REAL);
      NumericVector new_Incorrect(new_size, 0.0);

      std::copy(Exposure.begin(), Exposure.end(), new_Exposure.begin());
      std::copy(Infection.begin(), Infection.end(), new_Infection.begin());
      std::copy(Removal.begin(), Removal.end(), new_Removal.begin());
      std::copy(Incorrect.begin(), Incorrect.end(), new_Incorrect.begin());

      Exposure = new_Exposure;
      Infection = new_Infection;
      Removal = new_Removal;
      Incorrect = new_Incorrect;

      // no need to change next_free_index — it still points to the first unused slot
    }
  }
  
  // Remove the NA values from the vectors
  int valid_size = std::min(std::min(Exposure.size(), Infection.size()), std::min(Removal.size(), Incorrect.size()));
  NumericVector valid_Exposure(valid_size);
  NumericVector valid_Infection(valid_size);
  NumericVector valid_Removal(valid_size);
  NumericVector valid_Incorrect(valid_size);
  for (int i = 0; i < valid_size; ++i) {
    valid_Exposure[i] = Exposure[i];
    valid_Infection[i] = Infection[i];
    valid_Removal[i] = Removal[i];
    valid_Incorrect[i] = Incorrect[i];
  }
  
  // Create the DataFrame
  DataFrame result = DataFrame::create(
    Named("Exposure") = valid_Exposure,
    Named("Infection") = valid_Infection,
    Named("Removal") = valid_Removal,
    Named("Incorrect") = valid_Incorrect
  );
  
  return result;
}
