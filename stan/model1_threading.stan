// Model 1 with multi-threading to further parallelize sampling
// Uniform prior on shift parameter with limits of -17 and 27
// The first PCR within an infection is defined as Day 0.

functions {
  real partial_sum(int[] gstart_DAY_slice,
                   int start_idx, int end_idx,
                   int[] gend_DAY,
                   vector Y_DAY_imputed,
                   vector X_DAY, 
                   vector shift,
                   real alpha_mu,
                   vector intercept, 
                   vector b,
                   vector c,
                   vector sigma
                   ) {
    real logpd_sum = 0;
    for(g2 in start_idx:end_idx){
      int gstart_idx = g2 - start_idx + 1;
      int N_g = gend_DAY[g2] - gstart_DAY_slice[gstart_idx] + 1;
      vector[N_g] DAY_shifted1 = X_DAY[gstart_DAY_slice[gstart_idx]:gend_DAY[g2]] + shift[g2];
      vector[N_g] yhat1;
      yhat1 = intercept[g2] + b[g2] * DAY_shifted1 + c[g2] * DAY_shifted1 + 2 / alpha_mu * c[g2] * log1p_exp(-alpha_mu * DAY_shifted1); //see https://mc-stan.org/docs/2_21/functions-reference/composed-functions.html; old version: log(1 + exp(-alpha *  DAY_shifted1))
      logpd_sum += normal_lpdf(Y_DAY_imputed[gstart_DAY_slice[gstart_idx]:gend_DAY[g2]] | yhat1, sigma[g2]);
    }
    return logpd_sum;
  }
}


data {
  // Data for time course analysis
  int N_DAY;                                         // number of data points
  vector[N_DAY] Y_DAY;                               // outcome log10 viral load
  vector[N_DAY] X_DAY;                               // day of measurement in time series (day 0 is day of first PCR)
  int G;                                             // number of infections
  array[G] int gstart_DAY;                           // first day in time series for each infection; index for X_DAY and Y_DAY
  array[G] int gend_DAY;                             // last day in time series for each infection; index for X_DAY and Y_DAY
  int N_NegTests;                                    // total number negative test results
  array[N_NegTests] int idx_NegTests;                // position of negative test results in X_DAY and Y_DAY
  real L_median;                                     // upper limit for imputed viral loads for neg tests (based on median log10 vl of lowest 1% of positive PCRs)
  int K_PGH;                                         // number of covariates
  matrix[G,K_PGH] X_PGH;                             // covariates 
  int condition_on_data;                             // 1 for parameter estimation, 0 for prior predictive plots
  int grainsize;                                     // size of each data chunk processed in parallel
}

parameters {
  // grand means for key model parameters
  real log_slope_up_mu;
  real log_slope_down_mu;
  real log_intercept_mu;

  // covariates for key model parameters
  vector[K_PGH] betaPGH_intercept;
  vector[K_PGH] betaPGH_slope_up;
  vector[K_PGH] betaPGH_slope_down;

  // infection random effects for key model parameters
  real<lower=0> intercept_sigma;
  vector[G] intercept_raw; 
  real<lower=0> slope_up_sigma;
  vector[G] slope_up_raw;
  real<lower=0> slope_down_sigma;
  vector[G] slope_down_raw;
  real<lower=0> sigma_sigma;
  vector[G] sigma_raw;

  // smoothness parameter
  real<lower=0> alpha_mu;

  // shifts for infections
  vector<lower=-17,upper=27>[G] shift;

  // error variance
  real sigma_mu; // mean

  // imputed values for negative tests
  vector<lower=-200,upper=L_median>[N_NegTests] imp_neg;
}

transformed parameters {
  // error variance
  vector[G] sigma = exp(sigma_mu + sigma_sigma*sigma_raw);
  // initialize vectors for infection level parameters
  vector<lower=0>[G] intercept;
  vector<lower=0>[G] slope_up;
  vector<upper=0>[G] slope_down;
  vector[G] b;
  vector[G] c;
  { // calculate model parameters
     vector[G] log_intercept = 
                 log_intercept_mu + // intercept (fixed effect [FE], group level mean)
                 intercept_sigma * intercept_raw + // infection level random effects
                 X_PGH * betaPGH_intercept; // covariates (FE)
     vector[G] log_slope_up = 
                 log_slope_up_mu +
                 slope_up_sigma * slope_up_raw + 
                 X_PGH * betaPGH_slope_up; 
     vector[G] log_slope_down = 
                 log_slope_down_mu + 
                 slope_down_sigma * slope_down_raw + 
                 X_PGH * betaPGH_slope_down;
     // apply link function
     intercept = exp(log_intercept);
     slope_up = exp(log_slope_up);
     slope_down = -exp(log_slope_down);
  }

  b = (slope_up + slope_down) / 2;
  c = (slope_down - slope_up) / 2;
}

model {
  // DAY
  log_intercept_mu ~ normal(2.2,.15); 
  log_slope_up_mu ~ normal(.5,.3);
  log_slope_down_mu ~ normal(-.8,.5);

  intercept_sigma ~ gamma(4, 20);
  intercept_raw ~ std_normal();

  slope_up_sigma ~ gamma(5, 10);
  slope_up_raw ~ std_normal();

  slope_down_sigma ~ gamma(5, 10);
  slope_down_raw ~ std_normal();

  sigma_mu ~ std_normal();
  sigma_sigma ~ std_normal();
  sigma_raw ~ std_normal();

  betaPGH_intercept ~ normal(0,.125);
  betaPGH_slope_up ~ normal(0,.3);
  betaPGH_slope_down ~ normal(0,.25);

  alpha_mu ~ gamma(2, .25);

  if (condition_on_data == 1) {
    vector[N_DAY] Y_DAY_imputed = Y_DAY;
    for (j in 1:N_NegTests) {
      Y_DAY_imputed[idx_NegTests[j]] = imp_neg[j];
    }
  
    // Stan will choose data subset size automatically (for threading) when 
    // grainsize is set to 1.
    target += reduce_sum(partial_sum, gstart_DAY, grainsize,
                         gend_DAY, Y_DAY_imputed, X_DAY, shift, 
                         alpha_mu, intercept, b, c, sigma);
  }
}

generated quantities {
  vector<lower=0>[G] time2peak = intercept ./ slope_up;
  vector<lower=0>[G] time_from_peak = -intercept ./ slope_down;


  real slope_up_mu = exp(log_slope_up_mu);
  real slope_down_mu = -exp(log_slope_down_mu);
  real intercept_mu = exp(log_intercept_mu);

}