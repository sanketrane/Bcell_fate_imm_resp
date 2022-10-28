functions{
   // function that describes the changes in CAR+ counts in FO B cells
  real[] theta_spline(real[] solve_time,  real[] parms){
    //real M0 = exp(14); real nu = 0.01; real b0 = 18;
    int numdim = size(solve_time);
    real M0 = parms[3];
    real nu  = parms[4];
    real b0  = parms[5];
    real value[numdim];

    for (i in 1:numdim){
      value[i] = exp(M0) * (1 + exp(-nu * (solve_time[i] - b0)^2));
    }
    return value;
   }

   // function that contains the ODE equations to be used in ODE solver
   real[] ODE_sys(real time,  real[] y, real[] parms, real[] rdata, int[] idata) {
     // the array of parameters invloved in ode equations
     real alpha = parms[1];
     real delta  = parms[2];
     // the system of ODEs
     real dydt[1];
     dydt[1] = alpha * theta_spline(time, parms) - delta * y[1];
     return dydt;
   }

   real[,] solve_ODE_sys(real[] solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     real y_sol[numdim, 1];
     y_sol = integrate_ode_rk45(ODE_sys, init_cond, 1.0, solve_time, parms, {0.0}, {0});
     return y_sol;
   }
}

data{
  int<lower  = 1> numObs;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> numPred;
  real<lower = 0> solve_time[numObs];
  real<lower = 0> theta_obs[numObs];
  real<lower = 0> total_counts[numObs];
  real ts_pred[numPred];
  }

parameters{
  // parameters to sample with boundary conditions
  real<lower = 0> alpha;
  real<lower = 0> delta;
  real<lower = 0> M0;
  real<lower = 0> nu;
  real<lower = 0> b0;
  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  }


transformed parameters{
  real y_hat1[numObs, 1];     // declaring the array for ODE solution
  real theta_mean[numObs];
  real total_counts_mean[numObs];
  real parms[5];                  // declaring the array for parameters
  real init_cond[1];              // declaring the array for state variables
  real G0 = exp(11.7);            // transformed parameters for better/faster sampling

  // initial conditions and parameters
  init_cond[1] = G0;
  parms[1] = alpha;
  parms[2] = delta;
  parms[3] = M0;
  parms[4] = nu;
  parms[5] = b0;
  // solution of the system of ODEs for the predictor values
  y_hat1[1] = init_cond;
  y_hat1[2:] = solve_ODE_sys(solve_time[2:], init_cond, parms);

  for (i in 1:numObs){
    theta_mean[i] = theta_spline(solve_time, parms);
    total_counts_mean[i] = y_hat1[i, 1];
  }
}

model{
  // prior distribution for model parameters
  alpha ~ normal(0.01, 0.25);
  delta ~ normal(0.01, 0.25);
  M0 ~ normal(14, 1);
  nu ~ normal(0.01, 0.05);
  b0 ~ normal(18,1);
  sigma1 ~ normal(0, 2.5);
  sigma2 ~ normal(0, 2.5);
  // model fitting on to data
  log(theta_obs) ~ normal(log(theta_mean), sigma1);
  log(total_counts) ~ normal(log(total_counts_mean), sigma2);
}

generated quantities{
  // ODE predictions
  real y_hat_pred1[numPred, 1];
  // model predictions
  real y1_mean_pred[numPred];
  real y2_mean_pred[numPred];
  // model predictions with stdev
  real theta_pred[numPred];
  real totalcounts_pred[numPred];
  // Residuals
  vector[numObs] resid_d1;
  // log likelihoods
  vector[numObs] log_lik1;
  //ODE solution
  y_hat_pred1[1] = init_cond;
  y_hat_pred1[2:] = solve_ODE_sys(ts_pred[2:], init_cond, parms);
  y2_mean_pred = theta_spline(ts_pred, parms)
  for (i in 1:numPred){
    //CAR MZ
    y1_mean_pred[i] = y_hat_pred1[i, 1];
    totalcounts_pred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma2));
    theta_pred[i] = exp(normal_rng(log(y2_mean_pred[i]), sigma1));
  }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    resid_d1[n] = log(total_counts[n]) - log(total_counts_mean[n]);
    log_lik1[n] = normal_lpdf(log(total_counts[n]) | log(total_counts_mean[n]), sigma1);
  }
}
