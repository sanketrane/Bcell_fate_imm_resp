functions{
   // function that describes the changes in CAR+ counts in FO B cells
   real CAR_positive_FOB(real time){
     real F0 = exp(11.763019); real B0 = exp(4.717021); real n = 5.092933 ; real X = 7.121932 ; real q = 6.475884;
     real value = F0 + (B0 * time^n) * (1 - ((time^q)/((X^q) + (time^q))));
     return value;
    }

   // function that contains the ODE equations to be used in ODE solver
   real[] ODE(real time,  real[] y, real[] parms, real[] rdata, int[] idata) {
     // the array of parameters invloved in ode equations
     real alpha1 = parms[1];
     real lambda = parms[2];
     real mu     = parms[3];
     real delta  = parms[4];
     real alpha2 = parms[5];
     real nu = parms[6];
     real CAR_negative_FOB = exp(17.15);

     real alpha_GC = alpha2/(1 + exp(nu * (time)^2));

     // the system of ODEs
     real dydt[3];
     dydt[1] = alpha1 * (CAR_negative_FOB) - lambda * y[1];
     dydt[2] = alpha_GC * CAR_positive_FOB(time) + mu * y[3]  - delta * y[2];
     dydt[3] = alpha_GC * CAR_negative_FOB - mu * y[3] - delta * y[3];
     return dydt;
   }

   real[,] solve_ODE(real[] solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     real y_sol[numdim, 3];
     y_sol = integrate_ode_rk45(ODE, init_cond, 4.0, solve_time, parms, {0.0}, {0});
     return y_sol;
   }

   // functions for transformation of fractions in (0,a), where a >=1
   real logit_inverse(real x){
      real ans;
      ans = exp(x)/(1+exp(x));
      return ans;
    }
}

data{
  int<lower  = 1> numObs;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> n_shards;
  int<lower  = 1> numPred;
  int<lower  = 1> time_index[numObs];
  real<lower = 0> solve_time[n_shards];
  real<lower = 0> CAR_MZ_counts[numObs];
  real<lower = 0> GC_counts[numObs];
  real<lower = 0> CAR_GC_counts[numObs];
  real ts_pred[numPred];
  }

parameters{
  // parameters to sample with boundary conditions
  real<lower = 0> alpha1;
  real<lower = 0> alpha2;
  real<lower = 0> delta;
  real<lower = 0> lambda;
  real<lower = 0> mu;
  real<lower = 0> nu;

  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  }


transformed parameters{
  real y_hat[n_shards, 3];     // declaring the array for ODE solution
  real y1_mean[n_shards];      // predictions for dataset1 based on ODE solution
  real y2_mean[n_shards];      // predictions for dataset2 based on ODE solution
  real y3_mean[n_shards];      // predictions for dataset3 based on ODE solution

  real CAR_MZcounts_mean[numObs];
  real CAR_GCcounts_mean[numObs];
  real GC_counts_mean[numObs];

  real parms[6];                  // declaring the array for parameters
  real init_cond[3];              // declaring the array for state variables

  real G0 = exp(12.0);            // transformed parameters for better/faster sampling
  real CAR_MZ0 = exp(10.8);            // transformed parameters for better/faster sampling
  real CAR_GC0 = exp(11.1);              // transformed parameters for better/faster sampling


  // initial conditions and parameters
  init_cond[1] = CAR_MZ0;
  init_cond[2] = CAR_GC0;
  init_cond[3] = G0 - CAR_GC0;

  parms[1] = alpha1;
  parms[2] = lambda;
  parms[3] = mu;
  parms[4] = delta;
  parms[5] = alpha2;
  parms[6] = nu;

  y_hat[1] = init_cond;
  // solution of the system of ODEs for the predictor values
  y_hat[2:] = solve_ODE(solve_time[2:], init_cond, parms);

  for (i in 1:n_shards){
  // CAR fractions in MZ
    y1_mean[i] = y_hat[i, 1];
    // CAR fractions in GC
    y2_mean[i] = y_hat[i, 2];
    // total GC counts
    y3_mean[i] =  y_hat[i, 2] + y_hat[i, 3];
  }

  for (i in 1:numObs){
    CAR_MZcounts_mean[i] = y1_mean[time_index[i]];
    CAR_GCcounts_mean[i] = y2_mean[time_index[i]];
    GC_counts_mean[i] = y3_mean[time_index[i]];
  }
}

model{
  // prior distribution for model parameters
  alpha1 ~ normal(0.01, 0.5);
  lambda ~ normal(0.1, 0.5);
  mu ~ normal(0.1, 0.5);
  delta ~ normal(0.01, 0.5);
  alpha2 ~ normal(0.01, 0.5);
  nu ~ normal(0.01, 0.5);

  sigma1 ~ normal(0, 2.5);
  sigma2 ~ normal(0, 2.5);
  sigma3 ~ normal(0, 2.5);

  // model fitting on to data
  log(CAR_MZ_counts) ~ normal(log(CAR_MZcounts_mean), sigma1);
  log(CAR_GC_counts) ~ normal(log(CAR_GCcounts_mean), sigma2);
  log(GC_counts) ~ normal(log(GC_counts_mean), sigma3);
}

generated quantities{
  // ODE predictions
  real y_hat_pred[numPred, 3];
  // model predictions
  real y1_mean_pred[numPred]; real y2_mean_pred[numPred]; real y3_mean_pred[numPred];
  // model predictions with stdev
  real CAR_MZcounts_pred[numPred]; real CAR_GCcounts_pred[numPred]; real GCcounts_pred[numPred];
  // Residuals
  vector[numObs] resid_d1; vector[numObs] resid_d2; vector[numObs] resid_d3;
  // log likelihoods
  vector[numObs] log_lik1; vector[numObs] log_lik2; vector[numObs] log_lik3;

  //ODE solution
  y_hat_pred[1] = init_cond;
  y_hat_pred[2:] = solve_ODE(ts_pred[2:], init_cond, parms);

  for (i in 1:numPred){
    //MZ
    y1_mean_pred[i] = y_hat_pred[i, 1];
    CAR_MZcounts_pred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma1));

    //GC
    y2_mean_pred[i] = y_hat_pred[i, 2];
    CAR_GCcounts_pred[i] = exp(normal_rng(log(y2_mean_pred[i]), sigma2));

    y3_mean_pred[i] = y_hat_pred[i, 3] + y_hat_pred[i, 2];
    GCcounts_pred[i] = exp(normal_rng(log(y3_mean_pred[i]), sigma3));
  }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    resid_d1[n] = log(CAR_MZ_counts[n]) - log(CAR_MZcounts_mean[n]);
    resid_d2[n] = log(CAR_GC_counts[n]) - log(CAR_GCcounts_mean[n]);
    resid_d3[n] = log(GC_counts[n]) - log(GC_counts_mean[n]);

    log_lik1[n] = normal_lpdf(log(CAR_MZ_counts[n]) | log(CAR_MZcounts_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(log(CAR_GC_counts[n]) | log(CAR_GCcounts_mean[n]), sigma2);
    log_lik3[n] = normal_lpdf(log(GC_counts[n]) | log(GC_counts_mean[n]), sigma3);
  }
}
