functions{
  // function that describes the changes in CAR+ counts in FO B cells
  real CAR_positive_FOB(real time){
    real FO_mean = exp(17.4); real B0 = -5.7507999; real r = 0.4907642; real nu = 0.5807097;
    real value;
    int tau = 7;
    value = FO_mean * exp(B0) * exp(r * time)/(1 +  exp(nu * (time - tau)));
    return value;
   }

   // function that describes the changes in CAR- counts in FO B cells
   real CAR_negative_FOB(real time){
     real FO_mean = exp(17.4); real B0 = -5.7507999; real r = 0.4907642; real nu = 0.5807097;
     real value;
     int tau = 7;
     real phi = exp(B0) * exp(r *time)/(1 +  exp(nu * (time - tau)));
     value = FO_mean * (1 - phi);
     return value;
    }

   // function that contains the ODE equations to be used in ODE solver
   vector DDQode(real time,  vector y, vector parms) {

     // the array of parameters invloved in ode equations
     real alpha1  = parms[1];
     real alpha2  = parms[1];
     real lambda = parms[3];
     real mu     = parms[4];
     real delta  = parms[5];

     real M = exp(14.2);
     real F = exp(17.3);

     // the system of ODEs
     vector[3] dydt;

     dydt[1] = ((alpha1 * y[1])/M) * CAR_positive_FOB(time)  - lambda * y[1];

     dydt[2] = (alpha2/y[3]) * (CAR_positive_FOB(time) - F * y[2]) + mu * (1 - y[2]);
     dydt[3] = alpha2 * y[3] * F - delta * y[3];

     return dydt;
   }

   vector[] solve_ODE(real[] solve_time, vector init_cond, vector parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     array[3] vector[numdim] y_sol;
     y_sol = ode_rk45(DDQode, init_cond, 0.0, solve_time, parms);
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
  real<lower = 0> MZ_counts[numObs];
  real<lower = 0> MZ_fractions[numObs];
  real<lower = 0> GC_counts[numObs];
  real<lower = 0> GC_fractions[numObs];
  real ts_pred[numPred];
  }

parameters{
  // parameters to sample with boundary conditions
  real G0Log;
  real<lower = 0, upper=1> fM_0;
  real<lower = 0, upper=1> fG_0;
  real<lower = 0> alpha1Log;
  real<lower = 0> alpha2Log;
  real<lower = 0> delta;
  real<lower = 0> lambda;
  real<lower = 0> mu;

  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  }


transformed parameters{
  real y_hat[n_shards, 3];  // declaring the array for ODE solution
  real y1_mean[n_shards];      // predictions for dataset1 based on ODE solution
  real y2_mean[n_shards];      // predictions for dataset2 based on ODE solution
  real y3_mean[n_shards];      // predictions for dataset3 based on ODE solution

  real MZ_fractions_mean[numObs];
  real GC_counts_mean[numObs];
  real GC_fractions_mean[numObs];

  vector[5] parms;              // declaring the array for parameters
  vector[3] init_cond;          // declaring the array for state variables

  real G0    = exp(G0Log);      // transformed parameters for better/faster sampling
  real alpha1 = 10^(alpha1Log);      // transformed parameters for better/faster sampling
  real alpha2 = 10^(alpha2Log);      // transformed parameters for better/faster sampling

  // initial conditions and parameters
  init_cond[1] = fM_0;
  init_cond[2] = fG_0;
  init_cond[3] = G0;

  parms[1] = alpha1;
  parms[2] = alpha2;
  parms[3] = lambda;
  parms[4] = mu;
  parms[5] = delta;

  y_hat[1] = init_cond;
  // solution of the system of ODEs for the predictor values
  y_hat[2:] = solve_ODE(solve_time[2:], init_cond, parms);

  for (i in 1:n_shards){
  // CAR fractions in MZ
    y1_mean[i] = y_hat[i, 1];
    // CAR fractions in GC
    y2_mean[i] = y_hat[i, 2];
    // total GC counts
    y3_mean[i] = y_hat[i, 3];
  }

  for (i in 1:numObs){
    MZ_fractions_mean[i] = y1_mean[time_index[i]];
    GC_fractions_mean[i] = y2_mean[time_index[i]];
    GC_counts_mean[i] = y3_mean[time_index[i]];
  }
}

model{
  // prior distribution for model parameters
  alpha1Log ~ normal(-6, 2);
  alpha2Log ~ normal(-6, 2);
  lambda ~ normal(0.05, 0.2);
  mu ~ normal(0.1, 0.2);
  delta ~ normal(0.1, 0.2);

  G0Log ~ normal(11, 2);
  fM_0 ~ normal(0.01, 0.1);
  fG_0 ~ normal(0.3, 0.1);

  sigma1 ~ normal(0, 2);
  sigma2 ~ normal(0, 2);
  sigma3 ~ normal(0, 2);

  // model fitting on to data
  logit(MZ_fractions) ~ normal(logit(MZ_fractions_mean), sigma1);
  log(GC_counts) ~ normal(log(GC_counts_mean), sigma3);
  logit(GC_fractions) ~ normal(logit(GC_fractions_mean), sigma2);
}

generated quantities{
  // ODE predictions
  real y_hat_pred[numPred, 3];
  // model predictions
  real y1_mean_pred[numPred]; real y2_mean_pred[numPred]; real y3_mean_pred[numPred];
  // model predictions with stdev
  real MZfractions_pred[numPred]; real GCcounts_pred[numPred]; real GCfractions_pred[numPred];
  // log likelihoods
  vector[numObs] log_lik1; vector[numObs] log_lik2; vector[numObs] log_lik3;

  //ODE solution
  y_hat_pred[1] = init_cond;
  y_hat_pred[2:] = solve_ODE(ts_pred[2:], init_cond, parms);

  for (i in 1:numPred){
    //MZ
    y1_mean_pred[i] = y_hat_pred[i, 1];
    MZfractions_pred[i] = logit_inverse(normal_rng(logit(y1_mean_pred[i]), sigma1));

    //GC
    y2_mean_pred[i] = y_hat_pred[i, 2];
    GCfractions_pred[i] = logit_inverse(normal_rng(logit(y2_mean_pred[i]), sigma2));

    y3_mean_pred[i] = y_hat_pred[i, 3];
    GCcounts_pred[i] = exp(normal_rng(log(y3_mean_pred[i]), sigma3));

  }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(logit(MZ_fractions[n]) | logit(MZ_fractions_mean[n]), sigma1);
    log_lik3[n] = normal_lpdf(log(GC_counts[n]) | log(GC_counts_mean[n]), sigma3);
    log_lik2[n] = normal_lpdf(logit(GC_fractions[n]) | logit(GC_fractions_mean[n]), sigma2);
  }

}
