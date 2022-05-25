functions{
   // function that describes the changes in CAR+ counts in FO B cells
   real CAR_positive_FOB(real time){
     real F0 = exp(11.763019); real B0 = exp(4.717021); real n = 5.092933 ; real X = 7.121932 ; real q = 6.475884;
     real value = F0 + (B0 * time^n) * (1 - ((time^q)/((X^q) + (time^q))));
     return value;
    }

    real CAR_negative_MZB(real time){
     real M0 = exp(14); real nu = 0.01; real b0 = 18;
     real value = M0 * (1 + exp(-nu * (time - b0)^2));
     return value;
    }


   // function that contains the ODE equations to be used in ODE solver
   real[] ODE_sys1(real time,  real[] y, real[] parms, real[] rdata, int[] idata) {
     // the array of parameters invloved in ode equations
     real alpha1 = parms[1];
     real lambda_WT = parms[2];
     real mu     = parms[3];
     real delta  = parms[4];
     real alpha2 = parms[5];
     real nu     = parms[6];
     real alpha3 = parms[7];
     real CAR_negative_FOB = exp(17.15);
     real t0 = 4.0;

     real alpha_GC = alpha2/(1 + exp(nu * (time-t0)^2));
     // the system of ODEs
     real dydt[3];
     dydt[1] = alpha1 * y[2] + alpha3 * CAR_negative_MZB(time) - lambda_WT * y[1];
     dydt[2] = alpha_GC * CAR_positive_FOB(time) + mu * y[3]  - delta * y[2];
     dydt[3] = alpha_GC * CAR_negative_FOB - mu * y[3] - delta * y[3];
     return dydt;
   }

   // function that contains the ODE equations to be used in ODE solver
   real[] ODE_sys2(real time,  real[] y, real[] parms, real[] rdata, int[] idata) {
     // the array of parameters invloved in ode equations
     real alpha1 = parms[1];
     real lambda_WT = parms[2];
     real mu     = parms[3];
     real delta  = parms[4];
     real alpha2 = parms[5];
     real nu     = parms[6];
     real alpha3 = parms[7];
     real lambda_N2KO = parms[8];
     real CAR_negative_FOB = exp(17.15);
     real t0 = 4.0;

     real alpha_GC = alpha2/(1 + exp(nu * (time-t0)^2));
     // the system of ODEs
     real dydt[3];
     dydt[1] = alpha3 * CAR_negative_MZB(time) - lambda_N2KO * y[1];
     dydt[2] = alpha_GC * CAR_positive_FOB(time) + mu * y[3]  - delta * y[2];
     dydt[3] = alpha_GC * CAR_negative_FOB - mu * y[3] - delta * y[3];
     return dydt;
   }

   real[,] solve_ODE_sys1(real[] solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     real y_sol[numdim, 3];
     y_sol = integrate_ode_rk45(ODE_sys1, init_cond, 4.0, solve_time, parms, {0.0}, {0});
     return y_sol;
   }

   real[,] solve_ODE_sys2(real[] solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     real y_sol[numdim, 3];
     y_sol = integrate_ode_rk45(ODE_sys2, init_cond, 4.0, solve_time, parms, {0.0}, {0});
     return y_sol;
   }
}

data{
  int<lower  = 1> numObs1;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> numObs2;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> n_shards1;
  int<lower  = 1> n_shards2;
  int<lower  = 1> numPred;
  int<lower  = 1> time_index1[numObs1];
  int<lower  = 1> time_index2[numObs2];
  real<lower = 0> solve_time1[n_shards1];
  real<lower = 0> solve_time2[n_shards2];
  real<lower = 0> CAR_MZ_counts[numObs1];
  real<lower = 0> GC_counts[numObs1];
  real<lower = 0> CAR_GC_counts[numObs1];
  real<lower = 0> CAR_MZN2_counts[numObs2];
  real<lower = 0> GCN2_counts[numObs2];
  real<lower = 0> CAR_GCN2_counts[numObs2];
  real ts_pred[numPred];
  }

parameters{
  // parameters to sample with boundary conditions
  real<lower = 0> alpha1;
  real<lower = 0> alpha2;
  real<lower = 0> alpha3;
  real<lower = 0> delta;
  real<lower = 0> lambda_WT;
  real<lower = lambda_WT> lambda_N2KO;
  real<lower = 0> mu;
  real<lower = 0> nu;
  real<lower = 0> M0N2;

  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  real<lower = 0> sigma4;
  }


transformed parameters{
  real y_hat1[n_shards1, 3];     // declaring the array for ODE solution
  real y_hat2[n_shards2, 3];     // declaring the array for ODE solution
  real y1_mean[n_shards1];      // predictions for dataset1 based on ODE solution
  real y2_mean[n_shards1];      // predictions for dataset2 based on ODE solution
  real y3_mean[n_shards1];      // predictions for dataset3 based on ODE solution
  real y4_mean[n_shards2];      // predictions for dataset3 based on ODE solution
  real y5_mean[n_shards2];      // predictions for dataset1 based on ODE solution
  real y6_mean[n_shards2];      // predictions for dataset2 based on ODE solution

  real CAR_MZcounts_mean[numObs1];
  real CAR_GCcounts_mean[numObs1];
  real GC_counts_mean[numObs1];
  real CAR_MZN2counts_mean[numObs2];
  real CAR_GCN2counts_mean[numObs2];
  real GCN2_counts_mean[numObs2];

  real parms[8];                  // declaring the array for parameters
  real init_cond1[3];              // declaring the array for state variables
  real init_cond2[3];              // declaring the array for state variables

  real G0 = exp(11.7);            // transformed parameters for better/faster sampling
  real CAR_MZ0 = exp(10.6);            // transformed parameters for better/faster sampling
  real CAR_GC0 = exp(10.8);              // transformed parameters for better/faster sampling
  real CAR_MZ0N2k0 = exp(M0N2);              // transformed parameters for better/faster sampling


  // initial conditions and parameters
  init_cond1[1] = CAR_MZ0;
  init_cond1[2] = CAR_GC0;
  init_cond1[3] = G0 - CAR_GC0;

  init_cond2[1] = CAR_MZ0N2k0;
  init_cond2[2] = CAR_GC0;
  init_cond2[3] = G0 - CAR_GC0;

  parms[1] = alpha1;
  parms[2] = lambda_WT;
  parms[3] = mu;
  parms[4] = delta;
  parms[5] = alpha2;
  parms[6] = nu;
  parms[7] = alpha3;
  parms[8] = lambda_N2KO;

  y_hat1[1] = init_cond1;
  // solution of the system of ODEs for the predictor values
  y_hat1[2:] = solve_ODE_sys1(solve_time1[2:], init_cond1, parms);

  for (i in 1:n_shards1){
  // CAR counts in MZ
    y1_mean[i] = y_hat1[i, 1];
    // CAR counts in GC
    y2_mean[i] = y_hat1[i, 2];
    // total GC counts
    y3_mean[i] =  y_hat1[i, 2] + y_hat1[i, 3];
  }

  y_hat2[1] = init_cond2;
  // solution of the system of ODEs for the predictor values
  y_hat2[2:] = solve_ODE_sys2(solve_time2[2:], init_cond2, parms);

  for (i in 1:n_shards2){
  // CAR counts in MZ
    y4_mean[i] = y_hat2[i, 1];
    // CAR counts in GC
    y5_mean[i] = y_hat2[i, 2];
    // total GC counts
    y6_mean[i] =  y_hat2[i, 2] + y_hat2[i, 3];
  }


  for (i in 1:numObs1){
    CAR_MZcounts_mean[i] = y1_mean[time_index1[i]];
    CAR_GCcounts_mean[i] = y2_mean[time_index1[i]];
    GC_counts_mean[i] = y3_mean[time_index1[i]];
  }

  for (i in 1:numObs2){
    CAR_MZN2counts_mean[i] = y4_mean[time_index2[i]];
    CAR_GCN2counts_mean[i] = y2_mean[time_index2[i]];
    GCN2_counts_mean[i] = y3_mean[time_index2[i]];
  }
}

model{
  // prior distribution for model parameters
  alpha1 ~ normal(0.01, 0.5);
  lambda_WT ~ normal(0.1, 0.5);
  lambda_N2KO ~ normal(0.1, 0.5);
  mu ~ normal(0.1, 0.5);
  delta ~ normal(0.01, 0.5);
  alpha2 ~ normal(0.01, 0.5);
  nu ~ normal(0.01, 0.5);
  alpha3 ~ normal(0.01, 0.5);
  M0N2 ~ normal(8, 2);

  sigma1 ~ normal(0, 2.5);
  sigma2 ~ normal(0, 2.5);
  sigma3 ~ normal(0, 2.5);
  sigma4 ~ normal(0, 2.5);

  // model fitting on to data
  log(CAR_MZ_counts) ~ normal(log(CAR_MZcounts_mean), sigma1);
  log(CAR_GC_counts) ~ normal(log(CAR_GCcounts_mean), sigma2);
  log(GC_counts) ~ normal(log(GC_counts_mean), sigma3);
  log(CAR_MZN2_counts) ~ normal(log(CAR_MZN2counts_mean), sigma4);
  log(CAR_GCN2_counts) ~ normal(log(CAR_GCN2counts_mean), sigma2);
  log(GCN2_counts) ~ normal(log(GCN2_counts_mean), sigma3);
}

generated quantities{
  // ODE predictions
  real y_hat_pred1[numPred, 3];
  real y_hat_pred2[numPred, 3];
  // model predictions
  real y1_mean_pred[numPred]; real y2_mean_pred[numPred]; real y3_mean_pred[numPred];
  real y4_mean_pred[numPred]; real y5_mean_pred[numPred]; real y6_mean_pred[numPred];
  real alpha_pred[numPred];
  // model predictions with stdev
  real CAR_MZcounts_pred[numPred]; real CAR_GCcounts_pred[numPred]; real GCcounts_pred[numPred];
  real CAR_MZN2counts_pred[numPred]; real CAR_GCN2counts_pred[numPred]; real GCN2counts_pred[numPred];
  // Residuals
  vector[numObs1] resid_d1; vector[numObs1] resid_d2; vector[numObs1] resid_d3;
  vector[numObs2] resid_d4; vector[numObs2] resid_d5; vector[numObs2] resid_d6;
  // log likelihoods
  vector[numObs1] log_lik1; vector[numObs1] log_lik2; vector[numObs1] log_lik3;
  vector[numObs2] log_lik4;  vector[numObs2] log_lik5;  vector[numObs2] log_lik6;

  //ODE solution
  y_hat_pred1[1] = init_cond1;
  y_hat_pred1[2:] = solve_ODE_sys1(ts_pred[2:], init_cond1, parms);

  //ODE solution
  y_hat_pred2[1] = init_cond2;
  y_hat_pred2[2:] = solve_ODE_sys2(ts_pred[2:], init_cond2, parms);

  for (i in 1:numPred){
    //CAR MZ
    y1_mean_pred[i] = y_hat_pred1[i, 1];
    CAR_MZcounts_pred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma1));

    //CAR GC
    y2_mean_pred[i] = y_hat_pred1[i, 2];
    CAR_GCcounts_pred[i] = exp(normal_rng(log(y2_mean_pred[i]), sigma2));

    //total GC
    y3_mean_pred[i] = y_hat_pred1[i, 3] + y_hat_pred1[i, 2];
    GCcounts_pred[i] = exp(normal_rng(log(y3_mean_pred[i]), sigma3));

    //CAR MZ in N2KO
    y4_mean_pred[i] = y_hat_pred2[i, 1];
    CAR_MZN2counts_pred[i] = exp(normal_rng(log(y4_mean_pred[i]), sigma4));

    //CAR GC in N2KO
    y5_mean_pred[i] = y_hat_pred2[i, 2];
    CAR_GCN2counts_pred[i] = exp(normal_rng(log(y5_mean_pred[i]), sigma2));

    //total GC in N2KO
    y6_mean_pred[i] = y_hat_pred2[i, 3] + y_hat_pred2[i, 2];
    GCN2counts_pred[i] = exp(normal_rng(log(y6_mean_pred[i]), sigma3));

    //influx rate
    alpha_pred[i] = alpha2/(1 + exp(nu * (ts_pred[i] - 4.0)^2));
  }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs1) {
    resid_d1[n] = log(CAR_MZ_counts[n]) - log(CAR_MZcounts_mean[n]);
    resid_d2[n] = log(CAR_GC_counts[n]) - log(CAR_GCcounts_mean[n]);
    resid_d3[n] = log(GC_counts[n]) - log(GC_counts_mean[n]);

    log_lik1[n] = normal_lpdf(log(CAR_MZ_counts[n]) | log(CAR_MZcounts_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(log(CAR_GC_counts[n]) | log(CAR_GCcounts_mean[n]), sigma2);
    log_lik3[n] = normal_lpdf(log(GC_counts[n]) | log(GC_counts_mean[n]), sigma3);
  }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs2) {
    resid_d4[n] = log(CAR_MZN2_counts[n]) - log(CAR_MZN2counts_mean[n]);
    resid_d5[n] = log(CAR_GCN2_counts[n]) - log(CAR_GCN2counts_mean[n]);
    resid_d6[n] = log(GCN2_counts[n]) - log(GCN2_counts_mean[n]);

    log_lik4[n] = normal_lpdf(log(CAR_MZN2_counts[n]) | log(CAR_MZN2counts_mean[n]), sigma4);
    log_lik5[n] = normal_lpdf(log(CAR_GCN2_counts[n]) | log(CAR_GCN2counts_mean[n]), sigma2);
    log_lik6[n] = normal_lpdf(log(GCN2_counts[n]) | log(GCN2_counts_mean[n]), sigma3);
  }
}
