functions{
  // function that describes the changes in CAR+ counts in FO B cells
  real CAR_positive_FOB(real time){
    real F0 = exp(11.722278); real B0 = exp(4.475064); real n = 4.781548 ; real X = 6.943644 ; real q = 5;
    real value = F0 + (B0 * time^n) * (1 - ((time^q)/((X^q) + (time^q))));
    return value;
   }

   real CAR_negative_MZB(real time){
    real M0 = exp(14.06); real nu = 0.0033; real b0 = 20.58;
    real value = M0 * (1 + exp(-nu * (time - b0)^2));
    return value;
   }

   real Total_FoB(real time){
    real M0 = exp(16.7); real nu = 0.004; real b0 = 20;
    real value = M0 * (1 + exp(-nu * (time - b0)^2));
    return value;
   }

   // function that contains the ODE equations to be used in ODE solver
   real[] ODE_sys(real time,  real[] y, real[] parms, real[] rdata, int[] idata) {
     // the array of parameters invloved in ode equations
     real alpha = parms[1];
     real beta = parms[2];
     real mu = parms[3];
     real delta = parms[4];
     real lambda_WT = parms[5];
     real lambda_N2KO = parms[6];
     real nu = parms[7];

     real t0 = 4.0;
     real alpha_tau = alpha/(1 + exp(nu * (time-t0)^2));
     real beta_tau = beta/(1 + exp(nu * (time-t0)^2));

     // the system of ODEs
     real dydt[3];
     // CAR positive GCB cells in WT
     dydt[1] = alpha_tau * Total_FoB(time)  - delta * y[1];
     // CAR positive MZB cells in WT
     dydt[2] = mu * y[1] + beta_tau * CAR_negative_MZB(time) - lambda_WT * y[2];
     // CAR positive MZB cells in N2KO
     dydt[3] = beta_tau * CAR_negative_MZB(time) - lambda_N2KO * y[3];
     return dydt;
   }

   real[,] solve_ODE_sys(real[] solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     real y_sol[numdim, 3];
     y_sol = integrate_ode_rk45(ODE_sys, init_cond, 4.0, solve_time, parms, {0.0}, {0});
     return y_sol;
   }
}

data{
  int<lower  = 1> numObs1;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> numObs2;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> n_shards;
  int<lower  = 1> numPred;
  int<lower  = 1> time_index1[numObs1];
  int<lower  = 1> time_index2[numObs2];
  real<lower = 0> solve_time[n_shards];
  real<lower = 0> CAR_MZ_counts[numObs1];
  real<lower = 0> CAR_GC_counts[numObs1];
  real<lower = 0> CAR_MZN2_counts[numObs2];
  real<lower = 0> CAR_GCN2_counts[numObs2];
  real ts_pred[numPred];
  }

parameters{
  // parameters to sample with boundary conditions
  real<lower = 0> alpha;
  real<lower = 0> beta;
  real<lower = 0> mu;
  real<lower = 0> delta;
  real<lower = 0> lambda_WT;
  real<lower = 0> lambda_N2KO;
  real<lower = 0> nu;
  real<lower = 0> M0N2;

  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  }


transformed parameters{
  real y_hat[n_shards, 3];     // declaring the array for ODE solution

  real CAR_GCcounts_mean[numObs1];
  real CAR_MZcounts_mean[numObs1];
  real CAR_GCN2counts_mean[numObs2];
  real CAR_MZN2counts_mean[numObs2];

  real parms[7];                  // declaring the array for parameters
  real init_cond[3];              // declaring the array for state variables

  real CAR_GC0 = exp(11.5);              // transformed parameters for better/faster sampling
  real CAR_MZ0 = exp(10.8);              // transformed parameters for better/faster sampling
  real CAR_MZ0N2k0 = exp(M0N2);          // transformed parameters for better/faster sampling

  // initial conditions and parameters
  init_cond[1] = CAR_GC0;
  init_cond[2] = CAR_MZ0;
  init_cond[3] = CAR_MZ0N2k0;

  parms[1] = alpha;
  parms[2] = beta;
  parms[3] = mu;
  parms[4] = delta;
  parms[5] = lambda_WT;
  parms[6] = lambda_N2KO;
  parms[7] = nu;

  y_hat[1] = init_cond;
  // solution of the system of ODEs for the predictor values
  y_hat[2:] = solve_ODE_sys(solve_time[2:], init_cond, parms);

  for (i in 1:numObs1){
    CAR_GCcounts_mean[i] = y_hat[time_index1[i], 1];
    CAR_MZcounts_mean[i] = y_hat[time_index1[i], 2];
  }

  for (i in 1:numObs2){
    CAR_GCN2counts_mean[i] = y_hat[time_index2[i], 1];
    CAR_MZN2counts_mean[i] = y_hat[time_index2[i], 3];
  }
}

model{
  // prior distribution for model parameters
  alpha ~ normal(0.01, 0.5);
  beta ~ normal(0.01, 0.5);
  mu ~ normal(0.01, 0.5);
  nu ~ normal(0.01, 0.5);
  delta ~ normal(0.01, 0.5);
  lambda_WT ~ normal(0.01, 0.5);
  lambda_N2KO ~ normal(0.01, 0.5);
  M0N2 ~ normal(8, 1);

  sigma1 ~ normal(0, 2.5);
  sigma2 ~ normal(0, 2.5);
  sigma3 ~ normal(0, 2.5);

  // model fitting on to data
  log(CAR_GC_counts) ~ normal(log(CAR_GCcounts_mean), sigma1);
  log(CAR_MZ_counts) ~ normal(log(CAR_MZcounts_mean), sigma2);
  log(CAR_GCN2_counts) ~ normal(log(CAR_GCN2counts_mean), sigma1);
  log(CAR_MZN2_counts) ~ normal(log(CAR_MZN2counts_mean), sigma3);
}

generated quantities{
   // ODE predictions
   real y_hat_pred[numPred, 3];
   // variables for model predictions
   real y1_mean_pred[numPred]; real y2_mean_pred[numPred]; real y3_mean_pred[numPred]; real y4_mean_pred[numPred];
   real FOtoCARMZ_pred[numPred]; real MZtoCARMZ_pred[numPred]; real FOtoCARGC_pred[numPred];
   real alpha_pred[numPred];  real beta_pred[numPred];
   // variables for model predictions with stdev
   real CAR_MZcounts_pred[numPred]; real CAR_GCcounts_pred[numPred];
   real CAR_MZN2counts_pred[numPred]; real CAR_GCN2counts_pred[numPred];
   // Residuals
   vector[numObs1] resid_d1; vector[numObs1] resid_d2; vector[numObs2] resid_d3; vector[numObs2] resid_d4;
   // log likelihoods
   vector[numObs1] log_lik1; vector[numObs1] log_lik2; vector[numObs2] log_lik3; vector[numObs2] log_lik4;

   //ODE solution
   y_hat_pred[1] = init_cond;
   y_hat_pred[2:] = solve_ODE_sys(ts_pred[2:], init_cond, parms);

   // model predictions with stdev
   for (i in 1:numPred){
     //CAR GC
     y1_mean_pred[i] = y_hat_pred[i, 1];
     CAR_GCcounts_pred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma1));
     //CAR MZ
     y2_mean_pred[i] = y_hat_pred[i, 2];
     CAR_MZcounts_pred[i] = exp(normal_rng(log(y2_mean_pred[i]), sigma2));

     //CAR GC in N2KO
     y3_mean_pred[i] = y_hat_pred[i, 1];
     CAR_GCN2counts_pred[i] = exp(normal_rng(log(y3_mean_pred[i]), sigma1));
     //CAR MZ in N2KO
     y4_mean_pred[i] = y_hat_pred[i, 3];
     CAR_MZN2counts_pred[i] = exp(normal_rng(log(y4_mean_pred[i]), sigma3));

     // Influx into CAR MZ
     beta_pred[i] = mu/(1 + exp(nu *(ts_pred[i] - 4.0)^2));
     FOtoCARMZ_pred[i] = mu  * y1_mean_pred[i]/y2_mean_pred[i];
     MZtoCARMZ_pred[i] = beta_pred[i] * CAR_negative_MZB(ts_pred[i])/y2_mean_pred[i];
     // Influx into CAR GC
     alpha_pred[i] = alpha/(1 + exp(nu *(ts_pred[i] - 4.0)^2));
     FOtoCARGC_pred[i] = ((alpha/(1 + exp(nu * (ts_pred[i] - 4.0)^2))) * Total_FoB(ts_pred[i]))/y1_mean_pred[i];
   }

   // calculating the log predictive accuracy for each point
   for (n in 1:numObs1) {
     resid_d1[n] = log(CAR_GC_counts[n]) - log(CAR_GCcounts_mean[n]);
     resid_d2[n] = log(CAR_MZ_counts[n]) - log(CAR_MZcounts_mean[n]);

     log_lik1[n] = normal_lpdf(log(CAR_GC_counts[n]) | log(CAR_GCcounts_mean[n]), sigma1);
     log_lik2[n] = normal_lpdf(log(CAR_MZ_counts[n]) | log(CAR_MZcounts_mean[n]), sigma2);
   }

   // calculating the log predictive accuracy for each point
   for (n in 1:numObs2) {
     resid_d3[n] = log(CAR_GCN2_counts[n]) - log(CAR_GCN2counts_mean[n]);
     resid_d4[n] = log(CAR_MZN2_counts[n]) - log(CAR_MZN2counts_mean[n]);

     log_lik3[n] = normal_lpdf(log(CAR_GCN2_counts[n]) | log(CAR_GCN2counts_mean[n]), sigma1);
     log_lik4[n] = normal_lpdf(log(CAR_MZN2_counts[n]) | log(CAR_MZN2counts_mean[n]), sigma3);
   }
}
