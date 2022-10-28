functions{
  // function that describes the changes in source counts
  real theta_spline(real time){
    int t0 = 40;         // the earliest age at BMT
    real nu = 1.86e-03; real theta0 = 14.36;
      real theta;
       // the earliest age at BMT
      // nu is the rate of decline of source influx assuming that its proportional to the counts of source compartment
      //theta0 gives the initial counts of the source compartment
      theta = exp(theta0) * exp(-nu * (time-t0));
      return theta;
   }

    // function that contains the ODE equations to be used in ODE solver
    real[] shm(real time,  real[] k, real[] parms, real[] rdata, int[] idata) {
     // the array of parameters invloved in ode equations
     real psi    = parms[1];
     real delta = parms[2];

     // the system of ODEs
     real dkdt[1];

     //donor ki67hi
     dkdt[1] = psi * theta_spline(time) - delta * k[1];

     return dkdt;
   }

   real[] solve_shm(real solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     real value[1];

     if (solve_time <= 1.0){
       value = init_cond;
     } else {
       value = to_array_1d(integrate_ode_rk45(shm, init_cond, 1.0, rep_array(solve_time, 1), parms, {0.0}, {0}));
     }
     return value;
    }
}

data{
  int<lower  = 1> numObs;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> numPred;
  real<lower = 0> solve_time[n_shards];
  real<lower = 0> counts[numObs];
  }

transformed data{
  real y1[numObs];
  y1 = log(counts);                          // transforming cell counts of donor compartments to feed in to ODEs
  }

parameters{
  // parameters to sample with boundary conditions
  real y0_Log;
  real<lower = 0, upper=1> psi;
  real<lower = 0> delta;
  // stdev within individual datasets to be estimated
  real<lower = 0> sigma1;
  }

transformed parameters{
  real k_hat[n_shards];  // declaring the array for ODE solution
  real y1_mean[numObs];      // predictions for dataset1 based on ODE solution
  real y2_mean[numObs];      // predictions for dataset2 based on ODE solution
  real y3_mean[numObs];      // predictions for dataset3 based on ODE solution
  real y4_mean[numObs];      // predictions for dataset4 based on ODE solution

  real chi_counts_mean[n_shards];
  real donor_fractions_mean[n_shards];
  real host_ki_mean[n_shards];
  real donor_ki_mean[n_shards];

  real parms[4];              // declaring the array for parameters
  real init_cond[4];          // declaring the array for state variables

  real y0    = exp(y0_Log);      // transformed parameters for better/faster sampling

  // initial conditions and parameters
  // kapp0 = ki67 proportions within host at t0 estimated from the model fit
  // at t0 donor counts i.e Nd_0 = 0 --> Nh_0 = N0
  // also kappa_host_0 = kapp0
  init_cond[1] = 0.0;
  init_cond[2] = 0.0;
  init_cond[3] = y0 * kappa_0;
  init_cond[4] = y0 * (1 - kappa_0);

  parms[1] = psi;
  parms[2] = rho;
  parms[3] = delta;
  parms[4] = Beta;

  // solution of the system of ODEs for the predictor values
  k_hat = solve_ode_chi(solve_time, ageAtBMT, init_cond, parms);

  for (i in 1:n_shards){
    // total counts
    chi_counts_mean[i] = k_hat[4*i-3] + k_hat[4*i-2] + k_hat[4*i-1] + k_hat[4*i];

    // donor fractions normalised with chimerism in the source
    donor_fractions_mean[i] = (k_hat[4*i-3] + k_hat[4*i-2])/(chi_counts_mean[i] * Chi_T1(solve_time[i] - ageAtBMT[i]));

    // fractions of ki67 positive cells in the donor compartment
    donor_ki_mean[i] = k_hat[4*i-3]/(k_hat[4*i-3] + k_hat[4*i-2]);

    // fractions of ki67 positive cells in the host compartment
    host_ki_mean[i] = k_hat[4*i-1]/(k_hat[4*i-1] + k_hat[4*i]);
  }

  for (i in 1:numObs){
    // total counts
    y1_mean[i] = chi_counts_mean[time_index[i]];

    // donor fractions normalised with chimerism in the source
    y2_mean[i] = donor_fractions_mean[time_index[i]];

    // fractions of ki67 positive cells in the donor compartment
    y3_mean[i] = donor_ki_mean[time_index[i]];

    // fractions of ki67 positive cells in the host compartment
    y4_mean[i] = host_ki_mean[time_index[i]];
  }
}

model{
  // prior distribution for model parameters
  psi ~ normal(0.5, 0.25);
  y0_Log ~ normal(14, 1);
  kappa_0 ~  normal(0.1, 0.15);
  rho ~ normal(0.01, 0.2);
  delta ~ normal(0.04, 0.25);
  Beta ~ normal(3.5, 0.8);

  sigma1 ~ normal(0.6, 1);
  sigma2 ~ normal(3, 1);
  sigma3 ~ normal(0.3, 1);
  sigma4 ~ normal(0.9, 2);

  // model fitting on to data
  y1 ~ normal(log(y1_mean), sigma1);
  y2 ~ normal(logit_boundary_array(y2_mean), sigma2);
  y3 ~ normal(logit_boundary_array(y3_mean), sigma3);
  y4 ~ normal(logit_boundary_array(y4_mean), sigma4);
}

generated quantities{
  // ODE predictions for 3 different timecourses based on the age at BMT bins
  real k_hat_pred_age1[numPred*4];
  real k_hat_pred_age2[numPred*4];
  real k_hat_pred_age3[numPred*4];

  // model predictions for the datasets for 3 different timecourses based on the age at BMT bins
  real y1_mean_pred_age1[numPred]; real y2_mean_pred_age1[numPred]; real y3_mean_pred1[numPred]; real y4_mean_pred1[numPred];
  real y1_mean_pred_age2[numPred]; real y2_mean_pred_age2[numPred]; real y3_mean_pred2[numPred]; real y4_mean_pred2[numPred];
  real y1_mean_pred_age3[numPred]; real y2_mean_pred_age3[numPred]; real y3_mean_pred3[numPred]; real y4_mean_pred3[numPred];

  // model predictions for the datasets with stdev estimated from the fits for 3 different timecourses based on the age at BMT bins
  real countspred_age1[numPred]; real fdpred_age1[numPred]; real donor_kiprop_pred1[numPred]; real host_kiprop_pred1[numPred];
  real countspred_age2[numPred]; real fdpred_age2[numPred]; real donor_kiprop_pred2[numPred]; real host_kiprop_pred2[numPred];
  real countspred_age3[numPred]; real fdpred_age3[numPred]; real donor_kiprop_pred3[numPred]; real host_kiprop_pred3[numPred];

  // log likelihood for individual data sets and combined
  vector[numObs] log_lik1; vector[numObs] log_lik2; vector[numObs] log_lik3; vector[numObs] log_lik4;

  // parameters of interest
  real rho_inv = 1/rho;
  real delta_inv = 1/delta;

  real eps_donor = 0.99;

  //ODE solution for different age bins
  k_hat_pred_age1 = solve_ode_chi(ts_pred1, tb_pred1, init_cond, parms);
  k_hat_pred_age2 = solve_ode_chi(ts_pred2, tb_pred2, init_cond, parms);
  k_hat_pred_age3 = solve_ode_chi(ts_pred3, tb_pred3, init_cond, parms);

  // Total cell counts (donor + host) for different age bins

  for (i in 1:numPred){
    // age bin1
    y1_mean_pred_age1[i] = k_hat_pred_age1[4*i-3] + k_hat_pred_age1[4*i-2] + k_hat_pred_age1[4*i-1] + k_hat_pred_age1[4*i];
    countspred_age1[i] = exp(normal_rng(log(y1_mean_pred_age1[i]), sigma1));
    // age bin2
    y1_mean_pred_age2[i] = k_hat_pred_age2[4*i-3] + k_hat_pred_age2[4*i-2] + k_hat_pred_age2[4*i-1] + k_hat_pred_age2[4*i];
    countspred_age2[i] = exp(normal_rng(log(y1_mean_pred_age2[i]), sigma1));
    // age bin3
    y1_mean_pred_age3[i] = k_hat_pred_age3[4*i-3] + k_hat_pred_age3[4*i-2] + k_hat_pred_age3[4*i-1] + k_hat_pred_age3[4*i];
    countspred_age3[i] = exp(normal_rng(log(y1_mean_pred_age3[i]), sigma1));

    y2_mean_pred_age1[i] = (k_hat_pred_age1[4*i-3] + k_hat_pred_age1[4*i-2])/(y1_mean_pred_age1[i] * Chi_T1(ts_pred1[i] - tb_pred1[i]));
    fdpred_age1[i] = expit_boundary(normal_rng(logit_boundary(y2_mean_pred_age1[i]), sigma2));

    y2_mean_pred_age2[i] = (k_hat_pred_age2[4*i-3] + k_hat_pred_age2[4*i-2])/(y1_mean_pred_age2[i] * Chi_T1(ts_pred2[i] - tb_pred2[i]));
    fdpred_age2[i] = expit_boundary(normal_rng(logit_boundary(y2_mean_pred_age2[i]), sigma2));

    y2_mean_pred_age3[i] = (k_hat_pred_age3[4*i-3] + k_hat_pred_age3[4*i-2])/(y1_mean_pred_age3[i] * Chi_T1(ts_pred3[i] - tb_pred3[i]));
    fdpred_age3[i] = expit_boundary(normal_rng(logit_boundary(y2_mean_pred_age3[i]), sigma2));
  }


  // at t0 the ki67 proportions in donor of target = ki67 proportions in source donor compartemnt
  y3_mean_pred1[1] = eps_donor;
  donor_kiprop_pred1[1] = expit_boundary(normal_rng(logit_boundary(y3_mean_pred1[1]), sigma3));
  y3_mean_pred2[1] = eps_donor;
  donor_kiprop_pred2[1] = expit_boundary(normal_rng(logit_boundary(y3_mean_pred2[1]), sigma3));
  y3_mean_pred3[1] = eps_donor;
  donor_kiprop_pred3[1] = expit_boundary(normal_rng(logit_boundary(y3_mean_pred3[1]), sigma3));

  // ki67 proportions at t0 estimated from the model fit
  y4_mean_pred1[1] = kappa_0;
  host_kiprop_pred1[1] = expit_boundary(normal_rng(logit_boundary(y4_mean_pred1[1]), sigma4));
  y4_mean_pred2[1] = kappa_0;
  host_kiprop_pred2[1] = expit_boundary(normal_rng(logit_boundary(y4_mean_pred2[1]), sigma4));
  y4_mean_pred3[1] = kappa_0;
  host_kiprop_pred3[1] = expit_boundary(normal_rng(logit_boundary(y4_mean_pred3[1]), sigma4));

  for (i in 2:numPred){
    y3_mean_pred1[i] = k_hat_pred_age1[4*i-3]/(k_hat_pred_age1[4*i-3] + k_hat_pred_age1[4*i-2]);
    donor_kiprop_pred1[i] = expit_boundary(normal_rng(logit_boundary(y3_mean_pred1[i]), sigma3));

    y4_mean_pred1[i] = k_hat_pred_age1[4*i-1]/(k_hat_pred_age1[4*i-1] + k_hat_pred_age1[4*i]);
    host_kiprop_pred1[i] = expit_boundary(normal_rng(logit_boundary(y4_mean_pred1[i]), sigma4));

    y3_mean_pred2[i] = k_hat_pred_age2[4*i-3]/(k_hat_pred_age2[4*i-3] + k_hat_pred_age2[4*i-2]);
    donor_kiprop_pred2[i] = expit_boundary(normal_rng(logit_boundary(y3_mean_pred2[i]), sigma3));

    y4_mean_pred2[i] = k_hat_pred_age2[4*i-1]/(k_hat_pred_age2[4*i-1] + k_hat_pred_age2[4*i]);
    host_kiprop_pred2[i] = expit_boundary(normal_rng(logit_boundary(y4_mean_pred2[i]), sigma4));

    y3_mean_pred3[i] = k_hat_pred_age3[4*i-3]/(k_hat_pred_age3[4*i-3] + k_hat_pred_age3[4*i-2]);
    donor_kiprop_pred3[i] = expit_boundary(normal_rng(logit_boundary(y3_mean_pred3[i]), sigma3));

    y4_mean_pred3[i] = k_hat_pred_age3[4*i-1]/(k_hat_pred_age3[4*i-1] + k_hat_pred_age3[4*i]);
    host_kiprop_pred3[i] = expit_boundary(normal_rng(logit_boundary(y4_mean_pred3[i]), sigma4));
  }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(y1[n] | log(y1_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(y2[n] | logit_boundary(y2_mean[n]), sigma2);
    log_lik3[n] = normal_lpdf(y3[n] | logit_boundary(y3_mean[n]), sigma3);
    log_lik4[n] = normal_lpdf(y4[n] | logit_boundary(y4_mean[n]), sigma4);
  }
}
