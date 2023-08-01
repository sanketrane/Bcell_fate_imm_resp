functions{
  // function that describes the changes in CAR+ counts in FO B cells
  real CAR_positive_FOB(real time){
    real F0 = exp(11.722278); real B0 = exp(4.475064); real n = 4.781548 ; real X = 6.943644 ; real q = 5;
    real value = F0 + (B0 * time^n) * (1 - ((time^q)/((X^q) + (time^q))));
    return value;
   }

  real[] CAR_Positive_FOB_Compartments_ODE(real time, real[] y, real[] parms, real[] rdata, int[] idata) {
    real dydt[3];
    real alpha = parms[1];

    dydt[1] = CAR_positive_FOB(time) * alpha - alpha * y[1];
    dydt[2] = y[1] * alpha - alpha * y[2];
    dydt[3] = y[2] * alpha - alpha * y[3];

    return dydt;
  }
 
  real CAR_Positive_FOB_Solve(real solve_time, real alpha) {
    //int numdim = size(solve_time);
    real com_sol[1,3];
    real init_cond[3];

    real parms[1];
    parms[1] = alpha;

    init_cond[1] = 0.0;
    init_cond[2] = 0.0;
    init_cond[3] = 0.0;

    com_sol = integrate_ode_rk45(CAR_Positive_FOB_Compartments_ODE, init_cond, 0.0, rep_array(solve_time, 1), parms, {0.0}, {0});

    return com_sol[1,3];
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

   real pos_FOB(real time, real alpha){
    real F0 = exp(11.722278);
    if (time == 0){
      return F0;
    } else{
      return CAR_Positive_FOB_Solve(time, alpha);
    }
   }

   // function that contains the ODE equations to be used in ODE solver
   real[] ODE_sys(real time,  real[] y, real[] parms, real[] rdata, int[] idata) {
     // the array of parameters invloved in ode equations
     real alpha = parms[1];
     real beta = parms[2];
     real delta = parms[3];
     real lambda_WT = parms[4];
     real lambda_N2KO = parms[5];


     // the system of ODEs
     real dydt[3];
     // CAR positive GCB cells in WT
     dydt[1] = pos_FOB(time, alpha) - delta * y[1];
     // CAR positive MZB cells in WT
     dydt[2] = beta * CAR_negative_MZB(time) - lambda_WT * y[2];
     // CAR positive MZB cells in N2KO
     dydt[3] = beta * CAR_negative_MZB(time) - lambda_N2KO * y[3];
     return dydt;
   }

   real[,] solve_ODE_sys(real[] solve_time, real[] init_cond, real[] parms) {
     // solves the ode for each timepoint from t0
     int numdim = size(solve_time);
     real y_sol[numdim, 3];
     y_sol = integrate_ode_rk45(ODE_sys, init_cond, 0.0, solve_time, parms, {0.0}, {0});
     return y_sol;
   }
}
