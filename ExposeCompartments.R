library(rstan)

expose_stan_functions("Null_neutral_shift.stan")
time_seq <- seq(0.1, 30, length.out = 500)
parms <- numeric(5)           
init_cond <- numeric(6)        

CAR_GC0 <- exp(10.03)         
CAR_MZ0 <- exp(9.8)          
CAR_MZ0N2k0 <- exp(8.03)       

# initial conditions and parameters
init_cond[1] <- 0.0
init_cond[2] <- 0.0
init_cond[3] <- 0.0
init_cond[4] <- CAR_GC0
init_cond[5] <- CAR_MZ0
init_cond[6] <- CAR_MZ0N2k0

parms[1] <- 0.2
parms[2] <- 0.01
parms[3] <- 0.1
parms[4] <- 0.079
parms[5] <- 0.94

solution <- sapply(time_seq, solve_ODE_sys, init_cond, parms)
df <- as.data.frame(solution)
GC <- as.numeric(as.vector(df[4,]))
WTMZ <- as.numeric(as.vector(df[5,]))
N2KOMZ <- as.numeric(as.vector(df[6,]))

alldf <- data.frame(time_seq, GC, WTMZ, N2KOMZ)

ggplot(alldf) +
  geom_line(aes(x = time_seq, y = GC, color = "GC")) +
  geom_line(aes(x = time_seq, y = WTMZ, color = "WTMZ")) +
  geom_line(aes(x = time_seq, y = N2KOMZ, color = "N2KOMZ")) +
  labs(color = "Cell Types") +
  ylab("Cell Counts") +
  xlab("Time") +
  theme_minimal() + scale_y_log10()
# compartmentOutput <- sapply(time_seq, CAR_Positive_FOB_Solve, alpha = 0.1)
spline <- sapply(time_seq, CAR_positive_FOB)
  
FOBpop <- sapply(time_seq, pos_FOB, alpha = 2)
fodf <- as.data.frame(FOBpop)
ggplot(fodf)+
  geom_line(aes(x = time_seq, y=FOBpop, color = "Rep+FOB"))+
  geom_line(aes(x =time_seq, y=spline, color = "spline"))
  ylab("Cell Count")+
  xlab("Time")+
  theme_minimal()
