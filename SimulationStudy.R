#-----------------Simulation Study-----------
#Compare BART to other methods for causal inference
#Three levels of model misspecification for each of the 3 Simulations:
#No misspecification, Moderate misspecification, Strong Misspecification

#Simulation 1: Linear
#Simulation 2: Non-Linear
#Simulation 3: Non-Linear, Multiple Confounders


##------------------Simulation 1A------------
#Linear, Moderate Misspecifcation, No confounding

n=100
Sim1A_Data = data.frame()
Sim1A_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N = 1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + 0.3*X1 + 1.5*X2,
             sd = 1)
  
  Y <- T + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  Y1 <- grid[2] + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  Y2 <- grid[3] + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  Y3 <- grid[4] + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  Y4 <- grid[5] + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  truth <- cbind(Y0, Y1, Y2, Y3, Y4)
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, eY, T, Y))
  
  
  ###BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  ####GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  ####IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  ###Store results in row i of matrix
  Sim1A_Results[i,1] <- grid[1]
  Sim1A_Results[i,2] <- grid[2]
  Sim1A_Results[i,3] <- grid[3]
  Sim1A_Results[i,4] <- grid[4]
  Sim1A_Results[i,5] <- grid[5]
  
  Sim1A_Results[i,6] <- mean(Y0)
  Sim1A_Results[i,7] <- mean(Y1)
  Sim1A_Results[i,8] <- mean(Y2)
  Sim1A_Results[i,9] <- mean(Y3)
  Sim1A_Results[i,10] <- mean(Y4)
  
  Sim1A_Results[i,11] <- bart[1]
  Sim1A_Results[i,12] <- bart[2]
  Sim1A_Results[i,13] <- bart[3]
  Sim1A_Results[i,14] <- bart[4]
  Sim1A_Results[i,15] <- bart[5]
  
  Sim1A_Results[i,16] <- gam[1]
  Sim1A_Results[i,17] <- gam[2]
  Sim1A_Results[i,18] <- gam[3]
  Sim1A_Results[i,19] <- gam[4]
  Sim1A_Results[i,20] <- gam[5]
  
  Sim1A_Results[i,21] <- iptw[1]
  Sim1A_Results[i,22] <- iptw[2]
  Sim1A_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim1A_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim1A_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim1A_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #S#ave Data for Use in TMLE and GPS
  Sim1A_Data <- rbind(Sim1A_Data,data)
}
write.csv(Sim1A_Results,"Sim1A_Results.csv")
write.csv(Sim1A_Data,"Sim1A_Data.csv")




#------------------Simulation 1B------------
#Linear, Moderate Misspecifcation, with Confounding
n=100
Sim1B_Data = data.frame()
Sim1B_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N = 1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  C1 <- rnorm(N, 0, 1)
  
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + 0.3*X1 + 0.3*Z1 + 0.5*Z2 + 1.5*X2 + 0.5*C1,
             sd = 1)
  
  Y <- T + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.5*C1 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.5*C1 + eY
  Y1 <- grid[2] + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.5*C1 + eY
  Y2 <- grid[3] + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.5*C1 + eY
  Y3 <- grid[4] + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.5*C1 + eY
  Y4 <- grid[5] + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.5*C1 + eY
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, C1, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim1B_Results[i,1] <- grid[1]
  Sim1B_Results[i,2] <- grid[2]
  Sim1B_Results[i,3] <- grid[3]
  Sim1B_Results[i,4] <- grid[4]
  Sim1B_Results[i,5] <- grid[5]
  
  Sim1B_Results[i,6] <- mean(Y0)
  Sim1B_Results[i,7] <- mean(Y1)
  Sim1B_Results[i,8] <- mean(Y2)
  Sim1B_Results[i,9] <- mean(Y3)
  Sim1B_Results[i,10] <- mean(Y4)
  
  Sim1B_Results[i,11] <- bart[1]
  Sim1B_Results[i,12] <- bart[2]
  Sim1B_Results[i,13] <- bart[3]
  Sim1B_Results[i,14] <- bart[4]
  Sim1B_Results[i,15] <- bart[5]
  
  Sim1B_Results[i,16] <- gam[1]
  Sim1B_Results[i,17] <- gam[2]
  Sim1B_Results[i,18] <- gam[3]
  Sim1B_Results[i,19] <- gam[4]
  Sim1B_Results[i,20] <- gam[5]
  
  Sim1B_Results[i,21] <- iptw[1]
  Sim1B_Results[i,22] <- iptw[2]
  Sim1B_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim1B_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim1B_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim1B_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim1B_Data <- rbind(Sim1B_Data,data)
}
write.csv(Sim1B_Results,"Sim1B_Results.csv")
write.csv(Sim1B_Data,"Sim1B_Data.csv")






#------------------Simulation 1C------------
#Linear, Strong Misspecifcation, No Confounding
n=100
Sim1C_Data = data.frame()
Sim1C_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N = 1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  eY <- rnorm(N, 0, 1)
  
  T <- rnorm(N, mean = 3 + 0.7*Z1 + 1.5*Z2,
             sd = 1)
  
  Y <- T + Z1 + Z2 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] + Z1 + Z2 + eY
  Y1 <- grid[2] + Z1 + Z2 + eY
  Y2 <- grid[3] + Z1 + Z2 + eY
  Y3 <- grid[4] + Z1 + Z2 + eY
  Y4 <- grid[5] + Z1 + Z2 + eY
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1100)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim1C_Results[i,1] <- grid[1]
  Sim1C_Results[i,2] <- grid[2]
  Sim1C_Results[i,3] <- grid[3]
  Sim1C_Results[i,4] <- grid[4]
  Sim1C_Results[i,5] <- grid[5]
  
  Sim1C_Results[i,6] <- mean(Y0)
  Sim1C_Results[i,7] <- mean(Y1)
  Sim1C_Results[i,8] <- mean(Y2)
  Sim1C_Results[i,9] <- mean(Y3)
  Sim1C_Results[i,10] <- mean(Y4)
  
  Sim1C_Results[i,11] <- bart[1]
  Sim1C_Results[i,12] <- bart[2]
  Sim1C_Results[i,13] <- bart[3]
  Sim1C_Results[i,14] <- bart[4]
  Sim1C_Results[i,15] <- bart[5]
  
  Sim1C_Results[i,16] <- gam[1]
  Sim1C_Results[i,17] <- gam[2]
  Sim1C_Results[i,18] <- gam[3]
  Sim1C_Results[i,19] <- gam[4]
  Sim1C_Results[i,20] <- gam[5]
  
  Sim1C_Results[i,21] <- iptw[1]
  Sim1C_Results[i,22] <- iptw[2]
  Sim1C_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim1C_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim1C_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim1C_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim1C_Data <- rbind(Sim1C_Data,data)
}
write.csv(Sim1C_Results,"Sim1C_Results.csv")
write.csv(Sim1C_Data,"Sim1C_Data.csv")





#------------------Simulation 1D------------
#Linear, Strong Misspecifcation, with Confounding
n=100
Sim1D_Data = data.frame()
Sim1D_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N = 1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  C1 <- rnorm(N, 0, 1)
  eY <- rnorm(N, 0, 1)
  
  T <- rnorm(N, mean = 3 + 0.7*Z1 + 1.5*Z2 + 0.4*C1,
             sd = 1)
  
  Y <- T + Z1 + Z2 + 0.5*C1 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] + Z1 + Z2 + 0.5*C1 + eY
  Y1 <- grid[2] + Z1 + Z2 + 0.5*C1 + eY
  Y2 <- grid[3] + Z1 + Z2 + 0.5*C1 + eY
  Y3 <- grid[4] + Z1 + Z2 + 0.5*C1 + eY
  Y4 <- grid[5] + Z1 + Z2 + 0.5*C1 + eY
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, C1, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1300)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim1D_Results[i,1] <- grid[1]
  Sim1D_Results[i,2] <- grid[2]
  Sim1D_Results[i,3] <- grid[3]
  Sim1D_Results[i,4] <- grid[4]
  Sim1D_Results[i,5] <- grid[5]
  
  Sim1D_Results[i,6] <- mean(Y0)
  Sim1D_Results[i,7] <- mean(Y1)
  Sim1D_Results[i,8] <- mean(Y2)
  Sim1D_Results[i,9] <- mean(Y3)
  Sim1D_Results[i,10] <- mean(Y4)
  
  Sim1D_Results[i,11] <- bart[1]
  Sim1D_Results[i,12] <- bart[2]
  Sim1D_Results[i,13] <- bart[3]
  Sim1D_Results[i,14] <- bart[4]
  Sim1D_Results[i,15] <- bart[5]
  
  Sim1D_Results[i,16] <- gam[1]
  Sim1D_Results[i,17] <- gam[2]
  Sim1D_Results[i,18] <- gam[3]
  Sim1D_Results[i,19] <- gam[4]
  Sim1D_Results[i,20] <- gam[5]
  
  Sim1D_Results[i,21] <- iptw[1]
  Sim1D_Results[i,22] <- iptw[2]
  Sim1D_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim1D_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim1D_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim1D_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim1D_Data <- rbind(Sim1D_Data,data)
}
write.csv(Sim1D_Results,"Sim1D_Results.csv")
write.csv(Sim1D_Data,"Sim1D_Data.csv")




#------------------Simulation 1E------------
#Linear, No Misspecifcation, No confounding

n=100
Sim1E_Data = data.frame()
Sim1E_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N = 1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + 0.3*X1 + 1.5*X2 + X2*X2,
             sd = 1)
  
  Y <- T + X1 + X2 + X1*X1 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] + X1 + X2 + X1^2 + eY
  Y1 <- grid[2] + X1 + X2 + X1^2 + eY
  Y2 <- grid[3] + X1 + X2 + X1^2 + eY
  Y3 <- grid[4] + X1 + X2 + X1^2 + eY
  Y4 <- grid[5] + X1 + X2 + X1^2 + eY
  truth <- cbind(Y0, Y1, Y2, Y3, Y4)
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim1E_Results[i,1] <- grid[1]
  Sim1E_Results[i,2] <- grid[2]
  Sim1E_Results[i,3] <- grid[3]
  Sim1E_Results[i,4] <- grid[4]
  Sim1E_Results[i,5] <- grid[5]
  
  Sim1E_Results[i,6] <- mean(Y0)
  Sim1E_Results[i,7] <- mean(Y1)
  Sim1E_Results[i,8] <- mean(Y2)
  Sim1E_Results[i,9] <- mean(Y3)
  Sim1E_Results[i,10] <- mean(Y4)
  
  Sim1E_Results[i,11] <- bart[1]
  Sim1E_Results[i,12] <- bart[2]
  Sim1E_Results[i,13] <- bart[3]
  Sim1E_Results[i,14] <- bart[4]
  Sim1E_Results[i,15] <- bart[5]
  
  Sim1E_Results[i,16] <- gam[1]
  Sim1E_Results[i,17] <- gam[2]
  Sim1E_Results[i,18] <- gam[3]
  Sim1E_Results[i,19] <- gam[4]
  Sim1E_Results[i,20] <- gam[5]
  
  Sim1E_Results[i,21] <- iptw[1]
  Sim1E_Results[i,22] <- iptw[2]
  Sim1E_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim1E_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim1E_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim1E_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim1E_Data <- rbind(Sim1E_Data,data)
}
write.csv(Sim1E_Results,"Sim1E_Results.csv")
write.csv(Sim1E_Data,"Sim1E_Data.csv")





#------------------Simulation 1F------------
#Linear, No Misspecifcation, with Confounding
n=100
Sim1F_Data = data.frame()
Sim1F_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N = 1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  C1 <- rnorm(N, 0, 1)
  eY <- rnorm(N, 0, 1)
  
  T <- rnorm(N, mean = 3 + 0.7*X1 + 1.5*X2 + 0.4*C1,
             sd = 1)
  
  Y <- T + X1 + X2 + 0.5*C1 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] + X1 + X2 + 0.5*C1 + eY
  Y1 <- grid[2] + X1 + X2 + 0.5*C1 + eY
  Y2 <- grid[3] + X1 + X2 + 0.5*C1 + eY
  Y3 <- grid[4] + X1 + X2 + 0.5*C1 + eY
  Y4 <- grid[5] + X1 + X2 + 0.5*C1 + eY
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, C1, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1300)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim1F_Results[i,1] <- grid[1]
  Sim1F_Results[i,2] <- grid[2]
  Sim1F_Results[i,3] <- grid[3]
  Sim1F_Results[i,4] <- grid[4]
  Sim1F_Results[i,5] <- grid[5]
  
  Sim1F_Results[i,6] <- mean(Y0)
  Sim1F_Results[i,7] <- mean(Y1)
  Sim1F_Results[i,8] <- mean(Y2)
  Sim1F_Results[i,9] <- mean(Y3)
  Sim1F_Results[i,10] <- mean(Y4)
  
  Sim1F_Results[i,11] <- bart[1]
  Sim1F_Results[i,12] <- bart[2]
  Sim1F_Results[i,13] <- bart[3]
  Sim1F_Results[i,14] <- bart[4]
  Sim1F_Results[i,15] <- bart[5]
  
  Sim1F_Results[i,16] <- gam[1]
  Sim1F_Results[i,17] <- gam[2]
  Sim1F_Results[i,18] <- gam[3]
  Sim1F_Results[i,19] <- gam[4]
  Sim1F_Results[i,20] <- gam[5]
  
  Sim1F_Results[i,21] <- iptw[1]
  Sim1F_Results[i,22] <- iptw[2]
  Sim1F_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim1F_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim1F_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim1F_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim1F_Data <- rbind(Sim1F_Data,data)
}
write.csv(Sim1F_Results,"Sim1F_Results.csv")
write.csv(Sim1F_Data,"Sim1F_Data.csv")










#------------------Simulation 2A------------
#Non-Linear, Moderate Misspecifcation, No confounding

n=100
Sim2A_Data = data.frame()
Sim2A_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N=1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  C1 <- rnorm(N, 0, 1)
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + 0.5*X1 + 0.5*Z1 + 0.5*Z2 + 1.5*X2, sd = 1)
  
  Y <- T - 0.3*T^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.3*grid[1]^2 +  0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  Y1 <- grid[2] - 0.3*grid[2]^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  Y2 <- grid[3] - 0.3*grid[3]^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  Y3 <- grid[4] - 0.3*grid[4]^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  Y4 <- grid[5] - 0.3*grid[5]^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + eY
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim2A_Results[i,1] <- grid[1]
  Sim2A_Results[i,2] <- grid[2]
  Sim2A_Results[i,3] <- grid[3]
  Sim2A_Results[i,4] <- grid[4]
  Sim2A_Results[i,5] <- grid[5]
  
  Sim2A_Results[i,6] <- mean(Y0)
  Sim2A_Results[i,7] <- mean(Y1)
  Sim2A_Results[i,8] <- mean(Y2)
  Sim2A_Results[i,9] <- mean(Y3)
  Sim2A_Results[i,10] <- mean(Y4)
  
  Sim2A_Results[i,11] <- bart[1]
  Sim2A_Results[i,12] <- bart[2]
  Sim2A_Results[i,13] <- bart[3]
  Sim2A_Results[i,14] <- bart[4]
  Sim2A_Results[i,15] <- bart[5]
  
  Sim2A_Results[i,16] <- gam[1]
  Sim2A_Results[i,17] <- gam[2]
  Sim2A_Results[i,18] <- gam[3]
  Sim2A_Results[i,19] <- gam[4]
  Sim2A_Results[i,20] <- gam[5]
  
  Sim2A_Results[i,21] <- iptw[1]
  Sim2A_Results[i,22] <- iptw[2]
  Sim2A_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim2A_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim2A_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim2A_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim2A_Data <- rbind(Sim2A_Data,data)
}
write.csv(Sim2A_Results,"Sim2A_Results.csv")
write.csv(Sim2A_Data,"Sim2A_Data.csv")










#------------------Simulation 2B------------
#Non-Linear, Moderate Misspecifcation, with confounding

n=100
Sim2B_Data = data.frame()
Sim2B_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N=1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  C1 <- rnorm(N, 0, 1)
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + 0.5*X1 + 0.5*Z1 + 0.5*Z2 + 1.5*X2 + 0.6*C1, sd = 1)
  
  Y <- T - 0.3*T^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.4*C1 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.3*grid[1]^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.4*C1 + eY
  Y1 <- grid[2] - 0.3*grid[2]^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.4*C1 + eY
  Y2 <- grid[3] - 0.3*grid[3]^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.4*C1 + eY
  Y3 <- grid[4] - 0.3*grid[4]^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.4*C1 + eY
  Y4 <- grid[5] - 0.3*grid[5]^2 + 0.5*X1 + 0.5*Z1 + 0.5*X2 + 0.5*Z2 + 0.4*C1 + eY
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, C1, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim2B_Results[i,1] <- grid[1]
  Sim2B_Results[i,2] <- grid[2]
  Sim2B_Results[i,3] <- grid[3]
  Sim2B_Results[i,4] <- grid[4]
  Sim2B_Results[i,5] <- grid[5]
  
  Sim2B_Results[i,6] <- mean(Y0)
  Sim2B_Results[i,7] <- mean(Y1)
  Sim2B_Results[i,8] <- mean(Y2)
  Sim2B_Results[i,9] <- mean(Y3)
  Sim2B_Results[i,10] <- mean(Y4)
  
  Sim2B_Results[i,11] <- bart[1]
  Sim2B_Results[i,12] <- bart[2]
  Sim2B_Results[i,13] <- bart[3]
  Sim2B_Results[i,14] <- bart[4]
  Sim2B_Results[i,15] <- bart[5]
  
  Sim2B_Results[i,16] <- gam[1]
  Sim2B_Results[i,17] <- gam[2]
  Sim2B_Results[i,18] <- gam[3]
  Sim2B_Results[i,19] <- gam[4]
  Sim2B_Results[i,20] <- gam[5]
  
  Sim2B_Results[i,21] <- iptw[1]
  Sim2B_Results[i,22] <- iptw[2]
  Sim2B_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim2B_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim2B_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim2B_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim2B_Data <- rbind(Sim2B_Data,data)
}
write.csv(Sim2B_Results,"Sim2B_Results.csv")
write.csv(Sim2B_Data,"Sim2B_Data.csv")










#------------------Simulation 2C------------
#Non-Linear, Strong Misspecification, No confounding

n=100
Sim2C_Data = data.frame()
Sim2C_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N=1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 2 + Z1 + Z2, sd = 1.5)
  
  Y <- T - 0.5*T^2 + Z1 + Z2 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.3*grid[1]^2 + Z1 + Z2 + eY
  Y1 <- grid[2] - 0.3*grid[2]^2 + Z1 + Z2 + eY
  Y2 <- grid[3] - 0.3*grid[3]^2 + Z1 + Z2 + eY
  Y3 <- grid[4] - 0.3*grid[4]^2 + Z1 + Z2 + eY
  Y4 <- grid[5] - 0.3*grid[5]^2 + Z1 + Z2 + eY
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim2C_Results[i,1] <- grid[1]
  Sim2C_Results[i,2] <- grid[2]
  Sim2C_Results[i,3] <- grid[3]
  Sim2C_Results[i,4] <- grid[4]
  Sim2C_Results[i,5] <- grid[5]
  
  Sim2C_Results[i,6] <- mean(Y0)
  Sim2C_Results[i,7] <- mean(Y1)
  Sim2C_Results[i,8] <- mean(Y2)
  Sim2C_Results[i,9] <- mean(Y3)
  Sim2C_Results[i,10] <- mean(Y4)
  
  Sim2C_Results[i,11] <- bart[1]
  Sim2C_Results[i,12] <- bart[2]
  Sim2C_Results[i,13] <- bart[3]
  Sim2C_Results[i,14] <- bart[4]
  Sim2C_Results[i,15] <- bart[5]
  
  Sim2C_Results[i,16] <- gam[1]
  Sim2C_Results[i,17] <- gam[2]
  Sim2C_Results[i,18] <- gam[3]
  Sim2C_Results[i,19] <- gam[4]
  Sim2C_Results[i,20] <- gam[5]
  
  Sim2C_Results[i,21] <- iptw[1]
  Sim2C_Results[i,22] <- iptw[2]
  Sim2C_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim2C_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim2C_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim2C_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim2C_Data <- rbind(Sim2C_Data,data)
}
write.csv(Sim2C_Results,"Sim2C_Results.csv")
write.csv(Sim2C_Data,"Sim2C_Data.csv")






#------------------Simulation 2D------------
#Non-Linear, Strong Misspecifcation, with confounding

n=100
Sim2D_Data = data.frame()
Sim2D_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N=1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  C1 <- rnorm(N, 0, 1)
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + Z1 + Z2 + 0.6*C1, sd = 1)
  
  Y <- T - 0.3*T^2 + Z1 + Z2 + 0.4*C1 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.3*grid[1]^2 + Z1 + Z2 + 0.4*C1 + eY
  Y1 <- grid[2] - 0.3*grid[2]^2 + Z1 + Z2 + 0.4*C1 + eY
  Y2 <- grid[3] - 0.3*grid[3]^2 + Z1 + Z2 + 0.4*C1 + eY
  Y3 <- grid[4] - 0.3*grid[4]^2 + Z1 + Z2 + 0.4*C1 + eY
  Y4 <- grid[5] - 0.3*grid[5]^2 + Z1 + Z2 + 0.4*C1 + eY
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, C1, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim2D_Results[i,1] <- grid[1]
  Sim2D_Results[i,2] <- grid[2]
  Sim2D_Results[i,3] <- grid[3]
  Sim2D_Results[i,4] <- grid[4]
  Sim2D_Results[i,5] <- grid[5]
  
  Sim2D_Results[i,6] <- mean(Y0)
  Sim2D_Results[i,7] <- mean(Y1)
  Sim2D_Results[i,8] <- mean(Y2)
  Sim2D_Results[i,9] <- mean(Y3)
  Sim2D_Results[i,10] <- mean(Y4)
  
  Sim2D_Results[i,11] <- bart[1]
  Sim2D_Results[i,12] <- bart[2]
  Sim2D_Results[i,13] <- bart[3]
  Sim2D_Results[i,14] <- bart[4]
  Sim2D_Results[i,15] <- bart[5]
  
  Sim2D_Results[i,16] <- gam[1]
  Sim2D_Results[i,17] <- gam[2]
  Sim2D_Results[i,18] <- gam[3]
  Sim2D_Results[i,19] <- gam[4]
  Sim2D_Results[i,20] <- gam[5]
  
  Sim2D_Results[i,21] <- iptw[1]
  Sim2D_Results[i,22] <- iptw[2]
  Sim2D_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim2D_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim2D_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim2D_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim2D_Data <- rbind(Sim2D_Data,data)
}
write.csv(Sim2D_Results,"Sim2D_Results.csv")
write.csv(Sim2D_Data,"Sim2D_Data.csv")




#------------------Simulation 2E------------
#Non-Linear, No Misspecifcation, No confounding

n=100
Sim2E_Data = data.frame()
Sim2E_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N=1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  C1 <- rnorm(N, 0, 1)
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + 0.6*X1 + 0.4*X2, sd = 1)
  
  Y <- T - 0.3*T^2 + X1 + X2 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.3*grid[1]^2 + X1 + X2 + eY
  Y1 <- grid[2] - 0.3*grid[2]^2 + X1 + X2 + eY
  Y2 <- grid[3] - 0.3*grid[3]^2 + X1 + X2 + eY
  Y3 <- grid[4] - 0.3*grid[4]^2 + X1 + X2 + eY
  Y4 <- grid[5] - 0.3*grid[5]^2 + X1 + X2 + eY
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim2E_Results[i,1] <- grid[1]
  Sim2E_Results[i,2] <- grid[2]
  Sim2E_Results[i,3] <- grid[3]
  Sim2E_Results[i,4] <- grid[4]
  Sim2E_Results[i,5] <- grid[5]
  
  Sim2E_Results[i,6] <- mean(Y0)
  Sim2E_Results[i,7] <- mean(Y1)
  Sim2E_Results[i,8] <- mean(Y2)
  Sim2E_Results[i,9] <- mean(Y3)
  Sim2E_Results[i,10] <- mean(Y4)
  
  Sim2E_Results[i,11] <- bart[1]
  Sim2E_Results[i,12] <- bart[2]
  Sim2E_Results[i,13] <- bart[3]
  Sim2E_Results[i,14] <- bart[4]
  Sim2E_Results[i,15] <- bart[5]
  
  Sim2E_Results[i,16] <- gam[1]
  Sim2E_Results[i,17] <- gam[2]
  Sim2E_Results[i,18] <- gam[3]
  Sim2E_Results[i,19] <- gam[4]
  Sim2E_Results[i,20] <- gam[5]
  
  Sim2E_Results[i,21] <- iptw[1]
  Sim2E_Results[i,22] <- iptw[2]
  Sim2E_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim2E_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim2E_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim2E_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim2E_Data <- rbind(Sim2E_Data,data)
}
write.csv(Sim2E_Results,"Sim2E_Results.csv")
write.csv(Sim2E_Data,"Sim2E_Data.csv")



#------------------Simulation 2F------------
#Non-Linear, Strong Misspecifcation, with confounding

n=100
Sim2F_Data = data.frame()
Sim2F_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N=1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  C1 <- rnorm(N, 0, 1)
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + X1 + X2 + 0.6*C1, sd = 1)
  
  Y <- T - 0.3*T^2 + X1 + X2 + 0.4*C1 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.3*grid[1]^2 + X1 + X2 + 0.4*C1 + eY
  Y1 <- grid[2] - 0.3*grid[2]^2 + X1 + X2 + 0.4*C1 + eY
  Y2 <- grid[3] - 0.3*grid[3]^2 + X1 + X2 + 0.4*C1 + eY
  Y3 <- grid[4] - 0.3*grid[4]^2 + X1 + X2 + 0.4*C1 + eY
  Y4 <- grid[5] - 0.3*grid[5]^2 + X1 + X2 + 0.4*C1 + eY
  
  data <- data.frame(cbind(X1, X2, Z1, Z2, C1, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim2F_Results[i,1] <- grid[1]
  Sim2F_Results[i,2] <- grid[2]
  Sim2F_Results[i,3] <- grid[3]
  Sim2F_Results[i,4] <- grid[4]
  Sim2F_Results[i,5] <- grid[5]
  
  Sim2F_Results[i,6] <- mean(Y0)
  Sim2F_Results[i,7] <- mean(Y1)
  Sim2F_Results[i,8] <- mean(Y2)
  Sim2F_Results[i,9] <- mean(Y3)
  Sim2F_Results[i,10] <- mean(Y4)
  
  Sim2F_Results[i,11] <- bart[1]
  Sim2F_Results[i,12] <- bart[2]
  Sim2F_Results[i,13] <- bart[3]
  Sim2F_Results[i,14] <- bart[4]
  Sim2F_Results[i,15] <- bart[5]
  
  Sim2F_Results[i,16] <- gam[1]
  Sim2F_Results[i,17] <- gam[2]
  Sim2F_Results[i,18] <- gam[3]
  Sim2F_Results[i,19] <- gam[4]
  Sim2F_Results[i,20] <- gam[5]
  
  Sim2F_Results[i,21] <- iptw[1]
  Sim2F_Results[i,22] <- iptw[2]
  Sim2F_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim2F_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim2F_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim2F_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim2F_Data <- rbind(Sim2F_Data,data)
}
write.csv(Sim2F_Results,"Sim2F_Results.csv")
write.csv(Sim2F_Data,"Sim2F_Data.csv")










#---------------------Simulations 3-----------------
#-------Sim 3A----
#Moderate Misspecification, One Confounder
n=100
Sim3A_Data = data.frame()
Sim3A_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N = 1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  X3 <- rbinom(N, 1, 0.5)
  X4 <- rnorm(N, 0, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  Z3 <- (0.2*X1+X3)
  C1 <- rnorm(N, 0, 1)
  C2 <- rbinom(N, 1, 0.3)
  
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.6*C1, sd = 1)
  
  Y <- T - 0.5*T^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + T*X4 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.5*grid[1]^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + grid[1]*X4 + eY
  Y1 <- grid[2] - 0.5*grid[2]^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + grid[2]*X4 + eY
  Y2 <- grid[3] - 0.5*grid[3]^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + grid[3]*X4 + eY
  Y3 <- grid[4] - 0.5*grid[4]^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + grid[4]*X4 + eY
  Y4 <- grid[5] - 0.5*grid[5]^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + grid[5]*X4 + eY
  
  data <- data.frame(cbind(X1, X2, X3, X4, Z1, Z2, Z3, C1, C2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2 + X3 + X4,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2 + X3 + X4,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2 + X3 + X4,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim3A_Results[i,1] <- grid[1]
  Sim3A_Results[i,2] <- grid[2]
  Sim3A_Results[i,3] <- grid[3]
  Sim3A_Results[i,4] <- grid[4]
  Sim3A_Results[i,5] <- grid[5]
  
  Sim3A_Results[i,6] <- mean(Y0)
  Sim3A_Results[i,7] <- mean(Y1)
  Sim3A_Results[i,8] <- mean(Y2)
  Sim3A_Results[i,9] <- mean(Y3)
  Sim3A_Results[i,10] <- mean(Y4)
  
  Sim3A_Results[i,11] <- bart[1]
  Sim3A_Results[i,12] <- bart[2]
  Sim3A_Results[i,13] <- bart[3]
  Sim3A_Results[i,14] <- bart[4]
  Sim3A_Results[i,15] <- bart[5]
  
  Sim3A_Results[i,16] <- gam[1]
  Sim3A_Results[i,17] <- gam[2]
  Sim3A_Results[i,18] <- gam[3]
  Sim3A_Results[i,19] <- gam[4]
  Sim3A_Results[i,20] <- gam[5]
  
  Sim3A_Results[i,21] <- iptw[1]
  Sim3A_Results[i,22] <- iptw[2]
  Sim3A_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim3A_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim3A_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim3A_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim3A_Data <- rbind(Sim3A_Data,data)
}
write.csv(Sim3A_Results,"Sim3A_Results_1.csv")
write.csv(Sim3A_Data,"Sim3A_Data_1.csv")







#------------------Sim 3B-----
#Moderate Misspecification, Two Confounders
n=100
Sim3B_Data = data.frame()
Sim3B_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N = 1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  X3 <- rbinom(N, 1, 0.5)
  X4 <- rnorm(N, 0, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  Z3 <- (0.2*X1+X3)
  C1 <- rnorm(N, 0, 1)
  C2 <- rbinom(N, 1, 0.3)
  
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.6*C1 + 0.3*C2 + Z1*Z2, sd = 1)
  
  Y <- T - 0.5*T^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + 0.5*C2 + T*X4 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.5*grid[1]^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + 0.5*C2 + grid[1]*X4 + eY
  Y1 <- grid[2] - 0.5*grid[2]^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + 0.5*C2 + grid[2]*X4 + eY
  Y2 <- grid[3] - 0.5*grid[3]^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + 0.5*C2 + grid[3]*X4 + eY
  Y3 <- grid[4] - 0.5*grid[4]^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + 0.5*C2 + grid[4]*X4 + eY
  Y4 <- grid[5] - 0.5*grid[5]^2 + 0.5*X1 + 0.5*X2 + 0.5*X3 + 0.5*Z1 + 0.5*Z2 + 0.5*Z3 + 0.3*C1 + 0.5*C2 + grid[5]*X4 + eY
  
  data <- data.frame(cbind(X1, X2, X3, X4, Z1, Z2, Z3, C1, C2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2 + X3 + X4,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2 + X3 + X4,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2 + X3 + X4,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim3B_Results[i,1] <- grid[1]
  Sim3B_Results[i,2] <- grid[2]
  Sim3B_Results[i,3] <- grid[3]
  Sim3B_Results[i,4] <- grid[4]
  Sim3B_Results[i,5] <- grid[5]
  
  Sim3B_Results[i,6] <- mean(Y0)
  Sim3B_Results[i,7] <- mean(Y1)
  Sim3B_Results[i,8] <- mean(Y2)
  Sim3B_Results[i,9] <- mean(Y3)
  Sim3B_Results[i,10] <- mean(Y4)
  
  Sim3B_Results[i,11] <- bart[1]
  Sim3B_Results[i,12] <- bart[2]
  Sim3B_Results[i,13] <- bart[3]
  Sim3B_Results[i,14] <- bart[4]
  Sim3B_Results[i,15] <- bart[5]
  
  Sim3B_Results[i,16] <- gam[1]
  Sim3B_Results[i,17] <- gam[2]
  Sim3B_Results[i,18] <- gam[3]
  Sim3B_Results[i,19] <- gam[4]
  Sim3B_Results[i,20] <- gam[5]
  
  Sim3B_Results[i,21] <- iptw[1]
  Sim3B_Results[i,22] <- iptw[2]
  Sim3B_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim3B_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim3B_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim3B_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim3B_Data <- rbind(Sim3B_Data,data)
}
write.csv(Sim3B_Results,"Sim3B_Results_1.csv")
write.csv(Sim3B_Data,"Sim3B_Data_1.csv")







#------------------Simulation 3C------------
#Non-Linear, Strong Misspecification, One confounder

n=100
Sim3C_Data = data.frame()
Sim3C_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N=1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  X3 <- rbinom(N, 1, 0.5)
  X4 <- rnorm(N, 0, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  Z3 <- (0.2*X1+X3)
  C1 <- rnorm(N, 0, 1)
  C2 <- rbinom(N, 1, 0.3)
  
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + Z1 + Z2 + Z3 + 0.6*C1 + Z1*Z2, sd = 1)
  
  Y <- T - 0.5*T^2 + Z1 + Z2 + Z3 + 0.3*C1 + T*X4 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.5*grid[1]^2 + Z1 + Z2 + Z3 + 0.3*C1 + grid[1]*X4 + eY
  Y1 <- grid[2] - 0.5*grid[2]^2 + Z1 + Z2 + Z3 + 0.3*C1 + grid[2]*X4 + eY
  Y2 <- grid[3] - 0.5*grid[3]^2 + Z1 + Z2 + Z3 + 0.3*C1 + grid[3]*X4 + eY
  Y3 <- grid[4] - 0.5*grid[4]^2 + Z1 + Z2 + Z3 + 0.3*C1 + grid[4]*X4 + eY
  Y4 <- grid[5] - 0.5*grid[5]^2 + Z1 + Z2 + Z3 + 0.3*C1 + grid[5]*X4 + eY
  
  data <- data.frame(cbind(X1, X2, X3, X4, Z1, Z2, C1, C2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2 + X3 + X4,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2 + X3 + X4,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2 + X3 + X4,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim3C_Results[i,1] <- grid[1]
  Sim3C_Results[i,2] <- grid[2]
  Sim3C_Results[i,3] <- grid[3]
  Sim3C_Results[i,4] <- grid[4]
  Sim3C_Results[i,5] <- grid[5]
  
  Sim3C_Results[i,6] <- mean(Y0)
  Sim3C_Results[i,7] <- mean(Y1)
  Sim3C_Results[i,8] <- mean(Y2)
  Sim3C_Results[i,9] <- mean(Y3)
  Sim3C_Results[i,10] <- mean(Y4)
  
  Sim3C_Results[i,11] <- bart[1]
  Sim3C_Results[i,12] <- bart[2]
  Sim3C_Results[i,13] <- bart[3]
  Sim3C_Results[i,14] <- bart[4]
  Sim3C_Results[i,15] <- bart[5]
  
  Sim3C_Results[i,16] <- gam[1]
  Sim3C_Results[i,17] <- gam[2]
  Sim3C_Results[i,18] <- gam[3]
  Sim3C_Results[i,19] <- gam[4]
  Sim3C_Results[i,20] <- gam[5]
  
  Sim3C_Results[i,21] <- iptw[1]
  Sim3C_Results[i,22] <- iptw[2]
  Sim3C_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim3C_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim3C_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim3C_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim3C_Data <- rbind(Sim3C_Data,data)
}
write.csv(Sim3C_Results,"Sim3C_Results_1.csv")
write.csv(Sim3C_Data,"Sim3C_Data_1.csv")


#------------------Simulation 3D------------
#Non-Linear, Strong Misspecification, Two confounders

n=100
Sim3D_Data = data.frame()
Sim3D_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N=1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  X3 <- rbinom(N, 1, 0.5)
  X4 <- rnorm(N, 0, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  Z3 <- (0.2*X1+X3)
  C1 <- rnorm(N, 0, 1)
  C2 <- rbinom(N, 1, 0.3)
  
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + Z1 + Z2 + Z3 + 0.6*C1 + 0.3*C2 + Z1*Z2, sd = 1)
  
  Y <- T - 0.5*T^2 + Z1 + Z2 + Z3 + 0.3*C1 + 0.5*C2 + T*X4 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.5*grid[1]^2 + Z1 + Z2 + Z3 + 0.3*C1 + 0.5*C2 + grid[1]*X4 + eY
  Y1 <- grid[2] - 0.5*grid[2]^2 + Z1 + Z2 + Z3 + 0.3*C1 + 0.5*C2 + grid[2]*X4 + eY
  Y2 <- grid[3] - 0.5*grid[3]^2 + Z1 + Z2 + Z3 + 0.3*C1 + 0.5*C2 + grid[3]*X4 + eY
  Y3 <- grid[4] - 0.5*grid[4]^2 + Z1 + Z2 + Z3 + 0.3*C1 + 0.5*C2 + grid[4]*X4 + eY
  Y4 <- grid[5] - 0.5*grid[5]^2 + Z1 + Z2 + Z3 + 0.3*C1 + 0.5*C2 + grid[5]*X4 + eY
  
  data <- data.frame(cbind(X1, X2, X3, X4, Z1, Z2, C1, C2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2 + X3 + X4,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2 + X3 + X4,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2 + X3 + X4,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim3D_Results[i,1] <- grid[1]
  Sim3D_Results[i,2] <- grid[2]
  Sim3D_Results[i,3] <- grid[3]
  Sim3D_Results[i,4] <- grid[4]
  Sim3D_Results[i,5] <- grid[5]
  
  Sim3D_Results[i,6] <- mean(Y0)
  Sim3D_Results[i,7] <- mean(Y1)
  Sim3D_Results[i,8] <- mean(Y2)
  Sim3D_Results[i,9] <- mean(Y3)
  Sim3D_Results[i,10] <- mean(Y4)
  
  Sim3D_Results[i,11] <- bart[1]
  Sim3D_Results[i,12] <- bart[2]
  Sim3D_Results[i,13] <- bart[3]
  Sim3D_Results[i,14] <- bart[4]
  Sim3D_Results[i,15] <- bart[5]
  
  Sim3D_Results[i,16] <- gam[1]
  Sim3D_Results[i,17] <- gam[2]
  Sim3D_Results[i,18] <- gam[3]
  Sim3D_Results[i,19] <- gam[4]
  Sim3D_Results[i,20] <- gam[5]
  
  Sim3D_Results[i,21] <- iptw[1]
  Sim3D_Results[i,22] <- iptw[2]
  Sim3D_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim3D_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim3D_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim3D_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim3D_Data <- rbind(Sim3D_Data,data)
}
write.csv(Sim3D_Results,"Sim3D_Results_1.csv")
write.csv(Sim3D_Data,"Sim3D_Data_1.csv")








#-------Simulation 3E----
#No Misspecification, One Confounder
n=100
Sim3E_Data = data.frame()
Sim3E_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N = 1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  X3 <- rbinom(N, 1, 0.5)
  X4 <- rnorm(N, 0, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  Z3 <- (0.2*X1+X3)
  C1 <- rnorm(N, 0, 1)
  C2 <- rbinom(N, 1, 0.3)
  
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + 0.6*X1 + 0.4*X2 + 0.8*X3 + 0.3*X4 + 0.6*C1, sd = 1)
  
  Y <- T - 0.5*T^2 + X1 + X2 + X3 + 0.3*C1 + T*X4 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.5*grid[1]^2 + X1 + X2 + X3 + 0.3*C1 + grid[1]*X4 + eY
  Y1 <- grid[2] - 0.5*grid[2]^2 + X1 + X2 + X3 + 0.3*C1 + grid[2]*X4 + eY
  Y2 <- grid[3] - 0.5*grid[3]^2 + X1 + X2 + X3 + 0.3*C1 + grid[3]*X4 + eY
  Y3 <- grid[4] - 0.5*grid[4]^2 + X1 + X2 + X3 + 0.3*C1 + grid[4]*X4 + eY
  Y4 <- grid[5] - 0.5*grid[5]^2 + X1 + X2 + X3 + 0.3*C1 + grid[5]*X4 + eY
  
  data <- data.frame(cbind(X1, X2, X3, X4, Z1, Z2, Z3, C1, C2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2 + X3 + X4,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2 + X3 + X4,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2 + X3 + X4,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim3E_Results[i,1] <- grid[1]
  Sim3E_Results[i,2] <- grid[2]
  Sim3E_Results[i,3] <- grid[3]
  Sim3E_Results[i,4] <- grid[4]
  Sim3E_Results[i,5] <- grid[5]
  
  Sim3E_Results[i,6] <- mean(Y0)
  Sim3E_Results[i,7] <- mean(Y1)
  Sim3E_Results[i,8] <- mean(Y2)
  Sim3E_Results[i,9] <- mean(Y3)
  Sim3E_Results[i,10] <- mean(Y4)
  
  Sim3E_Results[i,11] <- bart[1]
  Sim3E_Results[i,12] <- bart[2]
  Sim3E_Results[i,13] <- bart[3]
  Sim3E_Results[i,14] <- bart[4]
  Sim3E_Results[i,15] <- bart[5]
  
  Sim3E_Results[i,16] <- gam[1]
  Sim3E_Results[i,17] <- gam[2]
  Sim3E_Results[i,18] <- gam[3]
  Sim3E_Results[i,19] <- gam[4]
  Sim3E_Results[i,20] <- gam[5]
  
  Sim3E_Results[i,21] <- iptw[1]
  Sim3E_Results[i,22] <- iptw[2]
  Sim3E_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim3E_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim3E_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim3E_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim3E_Data <- rbind(Sim3E_Data,data)
}
write.csv(Sim3E_Results,"Sim3E_Results.csv")
write.csv(Sim3E_Data,"Sim3E_Data.csv")








#-------Simulation 3F----
#No Misspecification, One Confounder
n=100
Sim3F_Data = data.frame()
Sim3F_Results <- data.frame(matrix(nrow = 100, ncol = 23))
for (i in 1:n ){
  N = 1000
  X1 <- rnorm(N, 0.5, 1)
  X2 <- rnorm(N, 0.5, 1)
  X3 <- rbinom(N, 1, 0.5)
  X4 <- rnorm(N, 0, 1)
  Z1 <- X2/exp(X1)
  Z2 <- log(abs(X2))
  Z3 <- (0.2*X1+X3)
  C1 <- rnorm(N, 0, 1)
  C2 <- rbinom(N, 1, 0.3)
  
  eY <- rnorm(N, 0, 1)
  
  
  T <- rnorm(N, mean = 3 + 0.6*X1 + 0.4*X2 + 0.8*X3 + 0.3*X4 + 0.6*C1 + 0.3*C2, sd = 1)
  
  Y <- T - 0.5*T^2 + X1 + X2 + X3 + 0.3*C1 +  0.5*C2 + T*X4 + eY
  
  grid <- as.vector(quantile(T, probs = seq(0, 1, 0.25)))
  Y0 <- grid[1] - 0.5*grid[1]^2 + X1 + X2 + X3 + 0.3*C1 + 0.5*C2 + grid[1]*X4 + eY
  Y1 <- grid[2] - 0.5*grid[2]^2 + X1 + X2 + X3 + 0.3*C1 + 0.5*C2 + grid[2]*X4 + eY
  Y2 <- grid[3] - 0.5*grid[3]^2 + X1 + X2 + X3 + 0.3*C1 + 0.5*C2 + grid[3]*X4 + eY
  Y3 <- grid[4] - 0.5*grid[4]^2 + X1 + X2 + X3 + 0.3*C1 + 0.5*C2 + grid[4]*X4 + eY
  Y4 <- grid[5] - 0.5*grid[5]^2 + X1 + X2 + X3 + 0.3*C1 + 0.5*C2 + grid[5]*X4 + eY
  
  data <- data.frame(cbind(X1, X2, X3, X4, Z1, Z2, Z3, C1, C2, eY, T, Y))
  
  
  #BART
  bart <-  bart_est(Y = Y,
                    treat = T,
                    outcome_formula = Y ~ T + X1 + X2 + X3 + X4,
                    data = data,
                    grid_val = grid,
                    ndpost=1000)
  
  bart <- as.vector(bart$param)
  
  #GAM
  gam <- gam_est(Y = Y,
                 treat = T,
                 treat_formula = T ~ X1 + X2 + X3 + X4,
                 data = data,
                 grid_val = grid,
                 treat_mod = "Normal")
  
  gam <- as.vector(gam$param)
  
  #IPTW
  iptw <- iptw_est(Y = Y,
                   treat = T,
                   treat_formula = T ~ X1 + X2 + X3 + X4,
                   numerator_formula = T ~ 1,
                   data = data,
                   degree = 1,
                   treat_mod = "Normal")
  
  iptw <- as.vector(iptw$param)
  
  #Store results in row i of matrix
  Sim3F_Results[i,1] <- grid[1]
  Sim3F_Results[i,2] <- grid[2]
  Sim3F_Results[i,3] <- grid[3]
  Sim3F_Results[i,4] <- grid[4]
  Sim3F_Results[i,5] <- grid[5]
  
  Sim3F_Results[i,6] <- mean(Y0)
  Sim3F_Results[i,7] <- mean(Y1)
  Sim3F_Results[i,8] <- mean(Y2)
  Sim3F_Results[i,9] <- mean(Y3)
  Sim3F_Results[i,10] <- mean(Y4)
  
  Sim3F_Results[i,11] <- bart[1]
  Sim3F_Results[i,12] <- bart[2]
  Sim3F_Results[i,13] <- bart[3]
  Sim3F_Results[i,14] <- bart[4]
  Sim3F_Results[i,15] <- bart[5]
  
  Sim3F_Results[i,16] <- gam[1]
  Sim3F_Results[i,17] <- gam[2]
  Sim3F_Results[i,18] <- gam[3]
  Sim3F_Results[i,19] <- gam[4]
  Sim3F_Results[i,20] <- gam[5]
  
  Sim3F_Results[i,21] <- iptw[1]
  Sim3F_Results[i,22] <- iptw[2]
  Sim3F_Results[i,23] <- iptw[1] + iptw[2]*grid[2]
  Sim3F_Results[i,24] <- iptw[1] + iptw[2]*grid[3]
  Sim3F_Results[i,25] <- iptw[1] + iptw[2]*grid[4]
  Sim3F_Results[i,26] <- iptw[1] + iptw[2]*grid[5]
  
  #Save Data for Use in TMLE and GPS
  Sim3F_Data <- rbind(Sim3F_Data,data)
}
write.csv(Sim3F_Results,"Sim3F_Results.csv")
write.csv(Sim3F_Data,"Sim3F_Data.csv")