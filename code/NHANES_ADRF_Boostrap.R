#----------------------Bootstrapping Standard Errors--------------
#Code to bootstrap standard errors (to obtain 95% CIs) for NHANES ADRF estimates

#------BART Function for Boostrap------
BART <- function(y, n, t, formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  grid <- quantile(d$t, probs = seq(0, 1, 1/20))
  estimate <- eval.parent(substitute(bart_est(Y = y,
                                              treat = t,
                                              outcome_formula = formula, 
                                              data=d,
                                              grid_val = grid,
                                              ndpost=n)))
  
  return(estimate$param)
}





#------------Total PFAS Bootstrapping-----------
#----Overall Total PFAS
PFAS_Boot <- boot(data=nhanes_final, 
                  statistic=bart,
                  R=100, 
                  y = eGFR, 
                  n=1000, 
                  t = PFAS_log,
                  formula=eGFR~PFAS_log+age+race+sex+bmi_log+diabetes+cotinine)

#----Male Total PFAS
PFAS_Male_Boot <- boot(data=nhanes_male, 
                       statistic=bart,
                       R=100, 
                       y = eGFR, 
                       n=1000, 
                       t = PFAS_log,
                       formula=eGFR~PFAS_log+age+race+bmi_log+diabetes+cotinine)

#----Female Total PFAS
PFOS_Female_Boot <- boot(data=nhanes_female, 
                         statistic=bart,
                         R=100, 
                         y = eGFR, 
                         n=1000, 
                         t = PFAS_log,
                         formula=eGFR~PFAS_log+age+race+bmi_log+diabetes+cotinine)




#------------PFOA Bootstrapping-----------
#----Overall PFOA
PFOA_Boot <- boot(data=nhanes_final, 
                  statistic=bart,
                  R=100, 
                  y = eGFR, 
                  n=1000, 
                  t = PFOA_log,
                  formula=eGFR~PFOA_log+age+race+sex+bmi_log+diabetes+cotinine)

#----Male PFOA
PFOA_Male_Boot <- boot(data=nhanes_male, 
                       statistic=bart,
                       R=100, 
                       y = eGFR, 
                       n=1000, 
                       t = PFOA_log,
                       formula=eGFR~PFOA_log+age+race+bmi_log+diabetes+cotinine)

#----Female PFOA
PFOA_Female_Boot <- boot(data=nhanes_female, 
                         statistic=bart,
                         R=100, 
                         y = eGFR, 
                         n=1000, 
                         t = PFOA_log,
                         formula=eGFR~PFOA_log+age+race+bmi_log+diabetes+cotinine)



#------------PFOS Bootstrapping-----------
#----Overall PFOS
PFOS_Boot <- boot(data=nhanes_final, 
                  statistic=bart,
                  R=100, 
                  y = eGFR, 
                  n=1000, 
                  t = PFOS_log,
                  formula=eGFR~PFOS_log+age+race+sex+bmi_log+diabetes+cotinine)

#----Male PFOS
PFOS_Male_Boot <- boot(data=nhanes_male, 
                       statistic=bart,
                       R=100, 
                       y = eGFR, 
                       n=1000, 
                       t = PFOS_log,
                       formula=eGFR~PFOS_log+age+race+bmi_log+diabetes+cotinine)

#----Female PFOS
PFOS_Female_Boot <- boot(data=nhanes_female, 
                         statistic=bart,
                         R=100, 
                         y = eGFR, 
                         n=1000, 
                         t = PFOS_log,
                         formula=eGFR~PFOS_log+age+race+bmi_log+diabetes+cotinine)




#------------PFHxS Bootstrapping-----------
#----Overall PFHxS
PFHS_Boot <- boot(data=nhanes_final, 
                  statistic=bart,
                  R=100, 
                  y = eGFR, 
                  n=1000, 
                  t = PFHS_log,
                  formula=eGFR~PFHS_log+age+race+sex+bmi_log+diabetes+cotinine)

#----Male PFHxS
PFHS_Male_Boot <- boot(data=nhanes_male, 
                       statistic=bart,
                       R=100, 
                       y = eGFR, 
                       n=1000, 
                       t = PFHS_log,
                       formula=eGFR~PFHS_log+age+race+bmi_log+diabetes+cotinine)

#----Female PFHxS
PFHS_Female_Boot <- boot(data=nhanes_female, 
                         statistic=bart,
                         R=100, 
                         y = eGFR, 
                         n=1000, 
                         t = PFHS_log,
                         formula=eGFR~PFHS_log+age+race+bmi_log+diabetes+cotinine)




#------------PFNA Bootstrapping-----------
#----Overall PFNA
PFNA_Boot <- boot(data=nhanes_final, 
                  statistic=bart,
                  R=100, 
                  y = eGFR, 
                  n=1000, 
                  t = PFNA_log,
                  formula=eGFR~PFNA_log+age+race+sex+bmi_log+diabetes+cotinine)

#----Male PFOS
PFNA_Male_Boot <- boot(data=nhanes_male, 
                       statistic=bart,
                       R=100, 
                       y = eGFR, 
                       n=1000, 
                       t = PFNA_log,
                       formula=eGFR~PFNA_log+age+race+bmi_log+diabetes+cotinine)

#----Female PFOS
PFNA_Female_Boot <- boot(data=nhanes_female, 
                         statistic=bart,
                         R=100, 
                         y = eGFR, 
                         n=1000, 
                         t = PFNA_log,
                         formula=eGFR~PFNA_log+age+race+bmi_log+diabetes+cotinine)
