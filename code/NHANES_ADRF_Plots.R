#---------------------Dose-Response Curves---------------
#Code for generating plots of ADRF curves

#Load Packages
library(tidyr)
library(causaldrf)
library(BayesTree)
library(ggplot2)

#Load Data
nhanes_final <- read.csv(file = 'nhanes_final.csv')
nhanes_female <- nhanes_final[nhanes_final$gender == 2,]
nhanes_male <- nhanes_final[nhanes_final$gender == 1,]



#--------------------Total PFAS Dose-Response Curves------------
#Total PFAS Grid of Treatment Values 
PFAS_Grid <- as.vector(quantile(nhanes_final$PFAS_log, probs = seq(0, 1, 1/20)))

#Total PFAS Effect Estimates
PFAS_est <- PFAS_Boot$t0

#Total PFAS Estimates Standard Errors
PFAS_se <- PFAS_Boot$se


#Confidence Intervals
PFAS_lower <- PFAS_est - 1.96*PFAS_se
PFAS_upper <- PFAS_est + 1.96*PFAS_se

#Create dataframe with all necessary columns
PFAS_graph <- as.data.frame(cbind(PFAS_grid, PFAS_est, PFAS_upper, PFAS_lower))

#Make Graph
ggplot(PFAS_graph, aes(x=PFAS_grid, y=PFAS_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFAS_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFAS_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  #ylim(50, 120) +
  xlab("ln(Total PFAS)") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#--------------Female - Total PFAS-------------
#Total PFAS Grid of Treatment Values 
PFAS_fem_grid <- as.vector(quantile(nhanes_female$PFAS_log, probs = seq(0, 1, 1/20)))


#BART Estimates
PFAS_fem_est <- PFAS_Female_Boot$t0

#Bootstrap Standard Errors
PFAS_fem_se <- PFAS_Female_Boot$se


#Confidence Intervals
PFAS_fem_lower <- PFAS_fem_est - 1.96*PFAS_fem_se
PFAS_fem_upper <- PFAS_fem_est + 1.96*PFAS_fem_se


#Create dataframe with all necessary columns
PFAS_fem_graph <- as.data.frame(cbind(PFAS_fem_est, PFAS_fem_grid, PFAS_fem_lower, PFAS_fem_upper))


#Make Graph
ggplot(PFAS_fem_graph, aes(x=PFAS_fem_grid, y=PFAS_fem_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFAS_fem_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFAS_fem_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  #ylim(50, 120) +
  xlab("ln(Total PFAS) - Females") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




#----------------Male - Total PFAS-------------
#Total PFAS Grid of Treatment Values 
PFAS_male_grid <- as.vector(quantile(nhanes_male$PFAS_log, probs = seq(0, 1, 1/20)))

#BART Estimates
PFAS_male_est <- PFAS_Male_Boot$t0

#Bootstrap Standard Errors
PFAS_male_se <- PFAS_Male_Boot$se


#Confidence Intervals
PFAS_male_lower <- PFAS_male_est - PFAS_male_se*1.96
PFAS_male_upper <- PFAS_male_est + PFAS_male_se*1.96

#Create dataframe with all necessary columns
PFAS_male_graph <- as.data.frame(cbind(PFAS_male_est, PFAS_male_grid, PFAS_male_se, PFAS_male_lower, PFAS_male_upper))

#Make Graph
ggplot(PFAS_male_graph, aes(x=PFAS_male_grid, y=PFAS_male_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFAS_male_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFAS_male_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  ylim(50, 120) +
  xlab("ln(Total PFAS) - Males") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))





#--------------------PFOA Dose-Respone Curves------------
#Total PFAS Grid of Treatment Values 
PFOA_Grid <- as.vector(quantile(nhanes_final$PFOA_log, probs = seq(0, 1, 1/20)))

#Total PFAS Effect Estimates
PFOA_est <- PFOA_Boot$t0

#Total PFAS Estimates Standard Errors
PFOA_se <- PFOA_Boot$se


#Confidence Intervals
PFOA_lower <- PFOA_est - 1.96*PFOA_se
PFOA_upper <- PFOA_est + 1.96*PFOA_se

#Create dataframe with all necessary columns
PFOA_graph <- as.data.frame(cbind(PFOA_grid, PFOA_est, PFOA_upper, PFOA_lower))

#Make Graph
ggplot(PFOA_graph, aes(x=PFOA_grid, y=PFOA_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOA_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOA_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  #ylim(50, 120) +
  xlab("ln(PFOA)") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#--------------Female - Total PFAS-------------
#Total PFAS Grid of Treatment Values 
PFOA_fem_grid <- as.vector(quantile(nhanes_female$PFOA_log, probs = seq(0, 1, 1/20)))


#BART Estimates
PFOA_fem_est <- PFOA_Female_Boot$t0

#Bootstrap Standard Errors
PFOA_fem_se <- PFOA_Female_Boot$se


#Confidence Intervals
PFOA_fem_lower <- PFOA_fem_est - 1.96*PFOA_fem_se
PFOA_fem_upper <- PFOA_fem_est + 1.96*PFOA_fem_se


#Create dataframe with all necessary columns
PFOA_fem_graph <- as.data.frame(cbind(PFOA_fem_est, PFOA_fem_grid, PFOA_fem_lower, PFOA_fem_upper))


#Make Graph
ggplot(PFOA_fem_graph, aes(x=PFOA_fem_grid, y=PFOA_fem_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOA_fem_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOA_fem_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  #ylim(50, 120) +
  xlab("ln(PFOA) - Females") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




#--------------Male - Total PFAS-------------
#Total PFAS Grid of Treatment Values 
PFOA_male_grid <- as.vector(quantile(nhanes_male$PFOA_log, probs = seq(0, 1, 1/20)))

#BART Estimates
PFOA_male_est <- PFOA_Male_Boot$t0

#Bootstrap Standard Errors
PFOA_male_se <- PFOA_Male_Boot$se


#Confidence Intervals
PFOA_male_lower <- PFOA_male_est - PFOA_male_se*1.96
PFOA_male_upper <- PFOA_male_est + PFOA_male_se*1.96

#Create dataframe with all necessary columns
PFOA_male_graph <- as.data.frame(cbind(PFOA_male_est, PFOA_male_grid, PFOA_male_se, PFOA_male_lower, PFOA_male_upper))

#Make Graph
ggplot(PFOA_male_graph, aes(x=PFOA_male_grid, y=PFOA_male_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOA_male_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOA_male_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  ylim(50, 120) +
  xlab("ln(PFOA) - Males") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))







#--------------------PFOS Dose-Respone Curves------------
#PFOS Grid of Treatment Values 
PFOS_Grid <- as.vector(quantile(nhanes_final$PFOS_log, probs = seq(0, 1, 1/20)))

#PFOS Effect Estimates
PFOS_est <- PFOS_Boot$t0

#PFOS Estimates Standard Errors
PFOS_se <- PFOS_Boot$se


#Confidence Intervals
PFOS_lower <- PFOS_est - 1.96*PFOS_se
PFOS_upper <- PFOS_est + 1.96*PFOS_se

#Create dataframe with all necessary columns
PFOS_graph <- as.data.frame(cbind(PFOS_grid, PFOS_est, PFOS_upper, PFOS_lower))

#Make Graph
ggplot(PFOS_graph, aes(x=PFOS_grid, y=PFOS_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOS_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOS_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  #ylim(50, 120) +
  xlab("ln(PFOS)") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#--------------Female - PFOS-------------
#PFOS Grid of Treatment Values 
PFOS_fem_grid <- as.vector(quantile(nhanes_female$PFOS_log, probs = seq(0, 1, 1/20)))


#BART Estimates
PFOS_fem_est <- PFOS_Female_Boot$t0

#Bootstrap Standard Errors
PFOS_fem_se <- PFOS_Female_Boot$se


#Confidence Intervals
PFOS_fem_lower <- PFOS_fem_est - 1.96*PFOS_fem_se
PFOS_fem_upper <- PFOS_fem_est + 1.96*PFOS_fem_se


#Create dataframe with all necessary columns
PFOS_fem_graph <- as.data.frame(cbind(PFOS_fem_est, PFOS_fem_grid, PFOS_fem_lower, PFOS_fem_upper))


#Make Graph
ggplot(PFOS_fem_graph, aes(x=PFOS_fem_grid, y=PFOS_fem_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOS_fem_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOS_fem_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  #ylim(50, 120) +
  xlab("ln(PFOS) - Females") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




#--------------Male - PFOS------------
#PFOS Grid of Treatment Values 
PFOS_male_grid <- as.vector(quantile(nhanes_male$PFOS_log, probs = seq(0, 1, 1/20)))

#BART Estimates
PFOS_male_est <- PFOS_Male_Boot$t0

#Bootstrap Standard Errors
PFOS_male_se <- PFOS_Male_Boot$se


#Confidence Intervals
PFOS_male_lower <- PFOS_male_est - PFOS_male_se*1.96
PFOS_male_upper <- PFOS_male_est + PFOS_male_se*1.96

#Create dataframe with all necessary columns
PFOA_male_graph <- as.data.frame(cbind(PFOS_male_est, PFOS_male_grid, PFOS_male_se, PFOS_male_lower, PFOS_male_upper))

#Make Graph
ggplot(PFOS_male_graph, aes(x=PFOS_male_grid, y=PFOS_male_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOS_male_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOS_male_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  ylim(50, 120) +
  xlab("ln(PFOS) - Males") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))










#--------------------PFHxS Dose-Respone Curves------------
#PFHS Grid of Treatment Values 
PFHS_Grid <- as.vector(quantile(nhanes_final$PFHS_log, probs = seq(0, 1, 1/20)))

#PFHS Effect Estimates
PFHS_est <- PFHS_Boot$t0

#PFHSEstimates Standard Errors
PFHS_se <- PFHS_Boot$se


#Confidence Intervals
PFHS_lower <- PFHS_est - 1.96*PFHS_se
PFHS_upper <- PFHS_est + 1.96*PFHS_se

#Create dataframe with all necessary columns
PFHS_graph <- as.data.frame(cbind(PFHS_grid, PFHS_est, PFHS_upper, PFHS_lower))

#Make Graph
ggplot(PFHS_graph, aes(x=PFHS_grid, y=PFHS_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOS_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOS_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  #ylim(50, 120) +
  xlab("ln(PFHxS)") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#--------------Female - PFHxS-------------
#PFHS Grid of Treatment Values 
PFHS_fem_grid <- as.vector(quantile(nhanes_female$PFHS_log, probs = seq(0, 1, 1/20)))


#BART Estimates
PFHS_fem_est <- PFHS_Female_Boot$t0

#Bootstrap Standard Errors
PFHS_fem_se <- PFHS_Female_Boot$se


#Confidence Intervals
PFHS_fem_lower <- PFHS_fem_est - 1.96*PFHS_fem_se
PFHS_fem_upper <- PFHS_fem_est + 1.96*PFHS_fem_se


#Create dataframe with all necessary columns
PFHS_fem_graph <- as.data.frame(cbind(PFHS_fem_est, PFHS_fem_grid, PFHS_fem_lower, PFHS_fem_upper))


#Make Graph
ggplot(PFHS_fem_graph, aes(x=PFHS_fem_grid, y=PFHS_fem_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOS_fem_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOS_fem_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  #ylim(50, 120) +   #comment out for unscaled plots
  xlab("ln(PFHxS) - Females") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




#--------------Male - PFHxS------------
#PFHS Grid of Treatment Values 
PFHS_male_grid <- as.vector(quantile(nhanes_male$PFHS_log, probs = seq(0, 1, 1/20)))

#BART Estimates
PFHS_male_est <- PFHS_Male_Boot$t0

#Bootstrap Standard Errors
PFHS_male_se <- PFHS_Male_Boot$se


#Confidence Intervals
PFHS_male_lower <- PFHS_male_est - PFHS_male_se*1.96
PFHS_male_upper <- PFHS_male_est + PFHS_male_se*1.96

#Create dataframe with all necessary columns
PFHS_male_graph <- as.data.frame(cbind(PFHS_male_est, PFHS_male_grid, PFHS_male_se, PFHS_male_lower, PFHS_male_upper))

#Make Graph
ggplot(PFHS_male_graph, aes(x=PFHS_male_grid, y=PFHS_male_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOS_male_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOS_male_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  ylim(50, 120) +
  xlab("ln(PFHS) - Males") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))







#--------------------PFNA Dose-Respone Curves------------
#PFHS Grid of Treatment Values 
PFHS_Grid <- as.vector(quantile(nhanes_final$PFHS_log, probs = seq(0, 1, 1/20)))

#PFHS Effect Estimates
PFHS_est <- PFHS_Boot$t0

#PFHSEstimates Standard Errors
PFHS_se <- PFHS_Boot$se


#Confidence Intervals
PFHS_lower <- PFHS_est - 1.96*PFHS_se
PFHS_upper <- PFHS_est + 1.96*PFHS_se

#Create dataframe with all necessary columns
PFHS_graph <- as.data.frame(cbind(PFHS_grid, PFHS_est, PFHS_upper, PFHS_lower))

#Make Graph
ggplot(PFHS_graph, aes(x=PFHS_grid, y=PFHS_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOS_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOS_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  #ylim(50, 120) +
  xlab("ln(PFHxS)") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#--------------Female - PFHxS-------------
#PFHS Grid of Treatment Values 
PFHS_fem_grid <- as.vector(quantile(nhanes_female$PFHS_log, probs = seq(0, 1, 1/20)))


#BART Estimates
PFHS_fem_est <- PFHS_Female_Boot$t0

#Bootstrap Standard Errors
PFHS_fem_se <- PFHS_Female_Boot$se


#Confidence Intervals
PFHS_fem_lower <- PFHS_fem_est - 1.96*PFHS_fem_se
PFHS_fem_upper <- PFHS_fem_est + 1.96*PFHS_fem_se


#Create dataframe with all necessary columns
PFHS_fem_graph <- as.data.frame(cbind(PFHS_fem_est, PFHS_fem_grid, PFHS_fem_lower, PFHS_fem_upper))


#Make Graph
ggplot(PFHS_fem_graph, aes(x=PFHS_fem_grid, y=PFHS_fem_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOS_fem_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOS_fem_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  #ylim(50, 120) +   #comment out for unscaled plots
  xlab("ln(PFHxS) - Females") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




#--------------Male - PFHxS------------
#PFHS Grid of Treatment Values 
PFHS_male_grid <- as.vector(quantile(nhanes_male$PFHS_log, probs = seq(0, 1, 1/20)))

#BART Estimates
PFHS_male_est <- PFHS_Male_Boot$t0

#Bootstrap Standard Errors
PFHS_male_se <- PFHS_Male_Boot$se


#Confidence Intervals
PFHS_male_lower <- PFHS_male_est - PFHS_male_se*1.96
PFHS_male_upper <- PFHS_male_est + PFHS_male_se*1.96

#Create dataframe with all necessary columns
PFHS_male_graph <- as.data.frame(cbind(PFHS_male_est, PFHS_male_grid, PFHS_male_se, PFHS_male_lower, PFHS_male_upper))

#Make Graph
ggplot(PFHS_male_graph, aes(x=PFHS_male_grid, y=PFHS_male_est)) + geom_smooth(method="loess", se=FALSE) + 
  geom_smooth(aes(y = PFOS_male_lower), se=FALSE, color = "gray", linetype=2, size=0.5) +
  geom_smooth(aes(y = PFOS_male_upper), se=FALSE, color = "gray", linetype=2, size=0.5) +
  ylim(50, 120) +
  xlab("ln(PFHS) - Males") + 
  ylab("eGFR (mL/min/1.73"~m^2*")") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))