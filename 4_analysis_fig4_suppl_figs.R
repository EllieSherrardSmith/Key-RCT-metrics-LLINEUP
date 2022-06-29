#####################################
##
## Figure 4 and supplementary figures

## Olyset Net, Olyset Plus, PermaNet 2.0, PermaNet 3.0
cols = c("red","darkorange","blue","darkgreen")

## Linear regressions for prevalence at cross-sectional surveys after net deployment

## Overall for the trial arm
## and at the cluster level
## and a GLM for 
## difference_obs_mod ~ net_type + time_of_survey + simulation
## compare parameters for nets do sim 4 and 5

## simulations for prevalence
scenario_0 = read.csv("simulation-summary/mod_sims_scenario_0.csv",header=TRUE)[,2:6]
scenario_1 = read.csv("simulation-summary/mod_sims_scenario_1.csv",header=TRUE)[,2:6]
scenario_1_1 = read.csv("simulation-summary/mod_sims_scenario_1_1.csv",header=TRUE)[,2:106]
scenario_1_2 = read.csv("simulation-summary/mod_sims_scenario_1_2.csv",header=TRUE)[,2:106]
scenario_1_3 = read.csv("simulation-summary/mod_sims_scenario_1_3.csv",header=TRUE)[,2:106]
scenario_1_4 = read.csv("simulation-summary/mod_sims_scenario_1_4.csv",header=TRUE)[,2:6]
scenario_2 = read.csv("simulation-summary/mod_sims_scenario_2.csv",header=TRUE)[,2:106]
scenario_3 = read.csv("simulation-summary/mod_sims_scenario_3.csv",header=TRUE)[,2:106]
scenario_4 = read.csv("simulation-summary/mod_sims_scenario_4.csv",header=TRUE)[,2:106]
scenario_5 = read.csv("simulation-summary/mod_sims_scenario_5.csv",header=TRUE)[,2:106]

# # Define target, here two prevalence measures:
# target <- c(test_data$Prevalence_baseline_2_10_yrs)
# # Time points at which to match target
# target_tt <- c(6*365+test_data$days_after_jan_2017-30)
# points(target ~ target_tt,col = cols2)
# 
# match_est = numeric(104)
# for(i in 1:104){
#   match_est[i] = as.numeric(scenario_3[target_tt[i],i+2])
# }
# plot(match_est ~ target)
# abline(a=0,b=1)

## input parameters for scenario
test_data0 = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_0.csv",header=TRUE)
test_data1 = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_1.csv",header=TRUE)
test_data1_1 = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_1-1.csv",header=TRUE)
test_data1_2 = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_1-2.csv",header=TRUE)
test_data1_3 = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_1-3.csv",header=TRUE)
test_data1_4 = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_1-4.csv",header=TRUE)
test_data2 = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_2.csv",header=TRUE)
test_data3 = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_3.csv",header=TRUE)
test_data4 = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_4.csv",header=TRUE)
test_data5 = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_5.csv",header=TRUE)

## re order to align with the order for the simulated runs (check these align)
test_data0 = test_data0[order(test_data0$arm),]
test_data1 = test_data1[order(test_data1$arm),]
test_data1_1 = test_data1_1[order(test_data1_1$arm,test_data1_1$cluster),]
test_data1_2 = test_data1_2[order(test_data1_2$arm,test_data1_2$cluster),]
test_data1_3 = test_data1_3[order(test_data1_3$arm,test_data1_3$cluster),]
test_data1_4 = test_data1_4[order(test_data1_4$arm),]
test_data2 = test_data2[order(test_data2$arm,test_data2$cluster),]
test_data3 = test_data3[order(test_data3$arm,test_data3$cluster),]
test_data4 = test_data4[order(test_data4$arm,test_data4$cluster),]
test_data5 = test_data5[order(test_data5$arm,test_data5$cluster),]


## looking across clusters 
## these are just for labelling so remain same for all
dat = read.csv("simulation-summary/scenario_2.csv",header=TRUE)
dat = dat[order(dat$arm, dat$cluster),]

clusters = c(rep(dat$cluster[dat$Net_Type == "Olyset Net"],4),
             rep(dat$cluster[dat$Net_Type == "Olyset Plus"],4),
             rep(dat$cluster[dat$Net_Type == "PermaNet 2.0"],4),
             rep(dat$cluster[dat$Net_Type == "PermaNet 3.0"],4))

months = c(rep(6,length(dat$Prevalence_6m[dat$Net_Type == "Olyset Net"])),
           rep(12,length(dat$Prevalence_12m[dat$Net_Type == "Olyset Net"])),
           rep(18,length(dat$Prevalence_18m[dat$Net_Type == "Olyset Net"])),
           rep(25,length(dat$Prevalence_25m[dat$Net_Type == "Olyset Net"])),
           
           rep(6,length(dat$Prevalence_6m[dat$Net_Type == "Olyset Plus"])),
           rep(12,length(dat$Prevalence_12m[dat$Net_Type == "Olyset Plus"])),
           rep(18,length(dat$Prevalence_18m[dat$Net_Type == "Olyset Plus"])),
           rep(25,length(dat$Prevalence_25m[dat$Net_Type == "Olyset Plus"])),
           
           rep(6,length(dat$Prevalence_6m[dat$Net_Type == "PermaNet 2.0"])),
           rep(12,length(dat$Prevalence_12m[dat$Net_Type == "PermaNet 2.0"])),
           rep(18,length(dat$Prevalence_18m[dat$Net_Type == "PermaNet 2.0"])),
           rep(25,length(dat$Prevalence_25m[dat$Net_Type == "PermaNet 2.0"])),
           
           rep(6,length(dat$Prevalence_6m[dat$Net_Type == "PermaNet 3.0"])),
           rep(12,length(dat$Prevalence_12m[dat$Net_Type == "PermaNet 3.0"])),
           rep(18,length(dat$Prevalence_18m[dat$Net_Type == "PermaNet 3.0"])),
           rep(25,length(dat$Prevalence_25m[dat$Net_Type == "PermaNet 3.0"])) )

## need to use the actual prev_observed so using test_data4
prev_obs = c(test_data4$Prevalence_6m[test_data4$Net_Type == "Olyset Net"],
              test_data4$Prevalence_12m[test_data4$Net_Type == "Olyset Net"],
              test_data4$Prevalence_18m[test_data4$Net_Type == "Olyset Net"],
              test_data4$Prevalence_25m[test_data4$Net_Type == "Olyset Net"],
              
              test_data4$Prevalence_6m[test_data4$Net_Type == "Olyset Plus"],
              test_data4$Prevalence_12m[test_data4$Net_Type == "Olyset Plus"],
              test_data4$Prevalence_18m[test_data4$Net_Type == "Olyset Plus"],
              test_data4$Prevalence_25m[test_data4$Net_Type == "Olyset Plus"],
              
              test_data4$Prevalence_6m[test_data4$Net_Type == "PermaNet 2.0"],
              test_data4$Prevalence_12m[test_data4$Net_Type == "PermaNet 2.0"],
              test_data4$Prevalence_18m[test_data4$Net_Type == "PermaNet 2.0"],
              test_data4$Prevalence_25m[test_data4$Net_Type == "PermaNet 2.0"],
              
              test_data4$Prevalence_6m[test_data4$Net_Type == "PermaNet 3.0"],
              test_data4$Prevalence_12m[test_data4$Net_Type == "PermaNet 3.0"],
              test_data4$Prevalence_18m[test_data4$Net_Type == "PermaNet 3.0"],
              test_data4$Prevalence_25m[test_data4$Net_Type == "PermaNet 3.0"])

prev_obs_means = c(mean(prev_obs[1:52],na.rm=TRUE),
                   mean(prev_obs[53:128],na.rm=TRUE),
                   mean(prev_obs[129:288],na.rm=TRUE),
                   mean(prev_obs[289:416],na.rm=TRUE))


## predicted by model
## scenario 0
prev_mod0 = c(scenario_0$arm1_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365*2+30))],
             scenario_0$arm2_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365*2+30))],
             scenario_0$arm3_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365*2+30))],
             scenario_0$arm4_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365*2+30))])

prev_mod_means0 = c(mean(prev_mod0[1:52],na.rm=TRUE),
                    mean(prev_mod0[53:128],na.rm=TRUE),
                    mean(prev_mod0[129:288],na.rm=TRUE),
                    mean(prev_mod0[289:416],na.rm=TRUE))


## scenario 1
prev_mod1 = c(scenario_1$arm1_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365*2+30))],
             scenario_1$arm2_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365*2+30))],
             scenario_1$arm3_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365*2+30))],
             scenario_1$arm4_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365+365/2),
                                    c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365*2+30))])

prev_mod_means1 = c(mean(prev_mod1[1:52],na.rm=TRUE),
                    mean(prev_mod1[53:128],na.rm=TRUE),
                    mean(prev_mod1[129:288],na.rm=TRUE),
                    mean(prev_mod1[289:416],na.rm=TRUE))

## scenario 1-1
prev_mod_ol1_1 = array(dim=c(4,length(which(test_data4$Net_Type == "Olyset Net"))))
olys = c(which(test_data4$Net_Type == "Olyset Net")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "Olyset Net")))){
  prev_mod_ol1_1[1,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365/2),olys[i]]
  prev_mod_ol1_1[2,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365),olys[i]]
  prev_mod_ol1_1[3,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),olys[i]]
  prev_mod_ol1_1[4,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),olys[i]] 
}

prev_mod_olP1_1 = array(dim=c(4,length(which(test_data4$Net_Type == "Olyset Plus"))))
olyP = c(which(test_data4$Net_Type == "Olyset Plus")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "Olyset Plus")))){
  prev_mod_olP1_1[1,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365/2),olyP[i]]
  prev_mod_olP1_1[2,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365),olyP[i]]
  prev_mod_olP1_1[3,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),olyP[i]]
  prev_mod_olP1_1[4,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),olyP[i]] 
}

prev_mod_Pe21_1 = array(dim=c(4,length(which(test_data4$Net_Type == "PermaNet 2.0"))))
Per2 = c(which(test_data4$Net_Type == "PermaNet 2.0")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "PermaNet 2.0")))){
  prev_mod_Pe21_1[1,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365/2),Per2[i]]
  prev_mod_Pe21_1[2,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365),Per2[i]]
  prev_mod_Pe21_1[3,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),Per2[i]]
  prev_mod_Pe21_1[4,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),Per2[i]] 
}

prev_mod_Pe31_1 = array(dim=c(4,length(which(test_data4$Net_Type == "PermaNet 3.0"))))
Per3 = c(which(test_data4$Net_Type == "PermaNet 3.0")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "PermaNet 3.0")))){
  prev_mod_Pe31_1[1,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365/2),Per3[i]]
  prev_mod_Pe31_1[2,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365),Per3[i]]
  prev_mod_Pe31_1[3,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),Per3[i]]
  prev_mod_Pe31_1[4,i] = scenario_1_1[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),Per3[i]] 
}

prev_mod1_1 = c(prev_mod_ol1_1[1,],prev_mod_ol1_1[2,],prev_mod_ol1_1[3,],prev_mod_ol1_1[4,],
             prev_mod_olP1_1[1,],prev_mod_olP1_1[2,],prev_mod_olP1_1[3,],prev_mod_olP1_1[4,],
             prev_mod_Pe21_1[1,],prev_mod_Pe21_1[2,],prev_mod_Pe21_1[3,],prev_mod_Pe21_1[4,],
             prev_mod_Pe31_1[1,],prev_mod_Pe31_1[2,],prev_mod_Pe31_1[3,],prev_mod_Pe31_1[4,])

prev_mod_means1_1 = c(mean(prev_mod1_1[1:52],na.rm=TRUE),
                   mean(prev_mod1_1[53:128],na.rm=TRUE),
                   mean(prev_mod1_1[129:288],na.rm=TRUE),
                   mean(prev_mod1_1[289:416],na.rm=TRUE))

## scenario 1-2

prev_mod_ol1_2 = array(dim=c(4,length(which(test_data4$Net_Type == "Olyset Net"))))
olys = c(which(test_data4$Net_Type == "Olyset Net")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "Olyset Net")))){
  prev_mod_ol1_2[1,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365/2),olys[i]]
  prev_mod_ol1_2[2,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365),olys[i]]
  prev_mod_ol1_2[3,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),olys[i]]
  prev_mod_ol1_2[4,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),olys[i]] 
}

prev_mod_olP1_2 = array(dim=c(4,length(which(test_data4$Net_Type == "Olyset Plus"))))
olyP = c(which(test_data4$Net_Type == "Olyset Plus")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "Olyset Plus")))){
  prev_mod_olP1_2[1,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365/2),olyP[i]]
  prev_mod_olP1_2[2,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365),olyP[i]]
  prev_mod_olP1_2[3,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),olyP[i]]
  prev_mod_olP1_2[4,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),olyP[i]] 
}

prev_mod_Pe21_2 = array(dim=c(4,length(which(test_data4$Net_Type == "PermaNet 2.0"))))
Per2 = c(which(test_data4$Net_Type == "PermaNet 2.0")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "PermaNet 2.0")))){
  prev_mod_Pe21_2[1,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365/2),Per2[i]]
  prev_mod_Pe21_2[2,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365),Per2[i]]
  prev_mod_Pe21_2[3,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),Per2[i]]
  prev_mod_Pe21_2[4,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),Per2[i]] 
}

prev_mod_Pe31_2 = array(dim=c(4,length(which(test_data4$Net_Type == "PermaNet 3.0"))))
Per3 = c(which(test_data4$Net_Type == "PermaNet 3.0")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "PermaNet 3.0")))){
  prev_mod_Pe31_2[1,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365/2),Per3[i]]
  prev_mod_Pe31_2[2,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365),Per3[i]]
  prev_mod_Pe31_2[3,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),Per3[i]]
  prev_mod_Pe31_2[4,i] = scenario_1_2[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),Per3[i]] 
}

prev_mod1_2 = c(prev_mod_ol1_2[1,],prev_mod_ol1_2[2,],prev_mod_ol1_2[3,],prev_mod_ol1_2[4,],
                prev_mod_olP1_2[1,],prev_mod_olP1_2[2,],prev_mod_olP1_2[3,],prev_mod_olP1_2[4,],
                prev_mod_Pe21_2[1,],prev_mod_Pe21_2[2,],prev_mod_Pe21_2[3,],prev_mod_Pe21_2[4,],
                prev_mod_Pe31_2[1,],prev_mod_Pe31_2[2,],prev_mod_Pe31_2[3,],prev_mod_Pe31_2[4,])

prev_mod_means1_2 = c(mean(prev_mod1_2[1:52],na.rm=TRUE),
                      mean(prev_mod1_2[53:128],na.rm=TRUE),
                      mean(prev_mod1_2[129:288],na.rm=TRUE),
                      mean(prev_mod1_2[289:416],na.rm=TRUE))

## scenario 1-3

prev_mod_ol1_3 = array(dim=c(4,length(which(test_data4$Net_Type == "Olyset Net"))))
olys = c(which(test_data4$Net_Type == "Olyset Net")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "Olyset Net")))){
  prev_mod_ol1_3[1,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365/2),olys[i]]
  prev_mod_ol1_3[2,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365),olys[i]]
  prev_mod_ol1_3[3,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),olys[i]]
  prev_mod_ol1_3[4,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),olys[i]] 
}

prev_mod_olP1_3 = array(dim=c(4,length(which(test_data4$Net_Type == "Olyset Plus"))))
olyP = c(which(test_data4$Net_Type == "Olyset Plus")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "Olyset Plus")))){
  prev_mod_olP1_3[1,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365/2),olyP[i]]
  prev_mod_olP1_3[2,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365),olyP[i]]
  prev_mod_olP1_3[3,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),olyP[i]]
  prev_mod_olP1_3[4,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),olyP[i]] 
}

prev_mod_Pe21_3 = array(dim=c(4,length(which(test_data4$Net_Type == "PermaNet 2.0"))))
Per2 = c(which(test_data4$Net_Type == "PermaNet 2.0")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "PermaNet 2.0")))){
  prev_mod_Pe21_3[1,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365/2),Per2[i]]
  prev_mod_Pe21_3[2,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365),Per2[i]]
  prev_mod_Pe21_3[3,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),Per2[i]]
  prev_mod_Pe21_3[4,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),Per2[i]] 
}

prev_mod_Pe31_3 = array(dim=c(4,length(which(test_data4$Net_Type == "PermaNet 3.0"))))
Per3 = c(which(test_data4$Net_Type == "PermaNet 3.0")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "PermaNet 3.0")))){
  prev_mod_Pe31_3[1,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365/2),Per3[i]]
  prev_mod_Pe31_3[2,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365),Per3[i]]
  prev_mod_Pe31_3[3,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),Per3[i]]
  prev_mod_Pe31_3[4,i] = scenario_1_3[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),Per3[i]] 
}

prev_mod1_3 = c(prev_mod_ol1_3[1,],prev_mod_ol1_3[2,],prev_mod_ol1_3[3,],prev_mod_ol1_3[4,],
                prev_mod_olP1_3[1,],prev_mod_olP1_3[2,],prev_mod_olP1_3[3,],prev_mod_olP1_3[4,],
                prev_mod_Pe21_3[1,],prev_mod_Pe21_3[2,],prev_mod_Pe21_3[3,],prev_mod_Pe21_3[4,],
                prev_mod_Pe31_3[1,],prev_mod_Pe31_3[2,],prev_mod_Pe31_3[3,],prev_mod_Pe31_3[4,])

prev_mod_means1_3 = c(mean(prev_mod1_3[1:52],na.rm=TRUE),
                      mean(prev_mod1_3[53:128],na.rm=TRUE),
                      mean(prev_mod1_3[129:288],na.rm=TRUE),
                      mean(prev_mod1_3[289:416],na.rm=TRUE))


## scenario 1-4

prev_mod1_4 = c(scenario_1_4$arm1_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365/2),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365+365/2),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Net"]+365*2+30))],
              scenario_1_4$arm2_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365/2),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365+365/2),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "Olyset Plus"]+365*2+30))],
              scenario_1_4$arm3_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365/2),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365+365/2),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 2.0"]+365*2+30))],
              scenario_1_4$arm4_prev[c(c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365/2),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365+365/2),
                                     c(6*365+test_data4$days_after_jan_2017[test_data4$Net_Type == "PermaNet 3.0"]+365*2+30))])

prev_mod_means1_4 = c(mean(prev_mod1_4[1:52],na.rm=TRUE),
                    mean(prev_mod1_4[53:128],na.rm=TRUE),
                    mean(prev_mod1_4[129:288],na.rm=TRUE),
                    mean(prev_mod1_4[289:416],na.rm=TRUE))


## scenario 2

prev_mod_ol2 = array(dim=c(4,length(which(test_data4$Net_Type == "Olyset Net"))))
olys = c(which(test_data4$Net_Type == "Olyset Net")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "Olyset Net")))){
  prev_mod_ol2[1,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365/2),olys[i]]
  prev_mod_ol2[2,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365),olys[i]]
  prev_mod_ol2[3,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),olys[i]]
  prev_mod_ol2[4,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),olys[i]] 
}

prev_mod_olP2 = array(dim=c(4,length(which(test_data4$Net_Type == "Olyset Plus"))))
olyP = c(which(test_data4$Net_Type == "Olyset Plus")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "Olyset Plus")))){
  prev_mod_olP2[1,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365/2),olyP[i]]
  prev_mod_olP2[2,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365),olyP[i]]
  prev_mod_olP2[3,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),olyP[i]]
  prev_mod_olP2[4,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),olyP[i]] 
}

prev_mod_Pe22 = array(dim=c(4,length(which(test_data4$Net_Type == "PermaNet 2.0"))))
Per2 = c(which(test_data4$Net_Type == "PermaNet 2.0")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "PermaNet 2.0")))){
  prev_mod_Pe22[1,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365/2),Per2[i]]
  prev_mod_Pe22[2,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365),Per2[i]]
  prev_mod_Pe22[3,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),Per2[i]]
  prev_mod_Pe22[4,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),Per2[i]] 
}

prev_mod_Pe32 = array(dim=c(4,length(which(test_data4$Net_Type == "PermaNet 3.0"))))
Per3 = c(which(test_data4$Net_Type == "PermaNet 3.0")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "PermaNet 3.0")))){
  prev_mod_Pe32[1,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365/2),Per3[i]]
  prev_mod_Pe32[2,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365),Per3[i]]
  prev_mod_Pe32[3,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),Per3[i]]
  prev_mod_Pe32[4,i] = scenario_2[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),Per3[i]] 
}

prev_mod2 = c(prev_mod_ol2[1,],prev_mod_ol2[2,],prev_mod_ol2[3,],prev_mod_ol2[4,],
                prev_mod_olP2[1,],prev_mod_olP2[2,],prev_mod_olP2[3,],prev_mod_olP2[4,],
                prev_mod_Pe22[1,],prev_mod_Pe22[2,],prev_mod_Pe22[3,],prev_mod_Pe22[4,],
                prev_mod_Pe32[1,],prev_mod_Pe32[2,],prev_mod_Pe32[3,],prev_mod_Pe32[4,])

prev_mod_means2 = c(mean(prev_mod2[1:52],na.rm=TRUE),
                      mean(prev_mod2[53:128],na.rm=TRUE),
                      mean(prev_mod2[129:288],na.rm=TRUE),
                      mean(prev_mod2[289:416],na.rm=TRUE))


## scenario 3

prev_mod_ol3 = array(dim=c(4,length(which(test_data3$Net_Type == "Olyset Net"))))
olys = c(which(test_data3$Net_Type == "Olyset Net")+1)
for(i in 1:length(c(which(test_data3$Net_Type == "Olyset Net")))){
  prev_mod_ol3[1,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365/2),olys[i]]
  prev_mod_ol3[2,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365),olys[i]]
  prev_mod_ol3[3,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365+365/2),olys[i]]
  prev_mod_ol3[4,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365*2+30),olys[i]] 
}

prev_mod_olP3 = array(dim=c(4,length(which(test_data3$Net_Type == "Olyset Plus"))))
olyP = c(which(test_data3$Net_Type == "Olyset Plus")+1)
for(i in 1:length(c(which(test_data3$Net_Type == "Olyset Plus")))){
  prev_mod_olP3[1,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365/2),olyP[i]]
  prev_mod_olP3[2,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365),olyP[i]]
  prev_mod_olP3[3,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365+365/2),olyP[i]]
  prev_mod_olP3[4,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365*2+30),olyP[i]] 
}

prev_mod_Pe23 = array(dim=c(4,length(which(test_data3$Net_Type == "PermaNet 2.0"))))
Per2 = c(which(test_data3$Net_Type == "PermaNet 2.0")+1)
for(i in 1:length(c(which(test_data3$Net_Type == "PermaNet 2.0")))){
  prev_mod_Pe23[1,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365/2),Per2[i]]
  prev_mod_Pe23[2,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365),Per2[i]]
  prev_mod_Pe23[3,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365+365/2),Per2[i]]
  prev_mod_Pe23[4,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365*2+30),Per2[i]] 
}

prev_mod_Pe33 = array(dim=c(4,length(which(test_data3$Net_Type == "PermaNet 3.0"))))
Per3 = c(which(test_data3$Net_Type == "PermaNet 3.0")+1)
for(i in 1:length(c(which(test_data3$Net_Type == "PermaNet 3.0")))){
  prev_mod_Pe33[1,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365/2),Per3[i]]
  prev_mod_Pe33[2,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365),Per3[i]]
  prev_mod_Pe33[3,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365+365/2),Per3[i]]
  prev_mod_Pe33[4,i] = scenario_3[c(6*365+test_data3$days_after_jan_2017[i]+365*2+30),Per3[i]] 
}

prev_mod3 = c(prev_mod_ol3[1,],prev_mod_ol3[2,],prev_mod_ol3[3,],prev_mod_ol3[4,],
                prev_mod_olP3[1,],prev_mod_olP3[2,],prev_mod_olP3[3,],prev_mod_olP3[4,],
                prev_mod_Pe23[1,],prev_mod_Pe23[2,],prev_mod_Pe23[3,],prev_mod_Pe23[4,],
                prev_mod_Pe33[1,],prev_mod_Pe33[2,],prev_mod_Pe33[3,],prev_mod_Pe33[4,])

prev_mod_means3 = c(mean(prev_mod3[1:52],na.rm=TRUE),
                      mean(prev_mod3[53:128],na.rm=TRUE),
                      mean(prev_mod3[129:288],na.rm=TRUE),
                      mean(prev_mod3[289:416],na.rm=TRUE))

## scenario 4

prev_mod_ol4 = array(dim=c(4,length(which(test_data4$Net_Type == "Olyset Net"))))
olys = c(which(test_data4$Net_Type == "Olyset Net")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "Olyset Net")))){
  prev_mod_ol4[1,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365/2),olys[i]]
  prev_mod_ol4[2,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365),olys[i]]
  prev_mod_ol4[3,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),olys[i]]
  prev_mod_ol4[4,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),olys[i]] 
}

prev_mod_olP4 = array(dim=c(4,length(which(test_data4$Net_Type == "Olyset Plus"))))
olyP = c(which(test_data4$Net_Type == "Olyset Plus")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "Olyset Plus")))){
  prev_mod_olP4[1,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365/2),olyP[i]]
  prev_mod_olP4[2,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365),olyP[i]]
  prev_mod_olP4[3,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),olyP[i]]
  prev_mod_olP4[4,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),olyP[i]] 
}

prev_mod_Pe24 = array(dim=c(4,length(which(test_data4$Net_Type == "PermaNet 2.0"))))
Per2 = c(which(test_data4$Net_Type == "PermaNet 2.0")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "PermaNet 2.0")))){
  prev_mod_Pe24[1,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365/2),Per2[i]]
  prev_mod_Pe24[2,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365),Per2[i]]
  prev_mod_Pe24[3,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),Per2[i]]
  prev_mod_Pe24[4,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),Per2[i]] 
}

prev_mod_Pe34 = array(dim=c(4,length(which(test_data4$Net_Type == "PermaNet 3.0"))))
Per3 = c(which(test_data4$Net_Type == "PermaNet 3.0")+1)
for(i in 1:length(c(which(test_data4$Net_Type == "PermaNet 3.0")))){
  prev_mod_Pe34[1,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365/2),Per3[i]]
  prev_mod_Pe34[2,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365),Per3[i]]
  prev_mod_Pe34[3,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365+365/2),Per3[i]]
  prev_mod_Pe34[4,i] = scenario_4[c(6*365+test_data4$days_after_jan_2017[i]+365*2+30),Per3[i]] 
}

prev_mod4 = c(prev_mod_ol4[1,],prev_mod_ol4[2,],prev_mod_ol4[3,],prev_mod_ol4[4,],
                prev_mod_olP4[1,],prev_mod_olP4[2,],prev_mod_olP4[3,],prev_mod_olP4[4,],
                prev_mod_Pe24[1,],prev_mod_Pe24[2,],prev_mod_Pe24[3,],prev_mod_Pe24[4,],
                prev_mod_Pe34[1,],prev_mod_Pe34[2,],prev_mod_Pe34[3,],prev_mod_Pe34[4,])

prev_mod_means4 = c(mean(prev_mod4[1:52],na.rm=TRUE),
                      mean(prev_mod4[53:128],na.rm=TRUE),
                      mean(prev_mod4[129:288],na.rm=TRUE),
                      mean(prev_mod4[289:416],na.rm=TRUE))

## scenario 5

prev_mod_ol5 = array(dim=c(4,length(which(test_data5$Net_Type == "Olyset Net"))))
olys = c(which(test_data5$Net_Type == "Olyset Net")+1)
for(i in 1:length(c(which(test_data5$Net_Type == "Olyset Net")))){
  prev_mod_ol5[1,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365/2),olys[i]]
  prev_mod_ol5[2,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365),olys[i]]
  prev_mod_ol5[3,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365+365/2),olys[i]]
  prev_mod_ol5[4,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365*2+30),olys[i]] 
}

prev_mod_olP5 = array(dim=c(4,length(which(test_data5$Net_Type == "Olyset Plus"))))
olyP = c(which(test_data5$Net_Type == "Olyset Plus")+1)
for(i in 1:length(c(which(test_data5$Net_Type == "Olyset Plus")))){
  prev_mod_olP5[1,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365/2),olyP[i]]
  prev_mod_olP5[2,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365),olyP[i]]
  prev_mod_olP5[3,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365+365/2),olyP[i]]
  prev_mod_olP5[4,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365*2+30),olyP[i]] 
}

prev_mod_Pe25 = array(dim=c(4,length(which(test_data5$Net_Type == "PermaNet 2.0"))))
Per2 = c(which(test_data5$Net_Type == "PermaNet 2.0")+1)
for(i in 1:length(c(which(test_data5$Net_Type == "PermaNet 2.0")))){
  prev_mod_Pe25[1,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365/2),Per2[i]]
  prev_mod_Pe25[2,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365),Per2[i]]
  prev_mod_Pe25[3,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365+365/2),Per2[i]]
  prev_mod_Pe25[4,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365*2+30),Per2[i]] 
}

prev_mod_Pe35 = array(dim=c(4,length(which(test_data5$Net_Type == "PermaNet 3.0"))))
Per3 = c(which(test_data5$Net_Type == "PermaNet 3.0")+1)
for(i in 1:length(c(which(test_data5$Net_Type == "PermaNet 3.0")))){
  prev_mod_Pe35[1,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365/2),Per3[i]]
  prev_mod_Pe35[2,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365),Per3[i]]
  prev_mod_Pe35[3,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365+365/2),Per3[i]]
  prev_mod_Pe35[4,i] = scenario_5[c(6*365+test_data5$days_after_jan_2017[i]+365*2+30),Per3[i]] 
}

prev_mod5 = c(prev_mod_ol5[1,],prev_mod_ol5[2,],prev_mod_ol5[3,],prev_mod_ol5[4,],
                prev_mod_olP5[1,],prev_mod_olP5[2,],prev_mod_olP5[3,],prev_mod_olP5[4,],
                prev_mod_Pe25[1,],prev_mod_Pe25[2,],prev_mod_Pe25[3,],prev_mod_Pe25[4,],
                prev_mod_Pe35[1,],prev_mod_Pe35[2,],prev_mod_Pe35[3,],prev_mod_Pe35[4,])

prev_mod_means5 = c(mean(prev_mod5[1:52],na.rm=TRUE),
                      mean(prev_mod5[53:128],na.rm=TRUE),
                      mean(prev_mod5[129:288],na.rm=TRUE),
                      mean(prev_mod5[289:416],na.rm=TRUE))

#################################
##
## Compare the means
## but ... e.g. 0 and 1 arms 1 and 3 should be the same
## they differ - which is due to the stochastic nature of the simulation
## and this is contributing to the conclusion... 

par(mfrow=c(1,1))
## use A for data without matched deployment months
plot(prev_obs ~ prev_mod0, ylim=c(0,0.25),xlim=c(0,0.25),
     main = "",
     ylab="Trial reported prevalence (%)",yaxt="n",
     xlab = "Model simulated prevalence (%)",xaxt="n",
     col="white",pch=19)
axis(1, at = seq(0,0.25,0.05), labels = seq(0,25,5))
axis(2, las = 2, at = seq(0,0.25,0.05), labels = seq(0,25,5))
abline(a=0,b=1,lty=2)

points(prev_obs_means ~ prev_mod_means0,col=cols,pch=1,cex=1.2)
points(prev_obs_means ~ prev_mod_means1,col=adegenet::transp(cols,0.4),pch=19,cex=1.2)
points(prev_obs_means ~ prev_mod_means1_1,col=adegenet::transp(cols,0.4),pch=15,cex=1.2)
points(prev_obs_means ~ prev_mod_means1_2,col=adegenet::transp(cols,0.4),pch=17,cex=1.2)
points(prev_obs_means ~ prev_mod_means1_3,col=adegenet::transp(cols,0.4),pch=8,cex=1.2)
points(prev_obs_means ~ prev_mod_means1_4,col=adegenet::transp(cols,0.4),pch=2,cex=1.2)
points(prev_obs_means ~ prev_mod_means2,col=adegenet::transp(cols,0.4),pch=5,cex=1.2)
points(prev_obs_means ~ prev_mod_means3,col=adegenet::transp(cols,0.4),pch=6,cex=1.2)
points(prev_obs_means ~ prev_mod_means4,col=adegenet::transp(cols,0.4),pch=3,cex=1.2)
points(prev_obs_means ~ prev_mod_means5,col=adegenet::transp(cols,0.4),pch=16,cex=1.2)

legend("topleft",title="Simulation",
       legend=c("0","1","1-1","1-2","1-3","1-4","2","3","4","5"),
       pch = c(1,19,15,17,8,2,5,6,3,16),
       ncol=3,bty="n")
  


## Draw the association for whichever is the best predictor

## absolute difference in prediction and observation
ad_0 = prev_obs - prev_mod0
ad_1 = prev_obs - prev_mod1
ad_1_1 = prev_obs - prev_mod1_1
ad_1_2 = prev_obs - prev_mod1_2
ad_1_3 = prev_obs - prev_mod1_3
ad_1_4 = prev_obs - prev_mod1_4
ad_2 = prev_obs - prev_mod2
ad_3 = prev_obs - prev_mod3
ad_4 = prev_obs - prev_mod4
ad_5 = prev_obs - prev_mod5

range(ad_0,na.rm=T)
range(ad_1,na.rm=T)
range(ad_1_1,na.rm=T)
range(ad_1_2,na.rm=T)
range(ad_1_3,na.rm=T)
range(ad_1_4,na.rm=T)
range(ad_2,na.rm=T)
range(ad_3,na.rm=T)
range(ad_4,na.rm=T)
range(ad_5,na.rm=T)

quantile(ad_0,c(0.25,0.75),na.rm=T)
quantile(ad_1,c(0.25,0.75),na.rm=T)
quantile(ad_1_1,c(0.25,0.75),na.rm=T)
quantile(ad_1_2,c(0.25,0.75),na.rm=T)
quantile(ad_1_3,c(0.25,0.75),na.rm=T)
quantile(ad_1_4,c(0.25,0.75),na.rm=T)
quantile(ad_2,c(0.25,0.75),na.rm=T)
quantile(ad_3,c(0.25,0.75),na.rm=T)
quantile(ad_4,c(0.25,0.75),na.rm=T)
quantile(ad_5,c(0.25,0.75),na.rm=T)

ad_0_means = prev_obs_means - prev_mod_means0
ad_1_means = prev_obs_means - prev_mod_means1
ad_1_1_means = prev_obs_means - prev_mod_means1_1
ad_1_2_means = prev_obs_means - prev_mod_means1_2
ad_1_3_means = prev_obs_means - prev_mod_means1_3
ad_1_4_means = prev_obs_means - prev_mod_means1_4
ad_2_means = prev_obs_means - prev_mod_means2
ad_3_means = prev_obs_means - prev_mod_means3
ad_4_means = prev_obs_means - prev_mod_means4
ad_5_means = prev_obs_means - prev_mod_means5

boxplot(ad_0,ad_1,ad_1_1,ad_1_2,ad_1_3,ad_1_4,
        ad_2,ad_3,ad_4,ad_5,
        xlab="Simulations",
        xaxt="n",yaxt="n",ylab="Absolute difference in simulation and empirical data")
axis(2, las=2, at=seq(-0.4,0.4,0.1),labels=seq(-0.4,0.4,0.1))
axis(1, at=1:10,labels=c("0","1","1-1","1-2","1-3","1-4","2","3","4","5"))
abline(h=0,lty=2)
points(ad_0[1:52] ~ sample(x=seq(0.65,0.8,length=20),replace = T,size = 52),pch=rep(c(1,8,19,17),each=13),col=adegenet::transp("red",0.4))
points(ad_0[53:128] ~ sample(x=seq(0.85,1,length=20),replace = T,size = 76),pch=rep(c(1,8,19,17),each=19),col=adegenet::transp("darkorange",0.4))
points(ad_0[129:288] ~ sample(x=seq(1.05,1.2,length=20),replace = T,size = 160),pch=rep(c(1,8,19,17),each=40),col=adegenet::transp("blue",0.4))
points(ad_0[289:416] ~ sample(x=seq(1.25,1.4,length=20),replace = T,size = 128),pch=rep(c(1,8,19,17),each=32),col=adegenet::transp("darkgreen",0.4))

points(ad_1[1:52] ~ sample(x=seq(1.65,1.8,length=20),replace = T,size = 52),pch=rep(c(1,8,19,17),each=13),col=adegenet::transp("red",0.4))
points(ad_1[53:128] ~ sample(x=seq(1.85,2,length=20),replace = T,size = 76),pch=rep(c(1,8,19,17),each=19),col=adegenet::transp("darkorange",0.4))
points(ad_1[129:288] ~ sample(x=seq(2.05,2.2,length=20),replace = T,size = 160),pch=rep(c(1,8,19,17),each=40),col=adegenet::transp("blue",0.4))
points(ad_1[289:416] ~ sample(x=seq(2.25,2.4,length=20),replace = T,size = 128),pch=rep(c(1,8,19,17),each=32),col=adegenet::transp("darkgreen",0.4))

points(ad_1_1[1:52] ~ sample(x=seq(2.65,2.8,length=20),replace = T,size = 52),pch=rep(c(1,8,19,17),each=13),col=adegenet::transp("red",0.4))
points(ad_1_1[53:128] ~ sample(x=seq(2.85,3,length=20),replace = T,size = 76),pch=rep(c(1,8,19,17),each=19),col=adegenet::transp("darkorange",0.4))
points(ad_1_1[129:288] ~ sample(x=seq(3.05,3.2,length=20),replace = T,size = 160),pch=rep(c(1,8,19,17),each=40),col=adegenet::transp("blue",0.4))
points(ad_1_1[289:416] ~ sample(x=seq(3.25,3.4,length=20),replace = T,size = 128),pch=rep(c(1,8,19,17),each=32),col=adegenet::transp("darkgreen",0.4))

points(ad_1_2[1:52] ~ sample(x=seq(3.65,3.8,length=20),replace = T,size = 52),pch=rep(c(1,8,19,17),each=13),col=adegenet::transp("red",0.4))
points(ad_1_2[53:128] ~ sample(x=seq(3.85,4,length=20),replace = T,size = 76),pch=rep(c(1,8,19,17),each=19),col=adegenet::transp("darkorange",0.4))
points(ad_1_2[129:288] ~ sample(x=seq(4.05,4.2,length=20),replace = T,size = 160),pch=rep(c(1,8,19,17),each=40),col=adegenet::transp("blue",0.4))
points(ad_1_2[289:416] ~ sample(x=seq(4.25,4.4,length=20),replace = T,size = 128),pch=rep(c(1,8,19,17),each=32),col=adegenet::transp("darkgreen",0.4))

points(ad_1_3[1:52] ~ sample(x=seq(4.65,4.8,length=20),replace = T,size = 52),pch=rep(c(1,8,19,17),each=13),col=adegenet::transp("red",0.4))
points(ad_1_3[53:128] ~ sample(x=seq(4.85,5,length=20),replace = T,size = 76),pch=rep(c(1,8,19,17),each=19),col=adegenet::transp("darkorange",0.4))
points(ad_1_3[129:288] ~ sample(x=seq(5.05,5.2,length=20),replace = T,size = 160),pch=rep(c(1,8,19,17),each=40),col=adegenet::transp("blue",0.4))
points(ad_1_3[289:416] ~ sample(x=seq(5.25,5.4,length=20),replace = T,size = 128),pch=rep(c(1,8,19,17),each=32),col=adegenet::transp("darkgreen",0.4))

points(ad_1_4[1:52] ~ sample(x=seq(5.65,5.8,length=20),replace = T,size = 52),pch=rep(c(1,8,19,17),each=13),col=adegenet::transp("red",0.4))
points(ad_1_4[53:128] ~ sample(x=seq(5.85,6,length=20),replace = T,size = 76),pch=rep(c(1,8,19,17),each=19),col=adegenet::transp("darkorange",0.4))
points(ad_1_4[129:288] ~ sample(x=seq(6.05,6.2,length=20),replace = T,size = 160),pch=rep(c(1,8,19,17),each=40),col=adegenet::transp("blue",0.4))
points(ad_1_4[289:416] ~ sample(x=seq(6.25,6.4,length=20),replace = T,size = 128),pch=rep(c(1,8,19,17),each=32),col=adegenet::transp("darkgreen",0.4))

points(ad_2[1:52] ~ sample(x=seq(6.65,6.8,length=20),replace = T,size = 52),pch=rep(c(1,8,19,17),each=13),col=adegenet::transp("red",0.4))
points(ad_2[53:128] ~ sample(x=seq(6.85,7,length=20),replace = T,size = 76),pch=rep(c(1,8,19,17),each=19),col=adegenet::transp("darkorange",0.4))
points(ad_2[129:288] ~ sample(x=seq(7.05,7.2,length=20),replace = T,size = 160),pch=rep(c(1,8,19,17),each=40),col=adegenet::transp("blue",0.4))
points(ad_2[289:416] ~ sample(x=seq(7.25,7.4,length=20),replace = T,size = 128),pch=rep(c(1,8,19,17),each=32),col=adegenet::transp("darkgreen",0.4))

points(ad_3[1:52] ~ sample(x=seq(7.65,7.8,length=20),replace = T,size = 52),pch=rep(c(1,8,19,17),each=13),col=adegenet::transp("red",0.4))
points(ad_3[53:128] ~ sample(x=seq(7.85,8,length=20),replace = T,size = 76),pch=rep(c(1,8,19,17),each=19),col=adegenet::transp("darkorange",0.4))
points(ad_3[129:288] ~ sample(x=seq(8.05,8.2,length=20),replace = T,size = 160),pch=rep(c(1,8,19,17),each=40),col=adegenet::transp("blue",0.4))
points(ad_3[289:416] ~ sample(x=seq(8.25,8.4,length=20),replace = T,size = 128),pch=rep(c(1,8,19,17),each=32),col=adegenet::transp("darkgreen",0.4))

points(ad_4[1:52] ~ sample(x=seq(8.65,8.8,length=20),replace = T,size = 52),pch=rep(c(1,8,19,17),each=13),col=adegenet::transp("red",0.4))
points(ad_4[53:128] ~ sample(x=seq(8.85,9,length=20),replace = T,size = 76),pch=rep(c(1,8,19,17),each=19),col=adegenet::transp("darkorange",0.4))
points(ad_4[129:288] ~ sample(x=seq(9.05,9.2,length=20),replace = T,size = 160),pch=rep(c(1,8,19,17),each=40),col=adegenet::transp("blue",0.4))
points(ad_4[289:416] ~ sample(x=seq(9.25,9.4,length=20),replace = T,size = 128),pch=rep(c(1,8,19,17),each=32),col=adegenet::transp("darkgreen",0.4))

points(ad_5[1:52] ~ sample(x=seq(9.65,9.8,length=20),replace = T,size = 52),pch=rep(c(1,8,19,17),each=13),col=adegenet::transp("red",0.4))
points(ad_5[53:128] ~ sample(x=seq(9.85,10,length=20),replace = T,size = 76),pch=rep(c(1,8,19,17),each=19),col=adegenet::transp("darkorange",0.4))
points(ad_5[129:288] ~ sample(x=seq(10.05,10.2,length=20),replace = T,size = 160),pch=rep(c(1,8,19,17),each=40),col=adegenet::transp("blue",0.4))
points(ad_5[289:416] ~ sample(x=seq(10.25,10.4,length=20),replace = T,size = 128),pch=rep(c(1,8,19,17),each=32),col=adegenet::transp("darkgreen",0.4))

##################################################
##
## Figure 4
##
##################################################

## a) e.g. of a cluster (from sim 4)
##    shapes for cross-sectional surveys
##    colour for net type
## b) prev ~ prev for sim 4 all clusters
##    shapes for cross-sectional surveys
##    colour for net type
##    Overlay overall estimates for trial arm level
## c) absolute difference - jitter points and show 95% uncertainty and median estimate

plot_4a = function(scen,
                   scen3,
                   scen4,
                   scen5,
                   td,tdB,net,a,colnet){
  plot(scen[,2] ~ scen$timestep,type="l",
       ylab = "Prevalence in children 6m to 10yrs (%)",yaxt="n",ylim=c(0,0.8),
       xlab = "Time in days",col="white",xlim=c(2000,3500)
  )
  axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
  abline(v=timestep,lty=2)
  
  td$target = c(td$Prevalence_baseline_2_10_yrs)
  td$target_tt = c(6*365+td$days_after_jan_2017-30)
  tdB$target = c(tdB$Prevalence_baseline_2_10_yrs)
  tdB$target_tt = c(6*365+tdB$days_after_jan_2017-30)
  
  select_a_run = c(which(td$Net_Type == net)+1)
  select_a_point = td$target[td$Net_Type == net]
  select_time = td$target_tt[td$Net_Type == net]
  
  select_a_runB = c(which(tdB$Net_Type == net)+1)
  select_a_pointB = tdB$target[tdB$Net_Type == net]
  select_timeB = tdB$target_tt[tdB$Net_Type == net]
  
  time6m = td$Prevalence_6m[td$Net_Type == net]
  time12m = td$Prevalence_12m[td$Net_Type == net]
  time18m = td$Prevalence_18m[td$Net_Type == net]
  time25m = td$Prevalence_25m[td$Net_Type == net]
  
  t6=c(6*365+tdB$days_after_jan_2017[tdB$Net_Type == net]+365/2)
  t12=c(6*365+tdB$days_after_jan_2017[tdB$Net_Type == net]+365)
  t18=c(6*365+tdB$days_after_jan_2017[tdB$Net_Type == net]+365+365/2)
  t25=c(6*365+tdB$days_after_jan_2017[tdB$Net_Type == net]+365+365+30)

  lines(scen[,select_a_run[a]] ~ scen$timestep,col=colnet)
  lines(scen3[,select_a_run[a]] ~ scen3$timestep,col=adegenet::transp(colnet,0.7),lty=2)
  lines(scen4[,select_a_run[a]] ~ scen4$timestep,col=adegenet::transp(colnet,0.4),lty=4)
  lines(scen5[,select_a_run[a]] ~ scen5$timestep,col=adegenet::transp(colnet,0.2),lty=6)
  
  points(select_a_point[a] ~ select_time[a],
         col = colnet,pch=19,cex=1.4)
  points(select_a_pointB[a] ~ select_timeB[a],
         col = "black",pch=8,cex=1.4)
  
  points(time6m[a]~
           t6[a],
         col=adegenet::transp(colnet,0.2),pch=19)
  points(time12m[a]~
           t12[a],
         col=adegenet::transp(colnet,0.2),pch=19)
  points(time18m[a]~
           t18[a],
         col=adegenet::transp(colnet,0.2),pch=19)
  points(time25m[a]~
           t25[a],
         col=adegenet::transp(colnet,0.2),pch=19)
  
}
# 4
par(mfrow=c(2,2))
plot_4a(scen=scenario_2,
        scen3=scenario_3,
        scen4=scenario_4,
        scen5=scenario_5,
        td=test_data2,tdB=test_data4,net="Olyset Net",a=4,colnet="darkred")

i=10
  plot_4a(scen=scenario_2,
          scen3=scenario_3,
          scen4=scenario_4,
          scen5=scenario_5,
          td=test_data2,tdB=test_data4,net="Olyset Plus",a=i,colnet="darkorange")

i=22
  plot_4a(scen=scenario_2,
          scen3=scenario_3,
          scen4=scenario_4,
          scen5=scenario_5,
          td=test_data2,
          tdB=test_data4,net="PermaNet 2.0",a=i,colnet="blue")

i=6
  plot_4a(scen=scenario_2,
        scen3=scenario_3,
        scen4=scenario_4,
        scen5=scenario_5,
        td=test_data2,
        tdB=test_data4,net="PermaNet 3.0",a=i,colnet="darkgreen")



#######################################################
##
## observed = a + b1*predicted + b2*monthsurvey + error
## a = 0
## 
#######################################################

DF = data.frame(cluster = clusters,
                prev_obs = prev_obs,
                months = months,
                prev_mod0 = prev_mod0,
                prev_mod1 = prev_mod1,
                prev_mod1_1 = prev_mod1_1,
                prev_mod1_2 = prev_mod1_2,
                prev_mod1_3 = prev_mod1_3,
                prev_mod1_4 = prev_mod1_4,
                prev_mod2 = prev_mod2,
                prev_mod3 = prev_mod3,
                prev_mod4 = prev_mod4,
                prev_mod5 = prev_mod5)

## Compare specific questions
## When predicting impact using parameters specifying pyr-only nets
## is there a difference in the estimate for PBOs?
  
DF$net_type = c(rep("Oly",52),rep("OlyPlus",76),rep("Per2",160),rep("Per3",128))
DF$logic_net = c(rep(0,52),rep(1,76),rep(0,160),rep(1,128))
DF$pos = round(1000*DF$prev_obs,0)
DF$neg = 1000 - DF$pos  

DF$ad_0 = prev_obs - prev_mod0
DF$ad_1 = prev_obs - prev_mod1
DF$ad_1_1 = prev_obs - prev_mod1_1
DF$ad_1_2 = prev_obs - prev_mod1_2
DF$ad_1_3 = prev_obs - prev_mod1_3
DF$ad_1_4 = prev_obs - prev_mod1_4
DF$ad_2 = prev_obs - prev_mod2
DF$ad_3 = prev_obs - prev_mod3
DF$ad_4 = prev_obs - prev_mod4
DF$ad_5 = prev_obs - prev_mod5

DF = DF[!(DF$cluster==78 | DF$cluster==27),]  

tapply(DF$ad_0,DF$months,median,na.rm=TRUE)
tapply(DF$ad_1,DF$months,median,na.rm=TRUE)
tapply(DF$ad_1_1,DF$months,median,na.rm=TRUE)
tapply(DF$ad_1_2,DF$months,median,na.rm=TRUE)
tapply(DF$ad_1_3,DF$months,median,na.rm=TRUE)
tapply(DF$ad_1_4,DF$months,median,na.rm=TRUE)
tapply(DF$ad_2,DF$months,median,na.rm=TRUE)
tapply(DF$ad_3,DF$months,median,na.rm=TRUE)
tapply(DF$ad_4,DF$months,median,na.rm=TRUE)
tapply(DF$ad_5,DF$months,median,na.rm=TRUE)

q_fun = function(dt){
  a = quantile(dt[DF$months == 6],c(0.025,0.975),na.rm=TRUE)
  b = quantile(dt[DF$months == 12],c(0.025,0.975),na.rm=TRUE)
  d = quantile(dt[DF$months == 18],c(0.025,0.975),na.rm=TRUE)
  f = quantile(dt[DF$months == 25],c(0.025,0.975),na.rm=TRUE)
  return(list(a,b,d,f))
}
q_fun(DF$ad_0)
q_fun(DF$ad_1)
q_fun(DF$ad_1_1)
q_fun(DF$ad_1_2)
q_fun(DF$ad_1_3)
q_fun(DF$ad_1_4)
q_fun(DF$ad_2)
q_fun(DF$ad_3)
q_fun(DF$ad_4)
q_fun(DF$ad_5)

quantile(DF$ad_0,c(0.5,0.025,0.975),na.rm=T)
quantile(DF$ad_1,c(0.5,0.025,0.975),na.rm=T)
quantile(DF$ad_1_1,c(0.5,0.025,0.975),na.rm=T)
quantile(DF$ad_1_2,c(0.5,0.025,0.975),na.rm=T)
quantile(DF$ad_1_3,c(0.5,0.025,0.975),na.rm=T)
quantile(DF$ad_1_4,c(0.5,0.025,0.975),na.rm=T)
quantile(DF$ad_2,c(0.5,0.025,0.975),na.rm=T)
quantile(DF$ad_3,c(0.5,0.025,0.975),na.rm=T)
quantile(DF$ad_4,c(0.5,0.025,0.975),na.rm=T)
quantile(DF$ad_5,c(0.5,0.025,0.975),na.rm=T)


simp0 = rstanarm::stan_glm(prev_obs ~ prev_mod0 + months,data=DF)
simp1 = rstanarm::stan_glm(prev_obs ~ prev_mod1 + months,data=DF)
simp11 = rstanarm::stan_glm(prev_obs ~ prev_mod1_1 + months,data=DF)
simp12 = rstanarm::stan_glm(prev_obs ~ prev_mod1_2 + months,data=DF)
simp13 = rstanarm::stan_glm(prev_obs ~ prev_mod1_3 + months,data=DF)
simp14 = rstanarm::stan_glm(prev_obs ~ prev_mod1_4 + months,data=DF)
simp2 = rstanarm::stan_glm(prev_obs ~ prev_mod2 + months,data=DF)
simp3 = rstanarm::stan_glm(prev_obs ~ prev_mod3 + months,data=DF)
simp4 = rstanarm::stan_glm(prev_obs ~ prev_mod4 + months,data=DF)
simp5 = rstanarm::stan_glm(prev_obs ~ prev_mod5 + months,data=DF)

## intercept estimate, gradient estimate
coef(simp0) 
coef(simp1) 
coef(simp11)
coef(simp12) 
coef(simp13) 
coef(simp14) 
coef(simp2) 
coef(simp3)
coef(simp4)
coef(simp5) 


sims0 = as.matrix(simp0)
sims1 = as.matrix(simp1)
sims11 = as.matrix(simp11)
sims12 = as.matrix(simp12)
sims13 = as.matrix(simp13)
sims14 = as.matrix(simp14)
sims2 = as.matrix(simp2)
sims3 = as.matrix(simp3)
sims4 = as.matrix(simp4)
sims5 = as.matrix(simp5)

## uncertainty on gradient estimate
quantile(sims0[,2],c(0.025,0.5,0.975));quantile(sims0[,3],c(0.025,0.5,0.975))
quantile(sims1[,2],c(0.025,0.5,0.975));quantile(sims1[,3],c(0.025,0.5,0.975))
quantile(sims11[,2],c(0.025,0.5,0.975));quantile(sims11[,3],c(0.025,0.5,0.975))
quantile(sims12[,2],c(0.025,0.5,0.975));quantile(sims12[,3],c(0.025,0.5,0.975))
quantile(sims13[,2],c(0.025,0.5,0.975));quantile(sims13[,3],c(0.025,0.5,0.975))
quantile(sims14[,2],c(0.025,0.5,0.975));quantile(sims14[,3],c(0.025,0.5,0.975))
quantile(sims2[,2],c(0.025,0.5,0.975));quantile(sims2[,3],c(0.025,0.5,0.975))
quantile(sims3[,2],c(0.025,0.5,0.975));quantile(sims3[,3],c(0.025,0.5,0.975))
quantile(sims4[,2],c(0.025,0.5,0.975));quantile(sims4[,3],c(0.025,0.5,0.975))
quantile(sims5[,2],c(0.025,0.5,0.975));quantile(sims5[,3],c(0.025,0.5,0.975))

# drop months
simp0 = rstanarm::stan_glm(prev_obs ~ prev_mod0,data=DF)
simp1 = rstanarm::stan_glm(prev_obs ~ prev_mod1,data=DF)
simp11 = rstanarm::stan_glm(prev_obs ~ prev_mod1_1,data=DF)
simp12 = rstanarm::stan_glm(prev_obs ~ prev_mod1_2,data=DF)
simp13 = rstanarm::stan_glm(prev_obs ~ prev_mod1_3,data=DF)
simp14 = rstanarm::stan_glm(prev_obs ~ prev_mod1_4,data=DF)
simp2 = rstanarm::stan_glm(prev_obs ~ prev_mod2,data=DF)
simp3 = rstanarm::stan_glm(prev_obs ~ prev_mod3,data=DF)
simp4 = rstanarm::stan_glm(prev_obs ~ prev_mod4,data=DF)
simp5 = rstanarm::stan_glm(prev_obs ~ prev_mod5,data=DF)

sims0 = as.matrix(simp0)
sims1 = as.matrix(simp1)
sims11 = as.matrix(simp11)
sims12 = as.matrix(simp12)
sims13 = as.matrix(simp13)
sims14 = as.matrix(simp14)
sims2 = as.matrix(simp2)
sims3 = as.matrix(simp3)
sims4 = as.matrix(simp4)
sims5 = as.matrix(simp5)

## intercept estimate, gradient estimate
coef(simp0) 
coef(simp1) 
coef(simp11)
coef(simp12) 
coef(simp13) 
coef(simp14) 
coef(simp2) 
coef(simp3)
coef(simp4)
coef(simp5) 

## uncertainty on gradient estimate
quantile(sims0[,2],c(0.025,0.5,0.975))
quantile(sims1[,2],c(0.025,0.5,0.975))
quantile(sims11[,2],c(0.025,0.5,0.975))
quantile(sims12[,2],c(0.025,0.5,0.975))
quantile(sims13[,2],c(0.025,0.5,0.975))
quantile(sims14[,2],c(0.025,0.5,0.975))
quantile(sims2[,2],c(0.025,0.5,0.975))
quantile(sims3[,2],c(0.025,0.5,0.975))
quantile(sims4[,2],c(0.025,0.5,0.975))
quantile(sims5[,2],c(0.025,0.5,0.975))

## fitting all together with a covariate for simulation
dt = expand.grid(prev_empirical = rep(prev_obs,10))
dt$prev_simulated = c(prev_mod0,prev_mod1,
                   prev_mod1_1,prev_mod1_2,prev_mod1_3,prev_mod1_4,
                   prev_mod2,prev_mod3,prev_mod4,prev_mod5)
dt$simulationflag = rep(1:10, each = 416)

simp = rstanarm::stan_glm(prev_empirical ~ prev_simulated + as.factor(simulationflag),data=dt)
sims = as.matrix(simp)
quantile(sims[,2],c(0.025,0.5,0.975))
quantile(sims[,3],c(0.025,0.5,0.975))
print(simp)

m1 = glm(prev_empirical ~ prev_simulated + as.factor(simulationflag),data=dt)
summary(m1)
summary.lm(m1)

## Linear Regressions - which model predicts trial data best?

summary.lm(lm(DF$prev_obs ~ DF$prev_mod0 + 0))
summary.lm(lm(DF$prev_obs ~ DF$prev_mod1 + 0))
summary.lm(lm(DF$prev_obs ~ DF$prev_mod1_1 + 0))
summary.lm(lm(DF$prev_obs ~ DF$prev_mod1_2 + 0))
summary.lm(lm(DF$prev_obs ~ DF$prev_mod1_3 + 0))
summary.lm(lm(DF$prev_obs ~ DF$prev_mod1_4 + 0))
summary.lm(lm(DF$prev_obs ~ DF$prev_mod2 + 0))
summary.lm(lm(DF$prev_obs ~ DF$prev_mod3 + 0))
summary.lm(lm(DF$prev_obs ~ DF$prev_mod4 + 0)) 
summary.lm(lm(DF$prev_obs ~ DF$prev_mod5 + 0))

#######################################################
##
## x = a + b1*X + b2*m + b3*C + error
##
## MODEL 1: 1 / logit(x) = (exp^x)/(1 - exp^x)
##
#######################################################



mod0 = rstanarm::stan_glm(cbind(pos,neg) ~ prev_mod0 + logic_net,
                          family = binomial(link = "logit"), 
                          data = DF)

print(mod0)
par(mfrow=c(2,2))
plot(DF$prev_obs ~ DF$prev_mod0,ylim=c(0,1),xlim=c(0,1),
     col = ifelse(DF$logic_net == 0, "blue","darkred"),
     pch=19,
     yaxt = "n", ylab = "Observed prevalence in 2 - 10 years (%)",
     xaxt = "n", xlab = "Simulated prevalence in 2 - 10 years (%)")
axis(1,at = seq(0,1,0.2), labels = seq(0,100,20))
axis(2,las=2,at = seq(0,1,0.2), labels = seq(0,100,20))
curve(rstanarm::invlogit(coef(mod0)[1] + coef(mod0)[2] * x), col="blue", add=TRUE)
curve(rstanarm::invlogit(coef(mod0)[1]+coef(mod0)[3] + coef(mod0)[2] * x), col="darkred", add=TRUE)
abline(a=0,b=1,lty=2,col="grey")

sims1 = as.matrix(mod0)
n_sims1 = nrow(sims1)
for(j in sample (n_sims1,100)){
  curve(rstanarm::invlogit(sims1[j,1] +sims1[j,2] * x), col=adegenet::transp("blue",0.1), add=TRUE)
  curve(rstanarm::invlogit(sims1[j,1]-sims1[j,3] + sims1[j,2] * x), col=adegenet::transp("darkred",0.1), add=TRUE)
  
}

legend("topleft",legend = c("Pyrethroid-only nets",
                            "Pyrethroid-PBO nets"),
       col=c("blue","darkred"),pch=19,bty="n")




mod1 = rstanarm::stan_glm(cbind(pos,neg) ~ prev_mod1 + logic_net,
                          family = binomial(link = "logit"), 
                          data = DF)

print(mod1)

plot(DF$prev_obs ~ DF$prev_mod1,ylim=c(0,1),xlim=c(0,1),
     col = ifelse(DF$logic_net == 0, "blue","darkred"),
     pch=19,
     yaxt = "n", ylab = "Observed prevalence in 2 - 10 years (%)",
     xaxt = "n", xlab = "Simulated prevalence in 2 - 10 years (%)")
axis(1,at = seq(0,1,0.2), labels = seq(0,100,20))
axis(2,las=2,at = seq(0,1,0.2), labels = seq(0,100,20))
curve(rstanarm::invlogit(coef(mod1)[1] + coef(mod1)[2] * x), col="blue", add=TRUE)
curve(rstanarm::invlogit(coef(mod1)[1]-coef(mod1)[3] + coef(mod1)[2] * x), col="darkred", add=TRUE)
abline(a=0,b=1,lty=2,col="grey")

sims1 = as.matrix(mod1)
n_sims1 = nrow(sims1)
for(j in sample (n_sims1,100)){
  curve(rstanarm::invlogit(sims1[j,1] +sims1[j,2] * x), col=adegenet::transp("blue",0.1), add=TRUE)
  curve(rstanarm::invlogit(sims1[j,1]-sims1[j,3] + sims1[j,2] * x), col=adegenet::transp("darkred",0.1), add=TRUE)
  
}

legend("topleft",legend = c("Pyrethroid-only nets",
                            "Pyrethroid-PBO nets"),
       col=c("blue","darkred"),pch=19,bty="n")


sim_comparison_func = function(prev_modelled,vec_plot,print_sim){
  mod1 = rstanarm::stan_glm(cbind(pos,neg) ~ prev_modelled + logic_net,
                            family = binomial(link = "logit"), 
                            data = DF)
  
  print(mod1)
  
  
  plot(DF$prev_obs ~ vec_plot,ylim=c(0,1),xlim=c(0,1),
       col = ifelse(DF$logic_net == 0, "grey","purple"),
       pch=19,
       main = print_sim,
       yaxt = "n", ylab = "Observed prevalence in 2 - 10 years (%)",
       xaxt = "n", xlab = "Simulated prevalence in 2 - 10 years (%)")
  axis(1,at = seq(0,1,0.2), labels = seq(0,100,20))
  axis(2,las=2,at = seq(0,1,0.2), labels = seq(0,100,20))
  curve(rstanarm::invlogit(coef(mod1)[1] + coef(mod1)[2] * x), col="grey", add=TRUE)
  curve(rstanarm::invlogit(coef(mod1)[1]+coef(mod1)[3] + coef(mod1)[2] * x), col="purple", add=TRUE)
  abline(a=0,b=1,lty=2,col="grey")
  
  sims1 = as.matrix(mod1)
  n_sims1 = nrow(sims1)
  for(j in sample (n_sims1,100)){
    curve(rstanarm::invlogit(sims1[j,1] +sims1[j,2] * x), col=adegenet::transp("grey",0.1), add=TRUE)
    curve(rstanarm::invlogit(sims1[j,1]+sims1[j,3] + sims1[j,2] * x), col=adegenet::transp("purple",0.1), add=TRUE)
    
  }
  
  legend("topleft",legend = c("Pyrethroid-only nets",
                              "Pyrethroid-PBO nets"),
         col=c("grey","purple"),pch=19,bty="n")

  
  boxplot(DF$prev_obs[DF$logic_net == 0],
          DF$prev_obs[DF$logic_net == 1],
          at=c(0.88,0.94),boxwex=.05,
          yaxt="n",xaxt="n",
          col=adegenet::transp(c("grey","purple"),0.4),add=TRUE)
  return(summary(mod1))
}


sim_comparison_func(prev_modelled = DF$prev_mod0,
                    vec_plot=DF$prev_mod0, print_sim = "Sim 0")

boxplot(DF$prev_obs[DF$net_type == "Oly"],
        DF$prev_obs[DF$net_type == "Per2"],
        DF$prev_obs[DF$net_type == "OlyPlus"],
        DF$prev_obs[DF$net_type == "Per3"],
        at=c(0.68,0.72,0.76,0.8),boxwex=.025,
        yaxt="n",xaxt="n",
        col=adegenet::transp(c("red","blue","orange","darkgreen"),0.4),add=TRUE)

sim_comparison_func(DF$prev_mod1,vec_plot=DF$prev_mod1, print_sim = "Sim 1")
sim_comparison_func(DF$prev_mod1_1,vec_plot=DF$prev_mod1_1, print_sim = "Sim 1-1")
sim_comparison_func(DF$prev_mod1_2,vec_plot=DF$prev_mod1_2, print_sim = "Sim 1-2")
sim_comparison_func(DF$prev_mod1_3,vec_plot=DF$prev_mod1_3, print_sim = "Sim 1-3")
sim_comparison_func(DF$prev_mod1_4,vec_plot=DF$prev_mod1_4, print_sim = "Sim 1-4")
sim_comparison_func(DF$prev_mod2,vec_plot=DF$prev_mod2, print_sim = "Sim 2")
sim_comparison_func(DF$prev_mod3,vec_plot=DF$prev_mod3, print_sim = "Sim 3")
sim_comparison_func(DF$prev_mod4,vec_plot=DF$prev_mod4, print_sim = "Sim 4")
sim_comparison_func(DF$prev_mod5,vec_plot=DF$prev_mod5, print_sim = "Sim 5")


##########################################
##
## lm plots
##
##########################################

lm_plot_fn = function(prevsim,sc,mn_sc_sim){
  plot(DF$prev_obs ~ prevsim, ylim=c(0,0.8),xlim=c(0,0.8),
       main = sc,
       ylab="Trial reported prevalence (%)",yaxt="n",
       xlab = "Model simulated prevalence (%)",xaxt="n",
       col=adegenet::transp(rep(cols,each=4),0.1),pch=19)
  axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
  axis(2, las = 2, at = seq(0,1,0.2), labels = seq(0,100,20))
  abline(a=0,b=1,lty=2)
  
  points(prev_obs_means ~ mn_sc_sim,col=cols,pch=19,cex=1.2)
}

lm_plot_fn(prevsim = DF$prev_mod0, sc = "Scenario 0", mn_sc_sim = prev_mod_means0)
boxplot(DF$prev_obs[DF$net_type == "Oly"],
        DF$prev_obs[DF$net_type == "Per2"],
        DF$prev_obs[DF$net_type == "OlyPlus"],
        DF$prev_obs[DF$net_type == "Per3"],
        at=c(0.68,0.72,0.76,0.8),boxwex=.025,
        yaxt="n",xaxt="n",
        col=adegenet::transp(c("red","blue","orange","darkgreen"),0.4),add=TRUE)

lm_plot_fn(prevsim = DF$prev_mod1, sc = "Scenario 1", mn_sc_sim = prev_mod_means1)
lm_plot_fn(prevsim = DF$prev_mod1_1, sc = "Scenario 1-1", mn_sc_sim = prev_mod_means1_1)
lm_plot_fn(prevsim = DF$prev_mod1_2, sc = "Scenario 1-2", mn_sc_sim = prev_mod_means1_2)
lm_plot_fn(prevsim = DF$prev_mod1_3, sc = "Scenario 1-3", mn_sc_sim = prev_mod_means1_3)
lm_plot_fn(prevsim = DF$prev_mod1_4, sc = "Scenario 1-4", mn_sc_sim = prev_mod_means1_4)
lm_plot_fn(prevsim = DF$prev_mod2, sc = "Scenario 2", mn_sc_sim = prev_mod_means2)
lm_plot_fn(prevsim = DF$prev_mod3, sc = "Scenario 3", mn_sc_sim = prev_mod_means3)
lm_plot_fn(prevsim = DF$prev_mod4, sc = "Scenario 4", mn_sc_sim = prev_mod_means4)
lm_plot_fn(prevsim = DF$prev_mod5, sc = "Scenario 5", mn_sc_sim = prev_mod_means5)
