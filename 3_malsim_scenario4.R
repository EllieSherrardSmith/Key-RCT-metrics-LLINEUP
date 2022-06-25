##############################################
##
## Validation simulation for the cluster analysis of LLINEUP data
##
## Stage 1 we want to match as best as possible the 
## trial data and see how closely we can recreate the
## cross-sectional survey results

## Functions from packages
# remotes::install_github("mrc-ide/cali") 
## just needed this once to load in cali

# remotes::install_github("mrc-ide/malariasimulation")

library(malariasimulation)
library(cali)

# library(ggplot2)


## Data
test_data = read.csv("simulation-summary/scenario_4.csv",header=TRUE)

##########################################
##
## Function to run all sims first for calibrations to baseline prevalence
## assuming this is 1 month prior to the trial initiation

## Adjusting this to include a 2023 distribution running to 2026
## so we can compare to the 2 year net distribution campaign
## create the parameters list for the malaria model simulation: 

cali_clusters_f = function(dat_row){
  year <- 365
  month <- 30
  sim_length <- 15 * year ## initially just til 2020 as need extra resistance data!
  ## This is spanning Jan 2014 - Dec 2025
  ## Assume that all places have received nets at start of 2014
  ## then relative to jan 2017 (see uganda_what_if...)
  ## then as noted for 2020
  
  
  ## and finally we wish to see the difference for these places running forward from 2020-2023 (1 or 3 years)
  human_population <- 6000
  
  ## This is calibrated to reflect the mosquito densities: raw data
  starting_EIR <- 1
  # dat_row = 1
  
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      # irs_correlation = 
      prevalence_rendering_min_ages = 2 * 365, ## Prev in 6 months to 14 years measured
      prevalence_rendering_max_ages = 10 * 365,
      
      model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
      
      ## match to estimates of seasonal trends from Pete Winskill GTS updates
      g0 = test_data$seasonal_a0[dat_row],
      g = c(test_data$seasonal_a1[dat_row], test_data$seasonal_a2[dat_row], test_data$seasonal_a3[dat_row]),
      h = c(test_data$seasonal_b1[dat_row], test_data$seasonal_b2[dat_row], test_data$seasonal_b3[dat_row]),
      
      Q0 = 0.72,
      individual_mosquitoes = FALSE ## Update next
    )
  )
  
  simparams <- set_equilibrium(simparams, starting_EIR)
  
  # set species
  ## match to trial data - gambiae sl to arabiensis: raw data
  simparams <- set_species(simparams,
                           species=list(gamb_params, arab_params, fun_params),
                           proportions=c(0.4,
                                         0.6, ## mosquitoes mostly became arabiensis dominated
                                         0))
  # set treatment
  ## an unknown and assuming from World Malaria Report country-wide averages
  simparams <- set_drugs(simparams, list(AL_params, SP_AQ_params, DHA_PQP_params))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,  ## In Uganda used for uncomplicated malaria
                                      ##https://www.cdc.gov/globalhealth/countries/uganda/pdf/uganda_factsheet.pdf
                                      ##Uganda population 2017 = 42862958
                                      ##ACT courses delivered 2017 = 27396300
                                      ## Estimated no of cases: 12140161
                                      ## Popn at risk: 44269584
                                      time=c(100),
                                      coverage=0.62) ## 27396300/44269584
  simparams <- set_clinical_treatment(simparams, 
                                      drug=2,  ## SP for severe... 
                                      time=c(100),
                                      coverage=0.07) ## estimated pregnancy
  
  
  ## Set up bed nets
  
  bednetparams <- simparams
  
  ## as done
  bednet_events = data.frame(
    timestep = c(0, 3, 6, 9) * year + c(0,
                                        ## raw data: net timing
                                        test_data$days_after_jan_2017[dat_row],
                                        test_data$days_after_jan_2017[dat_row],
                                        test_data$days_after_jan_2017[dat_row]),
    name=c("background",
           "background",
           "trial_nets",
           "2020_nets")
  )
  
  
  # Will run over a loop so this is for either net
  # these with be the pyr-ITNs/pob-ITNs depending on row of data
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    ## raw data: net use prior and immediately after the campaign
    coverages = c(test_data$LLIN_USE_BASELINE[dat_row], ## historic prior to RCT
                  test_data$LLIN_USE_BASELINE[dat_row], ## historic prior 
                  test_data$sleeping_under_net_6m[dat_row], ## RCT (on deployment) 
                  test_data$sleeping_under_net_6m[dat_row]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    ## net attrition from raw data
    retention = test_data$itn_leave_durMean[dat_row] * year, ## Keeping this as it was observed during RCT
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito) *note transposing so 
    ## raw data: resistance status
    dn0 = t(matrix(as.numeric(rep(c(test_data$dn0_med_2014[dat_row],	test_data$dn0_med_2017[dat_row]),each=6)), nrow=3, ncol=4)),
    rn =  t(matrix(as.numeric(rep(c(test_data$rn0_med_2014[dat_row], test_data$rn0_med_2017[dat_row]),each=6)), nrow=3, ncol=4)),
    rnm = matrix(c(.24, .24, .24, .24), nrow=4, ncol=3),
    gamman = as.numeric(rep(c(test_data$gamman_med_2014[dat_row],test_data$gamman_med_2017[dat_row]),each=2) * 365)
  )
  
  ##end
  ## assume the same people are getting nets each round
  correlationsb1 <- get_correlation_parameters(bednetparams_1)
  correlationsb1$inter_round_rho('bednets', 1)
  
  # Define target, here two prevalence measures:
  target <- c(test_data$Prevalence_baseline_2_10_yrs[dat_row])
  # Time points at which to match target
  target_tt <- c(6*365+test_data$days_after_jan_2017[dat_row]-30)
  
  ## Run the simulations
  output1 <- run_simulation(sim_length, bednetparams_1,correlationsb1)
  
  set.seed(123)
  out <- calibrate(parameters = bednetparams_1,
                   target = target,
                   target_tt = target_tt,
                   summary_function = summary_pfpr_2_10,
                   tolerance = 0.02, 
                   interval = c(1, 500))##upper bound needs to be high enough so negative differences are not returned in uniroot
  
  return(out$root)
  
}

# eir_est = numeric(nrow(test_data))
for(i in 1:nrow(test_data)){
  eir_est[i] = cali_clusters_f(i)
}
test_data$eir_est = eir_est

write.csv(test_data,"C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_4.csv")



plot(test_data$eir_est ~ test_data$Prevalence_baseline_2_10_yrs,
     cex = 1+log(1+test_data$LLIN_USE_BASELINE),col=adegenet::transp("grey",0.6),pch=19)


##################################
##
## Simulate 
test_data = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_4.csv",header=TRUE)

uganda_1_f = function(dat_row){
  year <- 365
  month <- 30
  sim_length <- 15 * year ## initially just til 2020 as need extra resistance data!
  ## This is spanning Jan 2014 - Dec 2025
  ## Assume that all places have received nets at start of 2014
  ## then relative to jan 2017 (see uganda_what_if...)
  ## then as noted for 2020
  
  
  ## and finally we wish to see the difference for these places running forward from 2020-2023 (1 or 3 years)
  human_population <- 6000
  
  ## This is calibrated to reflect the mosquito densities: raw data
  starting_EIR <- test_data$eir_est[dat_row]
  
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      # irs_correlation = 
      prevalence_rendering_min_ages = 2 * 365, ## Prev in 6 months to 14 years measured
      prevalence_rendering_max_ages = 10 * 365,
      
      model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
      
      ## match to estimates of seasonal trends from Pete Winskill GTS updates
      g0 = test_data$seasonal_a0[dat_row],
      g = c(test_data$seasonal_a1[dat_row], test_data$seasonal_a2[dat_row], test_data$seasonal_a3[dat_row]),
      h = c(test_data$seasonal_b1[dat_row], test_data$seasonal_b2[dat_row], test_data$seasonal_b3[dat_row]),
      
      Q0 = 0.72,
      individual_mosquitoes = FALSE ## Update next
    )
  )
  
  simparams <- set_equilibrium(simparams, starting_EIR)
  
  # set species
  ## match to trial data - gambiae sl to arabiensis: raw data
  simparams <- set_species(simparams,
                           species=list(gamb_params, arab_params, fun_params),
                           proportions=c(0.4,
                                         0.6, ## mosquitoes mostly became arabiensis dominated
                                         0))
  # set treatment
  ## an unknown and assuming from World Malaria Report country-wide averages
  simparams <- set_drugs(simparams, list(AL_params, SP_AQ_params, DHA_PQP_params))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,  ## In Uganda used for uncomplicated malaria
                                      ##https://www.cdc.gov/globalhealth/countries/uganda/pdf/uganda_factsheet.pdf
                                      ##Uganda population 2017 = 42862958
                                      ##ACT courses delivered 2017 = 27396300
                                      ## Estimated no of cases: 12140161
                                      ## Popn at risk: 44269584
                                      time=c(100),
                                      coverage=0.62) ## 27396300/44269584
  simparams <- set_clinical_treatment(simparams, 
                                      drug=2,  ## SP for severe... 
                                      time=c(100),
                                      coverage=0.07) ## estimated pregnancy
  
  
  ## Set up bed nets
  
  bednetparams <- simparams
  
  ## as done
  bednet_events = data.frame(
    timestep = c(0, 3, 6, 9) * year + c(0,
                                        ## raw data: net timing
                                        test_data$days_after_jan_2017[dat_row],
                                        test_data$days_after_jan_2017[dat_row],
                                        test_data$days_after_jan_2017[dat_row]),
    name=c("background",
           "background",
           "trial_nets",
           "2020_nets")
  )
  
  
  # Will run over a loop so this is for either net
  # these with be the pyr-ITNs/pob-ITNs depending on row of data
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    ## raw data: net use prior and immediately after the campaign
    coverages = c(test_data$LLIN_USE_BASELINE[dat_row], ## historic prior to RCT
                  test_data$LLIN_USE_BASELINE[dat_row], ## historic prior 
                  test_data$sleeping_under_net_6m[dat_row], ## RCT (on deployment) 
                  test_data$sleeping_under_net_6m[dat_row]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    ## net attrition from raw data
    retention = test_data$itn_leave_durMean[dat_row] * year, ## Keeping this as it was observed during RCT
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito) *note transposing so 
    ## raw data: resistance status
    dn0 = t(matrix(as.numeric(rep(c(test_data$dn0_med_2014[dat_row],	test_data$dn0_med_2017[dat_row]),each=6)), nrow=3, ncol=4)),
    rn =  t(matrix(as.numeric(rep(c(test_data$rn0_med_2014[dat_row], test_data$rn0_med_2017[dat_row]),each=6)), nrow=3, ncol=4)),
    rnm = matrix(c(.24, .24, .24, .24), nrow=4, ncol=3),
    gamman = as.numeric(rep(c(test_data$gamman_med_2014[dat_row],test_data$gamman_med_2017[dat_row]),each=2) * 365)
  )
  
  ##end
  ## assume the same people are getting nets each round
  correlationsb1 <- get_correlation_parameters(bednetparams_1)
  correlationsb1$inter_round_rho('bednets', 1)

  ## Run the simulations
  output1 <- run_simulation(sim_length, bednetparams_1,correlationsb1)
  output1$pv_730_3650 = output1$n_detect_730_3650/output1$n_730_3650
  
  return(data.frame(timestep = output1$timestep,
                    prev_730_3650 = output1$pv_730_3650))
}

mod_sims_scenario_4 = expand.grid(
  timestep = scenario_0$timestep
)
for(i in 1:nrow(test_data)){
  mod_sims_scenario_4[,i+1] = uganda_1_f(i)[,2]
}

for(i in 1:nrow(test_data)){
  colnames(mod_sims_scenario_4)[i+1] = paste0("cluster",test_data$cluster[i],"_prev")  
}


write.csv(mod_sims_scenario_4, "simulation-summary/mod_sims_scenario_4.csv")


########################################
##
## Starting here with presenting results
##
########################################

## simulations for prevalence
scenario_4 = read.csv("simulation-summary/mod_sims_scenario_4.csv",header=TRUE)[,2:106]
## input parameters for scenario
test_data = read.csv("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/simulation-summary/scenario_4.csv",header=TRUE)
## looking across clusters
dat = read.csv("simulation-summary/scenario_4.csv",header=TRUE)

plot(scenario_4$cluster1_prev ~ scenario_4$timestep,type="l",
     ylab = "Prevalence in children 6m to 10yrs (%)",yaxt="n",ylim=c(0,0.8),
     xlab = "Time in days",col="darkred",xlim=c(2300,3500)
)


cols2=ifelse(test_data$Net_Type == "Olyset Net", "red",
             ifelse(test_data$Net_Type == "Olyset Plus", "orange",
                    ifelse(test_data$Net_Type == "PermaNet 2.0", "darkblue","green")))

year=365
dat_row=1

# Define target, here two prevalence measures:
target <- c(test_data$Prevalence_baseline_2_10_yrs)
# Time points at which to match target
target_tt <- c(6*365+test_data$days_after_jan_2017-30)
points(target ~ target_tt,col = cols2)



##
## Need to check out 
## that the clusters match with the net types and 
## colours align 
par(mfrow=c(2,2))
plot(scenario_4[,2] ~ scenario_4$timestep,type="l",
     ylab = "Prevalence in children 6m to 10yrs (%)",yaxt="n",ylim=c(0,0.8),
     xlab = "Time in days",col="white",xlim=c(2300,3500)
)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
which(test_data$Net_Type == "Olyset Net")
for(i in c(which(test_data$Net_Type == "Olyset Net")+1)){
  lines(scenario_4[,i] ~ scenario_4$timestep,col="darkred")
}
unique(test_data$days_after_jan_2017)
timestep = 6 * year + c(unique(test_data$days_after_jan_2017))

abline(v=timestep,lty=2)

# Define target, here two prevalence measures:
target <- c(test_data$Prevalence_baseline_2_10_yrs)
# Time points at which to match target
target_tt <- c(6*365+test_data$days_after_jan_2017-30)
points(target[dat$Net_Type == "Olyset Net"] ~ target_tt[dat$Net_Type == "Olyset Net"],
       col = cols2[dat$Net_Type == "Olyset Net"])
points(dat$Prevalence_6m[dat$Net_Type == "Olyset Net"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "Olyset Net"]+365/2),
       col=adegenet::transp("red",0.2),pch=19)
points(dat$Prevalence_12m[dat$Net_Type == "Olyset Net"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "Olyset Net"]+365),
       col=adegenet::transp("red",0.2),pch=19)
points(dat$Prevalence_18m[dat$Net_Type == "Olyset Net"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "Olyset Net"]+365+365/2),
       col=adegenet::transp("red",0.2),pch=19)
points(dat$Prevalence_25m[dat$Net_Type == "Olyset Net"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "Olyset Net"]+365+365+30),
       col=adegenet::transp("red",0.2),pch=19)


plot(scenario_4[,2] ~ scenario_4$timestep,type="l",
     ylab = "Prevalence in children 6m to 10yrs (%)",yaxt="n",ylim=c(0,0.8),
     xlab = "Time in days",col="white",xlim=c(2300,3500)
)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v=timestep,lty=2)

for(i in c(which(test_data$Net_Type == "Olyset Plus")+1)){
  lines(scenario_4[,i] ~ scenario_4$timestep,col="orange")
}
points(target[dat$Net_Type == "Olyset Plus"] ~ target_tt[dat$Net_Type == "Olyset Plus"],
       col = cols2[dat$Net_Type == "Olyset Plus"])
points(dat$Prevalence_6m[dat$Net_Type == "Olyset Plus"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "Olyset Plus"]+365/2),
       col=adegenet::transp("darkorange",0.2),pch=19)
points(dat$Prevalence_12m[dat$Net_Type == "Olyset Plus"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "Olyset Plus"]+365),
       col=adegenet::transp("darkorange",0.2),pch=19)
points(dat$Prevalence_18m[dat$Net_Type == "Olyset Plus"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "Olyset Plus"]+365+365/2),
       col=adegenet::transp("darkorange",0.2),pch=19)
points(dat$Prevalence_25m[dat$Net_Type == "Olyset Plus"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "Olyset Plus"]+365+365+30),
       col=adegenet::transp("darkorange",0.2),pch=19)


plot(scenario_4[,2] ~ scenario_4$timestep,type="l",
     ylab = "Prevalence in children 6m to 10yrs (%)",yaxt="n",ylim=c(0,0.8),
     xlab = "Time in days",col="white",xlim=c(2300,3500)
)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v=timestep,lty=2)

for(i in c(which(test_data$Net_Type == "PermaNet 2.0")+1)){
  lines(scenario_4[,i] ~ scenario_4$timestep,col="blue")
}
points(target[dat$Net_Type == "PermaNet 2.0"] ~ target_tt[dat$Net_Type == "PermaNet 2.0"],
       col = cols2[dat$Net_Type == "PermaNet 2.0"])
points(dat$Prevalence_6m[dat$Net_Type == "PermaNet 2.0"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "PermaNet 2.0"]+365/2),
       col=adegenet::transp("darkblue",0.2),pch=19)
points(dat$Prevalence_12m[dat$Net_Type == "PermaNet 2.0"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "PermaNet 2.0"]+365),
       col=adegenet::transp("darkblue",0.2),pch=19)
points(dat$Prevalence_18m[dat$Net_Type == "PermaNet 2.0"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "PermaNet 2.0"]+365+365/2),
       col=adegenet::transp("darkblue",0.2),pch=19)
points(dat$Prevalence_25m[dat$Net_Type == "PermaNet 2.0"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "PermaNet 2.0"]+365+365+30),
       col=adegenet::transp("darkblue",0.2),pch=19)


plot(scenario_4$cluster1_prev ~ scenario_4$timestep,type="l",
     ylab = "Prevalence in children 6m to 10yrs (%)",yaxt="n",ylim=c(0,0.8),
     xlab = "Time in days",col="white",xlim=c(2300,3500)
)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v=timestep,lty=2)

for(i in c(which(test_data$Net_Type == "PermaNet 3.0")+1)){
  lines(scenario_4[,i] ~ scenario_4$timestep,col="darkgreen")
}
points(target[dat$Net_Type == "PermaNet 3.0"] ~ target_tt[dat$Net_Type == "PermaNet 3.0"],
       col = cols2[dat$Net_Type == "PermaNet 3.0"])
points(dat$Prevalence_6m[dat$Net_Type == "PermaNet 3.0"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "PermaNet 3.0"]+365/2),
       col=adegenet::transp("green",0.2),pch=19)
points(dat$Prevalence_12m[dat$Net_Type == "PermaNet 3.0"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "PermaNet 3.0"]+365),
       col=adegenet::transp("green",0.2),pch=19)
points(dat$Prevalence_18m[dat$Net_Type == "PermaNet 3.0"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "PermaNet 3.0"]+365+365/2),
       col=adegenet::transp("green",0.2),pch=19)
points(dat$Prevalence_25m[dat$Net_Type == "PermaNet 3.0"]~
         c(6*365+dat$days_after_jan_2017[dat$Net_Type == "PermaNet 3.0"]+365+365+30),
       col=adegenet::transp("green",0.2),pch=19)



#######
##
## points versus time points

prev_obs = c(dat$Prevalence_6m[dat$Net_Type == "Olyset Net"],
             dat$Prevalence_12m[dat$Net_Type == "Olyset Net"],
             dat$Prevalence_18m[dat$Net_Type == "Olyset Net"],
             dat$Prevalence_25m[dat$Net_Type == "Olyset Net"],
             
             dat$Prevalence_6m[dat$Net_Type == "Olyset Plus"],
             dat$Prevalence_12m[dat$Net_Type == "Olyset Plus"],
             dat$Prevalence_18m[dat$Net_Type == "Olyset Plus"],
             dat$Prevalence_25m[dat$Net_Type == "Olyset Plus"],
             
             dat$Prevalence_6m[dat$Net_Type == "PermaNet 2.0"],
             dat$Prevalence_12m[dat$Net_Type == "PermaNet 2.0"],
             dat$Prevalence_18m[dat$Net_Type == "PermaNet 2.0"],
             dat$Prevalence_25m[dat$Net_Type == "PermaNet 2.0"],
             
             dat$Prevalence_6m[dat$Net_Type == "PermaNet 3.0"],
             dat$Prevalence_12m[dat$Net_Type == "PermaNet 3.0"],
             dat$Prevalence_18m[dat$Net_Type == "PermaNet 3.0"],
             dat$Prevalence_25m[dat$Net_Type == "PermaNet 3.0"])

prev_mod_ol = array(dim=c(4,length(which(test_data$Net_Type == "Olyset Net"))))
olys = c(which(test_data$Net_Type == "Olyset Net")+1)
for(i in 1:length(c(which(test_data$Net_Type == "Olyset Net")))){
  prev_mod_ol[1,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365/2),olys[i]]
  prev_mod_ol[2,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365),olys[i]]
  prev_mod_ol[3,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365+365/2),olys[i]]
  prev_mod_ol[4,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365*2+30),olys[i]] 
}

prev_mod_olP = array(dim=c(4,length(which(test_data$Net_Type == "Olyset Plus"))))
olyP = c(which(test_data$Net_Type == "Olyset Plus")+1)
for(i in 1:length(c(which(test_data$Net_Type == "Olyset Plus")))){
  prev_mod_olP[1,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365/2),olyP[i]]
  prev_mod_olP[2,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365),olyP[i]]
  prev_mod_olP[3,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365+365/2),olyP[i]]
  prev_mod_olP[4,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365*2+30),olyP[i]] 
}

prev_mod_Pe2 = array(dim=c(4,length(which(test_data$Net_Type == "PermaNet 2.0"))))
Per2 = c(which(test_data$Net_Type == "PermaNet 2.0")+1)
for(i in 1:length(c(which(test_data$Net_Type == "PermaNet 2.0")))){
  prev_mod_Pe2[1,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365/2),Per2[i]]
  prev_mod_Pe2[2,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365),Per2[i]]
  prev_mod_Pe2[3,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365+365/2),Per2[i]]
  prev_mod_Pe2[4,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365*2+30),Per2[i]] 
}

prev_mod_Pe3 = array(dim=c(4,length(which(test_data$Net_Type == "PermaNet 3.0"))))
Per3 = c(which(test_data$Net_Type == "PermaNet 3.0")+1)
for(i in 1:length(c(which(test_data$Net_Type == "PermaNet 3.0")))){
  prev_mod_Pe3[1,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365/2),Per3[i]]
  prev_mod_Pe3[2,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365),Per3[i]]
  prev_mod_Pe3[3,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365+365/2),Per3[i]]
  prev_mod_Pe3[4,i] = scenario_4[c(6*365+dat$days_after_jan_2017[i]+365*2+30),Per3[i]] 
}

prev_mod = c(prev_mod_ol[1,],prev_mod_ol[2,],prev_mod_ol[3,],prev_mod_ol[4,],
             prev_mod_olP[1,],prev_mod_olP[2,],prev_mod_olP[3,],prev_mod_olP[4,],
             prev_mod_Pe2[1,],prev_mod_Pe2[2,],prev_mod_Pe2[3,],prev_mod_Pe2[4,],
             prev_mod_Pe3[1,],prev_mod_Pe3[2,],prev_mod_Pe3[3,],prev_mod_Pe3[4,])

plot(prev_obs ~ prev_mod, ylim=c(0,0.8),xlim=c(0,0.8),
     col=adegenet::transp(rep(cols,each=4),0.1),pch=19)
abline(a=0,b=1,lty=2)

prev_obs_means = c(mean(prev_obs[1:52],na.rm=TRUE),
                   mean(prev_obs[53:128],na.rm=TRUE),
                   mean(prev_obs[129:288],na.rm=TRUE),
                   mean(prev_obs[289:416],na.rm=TRUE))

prev_mod_means = c(mean(prev_mod[1:52],na.rm=TRUE),
                   mean(prev_mod[53:128],na.rm=TRUE),
                   mean(prev_mod[129:288],na.rm=TRUE),
                   mean(prev_mod[289:416],na.rm=TRUE))

points(prev_obs_means ~ prev_mod_means,col=cols,pch=19,cex=1.2)


