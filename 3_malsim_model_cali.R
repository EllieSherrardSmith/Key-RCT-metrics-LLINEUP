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

test_data = read.csv("raw data/data.csv",header=TRUE)


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
  human_population <- 10000
  
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

# eir_est = numeric(104)
for(i in 1:104){
  eir_est[i] = cali_clusters_f(i)
}

## 44 is higher!

test_data$eir_est = eir_est
write.csv(test_data,"C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/raw data/data_with_calibrated_eir.csv")

plot(test_data$eir_est ~ test_data$Prevalence_baseline_2_10_yrs,
     cex = 1+log(1+test_data$LLIN_USE_BASELINE),col=adegenet::transp("grey",0.6),pch=19)
points(test_data$eir_est[test_data$res_2 > 0.68] ~ 
         test_data$Prevalence_baseline_2_10_yrs[test_data$res_2 > 0.68],
       cex = 1+log(1+test_data$LLIN_USE_BASELINE),
       col="blue",pch=19)
# 
# xx = seq(0,1,length=20)
# alpha = -3.4
# beta =7
# yy = 1 / (1 + exp(-alpha - 
#                     beta*xx))
# xx
# yy_est = 100 * (1 / (1 + exp(-alpha - 
#                         beta*test_data$Prevalence_baseline_2_10_yrs)))




test_data$EIR = comparison_accept

dat_row = 2
year=365

place1 = malsim_actual_f(test_data,dat_row,eir=test_data$eir_est[dat_row])
plot(place1$prev_pyr ~ place1$timestep,xlim=c(5.5*365,8.5*365),ylim=c(0,1),type="l")
points(test_data$Prevalence_baseline_2_10_yrs[dat_row]~c(6*365+test_data$days_after_jan_2017[dat_row]-30),
       col="orange",pch=19)




timestep = c(0, 3, 6, 9) * year + c(0,
                                    ## raw data: net timing
                                    test_data$days_after_jan_2017[dat_row],
                                    test_data$days_after_jan_2017[dat_row],
                                    test_data$days_after_jan_2017[dat_row])
abline(v=timestep,lty=2)
points(test_data$Prevalence_6m[dat_row]~
         c(6*365+test_data$days_after_jan_2017[dat_row]+365/2),
       col="red",pch=19)
points(test_data$Prevalence_12m[dat_row]~
         c(6*365+test_data$days_after_jan_2017[dat_row]+365),
       col="red",pch=19)
points(test_data$Prevalence_18m[dat_row]~
         c(6*365+test_data$days_after_jan_2017[dat_row]+365+365/2),
       col="red",pch=19)
points(test_data$Prevalence_25m[dat_row]~
         c(6*365+test_data$days_after_jan_2017[dat_row]+365*2+30),
       col="red",pch=19)
