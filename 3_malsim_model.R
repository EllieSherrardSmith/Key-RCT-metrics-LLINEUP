##############################################
##
## Validation simulation for the cluster analysis of LLINEUP data
##
## Stage 1 we want to match as best as possible the 
## trial data and see how closely we can recreate the
## cross-sectional survey results

## Data

test_data = read.csv("")

## Function

library(malariasimulation)


##########################################
##
## Function to run all sims
# library(malariasimulation@feat/alt_resistance)

ellie_cluster_run_malsim <- function(run){
  rr = read.csv(paste0("E:/set_up_Uganda/RCT_modelling/uncertainty_runs/runs1.csv"),header=TRUE)
  test7b = read.csv(paste0("E:/set_up_Uganda/RCT_modelling/uncertainty_runs/runs",run,".csv"),header=TRUE)
  test7b$ITN_2020 = rr$ITN_2020
  data_store = list()
  for(i in 1:nrow(test7b)){
    # for(i in 1:2){
    data_store[[i]] = malsim_actual_f(test_data = test7b,
                                      dat_row = i)
  }
  ## Can do this for each replicated simulation set for uncertainty...
  saveRDS(data_store, file = paste0("RCT_modelling/uncertainty_runs/simulations_staedke_covid_2and3yrNets_",run,".rds"))
}

## Adjusting this to include a 2023 distribution running to 2026
## so we can compare to the 2 year net distribution campaign

malsim_actual_f = function(test_data,dat_row){
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
  starting_EIR <- 65
  
  
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
  # these with be the pyr-ITNs
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    ## raw data: net use prior and immediately after the campaign
    coverages = c(test_data$ITN_hist[dat_row], ## historic prior to RCT
                  test_data$ITN_hist[dat_row], ## historic prior 
                  test_data$ITN_2017[dat_row], ## RCT (on deployment) 
                  test_data$ITN_2017[dat_row]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    ## net attrition from raw data
    retention = test_data$itn_leave_dur[dat_row] * year, ## Keeping this as it was observed during RCT
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito) *note transposing so 
    ## raw data: resistance status
    dn0 = t(matrix(as.numeric(c(test_data$dITN.x[dat_row], test_data$dITN.x[dat_row], test_data$dITN.x[dat_row],
                                test_data$dITN.y[dat_row], test_data$dITN.y[dat_row], test_data$dITN.y[dat_row],
                                test_data$dITN[dat_row], test_data$dITN[dat_row], test_data$dITN[dat_row],
                                test_data$dITN[dat_row], test_data$dITN[dat_row], test_data$dITN[dat_row])), nrow=3, ncol=4)),
    rn = t(matrix(as.numeric(c(test_data$rITN.x[dat_row], test_data$rITN.x[dat_row], test_data$rITN.x[dat_row],
                               test_data$rITN.y[dat_row], test_data$rITN.y[dat_row], test_data$rITN.y[dat_row],
                               test_data$rITN[dat_row], test_data$rITN[dat_row], test_data$rITN[dat_row],
                               test_data$rITN[dat_row], test_data$rITN[dat_row], test_data$rITN[dat_row])), nrow=3, ncol=4)),
    rnm = matrix(c(.24, .24, .24, .24), nrow=5, ncol=3),
    gamman = as.numeric(c(test_data$hlf.x[dat_row], test_data$hlf.y[dat_row], test_data$hlf[dat_row], test_data$hlf[dat_row])) * 365
  )
  
  
  ## these will be the pyr-PBO ITNs
  bednetparams_2 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    ## raw data: net use prior and immediately after the campaign
    coverages = c(test_data$ITN_hist[dat_row], ## historic prior to RCT
                  test_data$ITN_hist[dat_row], ## historic prior 
                  test_data$ITN_2017[dat_row], ## RCT (on deployment) 
                  test_data$ITN_2017[dat_row]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    ## net attrition from raw data
    retention = test_data$itn_leave_dur[dat_row] * year, ## Keeping this as it was observed during RCT
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito) *note transposing so 
    ## raw data: resistance status
    dn0 = t(matrix(as.numeric(c(test_data$dITN.x[dat_row], test_data$dITN.x[dat_row], test_data$dITN.x[dat_row],
                                test_data$dITN.y[dat_row], test_data$dITN.y[dat_row], test_data$dITN.y[dat_row],
                                test_data$dITN[dat_row], test_data$dITN[dat_row], test_data$dITN[dat_row],
                                test_data$dITN[dat_row], test_data$dITN[dat_row], test_data$dITN[dat_row])), nrow=3, ncol=4)),
    rn = t(matrix(as.numeric(c(test_data$rITN.x[dat_row], test_data$rITN.x[dat_row], test_data$rITN.x[dat_row],
                               test_data$rITN.y[dat_row], test_data$rITN.y[dat_row], test_data$rITN.y[dat_row],
                               test_data$rITN[dat_row], test_data$rITN[dat_row], test_data$rITN[dat_row],
                               test_data$rITN[dat_row], test_data$rITN[dat_row], test_data$rITN[dat_row])), nrow=3, ncol=4)),
    rnm = matrix(c(.24, .24, .24, .24), nrow=5, ncol=3),
    gamman = as.numeric(c(test_data$hlf.x[dat_row], test_data$hlf.y[dat_row], test_data$hlf[dat_row], test_data$hlf[dat_row])) * 365
  )
  
  
  ##end
  ## assume the same people are getting nets each round
  correlationsb1 <- get_correlation_parameters(bednetparams_1)
  correlationsb2 <- get_correlation_parameters(bednetparams_2)
  correlationsb1$inter_round_rho('bednets', 1)
  correlationsb2$inter_round_rho('bednets', 1)
  
  
  
  ## Run the simulations
  output1_pyr <- run_simulation(sim_length, bednetparams_1,correlationsb1)
  output2_pbo <- run_simulation(sim_length, bednetparams_2,correlationsb1)
  
  return(data.frame(timestep = output1_matching_rct_and_2020$timestep,
                    
                    prev_pyr = output1_pyr$pv_730_3650,
                    prev_pbo = output2_pbo$pv_730_3650
                    
  ))
}