#############################################################
##
## Final script to determine parameters for net efficacy
##
##############################################################

###########################
##
## Stage 1 running each net type with just the generic probabilities
##         with each net type for the benefit in mortality 
##              a = pyrethroid only
##              b = pyr-PBO
##              c = pyr=chlor Interceptor G2 
##              d = IG2 excluding Burkina Faso data

## Stage 2 running net specific probabilities
##

## DECISION: to assume all nets are working equivalently
##           when it comes to expt hut associations
##           Use median estimates for these from 'a' all data
##           Include uncertainty for the mortality associations
##           That is bioassay and hut
##           And added mortality from new nets
#####################################

####################################################
##
## Fits look ok so now extract the uncertainti using the function 

## Now we need to add in the net parameter estimates 
## working from the function the estimate crude estimates and
## using the parameters determined in Table 2 main manuscipt

## Critical is to keep, for each fit
## the same row of parameters for that simualtion

library(rstan)

## keep rows associated from the bayes posterior draws
data_picker = sample(1:4000,size = 1000,replace=TRUE)

## b is when we use just pyrethroid and pyrethroid-PBO nets
## h is just permaNet 3
## j is just Olyset Plus

## draw from the posterior distribution with the respective inputs
## 
## THIS IS THE SAME FOR ALL COMBINATIONS
## ASSOCIATION 1 MORTALITY TO BIOASSAY
# setwd("Q:/RProjects/Mosquito-Net-Parameters/stan model outputs")
# setwd("C:/Users/esherrar/Documents/Rprojects/Mosquito-Net-Parameters/stan model outputs")
# setwd("stan model outputs")
getwd() ## confirm

# Global set
# ll_1 = readRDS("fit_ew_comb_log_logistic.rds")
ll_1 = readRDS("stan model outputs/log_logistic_fit.rds")
LL_fit <- extract(ll_1, permuted = TRUE)

#function shape:
x = seq(0,1,length=101) ## this is mosquito survival 

f_LOG_logistic <- function(x, b, a){
  surv = 1 - (1/(1+((1-x)/b)^(-a)))
  mort = (1-surv)
  
  return(mort)
}
## returns hut mortality

## keep uncertainty
LL_b_full <- LL_fit$b[data_picker]
LL_a_full <- LL_fit$a[data_picker]

## Just confirm nothing spurious here
hut_mort_LL <- f_LOG_logistic(x = seq(0,1,length=101), LL_b_full[5], LL_a_full[5])
plot(hut_mort_LL ~ seq(0,1,length=101),ylim=c(0,1),
     xlab="bioassay survival"
)




## When bioassay mortality is 0, hut mortality is 0.57 ish
## When bioassay mortality is 1. hut mortlity is 1

## 
## ASSOCIATION 2 BENEFIT OF PBO
## THIS REMAINS THE SAME FOR ALL COMBINATIONS
benefitpbo <- readRDS("stan model outputs/ento_pbo_benefit.rds")
pbo_bene <- extract(benefitpbo, permuted = TRUE)

fitbene_1 <- pbo_bene$alpha1[data_picker]
fitbene_2 <- pbo_bene$alpha2[data_picker]

f_LOG_logistic_with_benefit <- function(mort, alp1, alp2){
  mort_newnet = 1 / (1 + exp(-alp1 - alp2 * mort))
  # surv = (1-mort)
  return(mort_newnet)
}

hut_mort_LLPBO <- f_LOG_logistic_with_benefit(mort = seq(0,1,length=101), 
                                              alp1 = fitbene_1[5],
                                              alp2 = fitbene_2[5])

plot(hut_mort_LLPBO ~ seq(0,1,length=101),
     xlab = "hut_mort_LL",
     YLAB = "New net hut mortality")


## ONLY USING MEDIAN ESTIMATES FROM FULL FITS
## i.e. all data from WHO recommended nets
## and assuming these associations hold for all nets
## ASSOCIATION 3 MORTALITY TO DETERRENCE
# load fit - just using All nets
# setwd("E:/Mosquito net parameters/stan model outputs")
# getwd() ## confirm
# fit1_a <- readRDS("ento_deterrence_AllRec.rds")
# fit1_a_fit <- extract(fit1_a, permuted = TRUE)
fit1_a_fit <- readRDS("stan model outputs/April_2022_ento_deterrence_AllRec_extract.rds")

## All nets (WHO recommended as per Okumu & Finda 2021) 
fit1_a_c <- fit1_a_fit$c[c(data_picker)]
fit1_a_d <- fit1_a_fit$d[c(data_picker)]
fit1_a_e <- fit1_a_fit$e[c(data_picker)]

fit1_a_med_c = median(fit1_a_fit$c)
fit1_a_med_d = median(fit1_a_fit$d)
fit1_a_med_e = median(fit1_a_fit$e)

quantile(fit1_a_fit$c,c(0.1,0.9))
which(fit1_a_fit$c > 2.08 & fit1_a_fit$c < 2.087)
RN_e = 0.36 ## These are the original estimates frmo Nash et al 2021
RN_d = 0.49 ## but we now use only WHO recommended nets, unwashed, and
RN_c = 2.57 ## of any type (pyr-only, pyr-PBO, pyr-chlor).

mort = seq(0,1,length=101)
deterrence = fit1_a_med_e * (exp(fit1_a_med_d * (1 - exp(fit1_a_med_c * mort)) / fit1_a_med_c))
deterrence_low = quantile(fit1_a_fit$e,0.9) * (exp(quantile(fit1_a_fit$d,0.1) * (1 - exp(quantile(fit1_a_fit$c,0.1) * mort)) / quantile(fit1_a_fit$c,0.1)))
deterrence_upp = quantile(fit1_a_fit$e,0.1) * (exp(quantile(fit1_a_fit$d,0.9) * (1 - exp(quantile(fit1_a_fit$c,0.9) * mort)) / quantile(fit1_a_fit$c,0.9)))

deterrenceRN = RN_e * (exp(RN_d * (1 - exp(RN_c * mort)) / RN_c))

plot(deterrence ~ mort,xaxt="n",yaxt="n",
     xlab = "Hut mortality of any WHO Recommended net (%)",
     ylab = "Deterrence due to presence of that net (%)",
     ylim=c(0,1),xlim=c(0,1))
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
lines(deterrenceRN ~ mort)
polygon(c(mort,rev(mort)),c(deterrence_low,rev(deterrence_upp)),col=adegenet::transp("darkgreen",0.4),border=NA)

## 
## ASSOCIATION 4 MORTALITY TO feed
# fit3_a <- readRDS("ento_feeding_0washes_AllRec.rds")
# fit3_a_fit <- extract(fit3_a, permuted = TRUE)
fit3_a_fit <- readRDS("stan model outputs/April_2022_ento_feeding_0washes_AllRec_extract.rds")

fit3_a_f <- median(fit3_a_fit$a)
fit3_a_g <- median(fit3_a_fit$b)

feeding = (1 - (exp(fit3_a_g * (1 - exp(fit3_a_f * mort))/fit3_a_f)))
feeding_low = (1 - (exp(quantile(fit3_a_fit$b,0.1) * (1 - exp(quantile(fit3_a_fit$a,0.1) * mort))/quantile(fit3_a_fit$a,0.1))))
feeding_upp = (1 - (exp(quantile(fit3_a_fit$b,0.9) * (1 - exp(quantile(fit3_a_fit$a,0.9) * mort))/quantile(fit3_a_fit$a,0.9))))

plot(feeding ~ mort,xaxt="n",yaxt="n",
     xlab = "Hut mortality of any WHO Recommended net (%)",
     ylab = "Success feeding in presence of that net (%)",
     ylim=c(0,1),xlim=c(0,1))
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
polygon(c(mort,rev(mort)),c(feeding_low,rev(feeding_upp)),col=adegenet::transp("darkred",0.4),border=NA)


RN_f <- 4.66
RN_g <- 0.04

feedingRN = (1 - (exp(RN_g * (1 - exp(RN_f * mort))/RN_f)))
lines(feedingRN ~ mort)

## We are planning to update this in 2023 once data are available
## Half life we could stick to the original - but this is then
## assuming a logistic association between susc bioassay and hut surv...
## Half life from all data - with log-logistic
# fit_hl_a <- readRDS("ento_fit1_half_life_a.rds") ## All data
# 
# testa = rstan::extract(fit_hl_a)
# 
# fithalf_a_e <- testa$a[c(data_picker)]
# fithalf_b_e <- testa$b[c(data_picker)]
# 
# mort = seq(0,1,length=101)
# P_hl1 <- 1/(1 + exp(-(median(testa$a) + median(testa$b) * mort)))
# plot(P_hl1 ~ mort)


resistance_ITN_default_params_2_f = function(product, 
                                             # res, ## SURVIVAL IN SUSC BIOASSAY
                                             shape,## for half life with log-log or original logistic
                                             data_picker_rand){ ## random draw from posterior pred  
  
  ## The parameters included here are from the statistical analysis 
  ## following on from Nash et al 2021 
  
  ## PARAMETERS
  
  # give the uncertainty for the log-logistic function
  #Assay to hut mortality conversion - median estimates	
  param_b = LL_b_full[data_picker_rand] 
  param_a = LL_a_full[data_picker_rand] 
  
  # x = res
  
  f_LOG_logistic <- function(x, b, a){
    surv = 1 - (1/(1+((1-x)/b)^(-a)))
    mort = (1-surv)
    
    return(mort)
  }
  
  # mort = 1 - res ## added benefit is relative to mortality in pyr-only nets
  
  f_LOG_logistic_with_benefit <- function(mort, alp1, alp2){
    mort_newnet = 1 / (1 + exp(-alp1 - alp2 * mort))
    # surv = (1-mort)
    return(mort_newnet)
  }
  ## returns hut mortality
  
  hut_mort_LL <- f_LOG_logistic(x = seq(0,1,length=101),#** 
                                param_b, param_a)
  
  hut_mort_LLPBOtemp <- f_LOG_logistic_with_benefit(mort = hut_mort_LL, 
                                                    alp1 = fitbene_1[data_picker_rand],
                                                    alp2 = fitbene_2[data_picker_rand])
  

  #specify whichever net is used in the RCT
  #this will determine what the mortality is in the hut trial
  mort_huta = if(product==0) hut_mort_LL else if(product==1) hut_mort_LLPBOtemp 
  mort_hut = mort_huta
  
  hut_surv = 1 - mort_hut
  # ff, ff1, and ff2 all match up for the associations for net det and fed
  # when using the generic parameters as here
  
  #Decay in insecticide non-PBO net		
  mup =	-2.429#NEW 1.812767 ###hlf_mu# #mu_p ## sample(test$a,size=1) ##array(c(rep(-2.429,3),rep(-2.984025,3),rep(-1.866,3)),c(3,3)) ## ... gam.medians[9]
  rhop =-3.007 #NEW -2.581591 ###hlf_rho# #rho_p ##  sample(test$b,size=1)	##array(c(rep(-3.007,3),rep(-3.74,3),rep(-2.295,3)),c(3,3)) ## ... gam.medians[10]
  
  ## The maximum successful feeding probability per feeding attempt 
  ## (feeding and not dying) in the absence of interventions 
  kp0=0.699 ## derived from Lines et al 1987 and Curtis et al 1990 
  
  ## The half-life of the net relative to it's capacity to kill mosquitoes
  ## with the insecticide active ingredient (a pyrethroid) when there is 
  ## no resistance in mosquitoes. 
  net_halflife=2.64
  
  ## Now we work through the probability steps to determine the key input parameters for the model
  ## These probability relationships are determined by Rebecca Nash, Ben and Tom see email notes above
  
  fit1_a_med_c = median(fit1_a_fit$c)
  fit1_a_med_d = median(fit1_a_fit$d)
  fit1_a_med_e = median(fit1_a_fit$e)
  
  ## This is association with hut survival
  det_hut = fit1_a_med_e * (exp(fit1_a_med_d * (1 - exp(fit1_a_med_c * hut_surv)) / fit1_a_med_c))
  
  fit3_a_f <- median(fit3_a_fit$a)
  fit3_a_g <- median(fit3_a_fit$b)
  
  ## This is association with hut survival
  suc_hut = (1 - (exp(fit3_a_g * (1 - exp(fit3_a_f * hut_surv))/fit3_a_f)))
  
  ## This is association with hut survival
  rep_hut = (1 - suc_hut - mort_hut)
  
  xx = data.frame(hut_surv,mort_hut,suc_hut,rep_hut,det_hut)
  ## Combine to estimate the 3 key probable outcomes of feeding attempts
  ## Here we adjust for those mosquitoes not entering treated huts (determined by deterrence)
  n1n0 = 1-xx$det_hut
  kp1  = n1n0*xx$suc_hut
  lp1  = n1n0*xx$mort_hut
  jp1  = n1n0*xx$rep_hut+(1-n1n0)
  
  kp1 = ifelse(kp1 > kp0,kp0,kp1) ## Capping impact so max feeding is no bigger than assumed
  # # max feeding for no interventions (kp0 = 0.699, Griffin et al 2010)
  # ## (time = 0 time steps after net implementation)
  # 
  # kp0 = 1
  # 
  r_ITN0  = (1-kp1/kp0)*(jp1/(lp1+jp1))	#probability of repeating behaviour
  d_ITN0  = (1-kp1/kp0)*(lp1/(lp1+jp1))	#probability of dying with an encounter with ITN
  s_ITN0  = 1-d_ITN0-r_ITN0             #probability of successfully feeding (surviving and feeding)
  
  
  # plot(r_ITN0 ~ c(1-mort),ylim=c(0,1),xlim=c(0,1),xlab = "Susc bioassay survial",type="l",col="orange") 
  
  ## Repeat these to determin the maximum and minimum effects which combine to help determine ITN half life
  ## We will stick to pyr-params for half life and update in 2023 once new data are available
  mort_maxA   = if (shape=="log-log") f_LOG_logistic(x = 0,#** this is surv i.e. mort max when surv =0
                                                     param_b, param_a) else if (shape=="logistic") 1/(1+exp(param_a*param_b-param_a*(1)))
  
  
 
  mort_max = mort_maxA#if(product==0) mort_maxA else if(product==1) PBO_benefitA else if(product==2) G2_benefitA else if(product==3) G2_benefitB
  
  
  #{halflife}
  my_max_washes_a = mup +rhop*(mort_max-0.5)		
  # my_max_washes   = log(2)/(exp(my_max_washes_a)/(1+exp(my_max_washes_a)))
  
  my_max_washes   = log(2)/(1/(1+exp(-my_max_washes_a)))
  
  
  ## Uncertainty
  net_half_life_min = 2
  net_half_life_max = 3
  
  ## FOR NON PYR-ONLY NETS THIS WILL RETURN HIGH HALF LIFE
  ## WE RECOMMEND ONLY USING PY-ONLY HALF LIFE UNTIL WE CAN
  ## VALIDATE OTHER NETS 
  wash_decay_rate_a = mup +rhop*(mort_hut-0.5)
  # wash_decay_rate   = log(2)/(exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a)))
  wash_decay_rate   = log(2)/(1/(1+exp(-wash_decay_rate_a)))
  itn_half_life     = wash_decay_rate/my_max_washes*net_halflife
  itn_half_life_max     = wash_decay_rate/my_max_washes*net_half_life_max  
  itn_half_life_min     = wash_decay_rate/my_max_washes*net_half_life_min
  ##No need to re-adjust these anymore
  ##Final Parameter estimates for the transmission model
  ERG_d_ITN0 <- d_ITN0
  ERG_s_ITN0 <- s_ITN0
  ERG_r_ITN0 <- 1-ERG_d_ITN0-ERG_s_ITN0
  
  ## Print out these estimates to a data.frame as the function output
  uncertainty_resistance_params_nets = data.frame(ERG_d_ITN0,ERG_r_ITN0,itn_half_life,
                                                  itn_half_life_min,itn_half_life_max,
                                                  det_hut=xx$det_hut,suc_hut=xx$suc_hut,mort_hut=xx$mort_hut,rep_hut=xx$rep_hut,
                                                  n1n0,kp1,lp1,jp1,
                                                  bioassay_surv = round(seq(0,1,length=101),2))
  
  return(uncertainty_resistance_params_nets)
}

## 20 uncertainty draws
## load the orignial data
data_raw = read.csv("raw data/data_with_calibrated_eir.csv",header=TRUE)

## minimise key data needed for aligning uncertainty
align_dat = data.frame(cluster = data_raw$cluster,
                       net_type = data_raw$Net_Type,
                       res_1 = round(data_raw$res_1,2),
                       res_2 = round(data_raw$res_2,2))


rand = sample(1:1000,20,replace=FALSE)
for(i in 1:length(rand)){
  test = resistance_ITN_default_params_2_f(product = 0, ## PYRETHROID ONLY LLIN 
                                           # res = seq(0,1,length=101), ## SURVIVAL IN SUSC BIOASSAY
                                           shape = "logistic",## for half life with log-log or original logistic
                                           data_picker_rand = rand[i])[,c(1:3,14)] ## any number from 1 to 1000

  colnames(test)[4] = "res_1"
  SET1a = merge(align_dat,test,by="res_1",all.x = TRUE)	

  colnames(test)[4] = "res_2"
  SET1b = merge(align_dat,test,by="res_2",all.x = TRUE)	
  
  test2 = resistance_ITN_default_params_2_f(product = 1, ## PYRETHROID ONLY LLIN 
                                           # res = seq(0,1,length=101), ## SURVIVAL IN SUSC BIOASSAY
                                           shape = "logistic",## for half life with log-log or original logistic
                                           data_picker_rand = rand[i])[,c(1:3,14)]
  
  colnames(test2)[4] = "res_2"
  SET1c = merge(align_dat,test2,by="res_2",all.x = TRUE)	
  
  SET1all = SET1a
  SET1all$d_ITN0_2017 = ifelse(SET1all$net_type == "Olyset Net",SET1b$ERG_d_ITN0,
                               ifelse(SET1all$net_type == "PermaNet 2.0",SET1b$ERG_d_ITN0,
                                      SET1c$ERG_d_ITN0))
  SET1all$r_ITN0_2017 = ifelse(SET1all$net_type == "Olyset Net",SET1b$ERG_r_ITN0,
                               ifelse(SET1all$net_type == "PermaNet 2.0",SET1b$ERG_r_ITN0,
                                      SET1c$ERG_r_ITN0))
  
  SET1all$itn_half_2017 = ifelse(SET1all$net_type == "Olyset Net",SET1b$itn_half_life,
                               ifelse(SET1all$net_type == "PermaNet 2.0",SET1b$itn_half_life,
                                      SET1b$itn_half_life)) ## because we are assuming same half life for both nets
  
  write.csv(SET1all,paste0("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/raw data/default-net-params_",i,".csv"))  
}


#####################################################
##
## The above is the default using all nets to estimate params

#####################################################
##
## Next pyrethroid only nets

fit1_a <- readRDS("stan model outputs/LLINEUP_ento_deterrence_pyr.rds")
fit1_a_fit <- extract(fit1_a, permuted = TRUE)

fit1_a_med_c = median(fit1_a_fit$c)
fit1_a_med_d = median(fit1_a_fit$d)
fit1_a_med_e = median(fit1_a_fit$e)
## ASSOCIATION 4 MORTALITY TO feed
fit3_a <- readRDS("stan model outputs/LLINEUP_ento_feeding_0washes_pyr.rds")
fit3_a_fit <- extract(fit3_a, permuted = TRUE)

fit3_a_f <- median(fit3_a_fit$a)
fit3_a_g <- median(fit3_a_fit$b)

resistance_ITN_pyr_only_params_2_f = function(product, 
                                             # res, ## SURVIVAL IN SUSC BIOASSAY
                                             shape,## for half life with log-log or original logistic
                                             data_picker_rand){ ## random draw from posterior pred  
  
  ## The parameters included here are from the statistical analysis 
  ## following on from Nash et al 2021 
  
  ## PARAMETERS
  
  # give the uncertainty for the log-logistic function
  #Assay to hut mortality conversion - median estimates	
  param_b = LL_b_full[data_picker_rand] 
  param_a = LL_a_full[data_picker_rand] 
  
  # x = res
  
  f_LOG_logistic <- function(x, b, a){
    surv = 1 - (1/(1+((1-x)/b)^(-a)))
    mort = (1-surv)
    
    return(mort)
  }
  
  # mort = 1 - res ## added benefit is relative to mortality in pyr-only nets
  
  f_LOG_logistic_with_benefit <- function(mort, alp1, alp2){
    mort_newnet = 1 / (1 + exp(-alp1 - alp2 * mort))
    # surv = (1-mort)
    return(mort_newnet)
  }
  ## returns hut mortality
  
  hut_mort_LL <- f_LOG_logistic(x = seq(0,1,length=101),#** 
                                param_b, param_a)
  
  hut_mort_LLPBOtemp <- f_LOG_logistic_with_benefit(mort = hut_mort_LL, 
                                                    alp1 = fitbene_1[data_picker_rand],
                                                    alp2 = fitbene_2[data_picker_rand])
  
  
  #specify whichever net is used in the RCT
  #this will determine what the mortality is in the hut trial
  mort_huta = if(product==0) hut_mort_LL else if(product==1) hut_mort_LLPBOtemp 
  mort_hut = mort_huta
  
  hut_surv = 1 - mort_hut
  # ff, ff1, and ff2 all match up for the associations for net det and fed
  # when using the generic parameters as here
  
  #Decay in insecticide non-PBO net		
  mup =	-2.429#NEW 1.812767 ###hlf_mu# #mu_p ## sample(test$a,size=1) ##array(c(rep(-2.429,3),rep(-2.984025,3),rep(-1.866,3)),c(3,3)) ## ... gam.medians[9]
  rhop =-3.007 #NEW -2.581591 ###hlf_rho# #rho_p ##  sample(test$b,size=1)	##array(c(rep(-3.007,3),rep(-3.74,3),rep(-2.295,3)),c(3,3)) ## ... gam.medians[10]
  
  ## The maximum successful feeding probability per feeding attempt 
  ## (feeding and not dying) in the absence of interventions 
  kp0=0.699 ## derived from Lines et al 1987 and Curtis et al 1990 
  
  ## The half-life of the net relative to it's capacity to kill mosquitoes
  ## with the insecticide active ingredient (a pyrethroid) when there is 
  ## no resistance in mosquitoes. 
  net_halflife=2.64
  
  ## Now we work through the probability steps to determine the key input parameters for the model
  ## These probability relationships are determined by Rebecca Nash, Ben and Tom see email notes above
  
  fit1_a_med_c = median(fit1_a_fit$c)
  fit1_a_med_d = median(fit1_a_fit$d)
  fit1_a_med_e = median(fit1_a_fit$e)
  
  ## This is association with hut survival
  det_hut = fit1_a_med_e * (exp(fit1_a_med_d * (1 - exp(fit1_a_med_c * hut_surv)) / fit1_a_med_c))
  
  fit3_a_f <- median(fit3_a_fit$a)
  fit3_a_g <- median(fit3_a_fit$b)
  
  ## This is association with hut survival
  suc_hut = (1 - (exp(fit3_a_g * (1 - exp(fit3_a_f * hut_surv))/fit3_a_f)))
  
  ## This is association with hut survival
  rep_hut = (1 - suc_hut - mort_hut)
  
  xx = data.frame(hut_surv,mort_hut,suc_hut,rep_hut,det_hut)
  ## Combine to estimate the 3 key probable outcomes of feeding attempts
  ## Here we adjust for those mosquitoes not entering treated huts (determined by deterrence)
  n1n0 = 1-xx$det_hut
  kp1  = n1n0*xx$suc_hut
  lp1  = n1n0*xx$mort_hut
  jp1  = n1n0*xx$rep_hut+(1-n1n0)
  
  kp1 = ifelse(kp1 > kp0,kp0,kp1) ## Capping impact so max feeding is no bigger than assumed
  # # max feeding for no interventions (kp0 = 0.699, Griffin et al 2010)
  # ## (time = 0 time steps after net implementation)
  # 
  # kp0 = 1
  # 
  r_ITN0  = (1-kp1/kp0)*(jp1/(lp1+jp1))	#probability of repeating behaviour
  d_ITN0  = (1-kp1/kp0)*(lp1/(lp1+jp1))	#probability of dying with an encounter with ITN
  s_ITN0  = 1-d_ITN0-r_ITN0             #probability of successfully feeding (surviving and feeding)
  
  # plot(r_ITN0 ~ c(1-mort),ylim=c(0,1),xlim=c(0,1),xlab = "Susc bioassay survial",type="l",col="orange") 
  
  ## Repeat these to determin the maximum and minimum effects which combine to help determine ITN half life
  ## We will stick to pyr-params for half life and update in 2023 once new data are available
  mort_maxA   = if (shape=="log-log") f_LOG_logistic(x = 0,#** this is surv i.e. mort max when surv =0
                                                     param_b, param_a) else if (shape=="logistic") 1/(1+exp(param_a*param_b-param_a*(1)))
  
  
  mort_max = mort_maxA#if(product==0) mort_maxA else if(product==1) PBO_benefitA else if(product==2) G2_benefitA else if(product==3) G2_benefitB
  
  
  #{halflife}
  my_max_washes_a = mup +rhop*(mort_max-0.5)		
  # my_max_washes   = log(2)/(exp(my_max_washes_a)/(1+exp(my_max_washes_a)))
  
  my_max_washes   = log(2)/(1/(1+exp(-my_max_washes_a)))
  
  
  ## Uncertainty
  net_half_life_min = 2
  net_half_life_max = 3
  
  ## FOR NON PYR-ONLY NETS THIS WILL RETURN HIGH HALF LIFE
  ## WE RECOMMEND ONLY USING PY-ONLY HALF LIFE UNTIL WE CAN
  ## VALIDATE OTHER NETS 
  wash_decay_rate_a = mup +rhop*(mort_hut-0.5)
  # wash_decay_rate   = log(2)/(exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a)))
  wash_decay_rate   = log(2)/(1/(1+exp(-wash_decay_rate_a)))
  itn_half_life     = wash_decay_rate/my_max_washes*net_halflife
  itn_half_life_max     = wash_decay_rate/my_max_washes*net_half_life_max  
  itn_half_life_min     = wash_decay_rate/my_max_washes*net_half_life_min
  ##No need to re-adjust these anymore
  ##Final Parameter estimates for the transmission model
  ERG_d_ITN0 <- d_ITN0
  ERG_s_ITN0 <- s_ITN0
  ERG_r_ITN0 <- 1-ERG_d_ITN0-ERG_s_ITN0
  
  ## Print out these estimates to a data.frame as the function output
  uncertainty_resistance_params_nets = data.frame(ERG_d_ITN0,ERG_r_ITN0,itn_half_life,
                                                  itn_half_life_min,itn_half_life_max,
                                                  det_hut=xx$det_hut,suc_hut=xx$suc_hut,mort_hut=xx$mort_hut,rep_hut=xx$rep_hut,
                                                  n1n0,kp1,lp1,jp1,
                                                  bioassay_surv = round(seq(0,1,length=101),2))
  
  return(uncertainty_resistance_params_nets)
}

## 20 uncertainty draws
## load the orignial data
data_raw = read.csv("raw data/data_with_calibrated_eir.csv",header=TRUE)

## minimise key data needed for aligning uncertainty
align_dat = data.frame(cluster = data_raw$cluster,
                       net_type = data_raw$Net_Type,
                       res_1 = round(data_raw$res_1,2),
                       res_2 = round(data_raw$res_2,2))

# 
# rand = sample(1:1000,20,replace=FALSE)
# for(i in 1:length(rand)){
#   test = resistance_ITN_pyr_only_params_2_f(product = 0, ## PYRETHROID ONLY LLIN 
#                                            # res = seq(0,1,length=101), ## SURVIVAL IN SUSC BIOASSAY
#                                            shape = "logistic",## for half life with log-log or original logistic
#                                            data_picker_rand = rand[i])[,c(1:3,14)] ## any number from 1 to 1000
#   
#   colnames(test)[4] = "res_1"
#   SET1a = merge(align_dat,test,by="res_1",all.x = TRUE)	
#   
#   colnames(test)[4] = "res_2"
#   SET1b = merge(align_dat,test,by="res_2",all.x = TRUE)	
#   
#   test2 = resistance_ITN_default_params_2_f(product = 1, ## PYRETHROID ONLY LLIN 
#                                             # res = seq(0,1,length=101), ## SURVIVAL IN SUSC BIOASSAY
#                                             shape = "logistic",## for half life with log-log or original logistic
#                                             data_picker_rand = rand[i])[,c(1:3,14)]
#   
#   colnames(test2)[4] = "res_2"
#   SET1c = merge(align_dat,test2,by="res_2",all.x = TRUE)	
#   
#   SET1all = SET1a
#   SET1all$d_ITN0_2017 = ifelse(SET1all$net_type == "Olyset Net",SET1b$ERG_d_ITN0,
#                                ifelse(SET1all$net_type == "PermaNet 2.0",SET1b$ERG_d_ITN0,
#                                       SET1c$ERG_d_ITN0))
#   SET1all$r_ITN0_2017 = ifelse(SET1all$net_type == "Olyset Net",SET1b$ERG_r_ITN0,
#                                ifelse(SET1all$net_type == "PermaNet 2.0",SET1b$ERG_r_ITN0,
#                                       SET1c$ERG_r_ITN0))
#   
#   SET1all$itn_half_2017 = ifelse(SET1all$net_type == "Olyset Net",SET1b$itn_half_life,
#                                  ifelse(SET1all$net_type == "PermaNet 2.0",SET1b$itn_half_life,
#                                         SET1b$itn_half_life)) ## because we are assuming same half life for both nets
#   
#   # write.csv(SET1all,paste0("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/raw data/pbo-net-params_",i,".csv"))  
# }

###########################################################
##
## Next PBO nets

fit1_a <- readRDS("stan model outputs/LLINEUP_ento_deterrence_pbo.rds")
fit1_a_fit <- extract(fit1_a, permuted = TRUE)

fit1_a_med_c = median(fit1_a_fit$c)
fit1_a_med_d = median(fit1_a_fit$d)
fit1_a_med_e = median(fit1_a_fit$e)

fit3_a <- readRDS("stan model outputs/LLINEUP_ento_feeding_0washes_pbo.rds")
fit3_a_fit <- extract(fit3_a, permuted = TRUE)

fit3_a_f <- median(fit3_a_fit$a)
fit3_a_g <- median(fit3_a_fit$b)


resistance_ITN_pbo_params_2_f = function(product, 
                                             # res, ## SURVIVAL IN SUSC BIOASSAY
                                             shape,## for half life with log-log or original logistic
                                             data_picker_rand){ ## random draw from posterior pred  
  
  ## The parameters included here are from the statistical analysis 
  ## following on from Nash et al 2021 
  
  ## PARAMETERS
  
  # give the uncertainty for the log-logistic function
  #Assay to hut mortality conversion - median estimates	
  param_b = LL_b_full[data_picker_rand] 
  param_a = LL_a_full[data_picker_rand] 
  
  # x = res
  
  f_LOG_logistic <- function(x, b, a){
    surv = 1 - (1/(1+((1-x)/b)^(-a)))
    mort = (1-surv)
    
    return(mort)
  }
  
  # mort = 1 - res ## added benefit is relative to mortality in pyr-only nets
  
  f_LOG_logistic_with_benefit <- function(mort, alp1, alp2){
    mort_newnet = 1 / (1 + exp(-alp1 - alp2 * mort))
    # surv = (1-mort)
    return(mort_newnet)
  }
  ## returns hut mortality
  
  hut_mort_LL <- f_LOG_logistic(x = seq(0,1,length=101),#** 
                                param_b, param_a)
  
  hut_mort_LLPBOtemp <- f_LOG_logistic_with_benefit(mort = hut_mort_LL, 
                                                    alp1 = fitbene_1[data_picker_rand],
                                                    alp2 = fitbene_2[data_picker_rand])
  
  
  #specify whichever net is used in the RCT
  #this will determine what the mortality is in the hut trial
  mort_huta = if(product==0) hut_mort_LL else if(product==1) hut_mort_LLPBOtemp 
  mort_hut = mort_huta
  
  hut_surv = 1 - mort_hut
  # ff, ff1, and ff2 all match up for the associations for net det and fed
  # when using the generic parameters as here
  
  #Decay in insecticide non-PBO net		
  mup =	-2.429#NEW 1.812767 ###hlf_mu# #mu_p ## sample(test$a,size=1) ##array(c(rep(-2.429,3),rep(-2.984025,3),rep(-1.866,3)),c(3,3)) ## ... gam.medians[9]
  rhop =-3.007 #NEW -2.581591 ###hlf_rho# #rho_p ##  sample(test$b,size=1)	##array(c(rep(-3.007,3),rep(-3.74,3),rep(-2.295,3)),c(3,3)) ## ... gam.medians[10]
  
  ## The maximum successful feeding probability per feeding attempt 
  ## (feeding and not dying) in the absence of interventions 
  kp0=0.699 ## derived from Lines et al 1987 and Curtis et al 1990 
  
  ## The half-life of the net relative to it's capacity to kill mosquitoes
  ## with the insecticide active ingredient (a pyrethroid) when there is 
  ## no resistance in mosquitoes. 
  net_halflife=2.64
  
  ## Now we work through the probability steps to determine the key input parameters for the model
  ## These probability relationships are determined by Rebecca Nash, Ben and Tom see email notes above
  
  fit1_a_med_c = median(fit1_a_fit$c)
  fit1_a_med_d = median(fit1_a_fit$d)
  fit1_a_med_e = median(fit1_a_fit$e)
  
  ## This is association with hut survival
  det_hut = fit1_a_med_e * (exp(fit1_a_med_d * (1 - exp(fit1_a_med_c * hut_surv)) / fit1_a_med_c))
  
  fit3_a_f <- median(fit3_a_fit$a)
  fit3_a_g <- median(fit3_a_fit$b)
  
  ## This is association with hut survival
  suc_hut = (1 - (exp(fit3_a_g * (1 - exp(fit3_a_f * hut_surv))/fit3_a_f)))
  
  ## This is association with hut survival
  rep_hut = (1 - suc_hut - mort_hut)
  
  xx = data.frame(hut_surv,mort_hut,suc_hut,rep_hut,det_hut)
  ## Combine to estimate the 3 key probable outcomes of feeding attempts
  ## Here we adjust for those mosquitoes not entering treated huts (determined by deterrence)
  n1n0 = 1-xx$det_hut
  kp1  = n1n0*xx$suc_hut
  lp1  = n1n0*xx$mort_hut
  jp1  = n1n0*xx$rep_hut+(1-n1n0)
  
  kp1 = ifelse(kp1 > kp0,kp0,kp1) ## Capping impact so max feeding is no bigger than assumed
  # # max feeding for no interventions (kp0 = 0.699, Griffin et al 2010)
  # ## (time = 0 time steps after net implementation)
  # 
  # kp0 = 1
  # 
  r_ITN0  = (1-kp1/kp0)*(jp1/(lp1+jp1))	#probability of repeating behaviour
  d_ITN0  = (1-kp1/kp0)*(lp1/(lp1+jp1))	#probability of dying with an encounter with ITN
  s_ITN0  = 1-d_ITN0-r_ITN0             #probability of successfully feeding (surviving and feeding)
  
  # plot(r_ITN0 ~ c(1-mort),ylim=c(0,1),xlim=c(0,1),xlab = "Susc bioassay survial",type="l",col="orange") 
  
  ## Repeat these to determin the maximum and minimum effects which combine to help determine ITN half life
  ## We will stick to pyr-params for half life and update in 2023 once new data are available
  mort_maxA   = if (shape=="log-log") f_LOG_logistic(x = 0,#** this is surv i.e. mort max when surv =0
                                                     param_b, param_a) else if (shape=="logistic") 1/(1+exp(param_a*param_b-param_a*(1)))
  
  

  mort_max = mort_maxA#if(product==0) mort_maxA else if(product==1) PBO_benefitA else if(product==2) G2_benefitA else if(product==3) G2_benefitB
  
  
  #{halflife}
  my_max_washes_a = mup +rhop*(mort_max-0.5)		
  # my_max_washes   = log(2)/(exp(my_max_washes_a)/(1+exp(my_max_washes_a)))
  
  my_max_washes   = log(2)/(1/(1+exp(-my_max_washes_a)))
  
  
  ## Uncertainty
  net_half_life_min = 2
  net_half_life_max = 3
  
  ## FOR NON PYR-ONLY NETS THIS WILL RETURN HIGH HALF LIFE
  ## WE RECOMMEND ONLY USING PY-ONLY HALF LIFE UNTIL WE CAN
  ## VALIDATE OTHER NETS 
  wash_decay_rate_a = mup +rhop*(mort_hut-0.5)
  # wash_decay_rate   = log(2)/(exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a)))
  wash_decay_rate   = log(2)/(1/(1+exp(-wash_decay_rate_a)))
  itn_half_life     = wash_decay_rate/my_max_washes*net_halflife
  itn_half_life_max     = wash_decay_rate/my_max_washes*net_half_life_max  
  itn_half_life_min     = wash_decay_rate/my_max_washes*net_half_life_min
  ##No need to re-adjust these anymore
  ##Final Parameter estimates for the transmission model
  ERG_d_ITN0 <- d_ITN0
  ERG_s_ITN0 <- s_ITN0
  ERG_r_ITN0 <- 1-ERG_d_ITN0-ERG_s_ITN0
  
  ## Print out these estimates to a data.frame as the function output
  uncertainty_resistance_params_nets = data.frame(ERG_d_ITN0,ERG_r_ITN0,itn_half_life,
                                                  itn_half_life_min,itn_half_life_max,
                                                  det_hut=xx$det_hut,suc_hut=xx$suc_hut,mort_hut=xx$mort_hut,rep_hut=xx$rep_hut,
                                                  n1n0,kp1,lp1,jp1,
                                                  bioassay_surv = round(seq(0,1,length=101),2))
  
  return(uncertainty_resistance_params_nets)
}

## 20 uncertainty draws
## load the orignial data
data_raw = read.csv("raw data/data_with_calibrated_eir.csv",header=TRUE)

## minimise key data needed for aligning uncertainty
align_dat = data.frame(cluster = data_raw$cluster,
                       net_type = data_raw$Net_Type,
                       res_1 = round(data_raw$res_1,2),
                       res_2 = round(data_raw$res_2,2))


rand = sample(1:1000,20,replace=FALSE)
for(i in 1:length(rand)){
  test = resistance_ITN_pyr_only_params_2_f(product = 0, ## PYRETHROID ONLY LLIN 
                                           # res = seq(0,1,length=101), ## SURVIVAL IN SUSC BIOASSAY
                                           shape = "logistic",## for half life with log-log or original logistic
                                           data_picker_rand = rand[i])[,c(1:3,14)] ## any number from 1 to 1000
  
  colnames(test)[4] = "res_1"
  SET1a = merge(align_dat,test,by="res_1",all.x = TRUE)	
  
  colnames(test)[4] = "res_2"
  SET1b = merge(align_dat,test,by="res_2",all.x = TRUE)	
  
  test2 = resistance_ITN_pbo_params_2_f(product = 1, ## PYRETHROID PBO LLIN 
                                            # res = seq(0,1,length=101), ## SURVIVAL IN SUSC BIOASSAY
                                            shape = "logistic",## for half life with log-log or original logistic
                                            data_picker_rand = rand[i])[,c(1:3,14)]
  
  colnames(test2)[4] = "res_2"
  SET1c = merge(align_dat,test2,by="res_2",all.x = TRUE)	
  
  SET1all = SET1a
  SET1all$d_ITN0_2017 = ifelse(SET1all$net_type == "Olyset Net",SET1b$ERG_d_ITN0,
                               ifelse(SET1all$net_type == "PermaNet 2.0",SET1b$ERG_d_ITN0,
                                      SET1c$ERG_d_ITN0))
  SET1all$r_ITN0_2017 = ifelse(SET1all$net_type == "Olyset Net",SET1b$ERG_r_ITN0,
                               ifelse(SET1all$net_type == "PermaNet 2.0",SET1b$ERG_r_ITN0,
                                      SET1c$ERG_r_ITN0))
  
  SET1all$itn_half_2017 = ifelse(SET1all$net_type == "Olyset Net",SET1b$itn_half_life,
                                 ifelse(SET1all$net_type == "PermaNet 2.0",SET1b$itn_half_life,
                                        SET1b$itn_half_life)) ## because we are assuming same half life for both nets
  
  write.csv(SET1all,paste0("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/raw data/specific-net-params_",i,".csv"))  
}

store_dn1 = store_rn1 = store_hl1 = array(dim=c(104,20))
store_dn2 = store_rn2 = store_hl2 = array(dim=c(104,20))
for(i in 1:20){
  dataspec = read.csv(paste0("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP/raw data/specific-net-params_",i,".csv"),header=TRUE)
  dataspec = dataspec[order(dataspec$cluster), ]
  store_dn1[,i] = dataspec[,6]
  store_rn1[,i] = dataspec[,7]
  store_hl1[,i] = dataspec[,8]
  store_dn2[,i] = dataspec[,9]
  store_rn2[,i] = dataspec[,10]
  store_hl2[,i] = dataspec[,11]
  
}
spec_mean_params = data.frame(cluster = dataspec$cluster, 
                              dn1 = rowMeans(store_dn1),
                              rn1 = rowMeans(store_rn1),
                              hl1 = rowMeans(store_hl1),
                              dn2 = rowMeans(store_dn2),
                              rn2 = rowMeans(store_rn2),
                              hl2 = rowMeans(store_hl2))
write.csv(spec_mean_params,"raw data/specific-net-paramsMEANS.csv")




# test = resistance_ITN_default_params_2_f(product = 1, ## PYRETHROID ONLY LLIN 
#                                          # res = seq(0,1,length=101), ## SURVIVAL IN SUSC BIOASSAY
#                                          shape = "logistic",## for half life with log-log or original logistic
#                                          data_picker_rand = 45) ## any number from 1 to 1000
# head(test)



#################################
##
## Global estimates for all of us to use

vec = 1:1000
matrix_dn0 = matrix_rn0 = matrix_halflife = array(dim=c(nrow(test),length(vec)))
for(i in 1:1000){
  test = resistance_ITN_default_params_2_f(product = 0, ## PYRETHROID ONLY LLIN 
                                           # res = seq(0,1,length=101), ## SURVIVAL IN SUSC BIOASSAY
                                           shape = "logistic",## for half life with log-log or original logistic
                                           data_picker_rand = i) ## any number from 1 to 1000
  matrix_dn0[,i] = test$ERG_d_ITN0
  matrix_rn0[,i] = test$ERG_r_ITN0
  matrix_halflife[,i] = test$itn_half_life
  
}

dn0MEAN = rowMeans(matrix_dn0)
rn0MEAN = rowMeans(matrix_rn0)
halflifeMEAN = rowMeans(matrix_halflife)

dn0 = rn0 = gamman = array(dim=c(nrow(test),3))
for(j in 1:nrow(test)){
  dn0[j,] = c(as.numeric(quantile(matrix_dn0[j,],c(0.1,0.5,0.9))))
  rn0[j,] = c(as.numeric(quantile(matrix_rn0[j,],c(0.1,0.5,0.9))))
  gamman[j,] = c(as.numeric(quantile(matrix_halflife[j,],c(0.1,0.5,0.9))))
}

pyrethroidOnlyNets = data.frame(dn0_lo10 = dn0[,1],dn0_med = dn0[,2],dn0_up90 = dn0[,3],
                                rn0_lo10 = rn0[,1],rn0_med = rn0[,2],rn0_up90 = rn0[,3],
                                gamman_lo10 = gamman[,1],gamman_med = gamman[,2],gamman_up90 = gamman[,3])
pyrethroidOnlyNets$bioassay_mortality = seq(1,0,length=101)
head(pyrethroidOnlyNets)

## pyrethroid PBO

vec = 1:1000
matrix_dn0 = matrix_rn0 = matrix_halflife = array(dim=c(nrow(test),length(vec)))
for(i in 1:1000){
  test = resistance_ITN_default_params_2_f(product = 1, ## PYRETHROID ONLY LLIN 
                                           # res = seq(0,1,length=101), ## SURVIVAL IN SUSC BIOASSAY
                                           shape = "logistic",## for half life with log-log or original logistic
                                           data_picker_rand = i) ## any number from 1 to 1000
  matrix_dn0[,i] = test$ERG_d_ITN0
  matrix_rn0[,i] = test$ERG_r_ITN0
  matrix_halflife[,i] = test$itn_half_life
  
}

dn0MEAN = rowMeans(matrix_dn0)
rn0MEAN = rowMeans(matrix_rn0)
halflifeMEAN = rowMeans(matrix_halflife)

dn0 = rn0 = gamman = array(dim=c(nrow(test),3))
for(j in 1:nrow(test)){
  dn0[j,] = c(as.numeric(quantile(matrix_dn0[j,],c(0.1,0.5,0.9))))
  rn0[j,] = c(as.numeric(quantile(matrix_rn0[j,],c(0.1,0.5,0.9))))
  gamman[j,] = c(as.numeric(quantile(matrix_halflife[j,],c(0.1,0.5,0.9))))
}

pyrethroidPBONets = data.frame(dn0_lo10 = dn0[,1],dn0_med = dn0[,2],dn0_up90 = dn0[,3],
                               rn0_lo10 = rn0[,1],rn0_med = rn0[,2],rn0_up90 = rn0[,3],
                               gamman_lo10 = gamman[,1],gamman_med = gamman[,2],gamman_up90 = gamman[,3])
pyrethroidPBONets$bioassay_mortality = seq(1,0,length=101)
head(pyrethroidPBONets)

