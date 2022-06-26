# Plots with combined fits for ratio between total collected from control and ITN huts, successful blood feeding and exiting unfed vs survival in EHT

## Using the Mosquito-net-parameters Repository

## Part 1 for a supplementary material


library(rstan)
options(mc.cores=4)
rstan_options(auto_write = TRUE)
##############################################################################
##
## Figures for Prioritisation paper
## Oct 2019
##
###############################################################################


###################################################
##
## Consistent with the standard and PBO-nets defined in 
## Churcher et al 2016 / Nash et al 2021 / Sherrard-Smith et al 2021
##
#########################################################################


library(rstan)

# For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
# To avoid recompilation of unchanged Stan programs, we recommend calling
rstan_options(auto_write = TRUE)
# setwd("E:/Mosquito net parameters")
setwd("C:/Users/esherrar/Documents/Rprojects/Key-RCT-metrics-LLINEUP")

###########################################
##
## Part 1:
##

mt = read.csv("raw data/mortality_data_from_Nash2021.csv",header=TRUE)

## remove Ifakara huts
mt = subset(mt, mt$Hut != "Ifakara")
library(dplyr)
selectiona <- c("Olyset Net", "Interceptor", 
                "PermaNet 2.0","PermaNet", 
                "MiraNet",
                "DuraNet","Duranet",
                "Yahe",
                "MAGNet", "DawaPlus 2.0", "Yorkool",
                "Panda Net 2.0","PandaNet 2.0",
                "Royal Sentry")
mt <- mt %>% filter(Net%in%selectiona)

data_list_MT = list(S = nrow(mt),
                    N_b = mt$Bioassay_N_tested,
                    N_h = mt$N_total,
                    X_b = mt$Bioassay_N_tested - mt$Bioassay_N_dead,
                    X_h = mt$N_total - mt$N_dead,
                    nsite = length(unique(mt$Site)),
                    site = as.numeric(factor(mt$Site)),
                    S_test = 101,
                    theta_b_test = seq(0,1,length=101))

stan_base <- rstan::stan(file="R code/stan models/A1 Log-logistic Model.stan", 
                         data=data_list_MT, 
                         warmup=1000,
                         control = list(adapt_delta = 0.8,
                                        max_treedepth = 20),
                         iter=2000, chains=4)
base <- rstan::extract(stan_base)

median(base$b)
median(base$a)

x = seq(0,1,length=101)
surv = 1 - (1/(1+((1-x)/median(base$b))^(-median(base$a))))

surv_low = 1 - (1/(1+((1-x)/quantile(base$b,0.1))^(-quantile(base$a,0.1))))
surv_upp = 1 - (1/(1+((1-x)/quantile(base$b,0.9))^(-quantile(base$a,0.9))))

par(mfrow=c(2,2))
par(mar=c(5,6,3,3))

plot(surv ~ x,type="l",
     xlab = "Survival at discriminatory dose bioassay (%)",
     ylab = "Survival in experimental hut trial (%)",
     xaxt="n",yaxt="n",ylim=c(0,1),xlim=c(0,1),
     cex.lab=1.2,cex.axis=1.2)
axis(1,at=seq(0,1,by=0.2),labels=seq(0,100,by=20),cex.lab=1.2,cex.axis=1.2)
axis(2,las=2,at=seq(0,1,by=0.2),labels=seq(0,100,by=20),cex.lab=1.2,cex.axis=1.2)

dd = as.numeric(quantile(base$b,0.025))
ff = as.numeric(quantile(base$b,0.975))
which(base$b < dd+0.05 & base$b > dd)
which(base$b > ff & base$b < ff+0.05)
# for(i in c(305,2081)){
#   surv = 1 - (1/(1+((1-x)/base$b[i])^(-base$a[i])))
#   lines(surv ~ x,col="black")
# }
surv_upp = surv = 1 - (1/(1+((1-x)/base$b[which(base$b < dd+0.05 & base$b > dd)[1]])^(-median(base$a))))
surv_low = surv = 1 - (1/(1+((1-x)/base$b[which(base$b > ff & base$b < ff+0.05)[1]])^(-(median(base$a)))))
polygon(c(x,rev(x)),c(surv_low,rev(surv_upp)),col=adegenet::transp("yellow",0.4),border=NA)

points(c((mt$Bioassay_N_tested - mt$Bioassay_N_dead)/mt$Bioassay_N_tested),
       c((mt$N_total - mt$N_dead)/mt$N_total),col=adegenet::transp("orange",0.7),pch=19,
       cex = c(log(50*(mt$N_total/max(mt$N_total))+1)))

###########################################
##
## Part 2:
##


##
## Drawing the figure
##
# par(mfrow=c(1,1))
# par(mar=c(5,6,3,3))



############################
## 
## Add PBO line

#Global
# pbo_dat = read.csv("data/PPBO_Studies_v2.csv")
pbo_dat = read.csv("raw data/PPBO_Studies_v2.csv")
dim(pbo_dat)
#*MAKE SURE ORDERED SO THAT llin ARE NEXT TO PPBO nets
pbo_dat$Location = pbo_dat$Location
pbo_dat$pbo_deterred  = (pbo_dat$N_mosq_C - pbo_dat$N_mosq)/pbo_dat$N_mosq_C
pbo_dat$pbo_24hr_dead = (pbo_dat$N_dead_24/pbo_dat$N_mosq)
pbo_dat$pbo_mosq_survival = 1 - pbo_dat$pbo_24hr_dead
pbo_dat$pbo_exiting   = (pbo_dat$N_exited/pbo_dat$N_mosq) * (1 - pbo_dat$pbo_24hr_dead)
pbo_dat$pbo_successful = 1 - pbo_dat$pbo_24hr_dead - pbo_dat$pbo_exiting

pbo_dat$x = pbo_dat$Bioassay_mort_PERMETHRIN
pbo_dat$survival = 100 - pbo_dat$bioassay_mort_PERMETHRIN
pbo_dat$x2 = pbo_dat$bioassay_mort_DELTAMETHRIN
pbo_dat$survival2 = 100 - pbo_dat$bioassay_mort_DELTAMETHRIN




##Relationship 2 PBO benefit on top of Normal LN 
##LOGISTIC BENEFIT
d_t = c(pbo_dat$N_dead_24[pbo_dat$Intervention_type == "PPBO" & pbo_dat$N_washes == 0]) # number of mosquitoes dying IRS HUTS
n_t = c(pbo_dat$N_mosq[pbo_dat$Intervention_type == "PPBO" & pbo_dat$N_washes == 0]) # Number of mosquitoes entering IRS huts
x1 = c(pbo_dat$N_dead_24[pbo_dat$Intervention_type == "LLIN"  & pbo_dat$N_washes == 0])/
  c(pbo_dat$N_mosq[pbo_dat$Intervention_type == "LLIN"  & pbo_dat$N_washes == 0])

data_PBO = data.frame(d_t,n_t,x1)
data_PBO = data_PBO[complete.cases(data_PBO),]
data_PBO$Proportion_dead = data_PBO$d_t/data_PBO$n_t
x = seq(1,0,length=100)

data_list_PBO = list(N = nrow(data_PBO),
                     n_t = data_PBO$n_t,
                     d_t = data_PBO$d_t,
                     x = data_PBO$x1)

stan_base <- stan(file="R code/stan models/binomial_fit_nets.stan", 
                  data=data_list_PBO, 
                  warmup=1000,
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 20),
                  iter=2000, chains=4)
base <- rstan::extract(stan_base)

# saveRDS(stan_base,"ento_pbo_benefit.rds")
# traceplot(stan_base)
# stan_diag(stan_base,
#           information = c("sample","stepsize", "treedepth","divergence"),
#           chain = 0)
# stan_hist(stan_base, include = TRUE, unconstrain = FALSE,
# inc_warmup = FALSE)
# stan_rhat(stan_base,bins = 10)
mean(base$alpha2)
mean(base$alpha1)

# saveRDS(stan_base,"data/ento_pbo_benefit.RDS")

preds_pbo = array(dim=c(100,4000))
for(i in 1:4000){
  preds_pbo[,i] = 1 / (1 + exp(-base$alpha1[i] - base$alpha2[i]*x))
}
median_pbo_pred = 1 / (1 + exp(-quantile(base$alpha1,0.5) - quantile(base$alpha2,0.5)*x))
upper_pbo_pred = 1 / (1 + exp(-quantile(base$alpha1,0.975) - quantile(base$alpha2,0.975)*x))
lower_pbo_pred = 1 / (1 + exp(-quantile(base$alpha1,0.025) - quantile(base$alpha2,0.025)*x))


prop_killed_pbo = data_list_PBO$d_t/data_list_PBO$n_t


plot(prop_killed_pbo ~ data_list_PBO$x,
     ylab="Mosquito mortality pyrethroid-PBO ITN (%)",
     ylim=c(0,1),xlim=c(0,1),pch="",yaxt="n",xaxt="n",
     xlab = "Mosquito mortality pyrethroid-only LLIN (%)",cex.lab=1.2,cex.axis=1.2)
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.2)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.2)
abline(0,1,lty=2,col="blue")

points(prop_killed_pbo ~ data_list_PBO$x,cex=1.6,pch=19,col="grey")
polygon(c(x,rev(x)),c(upper_pbo_pred,rev(lower_pbo_pred)),col=adegenet::transp("grey",0.4),border=NA)


lines(median_pbo_pred ~ x,col="black",lwd=2)

##LOGISTIC BENEFIT - PermaNets & Olysets
# pbo_dat = read.csv("E:/ARCHIVE/Rebecca_Expt hut analysis/PPBO_Studies_v2.csv")
pbo_dat = read.csv("raw data/pbos.csv")
# pbo_dat = read.csv("Q:/RProjects/Mosquito-Net-Parameters/pbos.csv")
dim(pbo_dat)

##Relationship 2 PBO benefit on top of Normal LN 
##LOGISTIC BENEFIT - global
d_t = c(pbo_dat$N_dead_24_pbo[pbo_dat$N_washes_pbo == 0]) # number of mosquitoes dying IRS HUTS
n_t = c(pbo_dat$N_mosq_pbo[pbo_dat$N_washes_pbo == 0]) # Number of mosquitoes entering IRS huts
x1 = c(pbo_dat$N_dead_24[pbo_dat$N_washes == 0])/
  c(pbo_dat$N_mosq[pbo_dat$N_washes == 0])


data_PBO = data.frame(d_t,n_t,x1)
data_PBO = data_PBO[complete.cases(data_PBO),]
data_PBO$Proportion_dead = data_PBO$d_t/data_PBO$n_t
x = seq(1,0,length=100)

data_list_PBO = list(N = nrow(data_PBO),
                     n_t = data_PBO$n_t,
                     d_t = data_PBO$d_t,
                     x = data_PBO$x1)

stan_base <- stan(file="R code/stan models/binomial_fit_nets.stan", 
                  data=data_list_PBO, 
                  warmup=1000,
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 20),
                  iter=2000, chains=4)
base <- rstan::extract(stan_base)


preds_pbop = array(dim=c(100,4000))
for(i in 1:4000){
  preds_pbop[,i] = 1 / (1 + exp(-base$alpha1[i] - base$alpha2[i]*x))
}
median_pbo_predp = 1 / (1 + exp(-quantile(base$alpha1,0.5) - quantile(base$alpha2,0.5)*x))
upper_pbo_predp = 1 / (1 + exp(-quantile(base$alpha1,0.975) - quantile(base$alpha2,0.975)*x))
lower_pbo_predp = 1 / (1 + exp(-quantile(base$alpha1,0.025) - quantile(base$alpha2,0.025)*x))

prop_killed_pbop = data_list_PBO$d_t/data_list_PBO$n_t

points(prop_killed_pbop ~ data_list_PBO$x,cex=1.2,pch=8,col="black")
polygon(c(x,rev(x)),c(upper_pbo_predp,rev(lower_pbo_predp)),col=adegenet::transp("aquamarine3",0.4),border=NA)
lines(median_pbo_predp ~ x,col="darkgreen",lwd=2,lty=2)




##LOGISTIC BENEFIT - PermaNets & Olysets
# pbo_dat = read.csv("E:/ARCHIVE/Rebecca_Expt hut analysis/PPBO_Studies_v2.csv")
# pbo_dat = read.csv("Q:/RProjects/Mosquito-Net-Parameters/pbos.csv")
pbo_dat = read.csv("raw data/pbos.csv")
dim(pbo_dat)


##Relationship 2 PBO benefit on top of Normal LN 
##LOGISTIC BENEFIT - permaNet only
d_t = c(pbo_dat$N_dead_24_pbo[pbo_dat$N_washes_pbo == 0 & pbo_dat$Net_type_pbo == "PermaNet_3.0"]) # number of mosquitoes dying IRS HUTS
n_t = c(pbo_dat$N_mosq_pbo[pbo_dat$N_washes_pbo == 0& pbo_dat$Net_type_pbo == "PermaNet_3.0"]) # Number of mosquitoes entering IRS huts
x1 = c(pbo_dat$N_dead_24[pbo_dat$N_washes == 0 & pbo_dat$Net_type_pbo == "PermaNet_3.0"])/
  c(pbo_dat$N_mosq[pbo_dat$N_washes == 0 & pbo_dat$Net_type_pbo == "PermaNet_3.0"])


data_PBO = data.frame(d_t,n_t,x1)
data_PBO = data_PBO[complete.cases(data_PBO),]
data_PBO$Proportion_dead = data_PBO$d_t/data_PBO$n_t
x = seq(1,0,length=100)

data_list_PBO = list(N = nrow(data_PBO),
                     n_t = data_PBO$n_t,
                     d_t = data_PBO$d_t,
                     x = data_PBO$x1)

stan_base <- stan(file="R code/stan models/binomial_fit_nets.stan", 
                  data=data_list_PBO, 
                  warmup=1000,
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 20),
                  iter=2000, chains=4)
base <- rstan::extract(stan_base)


preds_pbopr = array(dim=c(100,4000))
for(i in 1:4000){
  preds_pbopr[,i] = 1 / (1 + exp(-base$alpha1[i] - base$alpha2[i]*x))
}
median_pbo_predpr = 1 / (1 + exp(-quantile(base$alpha1,0.5) - quantile(base$alpha2,0.5)*x))
upper_pbo_predpr = 1 / (1 + exp(-quantile(base$alpha1,0.975) - quantile(base$alpha2,0.975)*x))
lower_pbo_predpr = 1 / (1 + exp(-quantile(base$alpha1,0.025) - quantile(base$alpha2,0.025)*x))

prop_killed_pbopr = data_list_PBO$d_t/data_list_PBO$n_t

points(prop_killed_pbopr ~ data_list_PBO$x,cex=1.6,pch=17,col="red")
polygon(c(x,rev(x)),c(upper_pbo_predpr,rev(lower_pbo_predpr)),col=adegenet::transp("darkred",0.4),border=NA)
lines(median_pbo_predpr ~ x,col="red",lwd=2,lty=2)

legend("bottomright",legend=c("Any pyrethroid-PBO ITNs",
                              "RCT pyrethroid-PBO ITNs",
                              "PermaNet 3.0"),
       pch=c(19,8,17),col=c("grey","blue","red"),lty=c(1,2,3),bty="n")


################################
##
## Part 3: Probable outcomes of feeding attempts in huts
##
#################################

##################################################################
##
## Read in an updated file for ALL net studies 
## to work out the impacts from Interceptor G2, Pyrethroid only, or Pyr-PBO

## ** Still need info from Ben Amoah for the G2
# df <- read.csv("data/3AC_data_with_72_final_And_new_G2.csv") # mac
df <- read.csv("raw data/All Hut Studies April 2022 public.csv")
# Remove Ifakara and those without totals or number dead
library(dplyr)
df2 <- df %>%
  filter(N_dead_comb!="NA" & Totals_known!="0" & N_total!="NA" & Intervention!="Untreated" & Hut!="Ifakara")
# all the data

# # select the nets with 0 washes
selection0 <- c(0)
dfAll <- df2 %>% filter(No_washes%in%selection0)
# dfAll = df2
# unique(dfAll$Intervention)
## Pyrethroid only nets included are;
## (all WHO Recommended according to Okumu & Finda 2021 book chapter)
# Olyset Net
# Interceptor
# Royal Sentry     ## No data
# Royal Sentry 2.0 ## No data
# PermaNet 2.0
# DuraNet LLIN
# MiraNet      # No data
# MAGNet
# YaheLN       # No data
# DawaPlus 2.0 LLIN
# SafeNet      # No data
# Yorkool LN
# Panda Net2.0 LLIN      # No data

selectionb <- c("Olyset Net", "PermaNet 2.0") ## For pyrethroid-only ITN
dfpyr <- dfAll %>% filter(Intervention%in%selectionb)

## Pyrethroid-PBO nets included are;
## (all WHO Recommended according to Okumu & Finda 2021 book chapter)

# Olyset Plus
# PermaNet 3.0
# Veeralin
# Dawa Plus 3.0 (ALSO 4.0)
# Tsara Boost # No data
selectionc <- c("Olyset Plus", "PermaNet 3.0") ## For pyr-PBO ITN
dfpbo <- dfAll %>% filter(Intervention%in%selectionc)

add_pyr_baselines_f = function(df2){
  
  nsite <- max(df2$Site_number)
  site <- df2$Site_number
  
  for(i in 1:length(df2$R_ID)){
    df2$pyr_base_count[i] = round(mean(df2$N_total[df2$LLIN_Gen == 1 & df2$R_ID == df2$R_ID[i]],na.rm = TRUE),0)
    df2$pyr_base_killed[i] = round(mean(df2$N_dead_comb[df2$LLIN_Gen == 1 & df2$R_ID == df2$R_ID[i]],na.rm = TRUE),0)
    
  }
  
  df2$pyr_base_countN = ifelse(df2$LLIN_Gen == 1,df2$N_total,df2$pyr_base_count)  
  df2$pyr_base_killN = ifelse(df2$LLIN_Gen == 1,df2$N_dead_comb,df2$pyr_base_killed)
  
  df3 <- df2[complete.cases(df2$pyr_base_countN), ] 
  data2 <- df3[complete.cases(df3$pyr_base_killN), ] 
  
  return(data2)
}

dfPYRtest = add_pyr_baselines_f(dfpyr)
dfPBOtest = add_pyr_baselines_f(dfpbo)


setup_inputs_ALL = function(df2){
  # Check which studies/trials included in deterrence analysis
  det_IDs <- unique(df2$Orig_ID) # trial (n=34)
  det_refs <- unique(df2$Library_ref) # publication (n=28)
  
  
  # Prep data
  S <- as.numeric(length(df2$N_total))
  X_survive <- df2$N_total-df2$N_dead_comb
  N_caught <- df2$N_total
  nsite <- max(df2$Site_number)
  site <- df2$Site_number
  X_control <- df2$N_Control
  S_sim <- 300
  theta_sim <- seq(0, 1, length.out = S_sim)
  
  # Need to pass it the data:
  data_stan <- list(S=S, 
                    X_survive=X_survive,
                    N_caught=N_caught,
                    nsite=nsite,
                    site=site,
                    X_control=X_control,
                    S_sim=S_sim,
                    theta_sim=theta_sim)
  
  # df2 <- df2 %>%
  #   mutate(X_survive=N_total-N_dead_comb,
  #          frac_surv=X_survive / N_caught,
  #          change=(X_control-N_caught) / X_control)
  # 
  
  df2 <- df2 %>%
    mutate(X_survive=N_total-N_dead_comb,
           frac_surv=X_survive / N_caught,
           Perc_surv=frac_surv*100,
           change=(X_control-N_caught) / X_control,
           Perc_int_vs_control=change*100)
  
  return(data_stan)
  
}


setup_inputs_pyr = function(df2){
  # Check which studies/trials included in deterrence analysis
  det_IDs <- unique(df2$Orig_ID) # trial (n=34)
  det_refs <- unique(df2$Library_ref) # publication (n=28)
  
  
  # Prep data
  S <- as.numeric(length(df2$pyr_base_countN))
  X_survive <- df2$pyr_base_countN-df2$pyr_base_killN
  N_caught <- df2$pyr_base_countN
  nsite <- max(df2$Site_number)
  site <- df2$Site_number
  X_control <- df2$N_Control
  S_sim <- 300
  theta_sim <- seq(0, 1, length.out = S_sim)
  
  # Need to pass it the data:
  data_stan <- list(S=S, 
                    X_survive=X_survive,
                    N_caught=N_caught,
                    nsite=nsite,
                    site=site,
                    X_control=X_control,
                    S_sim=S_sim,
                    theta_sim=theta_sim)
  
  # df2 <- df2 %>%
  #   mutate(X_survive=N_total-N_dead_comb,
  #          frac_surv=X_survive / N_caught,
  #          change=(X_control-N_caught) / X_control)
  # 
  
  df2 <- df2 %>%
    mutate(X_survive=N_total-N_dead_comb,
           frac_surv=X_survive / N_caught,
           Perc_surv=frac_surv*100,
           change=(X_control-N_caught) / X_control,
           Perc_int_vs_control=change*100)
  
  return(data_stan)
  
}


# Compile model
# full_model<- stan_model("R code/stan models/Deter_model_NoRE_DETER_opt1.stan") # flex params
# 
# run1_f = function(df2,dataset){
#   data_stan = setup_inputs_ALL(df2)
#   
#   # Run model
#   fit_full <- sampling(full_model, data=data_stan, iter=2000, chains=4)
#   # params_a = extract(fit_full)
#   
#   # save fit
#   saveRDS(fit_full, paste0("stan model outputs/LLINEUP_ento_deterrence_",dataset,".rds"))
# }
# 
# run1_f(dfpyr,"pyr")
# run1_f(dfpbo,"pbo")
# 

fit1_b <- readRDS("stan model outputs/LLINEUP_ento_deterrence_pyr.rds")
fit1_c <- readRDS("stan model outputs/LLINEUP_ento_deterrence_pbo.rds")


print(fit1_b, pars=c("c","d","e", "kappa","P_control")) 
print(fit1_c, pars=c("c","d","e", "kappa","P_control")) 

# Extract ratio
plot_ratio_prob <- function(P_param, fit){
  P <- rstan::extract(fit, P_param)[[1]]
  P_med <- apply(P, 2, median)
  return(P_med*100)
}

P_ratio_b <- plot_ratio_prob("P_ratio", fit1_b)
P_ratio_c <- plot_ratio_prob("P_ratio", fit1_c)

create_output_f = function(fit_full){
  # Extract upper credible intervals
  plot_upper_funct <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_med <- apply(P, 2, median)
    P_upper <- apply(P, 2, function(x) quantile(x,0.025))
    return(P_upper)
  }
  
  P_ratio_upper <- plot_upper_funct("P_ratio", fit_full)
  P_det_upper <- plot_upper_funct("P_deter_sim", fit_full)
  P_inside_upper <- plot_upper_funct("P_inside_sim", fit_full)
  
  
  plot_med_funct <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_med <- apply(P, 2, median)
    return(P_med)
  }
  # f_plot_prob <- function(P_param, theta_sim, fit){
  #   P <- rstan::extract(fit_a, P_param)[[1]]
  #   P_med <- apply(P, 2, median)
  P_ratio_med <- plot_med_funct("P_ratio", fit_full)
  P_det_med <- plot_med_funct("P_deter_sim", fit_full)
  P_inside_med <- plot_med_funct("P_inside_sim", fit_full)
  
  
  # Extract lower credible intervals
  plot_lower_funct <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_med <- apply(P, 2, median)
    P_lower <- apply(P, 2, function(x) quantile(x,0.975))
    return(P_lower) 
  }
  
  P_ratio_lower <- plot_lower_funct("P_ratio", fit_full)
  P_det_lower <- plot_lower_funct("P_deter_sim", fit_full)
  P_inside_lower <- plot_lower_funct("P_inside_sim", fit_full)
  
  
  a_theta <- seq(0,1,length=300)
  my_surv <- a_theta*100
  
  cred_int_data <- data.frame(P_ratio_med, P_ratio_upper, P_ratio_lower,my_surv)
  det_cred <- data.frame(P_det_med, P_det_upper, P_det_lower,a_theta)
  inside_cred <- data.frame(P_inside_med, P_inside_upper, P_inside_lower,a_theta)
  
  return(list(cred_int_data,
              det_cred,
              inside_cred))  
}

mod_data_b = create_output_f(fit1_b)
mod_data_c = create_output_f(fit1_c)

S_sim <- 300
theta_sim <- seq(0, 1, length.out = S_sim)

plot(mod_data_b[[2]][,1] ~ theta_sim,type = "l",
     ylim=c(0,1),xlim=c(0,1),yaxt="n",
     col="grey45",
     main="",
     ylab="Probability of being deterred",
     xlab="Probability of hut survival",
     cex.lab=1.2,cex.axis=1.2)
axis(2,las=2,at=seq(0,1,0.20),labels=seq(0,1,0.2),
     cex.lab=1.2,cex.axis=1.2)
##    a. All data (all nets as per Nash et al)
polygon(c(theta_sim,rev(theta_sim)),c(mod_data_b[[2]][,2],rev(mod_data_b[[2]][,3])),border = NA,col=adegenet::transp("grey75",0.4))
polygon(c(theta_sim,rev(theta_sim)),c(mod_data_c[[2]][,2],rev(mod_data_c[[2]][,3])),border = NA,col=adegenet::transp("darkred",0.4))
lines(mod_data_b[[2]][,1] ~ theta_sim,lty=2,col="black",lwd=2)
lines(mod_data_c[[2]][,1] ~ theta_sim,lty=2,col="darkred",lwd=2)

legend("topright",legend=c("Pyrethroid-only ITNs",
                           "Pyrethroid-PBO ITNs"),
       pch=15, lty=2, col=adegenet::transp(c("grey75","darkred"),0.9),
       bty="n")


# extract P_control
extract_control <- function(P_param, fit_full){
  P <- rstan::extract(fit_full, P_param)[[1]]
  return(median(P))
}

P_control_b <- extract_control("P_control", fit1_b)
P_control_c <- extract_control("P_control", fit1_c)




# plot(mod_data_b[[3]][,1] ~ theta_sim,type = "l",
#      ylim=c(0,1),xlim=c(0,1),yaxt="n",
#      col="grey55",
#      ylab="Probability of being caught in control hut",
#      xlab="Probability of hut survival")
# axis(2,las=2,at=seq(0,1,0.20),labels=seq(0,1,0.2))
# ##    a. All data (all nets as per Nash et al)
# polygon(c(theta_sim,rev(theta_sim)),c(mod_data_b[[3]][,2],rev(mod_data_b[[3]][,3])),border = NA,col=adegenet::transp("grey75",0.4))
# lines(mod_data_b[[3]][,2] ~ theta_sim,lty=2,col="grey75")
# lines(mod_data_b[[3]][,3] ~ theta_sim,lty=2,col="grey75")
# abline(h=P_control_b,lwd=2,col="grey75")
# 
# ##    b. Pyr
# lines(mod_data_b[[3]][,1] ~ theta_sim,lty=1,col="black")
# lines(mod_data_b[[3]][,2] ~ theta_sim,lty=2,col="black")
# lines(mod_data_b[[3]][,3] ~ theta_sim,lty=2,col="black")
# abline(h=P_control_b,lwd=2,col="black")
# 
# ##    b. Pyr-PBO
# lines(mod_data_c[[3]][,1] ~ theta_sim,lty=1,col="purple")
# lines(mod_data_c[[3]][,2] ~ theta_sim,lty=2,col="purple")
# lines(mod_data_c[[3]][,3] ~ theta_sim,lty=2,col="purple")
# abline(h=P_control_c,lwd=2,col="purple")

########################
##
## Part 3 successful feeding
setup_inputs_fed = function(df2){
  
  
  fed_IDs <- unique(df2$Orig_ID) # 37 trials
  fed_refs <- unique(df2$Library_ref) # from 31 publications/reports
  fed_totalsunknown_ref <- df2$Library_ref[df2$Totals_known=="0"] # 10 without totals known (only gave prop fed)
  fed_totalsunknown_ID <- df2$Orig_ID[df2$Totals_known=="0"]
  fed_refstotalsunknown <- unique(fed_totalsunknown_ref) # Publication ref where totals unknown: 76  77  81  84 
  fed_IDstotalsunknown <- unique(fed_totalsunknown_ID) #  Trial IDs where totals unknown : 83  84  88  91 
  
  # Data for stan
  S <- as.numeric(length(df2$N_total))
  X_survive <- df2$N_total-df2$N_dead_comb
  N_caught <- df2$N_total
  X_fed <- df2$N_fed
  X_survfed <- as.integer(X_fed * (1-(df2$N_dead_comb/df2$N_total))) # No. fed * proportion survived
  Perc_survfed <- X_survfed/df2$N_total*100
  nsite <- max(df2$Site_number)
  site <- df2$Site_number
  S_sim <- 300
  theta_sim <- seq(0, 1, length.out = S_sim)
  
  # Need to pass it the data:
  data_stan <- list(S=S, 
                    X_survive=X_survive,
                    N_caught=N_caught,
                    X_succfed=X_survfed,
                    nsite=nsite,
                    site=site,
                    S_sim = S_sim,
                    theta_sim = theta_sim)
  
  return(data_stan)
}

# full_model<- stan_model("R code/stan models/Model_succfed_WithREs.stan")
# 
# run_f2 = function(df2,dataset){
#   data_stan = setup_inputs_fed(df2)
#   
#   # Fit model
#   fit_hut_ALL <- sampling(full_model, data=data_stan, iter=2000, chains=4)
#   
#   # save fit
#   saveRDS(fit_hut_ALL, paste0("stan model outputs/LLINEUP_ento_feeding_0washes_",dataset,".rds"))
# }
# 
# dfpyr = dfpyr[!is.na(dfpyr$N_fed),]
# dfpbo = dfpbo[!is.na(dfpbo$N_fed),]
# run_f2(dfpyr,"pyr")
# run_f2(dfpbo,"pbo")

fit3_b <- readRDS("stan model outputs/LLINEUP_ento_feeding_0washes_pyr.rds")
fit3_c <- readRDS("stan model outputs/LLINEUP_ento_feeding_0washes_pbo.rds")
print(fit3_b, pars=c("a","b")) 
print(fit3_c, pars=c("a","b")) 

create_output2_f = function(fit_hut_ALL){
  # Extract prob of succ fed
  plot_fed_prob <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_med <- apply(P, 2, median)
    return(P_med*100)
  }
  
  P_succfed <- plot_fed_prob("P_fed_sim", fit_hut_ALL)
  
  
  # Extract upper and lower cred intervals
  plot_upper_funct <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_upper <- apply(P, 2, function(x) quantile(x,0.025))
    return(P_upper)
  }
  
  plot_lower_funct <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_lower <- apply(P, 2, function(x) quantile(x,0.975))
    return(P_lower) 
  }
  
  P_fed_upper <- plot_upper_funct("P_fed_sim", fit_hut_ALL)*100
  P_fed_lower <- plot_lower_funct("P_fed_sim", fit_hut_ALL)*100
  
  
  a_theta <- seq(0,1,length=300)
  my_surv <- a_theta*100
  
  cred_int_data <- data.frame(P_succfed, P_fed_upper, P_fed_lower,my_surv)
  
  return(cred_int_data)  
}
cred_int_data_b = create_output2_f(fit3_b)
cred_int_data_c = create_output2_f(fit3_c)

fed.data <- df %>%
  filter(N_dead_comb!="NA" & N_fed!="NA" & Intervention!="Untreated" & Hut!="Ifakara") 

fed.data$pc_surv <- (fed.data$N_total-fed.data$N_dead_comb)/fed.data$N_total*100
fed.data$prop_surv_fed <- fed.data$pc_fed/100*fed.data$pc_surv/100
fed.data$suc.u <- fed.data$prop_surv_fed*fed.data$N_total
fed.data$pc_surv_fed <- 100* fed.data$suc.u/fed.data$N_total
fed.data$prop_dead <- fed.data$N_dead_comb/fed.data$N_total
fed.data$repelled<-(1-fed.data$prop_surv_fed-fed.data$prop_dead)*100  # proportion repelled * 100 = percentage repelled

fed.data$pt_cex = c(0.5 + fed.data$N_total/1000)

outline.hut <- ifelse(fed.data$Hut=="West","#00008b",
                      ifelse(fed.data$Hut=="East","#8b0000",
                             ifelse(fed.data$Hut=="Ifakara","#006400",NA)))

fed.datab <- fed.data %>% filter(Intervention%in%selectionb)
fed.datac <- fed.data %>% filter(Intervention%in%selectionc)

plot(fed.datab$pc_surv_fed ~ fed.datab$pc_surv, cex=fed.datab$pt_cex,
     col=adegenet::transp("white",0.4),
     pch=19,
     ylab = "Successfully blood fed (%)",
     xlab = "Mosquito survival in EHT (%)",
     yaxt="n",ylim=c(0,100),
     cex.lab=1.2,cex.axis=1.2)
axis(2,las=2,at=c(0,20,40,60,80,100),
     cex.lab=1.2,cex.axis=1.2)

points(fed.datab$pc_surv_fed ~ fed.datab$pc_surv,
       cex=fed.datab$pt_cex,
       col=adegenet::transp("grey",0.4),
       pch=19,ylim=c(0,100))
points(fed.datac$pc_surv_fed ~ fed.datac$pc_surv,
       cex=fed.datac$pt_cex,
       col=adegenet::transp("purple",0.4),
       pch=19,ylim=c(0,100))

lines(cred_int_data_b$P_succfed ~ cred_int_data_b$my_surv, col = "black",lwd=2)
polygon(c(cred_int_data_b$my_surv,rev(cred_int_data_b$my_surv)),
        c(cred_int_data_b$P_fed_upper,rev(cred_int_data_b$P_fed_lower)),
        border=NA,col=adegenet::transp("grey55",0.4))

lines(cred_int_data_c$P_succfed ~ cred_int_data_c$my_surv, col = "darkred",lwd=2)
polygon(c(cred_int_data_c$my_surv,rev(cred_int_data_c$my_surv)),
        c(cred_int_data_c$P_fed_upper,rev(cred_int_data_c$P_fed_lower)),
        border=NA,col=adegenet::transp("darkred",0.4))

