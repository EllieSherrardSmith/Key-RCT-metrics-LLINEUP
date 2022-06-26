## Part 1

## fitting the stan model for mosquito survival in susc bioassay 
## to predict the survival in experimental huts

library(rstan)

# For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
# To avoid recompilation of unchanged Stan programs, we recommend calling
rstan_options(auto_write = TRUE)

## data
mt = read.csv("data/mortality_data_from_Nash2021.csv",header=TRUE)

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

stan_base <- rstan::stan(file="R code/stan models/Bioassay/A1 Log-logistic Model.stan", 
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
plot(surv ~ x,type="l",
     xlab = "Survival at discriminatory dose bioassay (%)",
     ylab = "Survival in experimental hut trial (%)",
     xaxt="n",yaxt="n",ylim=c(0,1),xlim=c(0,1))
axis(1,at=seq(0,1,by=0.2),labels=seq(0,100,by=20))
axis(2,las=2,at=seq(0,1,by=0.2),labels=seq(0,100,by=20))

# quantile(base$a,0.025)
# quantile(base$b,0.975)
# for(i in c(2625,3995)){
#   surv = 1 - (1/(1+((1-x)/base$b[i])^(-base$a[i])))
#   lines(surv ~ x,col="black")
# }
surv_upp = surv = 1 - (1/(1+((1-x)/base$b[2625])^(-base$a[2625])))
surv_low = surv = 1 - (1/(1+((1-x)/base$b[3995])^(-base$a[3995])))
polygon(c(x,rev(x)),c(surv_low,rev(surv_upp)),col=adegenet::transp("yellow",0.4),border=NA)




points(c((mt$Bioassay_N_tested - mt$Bioassay_N_dead)/mt$Bioassay_N_tested),
       c((mt$N_total - mt$N_dead)/mt$N_total),col=adegenet::transp("orange",0.7),pch=19,
       cex = c(log(50*(mt$N_total/max(mt$N_total))+1)))

## saveRDS(stan_base,"stan model outputs/log_logistic_fit.RDS")
