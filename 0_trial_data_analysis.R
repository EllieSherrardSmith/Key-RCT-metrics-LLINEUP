## Read in data

data = read.csv("raw data/data.csv",header=TRUE)

## Additional vectors
data$net_type_brand = ifelse(data$Net_Type == "Olyset Net", 1,
                                    ifelse(data$Net_Type == "PermaNet 2.0", 2,
                                           ifelse(data$Net_Type == "Olyset Plus", 3, 4)))
data$net_type_class = ifelse(data$Net_Type == "Olyset Net", 1,
                                    ifelse(data$Net_Type == "PermaNet 2.0", 1,2))

data$efficacy_6m = (data$Prevalence_baseline_2_10_yrs - data$Prevalence_6m)/data$Prevalence_baseline_2_10_yrs
data$efficacy_12m = (data$Prevalence_baseline_2_10_yrs - data$Prevalence_12m)/data$Prevalence_baseline_2_10_yrs
data$efficacy_18m = (data$Prevalence_baseline_2_10_yrs - data$Prevalence_18m)/data$Prevalence_baseline_2_10_yrs
data$efficacy_25m = (data$Prevalence_baseline_2_10_yrs - data$Prevalence_25m)/data$Prevalence_baseline_2_10_yrs

## Figure 1 

## (A) 
## Prevalence in each trial arm with jittered points for clusters
## at baseline
par(mfrow=c(2,2))
par(mar=c(5,5,2,2))
prev_base = as.numeric(tapply(data$Prevalence_baseline_2_10_yrs,data$net_type_brand,mean))
boxplot(data$Prevalence_baseline_2_10_yrs~data$net_type_brand,ylim=c(0,1),yaxt="n",
        ylab="",
        xaxt="n",xlab="",bty=NA,frame=FALSE,
        col = adegenet::transp(c("red","blue","orange","darkgreen"),0.7))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(1,at=1:4,labels=c("Olyset Net","PermaNet 2.0","Olyset Plus","PermaNet 3.0"))
mtext(line = 3,side=2,text="Prevalence in children 2 to 10-years (%)")
pt_col = c("darkred","darkblue","darkorange","darkgreen")
rng = seq(0.7,4.3,length=4)
for(i in 1:4){
  points(data$Prevalence_baseline_2_10_yrs[data$net_type_brand == i]~
           c(rnorm(n = length(data$Prevalence_baseline_2_10_yrs[data$net_type_brand == i]),mean = i,sd=0.1)),
         col=pt_col[i],pch=19)
}

## (B - not using) 
## Prevalence in each trial arm with jittered points for clusters
## at 6,12,18,25 MONTHS
# boxplot(data$Prevalence_6m[data$net_type_brand == 1],
#         data$Prevalence_12m[data$net_type_brand == 1],
#         data$Prevalence_18m[data$net_type_brand == 1],
#         data$Prevalence_25m[data$net_type_brand == 1],
#         
#         data$Prevalence_6m[data$net_type_brand == 2],
#         data$Prevalence_12m[data$net_type_brand == 2],
#         data$Prevalence_18m[data$net_type_brand == 2],
#         data$Prevalence_25m[data$net_type_brand == 2],
#         
#         data$Prevalence_6m[data$net_type_brand == 3],
#         data$Prevalence_12m[data$net_type_brand == 3],
#         data$Prevalence_18m[data$net_type_brand == 3],
#         data$Prevalence_25m[data$net_type_brand == 3],
#         
#         data$Prevalence_6m[data$net_type_brand == 4],
#         data$Prevalence_12m[data$net_type_brand == 4],
#         data$Prevalence_18m[data$net_type_brand == 4],
#         data$Prevalence_25m[data$net_type_brand == 4],
#         ylim=c(0,1),yaxt="n",
#         ylab="Prevalence in children 2 to 10-years (%)",
#         xaxt="n",xlab="Month since deployment of net",bty=NA,frame=FALSE,
#         col = rep(adegenet::transp(c("red","blue","orange","darkgreen"),0.7),each=4))
# axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
# axis(1,at=1:16,labels=rep(c(6,12,18,25),4))
# 
# pt_col = c("darkred","darkblue","darkorange","darkgreen")
# rng = c(1,5,9,13)
# for(i in 1:4){
#   points(data$Prevalence_6m[data$net_type_brand == i]~
#            c(rnorm(n = length(data$Prevalence_6m[data$net_type_brand == i]),mean = rng[i],sd=0.1)),
#          col=pt_col[i],pch=19)
#   points(data$Prevalence_12m[data$net_type_brand == i]~
#            c(rnorm(n = length(data$Prevalence_6m[data$net_type_brand == i]),mean = rng[i]+1,sd=0.1)),
#          col=pt_col[i],pch=19)
#   points(data$Prevalence_18m[data$net_type_brand == i]~
#            c(rnorm(n = length(data$Prevalence_6m[data$net_type_brand == i]),mean = rng[i]+2,sd=0.1)),
#          col=pt_col[i],pch=19)
#   points(data$Prevalence_25m[data$net_type_brand == i]~
#            c(rnorm(n = length(data$Prevalence_6m[data$net_type_brand == i]),mean = rng[i]+3,sd=0.1)),
#          col=pt_col[i],pch=19)
# }

## (C)
par(mar=c(5,6,2,2))
boxplot(data$efficacy_6m[data$net_type_brand == 1],
        data$efficacy_12m[data$net_type_brand == 1],
        data$efficacy_18m[data$net_type_brand == 1],
        data$efficacy_25m[data$net_type_brand == 1],
        
        data$efficacy_6m[data$net_type_brand == 2],
        data$efficacy_12m[data$net_type_brand == 2],
        data$efficacy_18m[data$net_type_brand == 2],
        data$efficacy_25m[data$net_type_brand == 2],
        
        data$efficacy_6m[data$net_type_brand == 3],
        data$efficacy_12m[data$net_type_brand == 3],
        data$efficacy_18m[data$net_type_brand == 3],
        data$efficacy_25m[data$net_type_brand == 3],
        
        data$efficacy_6m[data$net_type_brand == 4],
        data$efficacy_12m[data$net_type_brand == 4],
        data$efficacy_18m[data$net_type_brand == 4],
        data$efficacy_25m[data$net_type_brand == 4],
        ylim=c(-1,1),yaxt="n",
        ylab="",
        xaxt="n",xlab="Month since deployment of net",bty=NA,frame=FALSE,
        col = rep(adegenet::transp(c("red","blue","orange","darkgreen"),0.7),each=4))
axis(2,las=2,at=seq(-0.6,1,0.2),labels=seq(-60,100,20))
axis(1,at=1:16,labels=rep(c(6,12,18,25),4))
mtext(side = 2,line = 4,text = "Relative reduction in prevalence since")
mtext(side = 2,line = 3,text = "baseline at each cross-sectional survey (%)")

pt_col = c("darkred","darkblue","darkorange","darkgreen")
rng = c(1,5,9,13)
for(i in 1:4){
  points(data$efficacy_6m[data$net_type_brand == i]~
           c(rnorm(n = length(data$efficacy_6m[data$net_type_brand == i]),mean = rng[i],sd=0.1)),
         col=pt_col[i],pch=19)
  points(data$efficacy_12m[data$net_type_brand == i]~
           c(rnorm(n = length(data$efficacy_12m[data$net_type_brand == i]),mean = rng[i]+1,sd=0.1)),
         col=pt_col[i],pch=19)
  points(data$efficacy_18m[data$net_type_brand == i]~
           c(rnorm(n = length(data$efficacy_18m[data$net_type_brand == i]),mean = rng[i]+2,sd=0.1)),
         col=pt_col[i],pch=19)
  points(data$efficacy_25m[data$net_type_brand == i]~
           c(rnorm(n = length(data$efficacy_25m[data$net_type_brand == i]),mean = rng[i]+3,sd=0.1)),
         col=pt_col[i],pch=19)
}


## (C) Efficacy for the net class

par(mar=c(5,6,2,2))
boxplot(data$efficacy_6m[data$net_type_brand == 1 | data$net_type_brand == 2],
        data$efficacy_12m[data$net_type_brand == 1 | data$net_type_brand == 2],
        data$efficacy_18m[data$net_type_brand == 1 | data$net_type_brand == 2],
        data$efficacy_25m[data$net_type_brand == 1 | data$net_type_brand == 2],
        
        data$efficacy_6m[data$net_type_brand == 3 | data$net_type_brand == 4],
        data$efficacy_12m[data$net_type_brand == 3 | data$net_type_brand == 4],
        data$efficacy_18m[data$net_type_brand == 3 | data$net_type_brand == 4],
        data$efficacy_25m[data$net_type_brand == 3 | data$net_type_brand == 4],
        
        ylim=c(-1,1),yaxt="n",
        ylab="",
        xaxt="n",xlab="Month since deployment of net",bty=NA,frame=FALSE,
        col = rep(adegenet::transp(c("lightblue","green"),0.7),each=4))
axis(2,las=2,at=seq(-0.6,1,0.2),labels=seq(-60,100,20))
axis(1,at=1:8,labels=rep(c(6,12,18,25),2))
mtext(side = 2,line = 4,text = "Relative reduction in prevalence since")
mtext(side = 2,line = 3,text = "baseline at each cross-sectional survey (%)")

pt_col = c("darkred","darkblue","darkorange","darkgreen")
rng2 = c(1,2,3,4)
points(data$efficacy_6m[data$net_type_brand == 1 | data$net_type_brand == 2]~
           c(rnorm(n = length(data$efficacy_6m[data$net_type_brand == 1 | data$net_type_brand == 2]),mean = rng2[1],sd=0.1)),
         col="blue",pch=19)
points(data$efficacy_12m[data$net_type_brand == 1 | data$net_type_brand == 2]~
           c(rnorm(n = length(data$efficacy_12m[data$net_type_brand == 1 | data$net_type_brand == 2]),mean = rng2[2],sd=0.1)),
         col="blue",pch=19)
points(data$efficacy_18m[data$net_type_brand == 1 | data$net_type_brand == 2]~
           c(rnorm(n = length(data$efficacy_18m[data$net_type_brand == 1 | data$net_type_brand == 2]),mean = rng2[3],sd=0.1)),
         col="blue",pch=19)
points(data$efficacy_25m[data$net_type_brand == 1 | data$net_type_brand == 2]~
           c(rnorm(n = length(data$efficacy_25m[data$net_type_brand == 1 | data$net_type_brand == 2]),mean = rng2[4],sd=0.1)),
         col="blue",pch=19)

points(data$efficacy_6m[data$net_type_brand == 3 | data$net_type_brand == 4]~
         c(rnorm(n = length(data$efficacy_6m[data$net_type_brand == 3 | data$net_type_brand == 4]),mean = 5,sd=0.1)),
       col="darkgreen",pch=19)
points(data$efficacy_12m[data$net_type_brand == 3 | data$net_type_brand == 4]~
         c(rnorm(n = length(data$efficacy_12m[data$net_type_brand == 3 | data$net_type_brand == 4]),mean = 6,sd=0.1)),
       col="darkgreen",pch=19)
points(data$efficacy_18m[data$net_type_brand == 3 | data$net_type_brand == 4]~
         c(rnorm(n = length(data$efficacy_18m[data$net_type_brand == 3 | data$net_type_brand == 4]),mean = 7,sd=0.1)),
       col="darkgreen",pch=19)
points(data$efficacy_25m[data$net_type_brand == 3 | data$net_type_brand == 4]~
         c(rnorm(n = length(data$efficacy_25m[data$net_type_brand == 3 | data$net_type_brand == 4]),mean = 8,sd=0.1)),
       col="darkgreen",pch=19)
abline(v=4.5,lty=2)
text(2.5,-0.7,"Pyrethroid-only ITNs")
text(6.5,-0.7,"Pyrethroid-PBO ITNs")

# (D) The prevalence map

####################################
##
## Uganda map

library(rgdal)
library(classInt)
library(dplyr)
library(rgeos)
library(rmapshaper)

world = readOGR(dsn="H:/GADM/version3.6/GADM36_1", layer = "gadm36_1_new", stringsAsFactors = FALSE)
afr_ad = world[which(world$CONTINENT == "Africa"),]

w0 = readOGR(dsn="H:/GADM/version3.6/GADM36_0", layer = "gadm36_0", stringsAsFactors = FALSE)
afr0_ad = w0[w0$NAME_0 %in% afr_ad$NAME_0, ]

rm(list = c("w0", "world"))

afr0 <- gSimplify(afr0_ad, tol=0.01, topologyPreserve=TRUE)

afr = SpatialPolygonsDataFrame(afr_ad, data = afr_ad@data)
afr0 = SpatialPolygonsDataFrame(afr0, data = afr0_ad@data)

afr@data <- plyr::join(afr@data, data, by="DIDE_CODE")
names(afr@data)

delc = c("Uganda")

afr = afr[which(afr$NAME_0 %in% delc),]


legend_names <- function(x, pretty=FALSE){
  if(pretty) x = prettyNum(x,big.mark=",")
  op = rep(NA, length(x)-1)
  op = paste(x[1:(length(x)-1)],"-",x[2:length(x)])
  return(op)
}


# png(file = paste0("Uganda_map.png"), width=260,height=260,units="mm",res=300)


p_break = c(-1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.8)
cols_temp = rev(RColorBrewer::brewer.pal(length(p_break), "YlGnBu"))
cols = c(RColorBrewer::brewer.pal(length(p_break), "YlGnBu"))
cols = c("white",cols_temp[1:9])
# cols[length(cols)+1] = "black"

par(mfrow=c(1,1), mar=c(0,0,2,0), oma=c(0,0,3,0), xpd=FALSE)


dat1 = afr@data[,14] ## CHECK COLUMN IS PREV AT BASELINE
dat1[is.na(dat1)] <- -1
class <- classIntervals(dat1, (length(p_break)-1), style="fixed",fixedBreaks=p_break)
colcode <- findColours(class, cols)
plot(afr, col=colcode, border="grey35")
plot(afr0, border="white", lwd=3, add=TRUE)
plot(afr0, lwd=1, add=TRUE)
# mtext("Lives saved MDA all-ages", 3, line = -3, font=2, cex=1.2)

par(xpd=NA)
options(scipen = 999)
legend("left", legend=c(legend_names(100*p_break[2:9], pretty=TRUE)), col=c(cols[2:9]), pch=15, pt.cex=2, cex=1.1,title="Baseline prevalence 2 to 10 yr olds (%)", bty="n")


#####################################################
##
## Figure 2 is the exploration of how variable the baseline
## metrics are to one another. 

## (a) clusters (map from LLINEUP TEAM)
## (b) densities (map from LLINEUP TEAM)
## (c) resistance (map)
## as above then

p_break = c(0,0.75,0.77,0.79,0.81,0.83,0.85,0.88,0.90)
cols_temp = rev(RColorBrewer::brewer.pal(length(p_break), "YlGnBu"))
cols = c(RColorBrewer::brewer.pal(length(p_break), "YlGnBu"))
cols = c("white",cols_temp[1:9])
# cols[length(cols)+1] = "black"

par(mfrow=c(1,1), mar=c(0,0,2,0), oma=c(0,0,3,0), xpd=FALSE)


dat1 = afr@data[,13] ## CHECK COLUMN IS PYRETHROID RESISTANCE
dat1[is.na(dat1)] <- 0
class <- classIntervals(dat1, (length(p_break)-1), style="fixed",fixedBreaks=p_break)
colcode <- findColours(class, cols)
plot(afr, col=colcode, border="grey35")
plot(afr0, border="white", lwd=3, add=TRUE)
plot(afr0, lwd=1, add=TRUE)
# mtext("Lives saved MDA all-ages", 3, line = -3, font=2, cex=1.2)

par(xpd=NA)
options(scipen = 999)
legend("left", legend=c(legend_names(100*p_break, pretty=TRUE)), col=c(cols), pch=15, pt.cex=2, cex=1.1,title="Estimated pyrethroid resistance at baseline (%)", bty="n")


## (d) baseline nets
## baseline

hist(data$X3_and_over_net_use_baseline,
     ylab="Frequency",main="",
     xlab="Proportion of people using nets at baseline",
     col = "darkblue",bty="n",
     breaks = 15,xlim=c(0,1),
     yaxt="n",xaxt="n")
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(2,las=2,at=c(0,5,10,15))

baseline3nets = subset(data,data$Net_Type != "PermaNet 2.0")
hist(baseline3nets$X3_and_over_net_use_baseline,
     ylab="Frequency",main="",
     xlab="Proportion of people using nets at baseline",
     col = "darkred",bty="n",
     breaks = 15,xlim=c(0,1),add=T,
     yaxt="n",xaxt="n")

baseline2nets = subset(data,data$Net_Type == "PermaNet 3.0" | data$Net_Type ==  "Olyset Plus")
hist(baseline2nets$X3_and_over_net_use_baseline,
     ylab="Frequency",main="",
     xlab="Proportion of people using nets at baseline",
     col = "darkgreen",bty="n",
     breaks = 15,xlim=c(0,1),add=T,
     yaxt="n",xaxt="n")

baseline1nets = subset(data,data$Net_Type ==  "Olyset Plus")
hist(baseline1nets$X3_and_over_net_use_baseline,
     ylab="Frequency",main="",
     xlab="Proportion of people using nets at baseline",
     col = "orange",bty="n",
     breaks = 15,xlim=c(0,1),add=T,
     yaxt="n",xaxt="n")

legend("topleft",legend = c("PermaNet 2.0","Olyset net",
                            "PermaNet 3.0","Olyset Plus"),
       col=c("darkblue","darkred","darkgreen","orange"),
       bty="n",pch=15)

## Can net type explain baseline net use?
summary(aov(data$X3_and_over_net_use_baseline ~ data$Net_Type))

## (e) time of deployment
## net_timing

data$net_colour = ifelse(data$Net_Type == "PermaNet 2.0", "darkblue",
                         ifelse(data$Net_Type == "PermaNet 3.0","darkgreen",
                                ifelse(data$Net_Type == "Olyset Net", "darkred","orange")))

tapply(data$net_colour,data$Distribution_Month,length)

data$net_delivery_order =  ifelse(data$Distribution_Month == "Mar-17", 1,
                                  ifelse(data$Distribution_Month == "May-17",2,
                                         ifelse(data$Distribution_Month == "Jul-17", 3,4)))


three_nets = subset(data,data$Net_Type != "PermaNet 2.0")
two_nets = subset(data,data$Net_Type == "PermaNet 3.0" | data$Net_Type ==  "Olyset Plus")
one_net = subset(data,data$Net_Type == "Olyset Plus")

all = as.numeric(tapply(data$Net_Type,data$net_delivery_order,length))
barplot(all,col="darkblue",yaxt="n",ylab="Number of clusters")
axis(2,las=2,at=c(0,10,20,30,40))
mtext("Month for net deployment",side=1,line=3)
axis(1,at=c(0.75,1.9,3.1,4.25),labels=c("Mar-17","May-17","Jul-17","Mar-18"))

three = as.numeric(tapply(three_nets$Net_Type,three_nets$net_delivery_order,length))
barplot(three,col="darkred",yaxt="n",add=T)
two = as.numeric(tapply(two_nets$Net_Type,two_nets$net_delivery_order,length))
barplot(two,col="darkgreen",yaxt="n",add=T)
one = as.numeric(tapply(one_net$Net_Type,one_net$net_delivery_order,length))
barplot(one,col="orange",yaxt="n",add=T)

legend("topleft",legend = c("PermaNet 2.0","Olyset net",
                            "PermaNet 3.0","Olyset Plus"),
       col=c("darkblue","darkred","darkgreen","orange"),
       bty="n",pch=15)

## (f) net use
###################################################
##
## 1 Site parameter set up

## Work out seasonal pattern relative to timing of net distributions

## Fit function to net use estimates for each of 104 clusters


library(rstan)
library(adegenet)

## Looping through the 104 clusters 
PARMS_USAGE = expand.grid(1:1000)
PARMS_USAGE2 = array(dim=c(2,104))
## Specify the times when data were observed
time_obs = c(0,6/12,12/12,18/12)

# Coverage 6-months: 71% ## From supplementary figures Staedke et al 2020
# Coverage 12-months: 63%
# Coverage 18-months: 49%
## Specify the data observed (including the uncertainties as this provides more data points for a better fit)
## where more specific data for different clusters are available this can be improved.

net_dataAll = data


for(i in 1:104){
  
  ## ORIGINAL USES NET_DATA ONLY AND 6-18 MONTHS
  standard_net_usage = c(0.95,as.numeric(net_dataAll[i,c(12,14,16)]))
  
  # standard_net_usage = c(0.71,0.63,0.49,
  #                        0.62,0.58,0.44, ##these are the estimated uncertainty ranges
  #                        0.72,0.68,0.54) ##these are the estimated uncertainty ranges
  
  ## Transform data for the model fitting
  y_standard_net_usage = log(standard_net_usage)
  
  ## Specify a time series to project outcomes across
  time_m = seq(0,3,0.01)
  
  ## Specify the exponential decay model in Rstan
  stanmodelcode <- "
  data {
  int<lower=0> N;
  int<lower=0> N2; //the size of the new_X matrix
  vector[N] y;
  vector[N] x;
  vector[N2] New_x;
  
  }
  parameters {
  real beta0;
  real beta1;
  real sigma;
  }
  transformed parameters {
  vector[N] m;
  m = beta0 + beta1 * x;
  } 
  model {
  // priors
  beta0 ~ cauchy(0, 10); 
  beta1 ~ cauchy(0, 10); 
  
  // likelihood
  y ~ normal(m, sigma);   
  }
  generated quantities {
  vector[N2] y_pred;
  y_pred = beta0 + beta1 * New_x; //the y values predicted by the model
  }
  "
  
  # Just need to run the below 1 time
  # stanDso <- stan_model(model_code = stanmodelcode) 
  
  ## Put data into a list
  dat_standard <- list(N = length(time_obs), 
                       N2 = length(seq(0,3,0.01)),
                       y = y_standard_net_usage, 
                       x = time_obs,
                       New_x = seq(0,3,0.01))
  
  ## Run the statistical model
  fit <- sampling(stanDso, 
                  data = dat_standard, ## specify data
                  iter = 4000,         ## speify iterations
                  warmup=2000)         ## specify warm up (default is 4 chains here)
  
  
  ## plotting the posterior distribution for the parameters
  # post_beta<-As.mcmc.list(fit,pars="beta0")
  # plot(post_beta)
  
  ## gradient is fit to the data for alpha
  ## standard_net_usage ~ exp(-alpha*time_obs)
  b0 <- extract(fit, 'beta0')
  b0<- unlist(b0, use.names=FALSE)
  b1 <- extract(fit, 'beta1')
  b1<- unlist(b1, use.names=FALSE)
  
  ## Back translate output estimates
  y_predicted = mean(b0) + mean(b1)*time_m; 
  y_predicted_stn_exp = exp(mean(b0)) * exp(mean(b1)*time_m)
  
  ## and select some uncertainty for the Bayesian estimates
  y_predicted_stn_exp_min = exp(quantile(b0,0.25)) * exp(quantile(b1,0.25)*time_m)
  y_predicted_stn_exp_max = exp(quantile(b0,0.75)) * exp(quantile(b1,0.75)*time_m)
  
  ## Pull out the parameter required to specify drop out from using ITNs in the transmission model 
  parms_usage = data.frame(itn_leave_dur_standardLLIN = b1[b1 >= quantile(b1,0.25) & b1 <= quantile(b1,0.75)]) ##
  parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN
  
  PARMS_USAGE[,i] = sample(parms_usage$itn_leave_dur_standardLLIN_bt,1000,replace=FALSE)
  PARMS_USAGE2[1,i] = mean(b0)
  PARMS_USAGE2[2,i] = mean(b1)
}

data$itn_leave_dur = as.numeric(colMeans(PARMS_USAGE))
data$itn_leave_dur = ifelse(data$itn_leave_dur > 15,15,data$itn_leave_dur)

data$itn_leave_durMean = -1/PARMS_USAGE2[2,]
data$itn_leave_durMean = ifelse(data$itn_leave_durMean > 15,15,data$itn_leave_durMean)
data$itn_leave_durMean = ifelse(data$itn_leave_durMean < 0,15,data$itn_leave_durMean)

data$b0 = PARMS_USAGE2[1,]
data$b1 = PARMS_USAGE2[2,]

write.csv(data, "raw data/dd_params_attrition.csv")

i = 2
y_predicted_stn_exp = exp(PARMS_USAGE2[1,i]) * exp(PARMS_USAGE2[2,i]*time_m)
plot(y_predicted_stn_exp ~ time_m,ylim=c(0,1))
points(c(0.95,as.numeric(net_dataAll[i,c(12,14,16)])) ~ time_obs,pch=19,col="orange")

standard_net_usage = c(0.95,as.numeric(data[i,c(12,14,16)]))

# standard_net_usage = c(0.71,0.63,0.49,
#                        0.62,0.58,0.44, ##these are the estimated uncertainty ranges
#                        0.72,0.68,0.54) ##these are the estimated uncertainty ranges

## Transform data for the model fitting
y_standard_net_usage = log(standard_net_usage)

## Specify a time series to project outcomes across
time_m = seq(0,3,0.01)

## Specify the exponential decay model in Rstan
stanmodelcode <- "
data {
int<lower=0> N;
int<lower=0> N2; //the size of the new_X matrix
vector[N] y;
vector[N] x;
vector[N2] New_x;

}
parameters {
real beta0;
real beta1;
real sigma;
}
transformed parameters {
vector[N] m;
m = beta0 + beta1 * x;
} 
model {
// priors
beta0 ~ cauchy(0, 10); 
beta1 ~ cauchy(0, 10); 

// likelihood
y ~ normal(m, sigma);   
}
generated quantities {
vector[N2] y_pred;
y_pred = beta0 + beta1 * New_x; //the y values predicted by the model
}
"
stanDso <- stan_model(model_code = stanmodelcode) 

## Put data into a list
dat_standard <- list(N = length(time_obs), 
                     N2 = length(seq(0,3,0.01)),
                     y = y_standard_net_usage, 
                     x = time_obs,
                     New_x = seq(0,3,0.01))

## Run the statistical model
fit <- sampling(stanDso, 
                data = dat_standard, ## specify data
                iter = 4000,         ## speify iterations
                warmup=2000)         ## specify warm up (default is 4 chains here)


## plotting the posterior distribution for the parameters
post_beta<-As.mcmc.list(fit,pars="beta0")
plot(post_beta)

## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)
b0 <- extract(fit, 'beta0')
b0<- unlist(b0, use.names=FALSE)
b1 <- extract(fit, 'beta1')
b1<- unlist(b1, use.names=FALSE)

## Back translate output estimates
y_predicted = mean(b0) + mean(b1)*time_m; 
y_predicted_stn_exp = exp(mean(b0)) * exp(mean(b1)*time_m)

## and select some uncertainty for the Bayesian estimates
y_predicted_stn_exp_min = exp(quantile(b0,0.25)) * exp(quantile(b1,0.25)*time_m)
y_predicted_stn_exp_max = exp(quantile(b0,0.75)) * exp(quantile(b1,0.75)*time_m)


###################################################
##
## Some fits are poorer than others
## Check

## Final output plotted to confirm the fit
par(mfrow = c(1,1))
plot(standard_net_usage[1:5] ~ ## plotting the mean estimates first
       time_obs[1:5],pch="",
     ylab="Households with at least one net per two occupants (%)", ## specify data used
     xlab="Time in years",yaxt="n",ylim=c(0,1),cex.lab=1.4,cex.axis=1.4,xlim=c(0,3))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.4,cex.axis=1.4)

# polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=NA,col=transp("grey","0.5"))
# lines(y_predicted_stn_exp ~ time_m,col="black",lty=2,lwd=1)

# ## add the range in the uncertainty from the observed data for the trial if available
# for(i in 1:4){
#   segments(x0=time_obs[i],x1=time_obs[i],
#            y0=standard_net_usage[i],y1=standard_net_usage[i+3],lty=1)
# }

for(i in 1:104){
  standard_net_usage = c(0.95,as.numeric(net_dataAll[i,c(12,14,16)]))
  points(standard_net_usage[1:4] ~ time_obs[1:4])
  
}

# ## Pull out the parameter required to specify drop out from using ITNs in the transmission model 
# parms_usage = data.frame(itn_leave_dur_standardLLIN = b1[b1 >= quantile(b1,0.2) & b1 <= quantile(b1,0.8)]) ##
# parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN

# PARMS_USAGE[,1] = sample(parms_usage$itn_leave_dur_standardLLIN_bt,1000,replace=FALSE)


for(i in c(1:54,56:90)){
  # Confirm this parameter is specifying what we expect
  D = median(PARMS_USAGE[,i])
  D_LOW = quantile(PARMS_USAGE[,i],0.01)
  D_UPP = quantile(PARMS_USAGE[,i],0.99)
  
  cover_at_start = 0.95
  cover_at_start_UPP = 1
  cover_at_start_LOW = 0.9
  
  ## The below replicate the way this parameter is included in the transmission model
  aa = cover_at_start * exp((-1/D)*time_m)
  aL = cover_at_start_LOW * exp((-1/D_LOW)*time_m)
  aU = cover_at_start_UPP * exp((-1/D_UPP)*time_m)
  
  lines(aa ~ time_m,col="blue")
  # lines(aL ~ time_m,col="blue",lty=3)
  # lines(aU ~ time_m,col="blue",lty=3)
  # 
}

##odd fits
# for(i in c(9,18,30,52,69,89,100,104)){
for(i in c(55)){
  # Confirm this parameter is specifying what we expect
  D = median(PARMS_USAGE[,i])
  D_LOW = quantile(PARMS_USAGE[,i],0.01)
  D_UPP = quantile(PARMS_USAGE[,i],0.99)
  
  cover_at_start = 0.95
  cover_at_start_UPP = 1
  cover_at_start_LOW = 0.9
  
  ## The below replicate the way this parameter is included in the transmission model
  aa = cover_at_start * exp((-1/D)*time_m)
  aL = cover_at_start_LOW * exp((-1/D_LOW)*time_m)
  aU = cover_at_start_UPP * exp((-1/D_UPP)*time_m)
  
  lines(aa ~ time_m,col="red")
  # lines(aL ~ time_m,col="blue",lty=3)
  # lines(aU ~ time_m,col="blue",lty=3)
  # 
}
oddies = c(9,18,30,52,65,66,76,69,83,89,100,104)
standard_net_usage_outliers = array(dim=c(4,length(oddies)))
for(i in 1:length(oddies)){
  standard_net_usage_outliers[,i] = c(0.95,as.numeric(net_data[oddies[i],c(4,6,8)]))
}


PARMS_USAGE2 = array(dim=c(1000,104))
for(i in 1:104){
  PARMS_USAGE2[,i] = ifelse(PARMS_USAGE[,i] > 15,15,PARMS_USAGE[,i])
  
}

PARMS_USAGE3 = array(dim=c(1000,104))
for(i in 1:104){
  PARMS_USAGE3[,i] = ifelse(PARMS_USAGE2[,i] < 0,1.5,PARMS_USAGE2[,i])
  
}
D = colMeans(PARMS_USAGE3)
hist(D,breaks = 50)

write.csv(PARMS_USAGE3)

######################
## Final output plotted to confirm the fit
par(mfrow = c(1,2))
plot(standard_net_usage[1:4] ~ ## plotting the mean estimates first
       time_obs[1:4],pch="",
     # ylab="Households with >= net per two occupants (%)", ## specify data used
     ylab="Estimated adherrence to net use (%)", ## specify data used
     xlab="Time in years since mass distribution",yaxt="n",ylim=c(0,1),cex.lab=1,cex.axis=1,xlim=c(0,3))
axis(2,las=2, at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1,cex.axis=1)

# polygon(c(time_m,rev(time_m)),c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=NA,col=transp("grey","0.5"))
# lines(y_predicted_stn_exp ~ time_m,col="black",lty=2,lwd=1)

# ## add the range in the uncertainty from the observed data for the trial if available
# for(i in 1:4){
#   segments(x0=time_obs[i],x1=time_obs[i],
#            y0=standard_net_usage[i],y1=standard_net_usage[i+3],lty=1)
# }
net_timing = read.csv("data/net_distribution.csv",header=TRUE)
## Additional vectors
net_timing$net_type_INPUTA = ifelse(net_timing$Net_Type == "Olyset Net", 1,
                                    ifelse(net_timing$Net_Type == "PermaNet 2.0", 2,
                                           ifelse(net_timing$Net_Type == "Olyset Plus", 3, 4)))
net_timing$net_type_INPUTB = ifelse(net_timing$Net_Type == "Olyset Net", 1,
                                    ifelse(net_timing$Net_Type == "PermaNet 2.0", 2,
                                           ifelse(net_timing$Net_Type == "Olyset Plus", 1, 2)))

cols_nets = c("blue","darkred","darkgreen","orange")
for(i in 1:104){
  standard_net_usage = c(0.95,as.numeric(net_data[i,c(4,6,8)]))
  points(standard_net_usage[1:4] ~ time_obs[1:4],
         col=cols_nets[net_timing$net_type_INPUTA[i]])
  
}

# ## Pull out the parameter required to specify drop out from using ITNs in the transmission model 
# parms_usage = data.frame(itn_leave_dur_standardLLIN = b1[b1 >= quantile(b1,0.2) & b1 <= quantile(b1,0.8)]) ##
# parms_usage$itn_leave_dur_standardLLIN_bt = -1/parms_usage$itn_leave_dur_standardLLIN

# PARMS_USAGE[,1] = sample(parms_usage$itn_leave_dur_standardLLIN_bt,1000,replace=FALSE)

D = colMeans(PARMS_USAGE3)
colours_nets = c(adegenet::transp("blue",0.2),"darkgreen",adegenet::transp("darkred"),"orange")
as.numeric(net_data$net_type)

colour_match = ifelse(as.numeric(net_data$net_type) == 1,colours_nets[1],
                      ifelse(as.numeric(net_data$net_type) == 2,colours_nets[2],
                             ifelse(as.numeric(net_data$net_type) == 3,colours_nets[3],colours_nets[4])))

for(i in 1:104){
  # Confirm this parameter is specifying what we expect
  D = median(PARMS_USAGE3[,i])
  D_LOW = quantile(PARMS_USAGE3[,i],0.01)
  D_UPP = quantile(PARMS_USAGE3[,i],0.99)
  
  cover_at_start = 0.95
  cover_at_start_UPP = 1
  cover_at_start_LOW = 0.9
  
  ## The below replicate the way this parameter is included in the transmission model
  aa = cover_at_start * exp((-1/D)*time_m)
  aL = cover_at_start_LOW * exp((-1/D_LOW)*time_m)
  aU = cover_at_start_UPP * exp((-1/D_UPP)*time_m)
  
  # lines(aa ~ time_m,col="blue")
  # lines(aL ~ time_m,col="blue",lty=3)
  # lines(aU ~ time_m,col="blue",lty=3)
  # polygon(c(time_m,rev(time_m)),
  #         c(aL,rev(aU)),col=adegenet::transp("blue",0.05),border = NA) 
  lines(aa ~ time_m,col=cols_nets[net_timing$net_type_INPUTA[i]])
  
}
for(i in 1:104){
  standard_net_usage = c(0.95,as.numeric(net_data[i,c(4,6,8)]))
  points(standard_net_usage[1:4] ~ time_obs[1:4],pch=19,col=colour_match[i])
  
}
D = colMeans(PARMS_USAGE3)
# colours_nets = c(adegenet::transp("orange",0.2),"orange",adegenet::transp("darkgreen"),"darkgreen")
as.numeric(net_data$net_type)

colour_match = ifelse(as.numeric(net_data$net_type) == 1,colours_nets[1],
                      ifelse(as.numeric(net_data$net_type) == 2,colours_nets[2],
                             ifelse(as.numeric(net_data$net_type) == 3,colours_nets[3],colours_nets[4])))
hist(D,breaks = 104,main="",yaxt="n",border=F,
     ylab = "Frequency",xlab="Adherence to net use (years)",col=colour_match)


abline(v=median(D),col="blue",lwd=2)
axis(2,las=2,at=c(0:7))
abline(v=quantile(D,0.9),col="blue",lwd=2,lty=2)
abline(v=quantile(D,0.1),col="blue",lwd=2,lty=2)


###########################
##
##
## Alternatives

hist(data$sleeping_under_net_6m,
     ylab="Frequency",main="",
     xlab="Proportion of people using nets",
     border = "darkblue",lty=1,
     col = adegenet::transp("darkblue",0.4),bty="n",
     breaks = 15,xlim=c(0,1),
     yaxt="n",xaxt="n")
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(2,las=2,at=c(0,5,10,15))

hist(data$sleeping_under_net_12m,
     ylab="",main="",
     xlab="",
     border = "blue",lty=2,
     col = adegenet::transp("blue",0.4),bty="n",
     breaks = 20,xlim=c(0,1),add=T,
     yaxt="n",xaxt="n")

hist(data$sleeping_under_net_18m,
     ylab="",main="",
     xlab="",
     border = "lightblue",lty=3,
     col = adegenet::transp("lightblue",0.4),bty="n",
     breaks = 25,xlim=c(0,1),add=T,
     yaxt="n",xaxt="n")

# hist(data$sleeping_under_net_25m,
#      ylab="",main="",
#      xlab="",
#      border = adegenet::transp("grey",0.1),lty=3,
#      col = adegenet::transp("lightblue",0.4),bty="n",
#      breaks = 25,xlim=c(0,1),add=T,
#      yaxt="n",xaxt="n")

legend("topleft",legend=c("After 6 months",
                          "After 12 months",
                          "After 18 months"),
       col=adegenet::transp(c("darkblue","blue","lightblue"),0.4),
       pch=15,lty=c(1,2,3),bty = "n")




NET_COUNT = as.numeric(tapply(data$Distribution_Month,data$Net_Type,length))
NET_NAMES = unique(sort(data$Net_Type))
pie(x=NET_COUNT, labels = NET_NAMES, radius = 1,
    col = c("darkred","orange","darkblue","darkgreen"),cex=1.5,
    main = "Proportion of clusters with each net", clockwise = T)
text(-0.4,0.4,"n = 32",col="white",cex=1.5)
text(0.3,0.7,"n = 13",col="white",cex=1.5)
text(0.6,0.1,"n = 19",col="black",cex=1.5)
text(-0.1,-0.5,"n = 40",col="white",cex=1.5)


