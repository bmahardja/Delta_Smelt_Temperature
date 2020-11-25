library(tidyverse)
library(lubridate)
library(readxl)
library(jagsUI)

#Path to local drive
root <- "~/GitHub/Delta_Smelt_Temperature"
setwd(root)

data_root<-file.path(root,"data_raw")
code_root <- file.path(root,"code")
output_root <- file.path(root,"output")

TNS <- read_excel(file.path(data_root, "Townet_Data_1959-2017.xlsx"),sheet="CatchPerTow",col_types=c("numeric","numeric","guess",rep("numeric",97)))
TNS$SampleDate<-as.Date(TNS$SampleDate,"%Y-%m-%d")
str(TNS)

TNS_clean <- TNS  %>% rename(Delta.Smelt='Delta Smelt',VolSampled=TowVolm3)
TNS_clean <- TNS_clean[complete.cases(TNS_clean[ , c("TemperatureTop","ConductivityTop","Secchi","VolSampled")]),]


TNS_clean$julianday<-yday(TNS_clean$SampleDate)

TNS_smelt<- TNS_clean %>% select(Year,Survey,StationCode,SampleDate,TowNumber,Index,Secchi,TemperatureTop,ConductivityTop,Delta.Smelt,julianday) %>% 
  mutate(Delta.Smelt=ifelse(Delta.Smelt>0,1,0))
TNS_smelt <- unique(TNS_smelt)

TNS_smelt<- spread(TNS_smelt,TowNumber,Delta.Smelt) %>% rename(Tow1='1',Tow2='2',Tow3='3')

plot(TNS_smelt$julianday,TNS_smelt$TemperatureTop)
#Remove apparent erroneous reading of ~35 C
TNS_smelt <- TNS_smelt %>% filter(TemperatureTop<31)

#Create simple model to adjust temperature by julian day, probably quadratic function is appropriate based on plot
TNS_smelt$julianday_quad<-(TNS_smelt$julianday)^2
model_temp<-lm(TemperatureTop~julianday+julianday_quad,data=TNS_smelt)
summary(model_temp)
#Not a good fit, but we can test it out anyways
TNS_smelt$PredictedTemp<-predict(model_temp,data=TNS_smelt)
TNS_smelt$Temp_anomaly<-TNS_smelt$PredictedTemp-TNS_smelt$TemperatureTop
plot(TNS_smelt$julianday,TNS_smelt$Temp_anomaly)

#Standardize covariates
TNS_smelt$Secchi_s<-(TNS_smelt$Secchi-mean(TNS_smelt$Secchi))/sd(TNS_smelt$Secchi)
TNS_smelt$ConductivityTop_s<-(TNS_smelt$ConductivityTop-mean(TNS_smelt$ConductivityTop))/sd(TNS_smelt$ConductivityTop)
TNS_smelt$TemperatureTop_s<-(TNS_smelt$TemperatureTop-mean(TNS_smelt$TemperatureTop))/sd(TNS_smelt$TemperatureTop)
TNS_smelt$julianday_s<-(TNS_smelt$julianday-mean(TNS_smelt$julianday))/sd(TNS_smelt$julianday)
TNS_smelt$Year_s<-as.numeric(factor(TNS_smelt$Year))
TNS_smelt$StationCode_s<-as.numeric(as.factor(TNS_smelt$StationCode))
TNS_smelt$YearStation_s<-TNS_smelt$Year_s*TNS_smelt$StationCode_s

unique(TNS_smelt$Year_s)
unique(TNS_smelt$StationCode_s)

#To use complete data for now (i.e., no Tow = NA)
TNS_smelt <- TNS_smelt %>% filter(complete.cases(Tow1,Tow2,Tow3))

#Set up data for JAGS
k= 3 # number tows per survey
nobs=as.numeric(nrow(TNS_smelt))
no.site= as.numeric(length(unique(TNS_smelt$StationCode_s))) # number sites
no.yr= as.numeric(length(unique(TNS_smelt$Year_s))) # number years
no.yrsite= no.site*no.yr # number site survey combinations
Y=as.matrix(TNS_smelt[,c("Tow1","Tow2","Tow3")])
dimnames(Y) <- NULL

# create JAGS data input 
#jags_data<-list(nobs=nobs, k=k,no.yr=no.yr,no.site=no.site, no.yrsite=no.yrsite, 
#            yr=TNS_smelt$Year_s,site=TNS_smelt$StationCode_s, yrsite=TNS_smelt$YearStation_s, detect=detect, 
#            covar.Y = dat$covar.Y,covar.X = dat$covar.X)

jags_data<-list(Y=Y, nobs=nobs ,no.yr=no.yr ,no.site=no.site, 
                k=k, julianday_s=TNS_smelt$julianday_s,Secchi_s=TNS_smelt$Secchi_s, TemperatureTop_s=TNS_smelt$TemperatureTop_s, 
                ConductivityTop_s=TNS_smelt$ConductivityTop_s,
                site=TNS_smelt$StationCode_s, yr=TNS_smelt$Year_s)

cat(file = "model.txt","
model {

# Priors

beta.occ.Secchi ~ dnorm(0, 1)
alpha.p ~ dnorm(0,1)
beta.p.Secchi ~ dnorm(0, 1)
beta.occ.TemperatureTop ~ dnorm(0, 1)
beta.occ.ConductivityTop ~ dnorm(0, 1)
beta.occ.julianday ~ dnorm(0, 1)

# Likelihood
 for (i in 1:nobs) { #start initial loop over the nobs sites
 # True state model for the partially observed true state
    z[i] ~ dbern(psi[i])		# True occupancy z at site i
    logit(psi[i]) <- beta.occ.julianday * julianday_s[i] + beta.occ.Secchi * Secchi_s[i] + beta.occ.TemperatureTop * TemperatureTop_s[i] + 
    beta.occ.ConductivityTop * ConductivityTop_s[i] + alpha.occ.Y[yr[i]] + alpha.occ.S[site[i]]

    for (j in 1:k) { # start a second loop over the k replicates
       # Observation model for the actual observations
       y[i,j] ~ dbern(eff.p[i,j])	# Detection-nondetection at i and j
       eff.p[i,j] <- z[i] * p[i,j]
       logit(p[i,j]) <- alpha.p + beta.p.Secchi * Secchi_s[i]

       # Computation of fit statistic (for Bayesian p-value)
       Presi[i,j] <- abs(y[i,j]-p[i,j])	 # Absolute residual
       y.new[i,j]~dbern(eff.p[i,j])
       Presi.new[i,j] <- abs(y.new[i,j]-p[i,j])
    
    # Random Effects
    # Year
    for(a in 1:no.yr){
    alpha.occ.Y[a] ~ dnorm(mu.int.Y,tau.int.Y)
    }
    # Site
    for(b in 1:no.site){
    alpha.occ.S[b] ~ dnorm(mu.int.S,tau.int.S)
    }
    #hyperparameters for random effects
    # Year psi
    mu.int.Y ~ dnorm(0, 0.001)
    tau.int.Y <- 1 / (sigma.int.Y * sigma.int.Y)
    sigma.int.Y ~ dunif(0, 10)

    # Site psi
    mu.int.S ~ dnorm(0, 0.001)
    tau.int.S <- 1 / (sigma.int.S * sigma.int.S)
    sigma.int.S ~ dunif(0, 10)
    }
 }

fit <- sum(Presi[,])# Discrepancy for actual data set
fit.new <- sum(Presi.new[,]) 		# Discrepancy for replicate data set

# Derived quantities
 occ.fs <- sum(z[])			# Number of occupied sites among 
}
")

# Inits function
zst <- apply(Y, 1, max)			# Good starting values for latent states essential !
inits <- function(){list(z = zst, beta.occ.julianday = runif(1, -5, 5), 
                         beta.occ.Secchi = runif(1, -5, 5),beta.occ.TemperatureTop = runif(1, -5, 5),
                         beta.occ.ConductivityTop = runif(1, -5, 5),mu.int.Y=0, sigma.int.Y =1,mu.int.S=0, sigma.int.S =1,
                         alpha.p = runif(1, -5, 5), beta.p.Secchi = runif(1, -5, 5))}

# Parameters to estimate
params <- c("mu.int.Y","sigma.int.Y","mu.int.S","sigma.int.S","beta.occ.julianday","beta.occ.Secchi","beta.occ.TemperatureTop",
            "beta.occ.ConductivityTop", "alpha.p", "beta.p.Secchi", "occ.fs", "fit", "fit.new")

# MCMC settings
na <- 1000  ;  nc <- 3  ;  ni <- 12000  ;  nb <- 2000  ;  nt <- 5

# Call JAGS, check convergence and summarize posteriors
out <- jags(jags_data, inits, params, "model.txt", n.adapt = na, n.thin = nt, n.chains = nc, 
            n.burnin = nb, n.iter = ni, parallel = TRUE)



par(mfrow = c(3, 3))  ;  traceplot(out)  ;  par(mfrow = c(1, 1))
print(out, dig = 3)














########################################################
# Peterson and Barajas 2018 code
########################################################


########################################################
# Create mock data
########################################################
k= 3 # number tows per survey
no.surv= 3 # number surveys per year
no.site= 5 # number sites
no.yr= 10 # number years
no.yrsite= no.site*no.yr # number site survey combinations
nobs= no.yrsite*no.surv # number of observations
# create data with year, station, and year by station unique identifiers
dat <- expand.grid(Yr.no=c(1:no.yr),Stn.no=c(1:no.site),Surv=c(1:no.surv))

dat$YrStn.no<-with(dat,as.numeric(as.factor(paste(Yr.no,Stn.no))))

# add covariates to dataframe 
dat$covar.X<-rnorm(nrow(dat))
dat$covar.Y<-rnorm(nrow(dat))

Y<-rnorm(nrow(dat))
# generate detection history assuming 3 tows 
detect<-cbind(round(runif(nrow(dat),0.6,3.4)),round(runif(nrow(dat),0.6,3.4)), round(runif(nrow(dat),0.6,3.4)))
###### End of mock data creation
# create JAGS data input 
jdata<-list(nobs=nobs, k=k,no.yr=no.yr,no.site=no.site, no.yrsite=no.yrsite, 
            yr=dat$Yr.no,site=dat$Stn.no, yrsite=dat$YrStn.no, detect=detect, 
            covar.Y = dat$covar.Y,covar.X = dat$covar.X)

# Parameters to monitor
params<-c("eta.Y.1.mu",
           "eta.Y.1.ss",
           "eta.S.1.ss",
           "eta.YS.1.ss",
           "eta.Y.2.mu",
           "eta.Y.2.ss",
           "eta.S.2.ss",
           "eta.YS.2.ss",
           "gam",
           "gam2",
           "beta.p",
           "delta.p",
           "abund.det")
# invoke JAGS function 
out.put<-jags(jdata, inits=inits,
  parameters.to.save =params, parallel = TRUE,
  model.file= "jag.model", n.chains = 3, n.thin = 1,
  n.iter = 5000,
  n.burnin = 1000, n.adapt=3000) summary(out.put)