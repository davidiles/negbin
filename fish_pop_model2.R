library(jagsUI)
library(pscl)

rm(list=ls())

###########################################
# Read in trawl data
###########################################
yc=read.csv("../../data/CPUE_ycs_trawl_withcobb_BLS_2017.csv",header=TRUE)
head(yc)

###########################################
# Reorganize
########################################### 
yc$logoffset=log(yc$numtrawlyr)
head(yc)
names(yc)

##fixed effects
yc$fyear=as.factor(yc$year_cap)
yc$flocation=as.factor(yc$location)
yc$fdepth=as.factor(yc$trawldepthft)

#Check out patterns in raw data
rawdata_plot = ggplot(data = yc, aes(x = year_cap, y = countinyear, col = fdepth, shape = fdepth)) +
  geom_jitter(width = 0.2, height = 0)+
  facet_grid(assnage~., scales = "free")+
  theme_bw()
print(rawdata_plot)

# Create dummy variables for fdepth
yc$fdepth60 = yc$fdepth90 = yc$fdepth120 = 0
yc$fdepth60[yc$fdepth == "60"] = 1
yc$fdepth90[yc$fdepth == "90"] = 1
yc$fdepth120[yc$fdepth == "120"] = 1

# Create empty array to store age-specific abundances each year
obs_year_age = table(yc[,c("year_cap","assnage")]) #number of observations per age class per year
years = min(yc$year_cap):max(yc$year_cap)
years = c(years[1] - 1, years)

### Array to store counts each year
# First dimension of array is years of study
# Second dimension of array is age classes
# Third dimension of array is repeated counts for that age in that year (currently, up to 17)
N_array = array(NA, dim = c(length(years), length(unique(yc$assnage)), max(obs_year_age)))

# Arrays to store covariates
yc$fdepth_fac = as.numeric(as.factor(yc$fdepth))
fdepth_array = N_array

for (y in unique(yc$year_cap)){
  for (age in unique(yc$assnage)){
    
    year_age_dat = subset(yc, year_cap == y & assnage == age)
    
    for (i in 1:nrow(year_age_dat)){
      array_dim1 = which(years == y)
      array_dim2 = which(sort(unique(yc$assnage)) == age)
      array_dim3 = i
      
      N_array[array_dim1, array_dim2, array_dim3] = year_age_dat$countinyear[i]
      fdepth_array[array_dim1, array_dim2, array_dim3] = year_age_dat$fdepth_fac[i]
      
      
    }
  }
}

#Provide fdpeth in years when no observations were recorded (this does not affect results)
#Analysis requires covariates for every observation (including missing observations)
fdepth_array[is.na(fdepth_array)] = 1 

#Check out observed counts for each year and age class
N_array[1,,] #Year 1 (2000); 17 obs for each age class
N_array[2,,] #Year 2 (2001); no data
N_array[3,,] #Year 3 (2002); 17 obs for each age class
N_array[4,,] #Year 4 (2003); no data
N_array[5,,] #Year 5 (2004); 16 obs for each age class
N_array[6,,] #Year 6 (2005); no data
N_array[7,,] #Year 7 (2006); 13 obs for each age class
N_array[8,,] #Year 8 (2007); no data
N_array[9,,] #Year 9 (2008); 17 obs for each age class
N_array[10,,] #Year 10 (2009); no data
N_array[11,,] #Year 11 (2010); 6 obs for each age class
N_array[12,,] #Year 12 (2011); no data
N_array[13,,] #Year 13 (2012); 13-15 obs for each age class
N_array[14,,] #Year 14 (2013); no data
N_array[15,,] #Year 15 (2014); 17 obs for each age class
N_array[16,,] #Year 16 (2015); 10 obs for each age class
N_array[17,,] #Year 17 (2016); 4 obs for each age class
N_array[18,,] #Year 18 (2017); 6 obs for each age class

# Specify model in JAGS language
sink("fish_pop_model.txt")
cat("
    model{
    
    #-----------------------------------------
    # Priors
    #-----------------------------------------
    
    # Reproduction
    reprod.2.intercept ~ dunif(0,100)
    log.reprod.2.intercept <- log(reprod.2.intercept)
    
    reprod.3.intercept ~ dunif(0,100)
    log.reprod.3.intercept <- log(reprod.3.intercept)
    
    reprod.4.intercept ~ dunif(0,100)
    log.reprod.4.intercept <- log(reprod.4.intercept)
    
    # Temporal variance in reproduction (not explained by covariates)
    log.reprod.2.sd ~ dunif(0,5)
    log.reprod.2.tau <- pow(log.reprod.2.sd,-2)
    log.reprod.3.sd ~ dunif(0,5)
    log.reprod.3.tau <- pow(log.reprod.3.sd,-2)
    log.reprod.4.sd ~ dunif(0,5)
    log.reprod.4.tau <- pow(log.reprod.4.sd,-2)
    
    # Survival
    surv.0.intercept ~ dunif(0,1)
    logit.surv.0.intercept <- log(surv.0.intercept/(1-surv.0.intercept))
    
    surv.1.intercept ~ dunif(0,1)
    logit.surv.1.intercept <- log(surv.1.intercept/(1-surv.1.intercept))
    
    surv.2.intercept ~ dunif(0,1)
    logit.surv.2.intercept <- log(surv.2.intercept/(1-surv.2.intercept))
    
    surv.3.intercept ~ dunif(0,1)
    logit.surv.3.intercept <- log(surv.3.intercept/(1-surv.3.intercept))
    
    # Temporal variance in survival (not explained by covariates)
    logit.surv.0.sd ~ dunif(0,10)
    logit.surv.0.tau <- pow(logit.surv.0.sd,-2)
    logit.surv.1.sd ~ dunif(0,10)
    logit.surv.1.tau <- pow(logit.surv.1.sd,-2)
    logit.surv.2.sd ~ dunif(0,10)
    logit.surv.2.tau <- pow(logit.surv.2.sd,-2)
    logit.surv.3.sd ~ dunif(0,10)
    logit.surv.3.tau <- pow(logit.surv.3.sd,-2)

    #Priors for observation effects
    #Effect of sampling depth on observed count
    beta.fdepth[1] <- 0
    beta.fdepth[2] ~ dnorm(0,0.01) #Link this to mu through log link
    beta.fdepth[3] ~ dnorm(0,0.01) #Link this to mu through log link
    beta.fdepth[4] ~ dnorm(0,0.01) #Link this to mu through log link
    
    # Each age has a different rate parameter (r) and a different zero-term (these affect variance in observed counts)
    for (age in 1:A){
      r[age] ~ dunif(0,50)
      zero.intercept.age[age] ~ dunif(0.001,0.999)
      zero.B0[age] <- log(zero.intercept.age[age]/(1-zero.intercept.age[age]))
      
    }
    
    # FIRST YEAR ABUNDANCE - STATE PROCESS
    mu.count[1,1] ~ dunif(0,10000)
    mu.count[1,2] ~ dunif(0,10000)
    mu.count[1,3] ~ dunif(0,10000)
    mu.count[1,4] ~ dunif(0,10000)
    mu.count[1,5] ~ dunif(0,10000)
    
    #-----------------------------------------
    # LIKELIHOOD
    #-----------------------------------------
    
    #Time-varying vital rates
    for (year in 1:Y){
    
      #Time-varying age-specific reproduction
      log.reprod.2[year] ~ dnorm(log.reprod.2.intercept, log.reprod.2.tau)
      reprod.2[year] <- exp(log.reprod.2[year])
      log.reprod.3[year] ~ dnorm(log.reprod.3.intercept, log.reprod.3.tau)
      reprod.3[year] <- exp(log.reprod.3[year])
      log.reprod.4[year] ~ dnorm(log.reprod.4.intercept, log.reprod.4.tau)
      reprod.4[year] <- exp(log.reprod.4[year])
      
      #Time-varying age-specific survival
      logit.surv.0[year] ~ dnorm(logit.surv.0.intercept, logit.surv.0.tau)
      surv.0[year] <- 1/(1+exp(logit.surv.0[year]))
      logit.surv.1[year] ~ dnorm(logit.surv.1.intercept, logit.surv.1.tau)
      surv.1[year] <- 1/(1+exp(logit.surv.1[year]))
      logit.surv.2[year] ~ dnorm(logit.surv.2.intercept, logit.surv.2.tau)
      surv.2[year] <- 1/(1+exp(logit.surv.2[year]))
      logit.surv.3[year] ~ dnorm(logit.surv.3.intercept, logit.surv.3.tau)
      surv.3[year] <- 1/(1+exp(logit.surv.3[year]))
      
    }
    
    # FIRST YEAR
    
    # Observations are contingent upon current mean (mu.count)
    for (obs in 1:O){
    
    N[1,1,obs] ~ dnegbin(p[1,1,obs],r[1]) #Age 0 individuals
    N[1,2,obs] ~ dnegbin(p[1,2,obs],r[2]) #Age 1 individuals
    N[1,3,obs] ~ dnegbin(p[1,3,obs],r[3]) #Age 2 individuals
    N[1,4,obs] ~ dnegbin(p[1,4,obs],r[4]) #Age 3 individuals
    N[1,5,obs] ~ dnegbin(p[1,5,obs],r[5]) #Age 4 individuals
    
    # Zero-Inflation for observations
    zero[1,1,obs] ~ dbern(prob[1,1,obs])
    prob[1,1,obs] <- 1/(1+exp(-logit.prob[1,1,obs]))
    logit.prob[1,1,obs] <-  zero.B0[1]
    
    zero[1,2,obs] ~ dbern(prob[1,2,obs])
    prob[1,2,obs] <- 1/(1+exp(-logit.prob[1,2,obs]))
    logit.prob[1,2,obs] <-  zero.B0[2]
    
    zero[1,3,obs] ~ dbern(prob[1,3,obs])
    prob[1,3,obs] <- 1/(1+exp(-logit.prob[1,3,obs]))
    logit.prob[1,3,obs] <-  zero.B0[3]
    
    zero[1,4,obs] ~ dbern(prob[1,4,obs])
    prob[1,4,obs] <- 1/(1+exp(-logit.prob[1,4,obs]))
    logit.prob[1,4,obs] <-  zero.B0[4]
    
    zero[1,5,obs] ~ dbern(prob[1,5,obs])
    prob[1,5,obs] <- 1/(1+exp(-logit.prob[1,5,obs]))
    logit.prob[1,5,obs] <-  zero.B0[5]
    
    #Adjust mu for effort covariates(fdepth)
    mu.count.adjusted.for.effort[1,1,obs] <- exp(log(mu.count[1,1]) + beta.fdepth[fdepth[1, 1, obs]])
    mu.count.adjusted.for.effort[1,2,obs] <- exp(log(mu.count[1,2]) + beta.fdepth[fdepth[1, 2, obs]])
    mu.count.adjusted.for.effort[1,3,obs] <- exp(log(mu.count[1,3]) + beta.fdepth[fdepth[1, 3, obs]])
    mu.count.adjusted.for.effort[1,4,obs] <- exp(log(mu.count[1,4]) + beta.fdepth[fdepth[1, 4, obs]])
    mu.count.adjusted.for.effort[1,5,obs] <- exp(log(mu.count[1,5]) + beta.fdepth[fdepth[1, 5, obs]])
    
    #Calculate the probability parameter (one of the two key components of the ZINB) based on mu, r, and zero
    p[1,1,obs] <- r[1]/(r[1]+(1-zero[1,1,obs])* mu.count.adjusted.for.effort[1,1,obs]) - 1e-10*zero[1,1,obs]
    p[1,2,obs] <- r[2]/(r[2]+(1-zero[1,2,obs])* mu.count.adjusted.for.effort[1,2,obs]) - 1e-10*zero[1,2,obs]
    p[1,3,obs] <- r[3]/(r[3]+(1-zero[1,3,obs])* mu.count.adjusted.for.effort[1,3,obs]) - 1e-10*zero[1,3,obs]
    p[1,4,obs] <- r[4]/(r[4]+(1-zero[1,4,obs])* mu.count.adjusted.for.effort[1,4,obs]) - 1e-10*zero[1,4,obs]
    p[1,5,obs] <- r[5]/(r[5]+(1-zero[1,5,obs])* mu.count.adjusted.for.effort[1,5,obs]) - 1e-10*zero[1,5,obs]
    
    }
    
    # SUBSEQUENT YEARS
    for (year in 2:Y){
    
    #Population Model
    mu.count[year,1] <- mu.count[year-1,3]*reprod.2[year-1] + mu.count[year-1,4]*reprod.3[year-1] + mu.count[year-1,5]*reprod.4[year-1]
    
    mu.count[year,2] <- mu.count[year-1,1]*surv.0[year-1]
    mu.count[year,3] <- mu.count[year-1,2]*surv.1[year-1]
    mu.count[year,4] <- mu.count[year-1,3]*surv.2[year-1]
    mu.count[year,5] <- mu.count[year-1,4]*surv.3[year-1]
    
    # OBSERVATION PROCESS (ZERO INFLATED NEGATIVE BINOMIAL)
    #Observations are contingent on current mean for each age class
    for (obs in 1:O){
    
    # Zero-Inflation for observations
    zero[year,1,obs] ~ dbern(prob[year,1,obs])
    prob[year,1,obs] <- 1/(1+exp(-logit.prob[year,1,obs]))
    logit.prob[year,1,obs] <-  zero.B0[1]
    
    zero[year,2,obs] ~ dbern(prob[year,2,obs])
    prob[year,2,obs] <- 1/(1+exp(-logit.prob[year,2,obs]))
    logit.prob[year,2,obs] <-  zero.B0[2]
    
    zero[year,3,obs] ~ dbern(prob[year,3,obs])
    prob[year,3,obs] <- 1/(1+exp(-logit.prob[year,3,obs]))
    logit.prob[year,3,obs] <-  zero.B0[3]
    
    zero[year,4,obs] ~ dbern(prob[year,4,obs])
    prob[year,4,obs] <- 1/(1+exp(-logit.prob[year,4,obs]))
    logit.prob[year,4,obs] <-  zero.B0[4]
    
    zero[year,5,obs] ~ dbern(prob[year,5,obs])
    prob[year,5,obs] <- 1/(1+exp(-logit.prob[year,5,obs]))
    logit.prob[year,5,obs] <-  zero.B0[5]
    
    #Adjust mu for effort covariates(fdepth)
    mu.count.adjusted.for.effort[year,1,obs] <- exp(log(mu.count[year,1]) + beta.fdepth[fdepth[year, 1, obs]])
    mu.count.adjusted.for.effort[year,2,obs] <- exp(log(mu.count[year,2]) + beta.fdepth[fdepth[year, 2, obs]])
    mu.count.adjusted.for.effort[year,3,obs] <- exp(log(mu.count[year,3]) + beta.fdepth[fdepth[year, 3, obs]])
    mu.count.adjusted.for.effort[year,4,obs] <- exp(log(mu.count[year,4]) + beta.fdepth[fdepth[year, 4, obs]])
    mu.count.adjusted.for.effort[year,5,obs] <- exp(log(mu.count[year,5]) + beta.fdepth[fdepth[year, 5, obs]])
    
    #Probability parameter (one of the two key components of the ZINB) based on mu, r, and zero
    p[year,1,obs] <- r[1]/(r[1]+(1-zero[year,1,obs])* mu.count.adjusted.for.effort[year,1,obs]) - 1e-10*zero[year,1,obs]
    p[year,2,obs] <- r[2]/(r[2]+(1-zero[year,2,obs])* mu.count.adjusted.for.effort[year,2,obs]) - 1e-10*zero[year,2,obs]
    p[year,3,obs] <- r[3]/(r[3]+(1-zero[year,3,obs])* mu.count.adjusted.for.effort[year,3,obs]) - 1e-10*zero[year,3,obs]
    p[year,4,obs] <- r[4]/(r[4]+(1-zero[year,4,obs])* mu.count.adjusted.for.effort[year,4,obs]) - 1e-10*zero[year,4,obs]
    p[year,5,obs] <- r[5]/(r[5]+(1-zero[year,5,obs])* mu.count.adjusted.for.effort[year,5,obs]) - 1e-10*zero[year,5,obs]
    
    N[year,1,obs] ~ dnegbin(p[year,1,obs],r[1])
    N[year,2,obs] ~ dnegbin(p[year,2,obs],r[2])
    N[year,3,obs] ~ dnegbin(p[year,3,obs],r[3])
    N[year,4,obs] ~ dnegbin(p[year,4,obs],r[4])
    N[year,5,obs] ~ dnegbin(p[year,5,obs],r[5])
    
    }
    }
    
    
    
    }
    ",fill = TRUE)
sink()

# Bundle data
jags.data <- list(N = N_array,
                  Y = dim(N_array)[1],
                  A = dim(N_array)[2],
                  O = dim(N_array)[3],
                  fdepth = fdepth_array)

# Initial values
inits <- function() list(zero = array(0,dim = dim(N_array)),
                         r = c(0.35, 1.2, 1.2, 0.55, 0.55))

# Parameters monitored
params <- c("r","zero.B0","beta.fdepth",
            
            "reprod.2.intercept","reprod.3.intercept","reprod.4.intercept",
            "log.reprod.2.sd","log.reprod.3.sd","log.reprod.4.sd",
            
            "surv.0.intercept","surv.1.intercept","surv.2.intercept","surv.3.intercept",
            "logit.surv.0.sd","logit.surv.1.sd","logit.surv.2.sd","logit.surv.3.sd",
            
            "reprod.2","reprod.3","reprod.4",
            "surv.0","surv.1","surv.2","surv.3","surv.4",
            "mu.count"
)

# MCMC settings
ni <- 500
nt <- 2
nb <- 100
nc <- 2

# Call JAGS from R
out <- jags(data = jags.data, 
            inits, 
            params, 
            model.file = "fish_pop_model.txt", 
            n.chains = nc, 
            n.thin = nt, 
            n.iter = ni, 
            n.burnin = nb)

out

mu.1.med = apply(out$sims.list$mu.count[,,1],2,function(x) quantile(x,0.500))
mu.1.lcl = apply(out$sims.list$mu.count[,,1],2,function(x) quantile(x,0.025))
mu.1.ucl = apply(out$sims.list$mu.count[,,1],2,function(x) quantile(x,0.975))

mu.2.med = apply(out$sims.list$mu.count[,,2],2,function(x) quantile(x,0.500))
mu.2.lcl = apply(out$sims.list$mu.count[,,2],2,function(x) quantile(x,0.025))
mu.2.ucl = apply(out$sims.list$mu.count[,,2],2,function(x) quantile(x,0.975))

mu.3.med = apply(out$sims.list$mu.count[,,3],2,function(x) quantile(x,0.500))
mu.3.lcl = apply(out$sims.list$mu.count[,,3],2,function(x) quantile(x,0.025))
mu.3.ucl = apply(out$sims.list$mu.count[,,3],2,function(x) quantile(x,0.975))

mu.4.med = apply(out$sims.list$mu.count[,,4],2,function(x) quantile(x,0.500))
mu.4.lcl = apply(out$sims.list$mu.count[,,4],2,function(x) quantile(x,0.025))
mu.4.ucl = apply(out$sims.list$mu.count[,,4],2,function(x) quantile(x,0.975))

mu.5.med = apply(out$sims.list$mu.count[,,5],2,function(x) quantile(x,0.500))
mu.5.lcl = apply(out$sims.list$mu.count[,,5],2,function(x) quantile(x,0.025))
mu.5.ucl = apply(out$sims.list$mu.count[,,5],2,function(x) quantile(x,0.975))


library(ggplot2)
p1 = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.1.lcl, ymax = mu.1.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.1.med)) +
  ylab("Rel. abundance")+
  theme_bw()

p2 = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.2.lcl, ymax = mu.2.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.2.med)) +
  ylab("Rel. abundance")+
  theme_bw()

p3 = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.3.lcl, ymax = mu.3.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.3.med)) +
  ylab("Rel. abundance")+
  theme_bw()

p4 = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.4.lcl, ymax = mu.4.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.4.med)) +
  ylab("Rel. abundance")+
  theme_bw()

p5 = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.5.lcl, ymax = mu.5.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.5.med)) +
  ylab("Rel. abundance")+
  theme_bw()

library(cowplot)

fullplot = plot_grid(p1,p2,p3,p4,p5, nrow=5)
print(fullplot)
