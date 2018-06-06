library(jagsUI)
library(pscl)
library(ggplot2)

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
  facet_grid(assnage~.)+
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
years = c(years, (max(years)+1):(max(years)+10))

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

cobble_year = aggregate(spawncob ~ year_hatch, FUN = mean, data = yc)
cobble_year = subset(cobble_year, year_hatch %in% years)
cobble_year = cobble_year[,2]
cobble_year = (cobble_year - mean(cobble_year))/sd(cobble_year)
cobble_year = c(cobble_year, rep(min(cobble_year),10)) #Use last cobble value for 2018
plot(cobble_year ~ years, type = "o")

# Specify model in JAGS language
sink("fish_pop_model.txt")
cat("
    model{
    
    #-----------------------------------------
    # Priors
    #-----------------------------------------
    
    #
    prob_detect_0 ~ dunif(0,1)
    
    # Reproduction
    reprod.2.intercept ~ dunif(0.001,50)
    reprod.3.intercept ~ dunif(0.001,50)
    reprod.4.intercept ~ dunif(0.001,50)
    
    log.reprod.2.intercept <- log(reprod.2.intercept)
    log.reprod.3.intercept <- log(reprod.3.intercept)
    log.reprod.4.intercept <- log(reprod.4.intercept)
    
    # Survival
    surv.0.intercept ~ dunif(0,1)
    surv.1.intercept ~ dunif(0,1)
    surv.2.intercept ~ dunif(0,1)
    surv.3.intercept ~ dunif(0,1)
    
    logit.surv.0.intercept <- log(surv.0.intercept/(1-surv.0.intercept))
    logit.surv.1.intercept <- log(surv.1.intercept/(1-surv.1.intercept))
    logit.surv.2.intercept <- log(surv.2.intercept/(1-surv.2.intercept))
    logit.surv.3.intercept <- log(surv.3.intercept/(1-surv.3.intercept))
    
    #effect of cobble
    surv.0.beta.cobble ~ dnorm(0,0.01)
    surv.1.beta.cobble ~ dnorm(0,0.01)
    surv.2.beta.cobble ~ dnorm(0,0.01)
    surv.3.beta.cobble ~ dnorm(0,0.01)
    
    reprod.2.beta.cobble ~ dnorm(0,0.01)
    reprod.3.beta.cobble ~ dnorm(0,0.01)
    reprod.4.beta.cobble ~ dnorm(0,0.01)
    
    #Priors for observation effects
    #Effect of sampling depth on observed count
    beta.fdepth[1] <- 0
    beta.fdepth[2] ~ dnorm(0,0.01) #Link this to mu through log link
    beta.fdepth[3] ~ dnorm(0,0.01) #Link this to mu through log link
    beta.fdepth[4] ~ dnorm(0,0.01) #Link this to mu through log link
    
    # Each age has a different rate parameter (r) and a different zero-term; together these affect variance in observed counts
    for (age in 1:A){
    r[age] ~ dunif(0,10)
    zero.intercept.age[age] ~ dunif(0.001,0.999)
    zero.B0[age] <- log(zero.intercept.age[age]/(1-zero.intercept.age[age]))
    }

    # FIRST YEAR ABUNDANCE - STATE PROCESS
    mu.count[1,1] <- mu.count[1,3]*reprod.2[1] + mu.count[1,4]*reprod.3[1] + mu.count[1,5]*reprod.3[1]
    
    mu.count[1,2] ~ dunif(0.01,500)
    mu.count[1,3] ~ dunif(0.01,500)
    mu.count[1,4] ~ dunif(0.01,500)
    mu.count[1,5] ~ dunif(0.01,500)
    
    #Only a fraction are actually available to be observed (relative detection probability)
    mu.obs[1,1] <- mu.count[1,1]*prob_detect_0
    
    sd.reprod ~ dunif(0,4)
    tau.reprod <- pow(sd.reprod,-2)
    
    sd.surv ~ dunif(0,5)
    tau.surv <- pow(sd.surv,-2)
    
    #-----------------------------------------
    # LIKELIHOOD
    #-----------------------------------------
    
    #Time-varying vital rates
    for (year in 1:Y){
    
    log.reprod.2[year] ~ dnorm(log.reprod.2.intercept + reprod.2.beta.cobble * cobble_year[year], tau.reprod)
    log.reprod.3[year] ~ dnorm(log.reprod.3.intercept + reprod.3.beta.cobble * cobble_year[year], tau.reprod)
    log.reprod.4[year] ~ dnorm(log.reprod.4.intercept + reprod.4.beta.cobble * cobble_year[year], tau.reprod)
    
    reprod.2[year] <- exp(log.reprod.2[year])
    reprod.3[year] <- exp(log.reprod.3[year])
    reprod.4[year] <- exp(log.reprod.4[year])
    
    logit.surv.0[year] ~ dnorm(logit.surv.0.intercept + surv.0.beta.cobble * cobble_year[year], tau.surv)
    logit.surv.1[year] ~ dnorm(logit.surv.1.intercept + surv.1.beta.cobble * cobble_year[year], tau.surv)
    logit.surv.2[year] ~ dnorm(logit.surv.2.intercept + surv.2.beta.cobble * cobble_year[year], tau.surv)
    logit.surv.3[year] ~ dnorm(logit.surv.3.intercept + surv.3.beta.cobble * cobble_year[year], tau.surv)
    
    surv.0[year] <- 1/(1+exp(-logit.surv.0[year]))
    surv.1[year] <- 1/(1+exp(-logit.surv.1[year]))
    surv.2[year] <- 1/(1+exp(-logit.surv.2[year]))
    surv.3[year] <- 1/(1+exp(-logit.surv.3[year]))
    
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
    mu.count.adjusted.for.effort[1,1,obs] <- exp(log(mu.obs[1,1]) + beta.fdepth[fdepth[1, 1, obs]])
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
    mu.count[year,1] <- mu.count[year,3]*reprod.2[year] + mu.count[year,4]*reprod.3[year] + mu.count[year,5]*reprod.3[year]
    
    #Only a fraction are actually available to be observed (relative detection probability)
    mu.obs[year,1] <- mu.count[year,1]*prob_detect_0
    
    mu.count[year,2] <- mu.count[year-1,1]*surv.1[year-1]
    mu.count[year,3] <- mu.count[year-1,2]*surv.1[year-1]
    mu.count[year,4] <- mu.count[year-1,3]*surv.2[year-1]
    mu.count[year,5] <- mu.count[year-1,4]*surv.2[year-1]
    
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
    mu.count.adjusted.for.effort[year,1,obs] <- exp(log(mu.obs[year,1]) + beta.fdepth[fdepth[year, 1, obs]])
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
    
for (year in 1:Y){
    mu.breeding[year] <- mu.count[year,3] + mu.count[year,4] + mu.count[year,5]
}

    
    
    }
    ",fill = TRUE)
sink()

# Bundle data
jags.data <- list(N = N_array,
                  Y = dim(N_array)[1],
                  A = dim(N_array)[2],
                  O = dim(N_array)[3],
                  fdepth = fdepth_array,
                  cobble_year = cobble_year)

# Initial values
#mu.count.init = matrix(NA,ncol = 5, nrow = length(years))
#mu.count.init[1,] = c(500,100,100,50,NA)
inits <- function() list(zero = array(0,dim = dim(N_array)))

# Parameters monitored
params <- c("r","zero.B0","beta.fdepth",
            
            "prob_detect_0",
            
            "surv.0.intercept","surv.1.intercept","surv.2.intercept","surv.3.intercept",
            "surv.0.beta.cobble","surv.1.beta.cobble","surv.2.beta.cobble","surv.3.beta.cobble",
            
            "reprod.2.intercept","reprod.3.intercept","reprod.4.intercept",
            "reprod.2.beta.cobble","reprod.3.beta.cobble","reprod.4.beta.cobble",
            
            "sd.surv","sd.reprod",
            
            "surv.0","surv.1","surv.2","surv.3",
            "reprod.2","reprod.3","reprod.4",
            "mu.count",
            "mu.breeding"
)

# MCMC settings
ni <- 5000
nt <- 2
nb <- 4000
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

p1 = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.1.lcl, ymax = mu.1.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.1.med)) +
  geom_jitter(data = subset(yc, assnage == 0), aes(x = year_cap, y = countinyear), width = 0.2, height = 0)+
  coord_cartesian(ylim=c(0,300))+
  ylab("Age 0")+
  theme_bw()

p2 = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.2.lcl, ymax = mu.2.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.2.med)) +
  geom_jitter(data = subset(yc, assnage == 1), aes(x = year_cap, y = countinyear), width = 0.2, height = 0)+
  coord_cartesian(ylim=c(0,200))+
  ylab("Age 1")+
  theme_bw()

p3 = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.3.lcl, ymax = mu.3.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.3.med)) +
  geom_jitter(data = subset(yc, assnage == 2), aes(x = year_cap, y = countinyear), width = 0.2, height = 0)+
  coord_cartesian(ylim=c(0,200))+
  ylab("Age 2")+
  theme_bw()

p4 = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.4.lcl, ymax = mu.4.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.4.med)) +
  coord_cartesian(ylim=c(0,75))+
  geom_jitter(data = subset(yc, assnage == 3), aes(x = year_cap, y = countinyear), width = 0.2, height = 0)+
  ylab("Age 3")+
  theme_bw()

p5 = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.5.lcl, ymax = mu.5.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.5.med)) +
  geom_jitter(data = subset(yc, assnage == 4), aes(x = year_cap, y = countinyear), width = 0.2, height = 0)+
  #coord_cartesian(ylim=c(0,1000))+
  ylab("Age 4")+
  theme_bw()

cobble_plot = ggplot() +
  geom_line(aes(x = years, y = cobble_year), col = "blue", lwd = 1.5) +
  ylab("Submerged Cobble")+
  theme_bw()

library(cowplot)

fullplot = plot_grid(p1,p2,p3,p4,p5,cobble_plot, nrow=3)
print(fullplot)



mu.breeding.med = apply(out$sims.list$mu.breeding,2,function(x) quantile(x,0.500))

mu.breeding.lcl = apply(out$sims.list$mu.count[,,3] + out$sims.list$mu.count[,,4] + out$sims.list$mu.count[,,5],2,function(x) quantile(x,0.025))
mu.breeding.ucl = apply(out$sims.list$mu.count[,,3] + out$sims.list$mu.count[,,4] + out$sims.list$mu.count[,,5],2,function(x) quantile(x,0.975))


#total counts
yc$breeder = 0
yc$breeder[yc$assnage >=2] = 1

breeder_counts = aggregate(countinyear~year_cap + trawlnum, data = subset(yc, breeder == 1), FUN = sum)

cobble_plot = ggplot() +
  geom_line(aes(x = years, y = cobble_year), col = "blue", lwd = 1.5) +
  ylab("Submerged Cobble")+
  theme_bw()

pb = ggplot() +
  geom_ribbon(aes(x = years, ymin = mu.breeding.lcl, ymax = mu.breeding.ucl), alpha = 0.2) +
  geom_line(aes(x = years, y = mu.breeding.med)) +
  ylab("Breeding Population")+
  coord_cartesian(ylim=c(0,150))+
  geom_point(aes(x = unique(yc$year_cap), y = rep(0,length(unique(yc$year_cap)))), shape = 2)+
  theme_bw()

fullplot = plot_grid(cobble_plot, pb, nrow=2)
print(fullplot)


cobble.pred = seq(min(cobble_year),max(cobble_year),length.out = 100)

reprod.2.pred.med = cobble.pred*NA
reprod.2.pred.lcl = cobble.pred*NA
reprod.2.pred.ucl = cobble.pred*NA

reprod.3.pred.med = cobble.pred*NA
reprod.3.pred.lcl = cobble.pred*NA
reprod.3.pred.ucl = cobble.pred*NA

surv.1.pred.med = cobble.pred*NA
surv.1.pred.lcl = cobble.pred*NA
surv.1.pred.ucl = cobble.pred*NA

surv.2.pred.med = cobble.pred*NA
surv.2.pred.lcl = cobble.pred*NA
surv.2.pred.ucl = cobble.pred*NA

for (i in 1:length(reprod.2.pred.med)){
  reprod.2.pred.i = exp(log(out$sims.list$reprod.2.intercept) + out$sims.list$reprod.2.beta.cobble*cobble.pred[i])
  reprod.2.pred.med[i] <- quantile(reprod.2.pred.i,0.5)
  reprod.2.pred.lcl[i] <- quantile(reprod.2.pred.i,0.025)
  reprod.2.pred.ucl[i] <- quantile(reprod.2.pred.i,0.975)
  
  reprod.3.pred.i = exp(log(out$sims.list$reprod.3.intercept) + out$sims.list$reprod.3.beta.cobble*cobble.pred[i])
  reprod.3.pred.med[i] <- quantile(reprod.3.pred.i,0.5)
  reprod.3.pred.lcl[i] <- quantile(reprod.3.pred.i,0.025)
  reprod.3.pred.ucl[i] <- quantile(reprod.3.pred.i,0.975)
  
  surv.1.pred.i = plogis(qlogis(out$sims.list$surv.1.intercept) + out$sims.list$surv.1.beta.cobble*cobble.pred[i])
  surv.1.pred.med[i] <- quantile(surv.1.pred.i,0.5)
  surv.1.pred.lcl[i] <- quantile(surv.1.pred.i,0.025)
  surv.1.pred.ucl[i] <- quantile(surv.1.pred.i,0.975)
  
  surv.2.pred.i = plogis(qlogis(out$sims.list$surv.2.intercept) + out$sims.list$surv.2.beta.cobble*cobble.pred[i])
  surv.2.pred.med[i] <- quantile(surv.2.pred.i,0.5)
  surv.2.pred.lcl[i] <- quantile(surv.2.pred.i,0.025)
  surv.2.pred.ucl[i] <- quantile(surv.2.pred.i,0.975)
  
}

reprod.plot = ggplot() + 
  
  geom_ribbon(aes(x = cobble.pred, ymin = reprod.2.pred.lcl, ymax = reprod.2.pred.ucl), alpha = 0.2)+
  geom_line(aes(x = cobble.pred, y = reprod.2.pred.med))+
  
  geom_ribbon(aes(x = cobble.pred, ymin = reprod.3.pred.lcl, ymax = reprod.3.pred.ucl), alpha = 0.2)+
  geom_line(aes(x = cobble.pred, y = reprod.3.pred.med))+
  
  
  theme_bw()
reprod.plot

surv.plot = ggplot() + 
  
  geom_ribbon(aes(x = cobble.pred, ymin = surv.1.pred.lcl, ymax = surv.1.pred.ucl), alpha = 0.2)+
  geom_line(aes(x = cobble.pred, y = surv.1.pred.med))+
  
  geom_ribbon(aes(x = cobble.pred, ymin = surv.2.pred.lcl, ymax = surv.2.pred.ucl), alpha = 0.2)+
  geom_line(aes(x = cobble.pred, y = surv.2.pred.med))+
  
  theme_bw()
surv.plot
