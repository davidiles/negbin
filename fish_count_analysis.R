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

plot(yc$avgtrawlcobsqkm ~ as.integer(as.character(yc$fyear)))

addi <- formula(yc$countinyear ~ assnage + spawncob + fdepth| assnage + spawncob + fdepth )
addzi<- zeroinfl(addi, offset=logoffset, dist = "negbin",link = "logit", data = yc)

summary(addzi)

# Create dummy variables for fdepth
yc$fdepth60 = yc$fdepth90 = yc$fdepth120 = 0
yc$fdepth60[yc$fdepth == "60"] = 1
yc$fdepth90[yc$fdepth == "90"] = 1
yc$fdepth120[yc$fdepth == "120"] = 1

# Specify model in JAGS language
sink("fish_count_model.txt")
cat("
    model{
    
    #------------------------------------------------
    #Zero-inflated Negative Binomial
    #------------------------------------------------
    
    #Priors
    r ~ dunif(0,50)
    
    count.intercept ~ dunif(0.001,200)
    count.B0 <- log(count.intercept)
    
    count.B0.offset <- 1 #Note that an offset is a predictor variable with coefficient set to 1
    
    count.B1 ~ dnorm(0,0.01) #assnage
    count.B2 ~ dnorm(0,0.01) #spawncob
    
    count.B3 ~ dnorm(0,0.01) #fdepth60
    count.B4 ~ dnorm(0,0.01) #fdepth90
    count.B5 ~ dnorm(0,0.01) #fdepth120
    
    zero.intercept ~ dunif(0.001,0.999)
    zero.B0 <- log(zero.intercept/(1-zero.intercept))
    
    zero.B1 ~ dnorm(0,0.01) #assnage
    zero.B2 ~ dnorm(0,0.01) #spawncob
    
    zero.B3 ~ dnorm(0,0.01) #fdepth60
    zero.B4 ~ dnorm(0,0.01) #fdepth90
    zero.B5 ~ dnorm(0,0.01) #fdepth120
    
    # Likelihood
    for(i in 1:N){
    
    y[i] ~ dnegbin(p[i],r)
    
    #Count Model
    p[i] <- r/(r+(1-zero[i])*mu.count[i]) - 1e-10*zero[i]
    mu.count[i] <- exp(log.mu.count[i])
    
    log.mu.count[i] <- count.B0 + count.B0.offset*logoffset[i] + count.B1*assnage[i] + count.B2*spawncob[i] + count.B3*fdepth60[i] + count.B4*fdepth90[i] + count.B5*fdepth120[i]
    
    # Zero-Inflation
    zero[i] ~ dbern(prob[i])
    prob[i] <- 1/(1+exp(-logit.prob[i]))
    logit.prob[i] <-  zero.B0 + zero.B1*assnage[i] + zero.B2*spawncob[i] + zero.B3*fdepth60[i] + zero.B4*fdepth90[i] + zero.B5*fdepth120[i]
    
    }
    }
    ",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = yc$countinyear,
                  N = length(yc$countinyear),
                  logoffset = yc$logoffset,
                  assnage = yc$assnage,
                  spawncob = yc$spawncob,
                  fdepth60 = yc$fdepth60,
                  fdepth90 = yc$fdepth90,
                  fdepth120 = yc$fdepth120
)

# Initial values
inits <- function() list(zero = rep(0,length(yc$countinyear)))

# Parameters monitored
params <- c("r","count.B0","count.B1","count.B2","count.B3","count.B4", "count.B5",
            "zero.B0","zero.B1","zero.B2","zero.B3","zero.B4","zero.B5")

# MCMC settings
ni <- 15000
nt <- 2
nb <- 10000
nc <- 3

# Call JAGS from R
out <- jags(data = jags.data, 
            inits, 
            params, 
            model.file = "fish_count_model.txt", 
            n.chains = nc, 
            n.thin = nt, 
            n.iter = ni, 
            n.burnin = nb)

out