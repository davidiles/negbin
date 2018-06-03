library(jagsUI)
rm(list=ls())

###########################################
# Zero-inflated negative binomial
###########################################
N = 3000

#Covariates
X = rnorm(N)
Y = rnorm(N)
Z = rnorm(N)

count.B0 = 3.218876 #equal to log(25)
count.B1 = 1
count.B2 = -1
count.B3 = 0.2

zero.B0 = -0.8472979 #equal to logit of 0.3 (30% chance of zero count)
zero.B1 = 0.5
zero.B2 = 0
zero.B3 = 0
zero.sd.resid = 0

#Count term (negative binomial) - affected by X and Y
count_mu = exp(count.B0 + count.B1*X + count.B2*Y + count.B3*Z)
count_r = 10
count = rnbinom(N,mu = count_mu,size = count_r)

#Zero term - affected by Y, Z, and other unmeasured factors
zero_mu = plogis(zero.B0 + zero.B1*X + zero.B2*Y + zero.B3*Z + rnorm(N,0,zero.sd.resid))
zero = rbinom(N,1,zero_mu)

#actual data
y = count * (1-zero)

#Now try to estimate parameters
library(pscl)
addi <- formula(y ~ X + Y + Z| X + Y + Z)
addzi<- zeroinfl(addi, 
                 #offset=logoffset,
                 dist = "negbin",link = "logit")
summary(addzi)

library(jagsUI)
# Specify model in JAGS language
sink("zinb.txt")
cat("
  model{
    
    #Priors
    r ~ dunif(0,50)
    
    count.intercept ~ dunif(0.001,200)
    count.B0 <- log(count.intercept)

    count.B1 ~ dnorm(0,0.01)
    count.B2 ~ dnorm(0,0.01)
    count.B3 ~ dnorm(0,0.01)
    
    zero.intercept ~ dunif(0.001,0.999)
    zero.B0 <- log(zero.intercept/(1-zero.intercept))

    zero.B1 ~ dnorm(0,0.01)
    zero.B2 ~ dnorm(0,0.01)
    zero.B3 ~ dnorm(0,0.01)
    
    # Likelihood
    for(i in 1:N){
      y[i] ~ dnegbin(p[i],r)
      p[i] <- r/(r+(1-zero[i])*mu.count[i]) - 1e-10*zero[i]
      mu.count[i] <- exp(log.mu.count[i])
      log.mu.count[i] <- count.B0 + count.B1*X[i] + count.B2*Y[i] + count.B3*Z[i]
      
      ## Zero-Inflation
      zero[i] ~ dbern(prob[i])
      prob[i] <- 1/(1+exp(-logit.prob[i]))
      logit.prob[i] <- zero.B0 + zero.B1*X[i] + zero.B2*Y[i] + zero.B3*Z[i]
    
    }
  }
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = y,
                  N = length(y),
                  X = X,
                  Y = Y,
                  Z = Z)

# Initial values
inits <- function() list(zero = rep(0,N))

# Parameters monitored
params <- c("r","count.B0","count.B1","count.B2","count.B3",
            "zero.B0","zero.B1","zero.B2","zero.B3")

# MCMC settings
ni <- 5000
nt <- 2
nb <- 4000
nc <- 3

# Call JAGS from R
out <- jags(data = jags.data, 
            inits, 
            params, 
            model.file = "zinb.txt", 
            n.chains = nc, 
            n.thin = nt, 
            n.iter = ni, 
            n.burnin = nb)

out

