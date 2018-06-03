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

# Create dummy variables for fdepth
yc$fdepth60 = yc$fdepth90 = yc$fdepth120 = 0
yc$fdepth60[yc$fdepth == "60"] = 1
yc$fdepth90[yc$fdepth == "90"] = 1
yc$fdepth120[yc$fdepth == "120"] = 1

# Create empty array to store age-specific abundances each year
obs_year_age = table(yc[,c("year_cap","assnage")]) #number of observations per age class per year
years = min(yc$year_cap):max(yc$year_cap)

### Array to store counts each year
# First dimension of array is years of study
# Second dimension of array is age classes
# Third dimension of array is repeated counts for that age in that year (currently, up to 17)
N_array = array(NA, dim = c(length(years), length(unique(yc$assnage)), max(obs_year_age)))

# Arrays to store covariates
fdepth_array = N_array

for (y in unique(yc$year_cap)){
  for (age in unique(yc$assnage)){
    
    year_age_dat = subset(yc, year_cap == y & assnage == age)
    
    for (i in 1:nrow(year_age_dat)){
      array_dim1 = which(years == y)
      array_dim2 = which(sort(unique(yc$assnage)) == age)
      array_dim3 = i
      
      N_array[array_dim1, array_dim2, array_dim3] = year_age_dat$countinyear[i]
      
      
    }
  }
}

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
    
    #------------------------------------------------
    #Zero-inflated Negative Binomial
    #------------------------------------------------
    
    #Priors
    # Each age has a different rate parameter (r), a different expected count, and a different zero-term
    for (age in 1:A){

      r[age] ~ dunif(0,50)
      count.intercept[age] ~ dunif(0.001,200)
      count.B0[age] <- log(count.intercept[age])

      zero.intercept.age[age] ~ dunif(0.001,0.999)
      zero.B0[age] <- log(zero.intercept.age[age]/(1-zero.intercept.age[age]))
    
    }

    # Likelihood
    for (year in 1:Y){
      for (age in 1:A){
        for (obs in 1:O){

          N[year,age,obs] ~ dnegbin(p[year,age,obs],r[age])

          #Count Model
              
              mu.count[year,age,obs] <- exp(log.mu.count[year,age,obs])
              log.mu.count[year,age,obs] <- count.B0[age]

              # Zero-Inflation
              zero[year,age,obs] ~ dbern(prob[year,age,obs])
              prob[year,age,obs] <- 1/(1+exp(-logit.prob[year,age,obs]))
              logit.prob[year,age,obs] <-  zero.B0[age]

              #Calculate the probability parameter (one of the two key components of the ZINB) based on mu, r, and zero
              p[year,age,obs] <- r[age]/(r[age]+(1-zero[year,age,obs])*mu.count[year,age,obs]) - 1e-10*zero[year,age,obs]

        }
      }
    
      
      
      }
    }
    ",fill = TRUE)
sink()

# Bundle data
jags.data <- list(N = N_array,
                  Y = dim(N_array)[1],
                  A = dim(N_array)[2],
                  O = dim(N_array)[3])

# Initial values
inits <- function() list(zero = array(0,dim = dim(N_array)))

# Parameters monitored
params <- c("r","count.B0","zero.B0")

# MCMC settings
ni <- 15000
nt <- 2
nb <- 10000
nc <- 3

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
