### Lifespan-analysis for Biologists ###

## Chapter 4.3: Survival analysis as a Bayesian ##

#Load packages and data
suppressMessages(library(dplyr))
suppressMessages(library(rstanarm))
data1 <- read.csv("https://github.com/ZajitschekTeam/lifespananalysis/raw/master/binder/data/expevol_male_flies.csv")
data1 <- data1 %>% mutate(across(where(is.integer), as.factor))
          
# Your first Bayesian Cox PH model
#  For now, we will assume that Gaussian errors are ok (this is the default)

stan_surv_model1 <- stan_surv(lifespan ~ assaydiet + (1|vial), data= data1, seed= 111)
plot(stan_surv_model1)
summary(stan_surv_model1)

# Try out the following on your local computer only
#library(shinystan)
#launch_shinystan(stan_surv_model1)  