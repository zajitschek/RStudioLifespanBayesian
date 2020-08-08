### Lifespan-analysis for Biologists ###

## Chapter 3.3: Bayesian with rstanarm::stan_lmer ##

#Load packages and data
suppressMessages(library(dplyr))
suppressMessages(library(rstanarm))
data1 <- read.csv("https://github.com/ZajitschekTeam/lifespananalysis/raw/master/binder/data/expevol_male_flies.csv")
data1 <- data1 %>% mutate(across(where(is.integer), as.factor))
          
# Your first Bayesian model
#  For now, we will assume that Gaussian errors are ok (this is the default)

stan_glmm_model1 <- stan_lmer(lifespan ~ assaydiet + (1|vial), data= data1, seed= 111)
plot(stan_glmm_model1)
summary(stan_glmm_model1)

summary(stan_glmm_model1, 
        pars = c("(Intercept)", "assaydiet3", "assaydiet4"),
        probs = c(0.025, 0.975),
        digits = 2)

# Try out the following on your local computer only
#library(shinystan)
#launch_shinystan(stan_glmm_model1)  