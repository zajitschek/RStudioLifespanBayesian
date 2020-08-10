### Lifespan-analysis for Biologists ###

## Chapter 4.3: Survival analysis as a Bayesian ##

## The following rastanarm code is largely based on the preprint:
## "Bayesian Survival Analysis Using the rstanarm R Package" by Brilleman et al. 2020
## https://arxiv.org/pdf/2002.09633.pdf

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

# Try out other available baseline hazard functions

stan_surv_model1_exp <- update(stan_surv_model1, basehaz = "exp")
stan_surv_model1_weibull <- update(stan_surv_model1, basehaz = "weibull")
stan_surv_model1_gompertz <- update(stan_surv_model1, basehaz = "gompertz")
stan_surv_model1_bspline <- update(stan_surv_model1, basehaz = "bs")  #B-splines
stan_surv_model1_mspline1 <- update(stan_surv_model1, basehaz = "ms")  #M-splines
stan_surv_model1_mspline2 <- update(stan_surv_model1, basehaz = "ms", basehaz_ops = list(df = 9))

# To plot the baseline, we use the plot function, and adjust it slightly to make it look nicer

plotfun <- function(model, title) {
 plot(model, plotfun = "basehaz") +
 coord_cartesian(ylim = c(0,0.4)) +
 labs(title = title) +
 theme(plot.title = element_text(hjust = 0.5))
 }

# Now let's create and save the plots

p_exp <- plotfun(stan_surv_model1_exp, "Exponential")
p_weibull <- plotfun(stan_surv_model1_weibull, "Weibull")
p_gompertz <- plotfun(stan_surv_model1_gompertz, "Gompertz")
p_bspline <- plotfun(stan_surv_model1_bspline, "B-splines with\ntwo internal knots"))
p_mspline1 <- plotfun(stan_surv_model1_mspline1, "M-splines with\ntwo internal knots"))
p_mspline2 <- plotfun(stan_surv_model1_mspline2, "M-splines with\nfive internal knots")

p_combined <- plot_grid(p_exp, p_weibull, p_gompertz, 
 p_bspline, p_mspline1, p_mspline2, ncol = 3)
p_combined

# Compare model fits

compare_models(loo(stan_surv_model1_exp),
 loo(stan_surv_model1_weibull), loo(stan_surv_model1_gompertz),
 loo(stan_surv_model1_bspline), loo(stan_surv_model1_mspline1),
 loo(stan_surv_model1_mspline2), nrow= 3)

# We can also predict survival, using the posterior distributions

nd <- data.frame(group = c("Low", "Standard", "High"))

ps <- posterior_survfit(stan_surv_model1,
 newdata = nd,
 times = 0,
 extrapolate = TRUE,
 control = list(edist = 70))  #set end of prediction (in time units)
head(ps)

panel_labels <- c('1' = "Low",
 '2' = "Standard",
 '3' = "High")
pps <- plot(ps) +
 facet_wrap(~ id, labeller = labeller(id = panel_labels))

# This also works for the predicted hazard or log hazard function 

ph <- posterior_survfit(stan_surv_model1, newdata = nd, type = "haz")
pl <- posterior_survfit(stan_surv_model1, newdata = nd, type = "loghaz")

# Plots

pph <- plot(ph) +
 facet_wrap(~ id, labeller = labeller(id = panel_labels))
ppl <- plot(pl) +
 facet_wrap(~ id, labeller = labeller(id = panel_labels))

pph

ppl



# Try out the following on your local computer only
#library(shinystan)
#launch_shinystan(stan_surv_model1)  