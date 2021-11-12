# load the data ---------------------------------------------------------------
rm(list=ls())
obs_Data <- read.csv("../data-poli_inclu/poli-inclu-study3-InclusionVsControl.csv")

source("wrapper-interventional.R")
MultipleMediators(outcome_name = "Y",
                  mediator_names = paste0("M",4:6), # subset as an example
                  treatment_name = "A",
                  covariate_names = c("Ideology","Gender","Age"),
                  obs_Data=obs_Data,
                  nboots = 500) # increase number of bootstrap draws in practice
