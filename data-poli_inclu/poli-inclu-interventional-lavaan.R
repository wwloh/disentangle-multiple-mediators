rm(list=ls())
# load the data ---------------------------------------------------------------
data <- read.csv("poli-inclu-study3-InclusionVsControl.csv")
# specify different models ----------------------------------------------------
models <- list()
# total effect
models[["total_effect"]] <- '
  Y ~ te*A + Ideology + Age + Gender
  Ideology ~~ 0*Age + 0*Gender
  Age ~~ 0*Gender
  A ~~ 0*Ideology + 0*Age + 0*Gender
  '
# parallel path model
models[["noMcovary"]] <- '
	Y ~ bA*A + b1*M1 + b2*M2 + b3*M3 + b4*M4 + b5*M5 + b6*M6 + Ideology + Age + Gender
	M1 ~ d1*A + Ideology + Age + Gender
	M2 ~ d2*A + Ideology + Age + Gender
	M3 ~ d3*A + Ideology + Age + Gender
	M4 ~ d4*A + Ideology + Age + Gender
	M5 ~ d5*A + Ideology + Age + Gender
	M6 ~ d6*A + Ideology + Age + Gender
	# fix covariates and treatment to be uncorrelated
	Ideology ~~ 0*Age + 0*Gender
	Age ~~ 0*Gender
	A ~~ 0*Ideology + 0*Age + 0*Gender
	# define interventional indirect effects via each mediator
	ie1 := b1*d1
	ie2 := b2*d2
	ie3 := b3*d3
	ie4 := b4*d4
	ie5 := b5*d5
	ie6 := b6*d6
  # define direct and joint indirect effects, and their sum
	de := bA
	ie_jt := ie1 + ie2 + ie3 + ie4 + ie5 + ie6
	de_ie_sum := de + ie_jt
	'	
# parallel path model but permit mediators to covary
models[["Mcovary"]] <- paste(
  models[["noMcovary"]],
  "M1 ~~ M2 + M3 + M4 + M5 + M6",
  "M2 ~~ M3 + M4 + M5 + M6",
  "M3 ~~ M4 + M5 + M6",
  "M4 ~~ M5 + M6",
  "M5 ~~ M6",
  sep=" \n ")

# fit models in laavaan -----------------------------------------------------
library("lavaan")
nboots <- 10000
res <- lapply(models, function(mm) {
  fit <- sem(mm, data = data, estimator = "ML", 
             se = "bootstrap", bootstrap = nboots, fixed.x=FALSE)
  summary(fit)
  fit.est <- parameterEstimates(fit)
  return(fit.est)
})
save(res,file="poli-inclu-study3-InclusionVsControl-interventional-lavaan.Rdata")

# summary output --------------------------------------------------------------
## total_effect --------------------------------------
# Optimization method                           NLMINB
# Number of free parameters                          9
# 
# Number of observations                           183
# 
# Estimator                                         ML
# Model Fit Test Statistic                       8.658
# Degrees of freedom                                 6
# P-value (Chi-square)                           0.194
# 
# Parameter Estimates:
#   
#   Standard Errors                            Bootstrap
# Number of requested bootstrap draws            10000
# Number of successful bootstrap draws           10000

## noMcovary -----------------------------------------
# Optimization method                           NLMINB
# Number of free parameters                         45
# 
# Number of observations                           183
# 
# Estimator                                         ML
# Model Fit Test Statistic                     259.495
# Degrees of freedom                                21
# P-value (Chi-square)                           0.000
# 
# Parameter Estimates:
#   
#   Standard Errors                            Bootstrap
# Number of requested bootstrap draws            10000
# Number of successful bootstrap draws            7129

## Mcovary --------------------------------------------
# Optimization method                           NLMINB
# Number of free parameters                         60
# 
# Number of observations                           183
# 
# Estimator                                         ML
# Model Fit Test Statistic                       8.658
# Degrees of freedom                                 6
# P-value (Chi-square)                           0.194
# 
# Parameter Estimates:
#   
#   Standard Errors                            Bootstrap
# Number of requested bootstrap draws            10000
# Number of successful bootstrap draws            9997

# format results --------------------------------------------------------------
load(file="poli-inclu-study3-InclusionVsControl-interventional-lavaan.Rdata")
library("data.table")
library("knitr")
lapply(res, function(fit.est) {
  res <- data.table(fit.est)
  setkey(res)
  kable(res[grep(pattern="e",label), .(label, est, se, ci.lower, ci.upper)], 
        digits=4)
})
# $total_effect
# 
# 
# |label |     est|     se| ci.lower| ci.upper|
#   |:-----|-------:|------:|--------:|--------:|
#   |te    | -0.0756| 0.0301|  -0.1342|   -0.016|
#   
#   $noMcovary
# 
# 
# |label     |     est|     se| ci.lower| ci.upper|
#   |:---------|-------:|------:|--------:|--------:|
#   |de        |  0.0394| 0.0220|  -0.0039|   0.0817|
#   |de_ie_sum | -0.0756| 0.0304|  -0.1348|  -0.0173|
#   |ie1       | -0.0028| 0.0071|  -0.0157|   0.0128|
#   |ie2       |  0.0000| 0.0022|  -0.0045|   0.0054|
#   |ie3       | -0.0045| 0.0053|  -0.0168|   0.0039|
#   |ie4       |  0.0017| 0.0037|  -0.0054|   0.0099|
#   |ie5       | -0.0121| 0.0065|  -0.0269|  -0.0017|
#   |ie6       | -0.0973| 0.0238|  -0.1461|  -0.0535|
#   |ie_jt     | -0.1150| 0.0279|  -0.1701|  -0.0626|
#   
#   $Mcovary
# 
# 
# |label     |     est|     se| ci.lower| ci.upper|
#   |:---------|-------:|------:|--------:|--------:|
#   |de        |  0.0394| 0.0220|  -0.0033|   0.0826|
#   |de_ie_sum | -0.0756| 0.0303|  -0.1354|  -0.0162|
#   |ie1       | -0.0028| 0.0070|  -0.0157|   0.0125|
#   |ie2       |  0.0000| 0.0022|  -0.0045|   0.0050|
#   |ie3       | -0.0045| 0.0053|  -0.0172|   0.0036|
#   |ie4       |  0.0017| 0.0037|  -0.0056|   0.0099|
#   |ie5       | -0.0121| 0.0066|  -0.0274|  -0.0016|
#   |ie6       | -0.0973| 0.0237|  -0.1464|  -0.0541|
#   |ie_jt     | -0.1150| 0.0277|  -0.1713|  -0.0638|

# difference between total effect and sum of direct and joint indirect effects
abs(data.table(res$total_effect)[label=="te", est]-
      data.table(res$Mcovary)[label=="de_ie_sum", est])
# [1] 1.033235e-08