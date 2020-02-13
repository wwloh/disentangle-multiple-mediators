rm(list=ls())
libraries_check <- c("lavaan","data.table")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

a01 <- 1; a1 <- 0.5; aU1 <- 1.5*a1
a02 <- 1; a2 <- 1; aU2 <- 1.5*a2; e12 <- -3
s1 <- 1; s2 <- 1; sY <- 1
b0 <- 1; bA <- 1; b1 <- 1; b2 <- 1

OneData <- function(n=400) {
  A <- rbinom(n, size = 1, prob = 0.5)
  U <- rnorm(n, mean = 1)
  M1 <- rnorm(n, mean = a01 + a1*A + aU1*U, sd = s1)
  M2 <- rnorm(n, mean = a02 + a2*A + aU2*U + e12*M1, sd = s2)
  Y <- rnorm(n, mean = b0 + bA*A + b1*M1 + b2*M2, sd = sY)
  return(data.frame(id = 1:n, A, M1, M2, Y))
}

true_effects <- list()
true_effects[["ie1"]] <- b1*a1
true_effects[["ie2"]] <- b2*a2+(b2*e12*a1)
true_effects[["pe1"]] <- b1*a1
true_effects[["pe2"]] <- b2*a2
true_effects[["pe12"]] <- b2*e12*a1

OnePAestimator <- function() {
  data=OneData()
  models <- list()
  # (1) correctly specified path model with U unobserved 
  models[["correct"]] <- '
    M1 ~ a1*A
    M2 ~ a2*A + e12*M1
    Y ~ bA*A + b1*M1 + b2*M2
    pe1 := b1*a1
    pe2 := b2*a2
    pe12 := b2*e12*a1
    ie1 := pe1
    ie2 := pe2+pe12
  '
  # (2) parallel path model assuming M1,M2 do not affect each other
  models[["parallel"]] <- '
    M1 ~ a1*A
    M2 ~ a2*A
    Y ~ bA*A + b1*M1 + b2*M2
    M1 ~~ M2
    ie1 := b1*a1
    ie2 := b2*a2
  '
  # (3) parallel path model but M1,M2 not allowed to covary
  models[["uncorM"]] <- '
    M1 ~ a1*A
    M2 ~ a2*A
    Y ~ bA*A + b1*M1 + b2*M2
    M1 ~~ 0*M2
    ie1 := b1*a1
    ie2 := b2*a2
  '
  
  res <- lapply(1:length(models), function(mm) {
    stm <- tryCatch(system.time(
      fit <- sem(model=models[[mm]], data=data, estimator="ML", se="none")
    )[3], error=function(cond) return(NA))
    # check convergence criteria
    if (is.na(stm)) {
      # error returned from model fit
      converge.crit <- FALSE
    } else {
      converge.crit <- 
        # convergence OK
        lavInspect(fit,what="converged") & 
        # no NA estimates
        all(unlist(lapply(lavInspect(fit,what="est"),
                          function(x) all(!is.na(x)) )))
    }
    if (converge.crit) {
      return(fit)
    } else {
      return(NULL)
    }
  })
  names(res) <- names(models)
  return(res)
}

nsims <- 10000
ptm <- proc.time()[3]
simres <- lapply(1:nsims, function(x) OnePAestimator())
proc.time()[3]-ptm
# 38614.81 for 10000 sims

save(simres,file="disentangle-sim1-2M-gof.Rdata")
q()

# results ---------------------------------------------------------------------
# keep only datasets where all fitted models converged
load("disentangle-sim1-2M-gof.Rdata")
null_res <- do.call(rbind,lapply(simres, function(x) 
  unlist(lapply(x, is.null))))
todrop <- rowSums(null_res)
sim_results_list <- simres[todrop==0]
# diagnostic: datasets that did not converge
sim_results_NAs <- null_res[todrop>0,]
## number of NA's
colSums(sim_results_NAs)
colSums(sim_results_NAs)/nrow(null_res)
rm(simres)

# parameter estimates
res <- rbindlist(lapply(1:length(sim_results_list), function(ss) {
  sim <- sim_results_list[[ss]]
  rbindlist(lapply(1:length(sim), function(mm) {
    fit <- sim[[mm]]
    fit.est <- data.table(parameterEstimates(fit))
    setkey(fit.est)
    fit.est <- fit.est[grepl("pe",label) | grepl("ie",label), list(label,est)]
    fit.est <- cbind("sim"=ss,"fit"=names(sim)[mm],fit.est)
    return(fit.est)
  }))
}))
setkey(res)

# number of sims & MC error
(nsims <- unique(res[,.N,by=list(fit,label)][,N]))
(mc_err <- sqrt((1/nsims)*(1-1/nsims)))

# empirical mean and std err
res.mean <- res[,lapply(.SD, mean),.SDcols="est",by=list(fit,label)]
res.ese <- res[,lapply(.SD, sd),.SDcols="est",by=list(fit,label)]
setnames(res.ese,"est","ese")
setkey(res.mean)
setkey(res.ese)
res.summ <- merge(res.mean,res.ese)
setkey(res.summ)
for (eff in 1:length(true_effects)) {
  res.summ[label==names(true_effects)[eff], "true" := true_effects[[eff]]]
}

(meths <- res.summ[,unique(fit)])
res.table <- res.summ[fit==meths[1], list(
  label,true,paste0(round(est,2), " (", round(ese,2), ")"))]
setnames(res.table, "V3", meths[1])
setkey(res.table)
for (ff in 2:length(meths)) {
  ff.out <- res.summ[fit==meths[ff], list(
    label,true,paste0(round(est,2), " (", round(ese,2), ")"))]
  setnames(ff.out, "V3", meths[ff])
  setkey(ff.out)
  res.table <- merge(res.table,ff.out,all=TRUE)
  setkey(res.table)
}
library("xtable")
xtable(res.table)

## point estimates of interventional effects in correct vs. parallel models
ie1.correct_parallel <- res[,.SD[fit=="correct" & label=="ie1", est]-
                              .SD[fit=="parallel" & label=="ie1", est],by=sim][,V1]
ie2.correct_parallel <- res[,.SD[fit=="correct" & label=="ie2", est]-
                              .SD[fit=="parallel" & label=="ie2", est],by=sim][,V1]
max(abs(ie1.correct_parallel))
max(abs(ie2.correct_parallel))
pdf("plot-sim1-2M-gof-correct_parallel.pdf",width = 6,height = 3)
par(mfrow=c(1,2))
hist(ie1.correct_parallel,main="IE1",probability = TRUE,
     xlab="Correct vs. parallel")
abline(v=0,lwd=2)
hist(ie2.correct_parallel,main="IE2",probability = TRUE,
     xlab="Correct vs. parallel")
abline(v=0,lwd=2)
dev.off()

## point estimates of interventional effects in parallel models 
## M1,M2 allowed to covary vs. not allowed
ie1.uncorM_parallel <- res[,.SD[fit=="uncorM" & label=="ie1", est]-
                             .SD[fit=="parallel" & label=="ie1", est],by=sim][,V1]
ie2.uncorM_parallel <- res[,.SD[fit=="uncorM" & label=="ie2", est]-
                             .SD[fit=="parallel" & label=="ie2", est],by=sim][,V1]
max(abs(ie1.uncorM_parallel))
max(abs(ie2.uncorM_parallel))
pdf("plot-sim1-2M-gof-uncorM_parallel.pdf",width = 6,height = 3)
par(mfrow=c(1,2))
hist(ie1.uncorM_parallel,main="IE1",probability = TRUE,
     xlab="Uncorrelated vs. correlated")
abline(v=0,lwd=2)
hist(ie2.uncorM_parallel,main="IE2",probability = TRUE,
     xlab="Uncorrelated vs. correlated")
abline(v=0,lwd=2)
dev.off()

# goodness-of-fit statistics ==================================================
res_fits <- rbindlist(lapply(1:length(sim_results_list), function(ss) {
  sim <- sim_results_list[[ss]]
  sim_fit <- do.call(rbind,lapply(sim, fitmeasures))
  data.table("sim"=ss,"fit"=rownames(sim_fit),sim_fit)
}))
setkey(res_fits)

# saturated models
res_fits[fit %in% meths[1:3], unique(df), by=sim][,unique(V1)]
res_fits[fit %in% meths[1:3], unique(rmsea), by=sim][,unique(V1)]
res_fits[fit %in% meths[1:3], unique(cfi), by=sim][,unique(V1)]
res_fits[fit %in% meths[1:3], unique(tli), by=sim][,unique(V1)]
# identical GOF for the saturated models
res_fits[fit %in% meths[1:3], diff(range(aic)), by=sim][,max(V1)]
res_fits[fit %in% meths[1:3], diff(range(bic)), by=sim][,max(V1)]
res_fits[fit %in% meths[1:3], diff(range(chisq)), by=sim][,max(V1)]
res_fits[fit %in% meths[1:3], diff(range(baseline.chisq)), by=sim][,max(V1)]
res_fits[fit %in% meths[1:3], diff(range(logl)), by=sim][,max(V1)]
