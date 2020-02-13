rm(list=ls())
libraries_check <- c("lavaan","data.table")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)
sessionInfo()

a1 <- 1; aC1 <- 1; aU1 <- 1
e12 <- 1; aC2 <- 1; aU2 <- 1
e13 <- .5; aC3 <- 1; aU3 <- 1
aC4 <- 1; aU4 <- 1
b2 <- 1; b3 <- 2; b4 <- 1; aCY <- 1;
s1 <- 1; s2 <- 1; s3 <- 1; s4 <- 1; sY <- 1

OneData <- function(n=400) {
  U <- rnorm(n)
  C <- rnorm(n)
  A <- rbinom(n, size = 1, prob = .5)
  M1 <- rnorm(n, mean = a1*A + aC1*C + aU1*U, sd = s1)
  M2 <- rnorm(n, mean = e12*M1 + aC2*C + aU2*U, sd = s2)
  M3 <- rnorm(n, mean = e13*M1 + aC3*C + aU3*U, sd = s3)
  M4 <- rnorm(n, mean = aC4*C + aU4*U, sd = s4)
  Y <- rnorm(n, mean = b2*M2 + b3*M3 + b4*M4 + aCY*C, sd = sY)
  return(data.frame(id = 1:n, A, C, M1, M2, M3, M4, Y))
}

true_effects <- list()
true_effects[["pe12"]] <- b2*e12*a1
true_effects[["pe13"]] <- b3*e13*a1
true_effects[["pe14"]] <- 0
true_effects[["pe124"]] <- 0
true_effects[["pe134"]] <- 0
true_effects[["ie1"]] <- 0
true_effects[["ie2"]] <- true_effects[["pe12"]]
true_effects[["ie3"]] <- true_effects[["pe13"]]
true_effects[["ie4"]] <- 0

OnePAestimator <- function() {
  data=OneData()
  models <- list()
  # (1) posited model
  models[["posited"]] <- '
    M1 ~ a1*A + C
    M2 ~ e12*M1 + C
    M3 ~ e13*M1 + C
    M4 ~ e14*M1 + e24*M2 + e34*M3 + C
    Y ~ b2*M2 + b3*M3 + b4*M4 + C
    pe12 := b2*e12*a1
    pe13 := b3*e13*a1
    pe14 := b4*e14*a1
    pe124 := b4*e24*e12*a1
    pe134 := b4*e34*e13*a1
    ie2 := pe12
    ie3 := pe13
    ie4 := pe14 + pe124 + pe134
  '
  # (2) parallel mediator model
  models[["parallel"]] <- '
    M1 ~ a1*A + C
    M2 ~ a2*A + C
    M3 ~ a3*A + C
    M4 ~ a4*A + C
    Y ~ bA*A + b1*M1 + b2*M2 + b3*M3 + b4*M4 + C
    ie1 := b1*a1
    ie2 := b2*a2
    ie3 := b3*a3
    ie4 := b4*a4
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
ptm=proc.time()[3]
simres <- lapply(1:nsims, function(x) OnePAestimator())
proc.time()[3]-ptm
# 21203.75

save(simres,file="disentangle-sim4-4M-gof.Rdata")
q()

# results ---------------------------------------------------------------------
# keep only datasets where all fitted models converged
load("disentangle-sim4-4M-gof.Rdata")
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

# goodness-of-fit statistics ==================================================
res_fits <- do.call(rbind,lapply(1:length(sim_results_list), function(ss) {
  sim <- sim_results_list[[ss]]
  sim_fit <- data.table(do.call(rbind,lapply(sim, fitmeasures)),
                        keep.rownames = TRUE)
  setkey(sim_fit)
  c("sim"=ss,
    ## compare models with same df
    "posit.vs.para"=unlist(sim_fit[rn=="posited",list(aic,bic,tli,cfi)])-
      unlist(sim_fit[rn=="parallel",list(aic,bic,tli,cfi)]),
    "posit.srmr"=sim_fit[rn=="posited",srmr],
    "para.srmr"=sim_fit[rn=="parallel",srmr],
    "posit.rmsea"=sim_fit[rn=="posited",rmsea],
    "para.rmsea"=sim_fit[rn=="parallel",rmsea]
  )
}))
res_fits <- data.table(res_fits)
setkey(res_fits)
## probability that posited model has better fit than parallel model
mean(res_fits$posit.vs.para.aic<=0)
mean(res_fits$posit.vs.para.bic<=0)
mean(res_fits$posit.vs.para.tli>=0)
mean(res_fits$posit.vs.para.cfi>=0)
mean(res_fits$posit.srmr<.08)
mean(res_fits$para.srmr<.08)

# specific path effect estimates ==============================================
rm(res)
res <- rbindlist(lapply(1:length(sim_results_list), function(ss) {
  sim <- sim_results_list[[ss]]
  rbindlist(lapply(1:length(sim), function(mm) {
    fit <- sim[[mm]]
    fit.est <- data.table(parameterEstimates(fit))
    setkey(fit.est)
    fit.est <- fit.est[(label!="") & !grepl("pe",label) & !grepl("ie",label),
                       list(label,est)]
    if (mm==1) {
      fit.est <- rbind(fit.est,data.frame(
        "label"=c("a2","a3","a4"),"est"=rep(0.0,3)))
      fit.est[label=="a2", 
              est := fit.est[label=="a1",est]*fit.est[label=="e12",est]]
      fit.est[label=="a3", 
              est := fit.est[label=="a1",est]*fit.est[label=="e13",est]]
      fit.est[label=="a4", est := fit.est[label=="a1",est]*
                (fit.est[label=="e14",est]+
                   fit.est[label=="e12",est]*fit.est[label=="e24",est]+
                   fit.est[label=="e13",est]*fit.est[label=="e34",est])]
    }
    
    fit.est <- cbind("sim"=ss,"fit"=names(sim)[mm],fit.est)
    return(fit.est)
  }))
}))
setkey(res)

# empirical mean and std err
res.mean <- res[,lapply(.SD, mean),.SDcols="est",by=list(fit,label)]
res.ese <- res[,lapply(.SD, sd),.SDcols="est",by=list(fit,label)]
setnames(res.ese,"est","ese")
setkey(res.mean)
setkey(res.ese)
res.summ <- merge(res.mean,res.ese)
setkey(res.summ)
res.summ[label=="a1", true := a1]
res.summ[label=="a2", true := a1*e12]
res.summ[label=="a3", true := a1*e13]
res.summ[label=="a4", true := 0]
res.summ[label=="b1", true := 0]
res.summ[label=="b2", true := b2]
res.summ[label=="b3", true := b3]
res.summ[label=="b4", true := b4]
res.summ[label=="bA", true := 0]
res.summ[label=="e12", true := e12]
res.summ[label=="e13", true := e13]
res.summ[label=="e14", true := 0]
res.summ[label=="e24", true := 0]
res.summ[label=="e34", true := 0]
res.summ
