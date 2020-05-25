rm(list=ls())

OneData <- function(
  n, # sample size
  UA, # unobserved confounder of mediators affected by treatment
  Mint # M-M interactions in outcome model
  ) {
  
  # coefficients in mediator and outcome models
  a1 <- 12.8
  a2 <- 1.5
  a3 <- 1.92
  a04 <- 1
  b2 <- 0.08; b3 <- 0.08; b4 <- 0.19; bC <- 1;
  sM <- 1; sY <- 1
  aC <- rep(1,4)
  aU <- rep(1,4)
  
  # unobserved confounder affected by treatment
  aUA <- rep(UA,4)
  
  # mediator-mediator interactions in outcome model
  if (Mint==TRUE) {
    b23 <- -(b2+b3)
  } else {
    b23 <- 0
  }

  U <- rnorm(n)
  C <- rnorm(n)
  A <- rbinom(n, size = 1, prob = .5)
  M1 <- a1*A + aC[1]*C + aU[1]*U + aUA[1]*(U*A) + rnorm(n,sd=sM)
  M2 <- a2*A + aC[2]*C + aU[2]*U + aUA[2]*(U*A) + rnorm(n,sd=sM)
  M3 <- a3*A + aC[3]*C + aU[3]*U + aUA[3]*(U*A) + rnorm(n,sd=sM)
  M4 <- a04 + aC[4]*C + aU[4]*U + aUA[4]*(U*A) + rnorm(n,sd=sM)
  Y <- rnorm(n, mean = b2*M2 + b3*M3 + b4*M4 + bC*C + b23*M2*M3, sd = sY)
  Data <- data.frame(id = 1:n, A, C, M1, M2, M3, M4, Y)
  rm(A, C, M1, M2, M3, M4, Y, U)
  
  true_effects <- NULL
  true_effects["ie1"] <- 0
  true_effects["ie2"] <- a2*b2
  true_effects["ie3"] <- a3*(b3+b23*a2)
  true_effects["ie4"] <- 0
  true_effects["iemu"] <- b23*(aU[2]*aUA[3]+aU[3]*aUA[2]+aUA[2]*aUA[3])
  true_effects["de"] <- 0
  true_effects["pe12_"] <- 0
  true_effects["pe13_"] <- 0
  true_effects["pe14_"] <- 0
  true_effects["pe124"] <- 0
  true_effects["pe134"] <- 0
  
  # marginal mediator models
  fitM <- lapply(1:4, function(s) {
    lm(as.formula(paste0("M",s,"~A+C")), data=Data)
  })
  # mean mediators under treatment and control
  Mmean_A <- lapply(fitM, function(fitMs) {
    c("A0"=as.numeric(coef(fitMs)["(Intercept)"]+coef(fitMs)["C"]*mean(Data$C)),
      "A1"=as.numeric(coef(fitMs)["(Intercept)"]+coef(fitMs)["C"]*mean(Data$C)+
                        coef(fitMs)["A"]))
  })
  
  res <- NULL
  # (1) outcome model with all mediator-mediator interactions
  fitY <- lm(Y~A+M1+M2+M3+M4+C+M1:M2+M1:M3+M1:M4+M2:M3+M2:M4+M3:M4, data=Data)
  
  res[["ie1.int"]] <- 
    (coef(fitY)["M1"]+
       coef(fitY)["M1:M2"]*Mmean_A[[2]]["A0"]+
       coef(fitY)["M1:M3"]*Mmean_A[[3]]["A0"]+
       coef(fitY)["M1:M4"]*Mmean_A[[4]]["A0"])*coef(fitM[[1]])["A"]
  
  res[["ie2.int"]] <- 
    (coef(fitY)["M2"]+
       coef(fitY)["M1:M2"]*Mmean_A[[1]]["A1"]+
       coef(fitY)["M2:M3"]*Mmean_A[[3]]["A0"]+
       coef(fitY)["M2:M4"]*Mmean_A[[4]]["A0"])*coef(fitM[[2]])["A"]
  
  res[["ie3.int"]] <- 
    (coef(fitY)["M3"]+
       coef(fitY)["M1:M3"]*Mmean_A[[1]]["A1"]+
       coef(fitY)["M2:M3"]*Mmean_A[[2]]["A1"]+
       coef(fitY)["M3:M4"]*Mmean_A[[4]]["A0"])*coef(fitM[[3]])["A"]
  
  res[["ie4.int"]] <- 
    (coef(fitY)["M4"]+
       coef(fitY)["M1:M4"]*Mmean_A[[1]]["A1"]+
       coef(fitY)["M2:M4"]*Mmean_A[[3]]["A1"]+
       coef(fitY)["M3:M4"]*Mmean_A[[3]]["A1"])*coef(fitM[[4]])["A"]
  
  SigmaM_A <- cov(Data[Data$A==1,paste0("M",1:4)])-
    cov(Data[Data$A==0,paste0("M",1:4)])
  
  ie.mu <- 0
  for (k in 1:4) {
    for (l in k:4) {
      if (k<l) {
        ie.mu <- ie.mu + coef(fitY)[paste0("M",k,":","M",l)]*SigmaM_A[k,l]
      }
    }
  }
  res[["iemu.int"]] <- ie.mu
  res[["de.int"]] <- coef(fitY)["A"]
  
  # (2) parallel mediator model
  fitY <- lm(Y~A+M1+M2+M3+M4+C, data=Data)
  for (s in 1:4) {
    res[[paste0("ie",s,".parallel")]] <- 
      coef(fitY)[paste0("M",s)]*coef(fitM[[s]])["A"]
  }
  res[["de.parallel"]] <- coef(fitY)["A"]
  
  if (UA==0 && Mint==FALSE) {
    # setting for study 1
    # (3) posited model
    fitM1 <- lm(M1 ~ A+C, data=Data)
    fitM2 <- lm(M2 ~ M1+C, data=Data)
    fitM3 <- lm(M3 ~ M1+C, data=Data)
    fitM4 <- lm(M4 ~ M1+M2+M3+C, data=Data)
    fitY <- lm(Y~M2+M3+M4+C, data=Data)
    
    res["pe12_"] <- coef(fitY)["M2"]*coef(fitM2)["M1"]*coef(fitM1)["A"]
    res["pe13_"] <- coef(fitY)["M3"]*coef(fitM3)["M1"]*coef(fitM1)["A"]
    res["pe14_"] <- coef(fitY)["M4"]*coef(fitM4)["M1"]*coef(fitM1)["A"]
    res["pe124"] <- coef(fitY)["M4"]*coef(fitM4)["M2"]*
      coef(fitM2)["M1"]*coef(fitM1)["A"]
    res["pe134"] <- coef(fitY)["M4"]*coef(fitM4)["M3"]*
      coef(fitM3)["M1"]*coef(fitM1)["A"]  
  } else {
    res["pe12_"] <- res["pe13_"] <- res["pe14_"] <- 
      res["pe124"] <- res["pe134"] <- NA
  }
  
  return(c("true"=true_effects,res))
}

# simulation settings
simsettings <- expand.grid(
  "UA"=0:1,
  "Mint"=c(FALSE,TRUE),
  "n"=c(50,200,1000))
simsettings <- simsettings[(simsettings$UA<=simsettings$Mint),]

nsims <- 10000
simres <- list()
for (ss in 1:nrow(simsettings)) {
  ptm <- proc.time()[3]
  res <- replicate(nsims, OneData(n=simsettings[ss,"n"], 
                                  UA=simsettings[ss,"UA"], 
                                  Mint=simsettings[ss,"Mint"]),
                   simplify="matrix")
  simres[[ss]] <- cbind(simsettings[ss,],t(res),row.names=NULL)
  rm(res)
  
  filename <- paste(c(rbind(names(simsettings),
                            unlist(simsettings[ss,]))),collapse="_")
  cat(filename,"; time taken (mins) =", round((proc.time()[3]-ptm)/60),"\n")
}

save(simres,file="disentangle-sim4-4M-mutual.Rdata")

# results =====================================================================
library("xtable")
load(file="disentangle-sim4-4M-mutual.Rdata")
res <- do.call(rbind,simres)

## study 1 
effects <- paste0("pe",c("12_","13_","14_","124","134"))
res_pe <- NULL
for (eff in effects) {
  res_eff <- res[res$UA==0 & res$Mint==FALSE,
                 c("n",grep(pattern=eff,x=names(res),value=TRUE))]
  res_eff <- do.call(rbind,by(res_eff,res_eff$n,FUN=function(x) 
    c("est"=colMeans(x),"ese"=apply(x,2,sd))))
  # true value
  true_eff <- unique(res_eff[,paste0("est.true.",eff)])
  # flatten the results to a single row
  res_eff <- data.frame(matrix(
    t(res_eff[,paste0(c("est.","ese."),eff)]),nrow=1))
  res_eff <- cbind(eff,true_eff,res_eff)
  res_pe <- rbind(res_pe,res_eff)
  rm(res_eff)
}
xtable(res_pe)

effects <- c(paste0("ie",1:4),"de")
res_ie <- NULL
for (eff in effects) {
  res_eff <- res[res$UA==0 & res$Mint==FALSE,
                 c("n",
                   grep(pattern=paste0("true.",eff),x=names(res),value=TRUE),
                   grep(pattern=paste0(eff,".parallel"),x=names(res),value=TRUE)
                   )]
  res_eff <- do.call(rbind,by(res_eff,res_eff$n,FUN=function(x) 
    c("est"=colMeans(x),"ese"=apply(x,2,sd))))
  # true value
  true_eff <- unique(res_eff[,paste0("est.true.",eff)])
  # flatten the results to a single row
  res_eff <- data.frame(matrix(
    t(res_eff[,grep(pattern="parallel",x=colnames(res_eff),value=TRUE)]),nrow=1))
  res_eff <- cbind(eff,true_eff,res_eff)
  res_ie <- rbind(res_ie,res_eff)
  rm(res_eff)
}
xtable(res_ie)

## study 2
effects <- c(paste0("ie",c(1:4,"mu")),"de")
res_ie <- NULL
res_ie.parallel <- NULL
for (eff in effects) {
  res_eff <- res[res$Mint==TRUE,
                 c("UA","Mint","n",grep(pattern=eff,x=names(res),value=TRUE))]
  res_eff <- do.call(rbind,by(res_eff,list(res_eff$n,res_eff$UA,res_eff$Mint),
                              FUN=function(x) 
                                c("est"=colMeans(x),"ese"=apply(x,2,sd))))
  # different settings 
  res_settings <- unique(res_eff[,c("est.UA","est.Mint")])
  res_ie_eff <- NULL
  res_ie_eff.parallel <- NULL
  for (i in 1:nrow(res_settings)) {
    res_eff_i <- res_eff[res_eff[,"est.UA"]==res_settings[i,"est.UA"] &
                           res_eff[,"est.Mint"]==res_settings[i,"est.Mint"],]
    # true value
    true_eff <- unique(res_eff_i[,paste0("est.true.",eff)])
    # flatten the results to a single row
    res_ie_eff[[i]] <- data.frame(matrix(
      t(res_eff_i[,grep(pattern="[.]int",x=colnames(res_eff),value=TRUE)]),
      nrow=1))
    res_ie_eff.parallel[[i]] <- data.frame(matrix(
      t(res_eff_i[,grep(pattern="[.]parallel",x=colnames(res_eff),value=TRUE)]),
      nrow=1))
    res_ie_eff[[i]] <- cbind(eff,true_eff,res_ie_eff[[i]])
    res_ie_eff.parallel[[i]] <- cbind(eff,true_eff,res_ie_eff.parallel[[i]])
  }
  res_ie_eff <- do.call(rbind,res_ie_eff)
  res_ie_eff.parallel <- do.call(rbind,res_ie_eff.parallel)
  
  res_ie <- rbind(res_ie,res_ie_eff)
  if (eff != "iemu") {
    res_ie.parallel <- rbind(res_ie.parallel,res_ie_eff.parallel)  
  }
}
library("xtable")
xtable(res_ie)
xtable(res_ie.parallel)
