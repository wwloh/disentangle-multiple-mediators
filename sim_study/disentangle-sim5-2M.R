rm(list=ls())

OneData <- function(
  n # sample size
  ) {
  
  # coefficients in mediator and outcome models
  a1 <- -0.09
  a2 <- 0
  e12 <- 1
  b1 <- 0
  b2 <- -0.66
  bA <- 0
  bC <- 1
  sM <- 1; sY <- 1
  aC <- rep(1,2)
  aU <- rep(1,2)
  
  U <- rnorm(n,mean=1)
  C <- rnorm(n,mean=1)
  A <- rbinom(n, size = 1, prob = .5)
  M1 <- a1*A + aC[1]*C + aU[1]*U + rnorm(n,sd=sM)
  M2 <- a2*A + e12*M1 + aC[2]*C + aU[2]*U + rnorm(n,sd=sM)
  Y <- bA*A + b1*M1 + b2*M2 + bC*C + rnorm(n,sd=sY)
  Data <- data.frame(id = 1:n, A, C, M1, M2, Y)
  rm(A, C, M1, M2, Y, U)
  
  true_effects <- NULL
  true_effects["ie1"] <- a1*b1
  true_effects["ie2"] <- (a2+e12*a1)*b2
  true_effects["de"] <- bA
  
  # marginal mediator models
  fitM <- lapply(1:2, function(s) {
    lm(as.formula(paste0("M",s,"~A+C")), data=Data)
  })
  
  res <- NULL
  # parallel path model
  fitY <- lm(Y~A+M1+M2+C, data=Data)
  for (s in 1:2) {
    res[[paste0("ie",s,".parallel")]] <- 
      coef(fitY)[paste0("M",s)]*coef(fitM[[s]])["A"]
  }
  res[["de.parallel"]] <- coef(fitY)["A"]
  
  return(unlist(c("true"=true_effects,res)))
}

# simulation settings
simsettings <- expand.grid(
  "n"=c(50,200,1000))

nsims <- 10000
simres <- list()
for (ss in 1:nrow(simsettings)) {
  ptm <- proc.time()[3]
  res <- replicate(nsims, OneData(n=simsettings[ss,"n"]),simplify="matrix")
  simres[[ss]] <- cbind("n"=simsettings[ss,],t(res),row.names=NULL)
  rm(res)
  
  filename <- paste(c(rbind(names(simsettings),
                            unlist(simsettings[ss,]))),collapse="_")
  cat(filename,"; time taken (mins) =", round((proc.time()[3]-ptm)/60),"\n")
}

save(simres,file="disentangle-sim5-2M.Rdata")

# results =====================================================================
library("xtable")
load(file="disentangle-sim5-2M.Rdata")
res <- do.call(rbind,simres)

effects <- c(paste0("ie",c(1:2)),"de")
res_ie.parallel <- NULL
for (eff in effects) {
  res_eff <- data.frame(res[,c("n",grep(pattern=eff,x=colnames(res),value=TRUE))])
  res_eff <- by(res_eff,res_eff$n, function(x) c("est"=colMeans(x),
                                                 "ese"=apply(x,2,sd)))
  res_eff <- do.call(rbind,res_eff)
  # true value
  true_eff <- res_eff[,grep("true",colnames(res_eff),value=TRUE)][1]
  # flatten the results to a single row
  res_ie.parallel[[eff]] <- cbind(true_eff,data.frame(matrix(
    t(res_eff[,grep(pattern="[.]parallel",x=colnames(res_eff),value=TRUE)]),
  nrow=1)))
}
xtable(do.call(rbind,res_ie.parallel))
