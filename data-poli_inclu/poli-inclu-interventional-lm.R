rm(list=ls())
# specify different models ----------------------------------------------------
models <- list()
# total effect
models[["total_effect"]] <- as.formula(Y ~ A + Ideology + Age + Gender)
# mediator models
models_M <- sapply(1:6, function(s) 
  as.formula(paste0("M",s, " ~ A + Ideology + Age + Gender")))
names(models_M) <- paste0("M",1:6)
# outcome model with interactions
intM <- NULL
for (s in 1:6) {
  for (p in (s+1):6) {
    if (p > s & p <= 6) {
      intM <- c(intM,paste(paste0("M",s),paste0("M",p),sep=":"))  
    }
  }    
}
length(intM) == 6*(6-1)/2

models[["Y_Mint"]] <- as.formula(
  paste0("Y ~ A +", paste0("M",1:6,collapse="+"), 
         "+", paste0("A:M",1:6,collapse="+"),
         "+", paste(intM,collapse="+"),
         "+", paste(paste0("A:",intM),collapse="+"),
         "+ Ideology + Age + Gender"))

# consider all permutations for sensitivity analysis --------------------------
## appropriate only for small number of mediators
allperms <- expand.grid(lapply(1:6, function (x) 6:1))
allperms <- allperms[apply(allperms, 1, function(x) 
  length(unique(x))==length(x)),]
rownames(allperms) <- NULL
colnames(allperms) <- paste0("Mp",1:6)
head(allperms)

## alternative for larger number of mediators
if (FALSE) {
  # install.packages("gtools")
  library("gtools")
  allperms.gtools <- gtools::permutations(n=6,r=6)
  head(allperms.gtools)
  nrow(allperms.gtools)
}
  
# load the data ---------------------------------------------------------------
obs_Data <- read.csv("poli-inclu-study3-InclusionVsControl.csv")

# fit models ------------------------------------------------------------------
OneEst <- function(boot=FALSE) {
  
  if (boot) {
    # bootstrap sample
    Data <- obs_Data[sample(nrow(obs_Data),replace=TRUE),]  
  } else {
    Data <- obs_Data
  }
  
  res <- NULL
  # total effect
  fit_te <- lm(models$total_effect,data=Data)
  res["te"] <- coef(fit_te)["A"]
  rm(fit_te)
  
  # marginal effects for each mediator
  deltas <- NULL
  for (s in 1:6) {
    fit_M <- lm(models_M[[s]],data=Data)
    deltas[s] <- coef(fit_M)["A"]
    rm(fit_M)
  }
  names(deltas) <- names(models_M)
  
  # mediator-mediator interactions requires considering different pemutations
  res.allperms <- NULL
  for (mo in 1:nrow(allperms)) {
    neworder <- as.integer(allperms[mo,])
    Data.mo <- Data[,paste0("M",1:6)][,neworder]
    # head(Data.mo)
    colnames(Data.mo) <- paste0("M",1:6)
    Data.mo <- cbind(Data[,grep("M",names(Data),invert=TRUE,value=TRUE)],
                     Data.mo)
    
    # observed mean mediator value under each treatment
    M_A.hat <- NULL
    for (s in 1:6) {
      Mmean_A <- as.numeric(by(Data.mo[,paste0("M",s)],Data.mo$A,mean))
      names(Mmean_A) <- paste0("A",0:1)
      M_A.hat[[s]] <- Mmean_A
      rm(Mmean_A)
    }
    names(M_A.hat) <- names(models_M)
    
    # outcome model with mediator-mediator interactions
    fitY <- lm(models$Y_Mint,data=Data.mo)
    betas <- coef(fitY)[grep("M",names(coef(fitY)),value=TRUE)]
    betas.int <- NULL
    for (s in 1:6) {
      # regression coefficients for each mediator
      betas_s <- betas[grep(s,names(betas),value=TRUE)]
      # main effect
      betas_s.main <- betas_s[paste0(c("","A:"),paste0("M",s))]
      # mediator-mediator interaction terms
      betas_s.int <- betas_s[!(names(betas_s) %in% names(betas_s.main))]
      ## note order of coefficients
      # moderated effect of each mediator on outcome
      # group-specific mean values of other mediators
      M_nots_means <- NULL
      if (s>1) {
        M_nots_means <- c(
          M_nots_means,
          unlist(lapply(M_A.hat[paste0("M",(1:(s-1)))],"[","A0")))
      }
      if (s<6) {
        M_nots_means <- c(
          M_nots_means,
          unlist(lapply(M_A.hat[paste0("M",((s+1):6))],"[","A1")))
      }
      betas.int[s] <- sum(betas_s.main) + 
        (rep(M_nots_means,times=2) %*% betas_s.int)[,1]
    }
    # original label of mediators
    names(betas.int) <- paste0("ie",neworder,".Mint")
    betas.int <- betas.int[order(names(betas.int))]
    res.allperms[[mo]] <- betas.int*deltas
    rm(Data.mo)
    
    # only under original ordering
    if (all(neworder==(1:6))) {
      # observed mean mediator-mediator values under control
      MM_A.hat <- unlist(lapply(strsplit(intM,split=":"), function(mm) {
        mean(apply(Data[Data$A==0,mm],1,prod))
      }))
      names(MM_A.hat) <- intM
      # direct effect
      res["de.Mint"] <- coef(fitY)["A"] +
        sum(coef(fitY)[paste0("A:",names(M_A.hat))]*
              unlist(lapply(M_A.hat,"[","A0"))) +
        sum(coef(fitY)[paste0("A:",names(MM_A.hat))]*MM_A.hat)
        
      # mutual indirect effects
      covM <- cov(Data[Data$A==1,grep("M",colnames(Data))]) - 
        cov(Data[Data$A==0,grep("M",colnames(Data))])
      res.mu <- NULL
      for (s in 1:length(intM)) {
        idx <- as.integer(lapply(strsplit(strsplit(intM[s],":")[[1]],"M"),"[",2))
        res.mu[s] <- sum(betas[grep(intM[s],names(betas))])*covM[idx[1],idx[2]]
      }
      names(res.mu) <- paste0("mu.",intM)
      res.mu <- c(res.mu,"mu.sum"=sum(res.mu))
      res <- c(res,res.mu)
    }
    rm(fitY,betas)
  }
  res.allperms <- cbind("perm"=1:nrow(allperms),do.call(rbind,res.allperms))
  return(list(res,res.allperms))
}
# estimates based on observed data
res.obs <- OneEst(boot=FALSE)
# bootstrap samples
if (file.exists(file="interventional-lm-6M-allperms.Rdata")) {
  load(file="interventional-lm-6M-allperms.Rdata")
} else {
  res.boot <- replicate(1000, expr=OneEst(boot=TRUE))
  save(res.boot,file="interventional-lm-6M-allperms.Rdata")
}

# total, direct, mutual
res.obs.te <- res.obs[[1]]
## bootstrap SE and 95% CIs
res.boot.te <- data.frame(do.call(rbind,res.boot[1,]))
res.ci.te <- data.frame(
  cbind(res.obs.te,
        t(apply(res.boot.te, 2, function(x)
          c(sd(x),quantile(x,probs=c(.027,.975)))))))
xtable(t(apply(res.ci.te[c("te","de","mu.sum"),],1,function(x) {
  x <- format(round(x, 3), nsmall = 3)
  c(x[1:2],paste0("(",x[3],", ",x[4],")"))
})))

# moderated indirect effects
res.obs.ie <- res.obs[[2]]
## range of point estimates for different decompositions
res.ie.range <- t(apply(res.obs.ie,2,range)[,-1])

## bootstrap SE and 95% CIs
res.ie <- data.frame(do.call(rbind,res.boot[2,]))
res.ie.byperm <- by(res.ie,INDICES=res.ie$perm,function(x) {
  # observed estimate
  x.obs <- res.obs.ie[res.obs.ie[,"perm"]==unique(x$perm),]
  cbind("est"=x.obs,"se"=apply(x,2,sd),
        t(apply(x,2,quantile,probs=c(.025,.975))))
})

#### original permutation
xtable(t(apply(res.ie.byperm[[1]],1,function(x) {
  x <- format(round(x, 3), nsmall = 3)
  c(x[1:2],paste0("(",x[3],", ",x[4],")"))
})))

## range of 95% CIs for different decompositions 
res.ie.range.boot <- t(apply(do.call(rbind,lapply(res.ie.byperm, function(x) 
  c("lower"=x[-1,"2.5%"],"upper"=x[-1,"97.5%"]))),2,range))

library(xtable)
xtable(cbind(res.ie.range,
             res.ie.range.boot[1:6,1],
             res.ie.range.boot[7:12,2]),
       digits=3)
