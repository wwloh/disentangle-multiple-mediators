MultipleMediators <- function(
  outcome_name,
  mediator_names,
  treatment_name,
  covariate_names=NULL,
  obs_Data, 
  nboots) {
  
  var_names <- colnames(obs_Data)
  # replace outcome name
  var_names <- gsub(pattern=outcome_name,replacement="Y",x=var_names)
  # replace treatment name
  var_names <- gsub(pattern=treatment_name,replacement="A",x=var_names)
  # number of mediators  
  p <- length(mediator_names)
  for (s in 1:p) {
    # replace mediator names in turn
    var_names <- gsub(pattern=mediator_names[s],
                      replacement=paste0("M",s),
                      x=var_names)
  }
  colnames(obs_Data) <- var_names
  l_names <- covariate_names # confounders
  
  # specify different models ==================================================
  models <- list()
  # total effect
  models[["total_effect"]] <- as.formula(
    paste0("Y~A+",paste(l_names,collapse="+")))
  # mediator models
  models_M <- sapply(1:p, function(s) 
    as.formula(paste0("M",s,"~A+",paste(l_names,collapse="+"))))
  names(models_M) <- paste0("M",1:p)
  # outcome model with interactions
  intM <- NULL
  for (s in 1:p) {
    for (q in (s+1):p) {
      if (q > s & q <= p) {
        intM <- c(intM,paste(paste0("M",s),paste0("M",q),sep=":"))  
      }
    }
  }
  if (length(intM) != choose(p,2)) return(NULL)
  
  models[["Y_Mint"]] <- as.formula(
    paste0("Y ~ A +", paste0("M",1:p,collapse="+"), 
           "+", paste0("A:M",1:p,collapse="+"),
           "+", paste(intM,collapse="+"),
           "+", paste(paste0("A:",intM),collapse="+"),
           "+", paste(l_names,collapse="+")))
  
  # consider all permutations of mediator ordering for sensitivity analysis ===
  ## appropriate only for small number of mediators
  allperms <- expand.grid(lapply(1:p, function (x) p:1))
  allperms <- allperms[apply(allperms, 1, function(x) 
    length(unique(x))==length(x)),]
  row.names(allperms) <- colnames(allperms) <- NULL
  
  # fit models ================================================================
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
    for (s in 1:p) {
      fit_M <- lm(models_M[[s]],data=Data)
      deltas[s] <- coef(fit_M)["A"]
      rm(fit_M)
    }
    names(deltas) <- names(models_M)
    
    # mediator-mediator interactions requires considering different pemutations
    res.allperms <- NULL
    for (mo in 1:nrow(allperms)) {
      neworder <- as.integer(allperms[mo,])
      Data.mo <- Data[,paste0("M",1:p)][,neworder]
      # head(Data.mo)
      colnames(Data.mo) <- paste0("M",1:p)
      Data.mo <- cbind(Data[,grep("M",names(Data),invert=TRUE,value=TRUE)],
                       Data.mo)
      
      # observed mean mediator value under each treatment
      M_A.hat <- NULL
      for (s in 1:p) {
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
      for (s in 1:p) {
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
        if (s<p) {
          M_nots_means <- c(
            M_nots_means,
            unlist(lapply(M_A.hat[paste0("M",((s+1):p))],"[","A1")))
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
      if (all(neworder==(1:p))) {
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
  res.boot <- replicate(nboots, expr=OneEst(boot=TRUE))
  # for storing formatted results
  res_out <- list()
  
  # total, direct, mutual
  res.obs.te <- res.obs[[1]]
  ## bootstrap SE and 95% CIs
  res.boot.te <- data.frame(do.call(rbind,res.boot[1,]))
  res.ci.te <- data.frame(
    cbind(res.obs.te,
          t(apply(res.boot.te, 2, function(x)
            c(sd(x),quantile(x,probs=c(.027,.975)))))))
  colnames(res.ci.te) <- c("est","se","2.5%","97.5%")
  res_out[["te"]] <- res.ci.te[c("te","de","mu.sum"),]
  row.names(res_out[["te"]]) <- c("Total","Direct","Mutual Dependence")
  
  # moderated indirect effects
  res.obs.ie <- res.obs[[2]]
  
  ## bootstrap SE and 95% CIs
  res.ie <- data.frame(do.call(rbind,res.boot[2,]))
  res.ie.byperm <- by(res.ie,INDICES=res.ie$perm,function(x) {
    # observed estimate
    x.obs <- res.obs.ie[res.obs.ie[,"perm"]==unique(x$perm),]
    cbind("est"=x.obs,"se"=apply(x,2,sd),
          t(apply(x,2,quantile,probs=c(.025,.975))))
  })
  names(res.ie.byperm) <- paste0("perm.",names(res.ie.byperm))
  res_out[["ie_allperms"]] <- lapply(res.ie.byperm,function(x) {
    x <- x[-1,]
    row.names(x) <- mediator_names
    return(x)
  })
  
  ## range of point estimates and 95% CIs for different decompositions 
  res_out[["ie"]] <- lapply(mediator_names, function(m) {
    m_out <- apply(do.call(rbind,lapply(res_out[["ie_allperms"]], function(x) 
      x[m,c("est","2.5%","97.5%")])),2,range)
    row.names(m_out) <- c("min","max")
    return(m_out)
  })
  names(res_out[["ie"]]) <- mediator_names
  
  return(res_out)
}
