utils::globalVariables("estpx")
utils::globalVariables("newx")


##### Linear #####

### Basic ###
WEE_cts <- function(formula, D, data, pd_pop, boot=0, ...) {
  mf = model.frame(formula, data=data)
  y = model.response(mf, "numeric")
  namesx = all.vars(formula)[-1]
  xx = model.matrix(attr(mf, "terms"), data=mf)
  x = matrix(xx[, -1], nrow=nrow(data), byrow=F); colnames(x)=namesx
  temp.data = data.frame(cbind(D, y, x))
  f = update(formula, y  ~ . )
  n1 = sum(D == 1); n0 = sum(D == 0)      
  
  # compute the weight p(D|X)
  gamma = coef(glm(D~x, family="binomial"))
  PO = function(gamma0){
    gamma[1] = gamma0
    (mean( exp(xx%*%gamma)/(1+exp(xx%*%gamma)))  -pd_pop)^2
  }
  gamma[1] = suppressWarnings(optim(gamma[1],PO)$par)
  temp.data$estpx = exp(xx%*%gamma)/(1+exp(xx%*%gamma))
  
  # estimate py in cases and controls separately  
  pyD1 = lm(f, data=temp.data[which(temp.data$D == 1),]) # fit the case
  pyD0 = lm(f, data=temp.data[which(temp.data$D == 0),]) # fit the control
  py1 = predict(pyD0, newdata=data.frame(temp.data[which(temp.data$D == 1),])) # generate pseudo control
  py0 = predict(pyD1, newdata=data.frame(temp.data[which(temp.data$D == 0),])) # generate pseudo case
  data1 = cbind( rep(0, n1), py1, temp.data[which(D == 1), c(namesx, "estpx")] ) 
  data0 = cbind( rep(1, n0), py0, temp.data[which(D == 0), c(namesx, "estpx")] )
  colnames(data1)[1:2] = c("D", "y")
  colnames(data0)[1:2] = c("D", "y")
  alldat = rbind(temp.data, data1, data0)
  alldat$estpx[which(alldat$D==0)] = 1- alldat$estpx[which(alldat$D == 0)]
  
  # the point estmate
  # with optional arguments in lm
  dots <- list(...)
  arg_cur = modifyList(list(f = f, data = alldat, weights = alldat$estpx),dots)
  coef=do.call(lm,arg_cur)$coef
  
  # bootstrap SE	
  if(boot == 0){
    TAB = list(coefficients = coef)
    TAB
  } else 
  { 
    sample_cases = temp.data[which(temp.data$D == 1),]
    sample_controls = temp.data[which(temp.data$D == 0),]
    
    bootcoef = NULL	
    
    nCores=detectCores()
    registerDoParallel(nCores - 1)
  
    bootcoef = foreach ( iboot=1:boot, .combine=rbind ) %dopar% {
      
      boot_cases_sample = sample_cases[sample(n1, n1, replace=T),]
      boot_controls_sample = sample_controls[sample(n0, n0, replace=T),]
      bootsample = rbind(boot_cases_sample,boot_controls_sample)
      
      bootmf = model.frame(f, data=bootsample); 
      bootxx = model.matrix(attr(bootmf, "terms"), data=bootmf)
      
      # compute the weight p(D|X)
      gamma = coef(glm(update(formula, D  ~ . ), family="binomial", data=bootsample))
      PO = function(gamma0){
        gamma[1]=gamma0
        (mean( exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma)))  -pd_pop)^2
      }
      gamma[1] = suppressWarnings(optim(gamma[1],PO)$par)
      bootsample$estpx = exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma))
      
      pyD1 = lm(f, data=boot_cases_sample)
      pyD0 = lm(f, data= boot_controls_sample)
      py1 = predict(pyD0, newdata=data.frame(boot_cases_sample)) 
      py0 = predict(pyD1, newdata=data.frame(boot_controls_sample))
      data1 = cbind( rep(0, n1), py1, boot_cases_sample[, c(namesx,"estpx")] ) 
      data0 = cbind( rep(1, n0), py0, boot_controls_sample[, c(namesx,"estpx")] )
      colnames(data1)[1:2] = c("D",  "y")
      colnames(data0)[1:2] = c("D", "y")
      
      alldat = rbind(bootsample, data1, data0)
      alldat$estpx[which(alldat$D==0)] = 1- alldat$estpx[which(alldat$D==0)]
      bootcoef = lm(f, data=alldat, weights=estpx)$coef
      
    }
    
    var = apply(bootcoef, 2, var)
    cov = cov(bootcoef)
    chisq = coef^2/var
    pvalue = pchisq(chisq, df=1, lower.tail=F)
    
    TAB = list(coefficients=coef, StdErr=sqrt(var), Chisq = chisq, p.value=pvalue, Covariance=cov)
    TAB    
  }
}

### Make it Generic ###
WEE.linear <- function(formula, D, data, pd_pop, boot=0, ...) UseMethod("WEE.linear")

WEE.linear <- function(formula, D, data, pd_pop, boot=0, ...){
  D = as.numeric(D)
  data = as.data.frame(data)
  pd_pop = as.numeric(pd_pop)
  
  est = WEE_cts(formula, D, data, pd_pop, boot=boot)
  est$call = match.call()
  class(est) = "WEE.linear"
  est
}

print.WEE.linear <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

### Basic summary ###
summary.WEE.linear <- function(object, ...) {
  TAB = cbind(Estimates=object$coefficients, StdErr=object$StdErr, Chisq=object$Chisq, p.value=object$p.value)

  if (length(object) == 6){
    res = list(call=object$call,
               coefficients=TAB,
               var.cov=object$Covariance)
  } else
  {
    res = list(call=object$call,
               coefficients=TAB)
  }
  
  class(res) = "summary.WEE.linear"
  res
}

print.summary.WEE.linear <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat("Coefficients:\n")
  if (ncol(x$coefficients) > 1) {
    printCoefmat(x$coefficients,P.values = TRUE, has.Pvalue = TRUE)
  } 
  else{
    print(as.data.frame(x$coefficients))
  }
  cat("\n")
  
  if (length(x) == 3){
    cat("Covariance:\n")
    print(as.data.frame(x$var.cov))
  }
}

### Predict based on the model fitted ###
predict.WEE.linear <- function(object, newx, ...){
  len = length(object$coefficients) - 1
  newx = matrix(as.numeric(as.matrix(newx)),ncol = len)
  colnames (newx) = names(object$coefficients)[-1]
  
  prediction =  cbind(1, newx) %*% as.matrix(object$coefficients)
  colnames(prediction) = "Predicted Values"
  
  if (length(object) == 6) {
    se = sqrt( diag( cbind(1, newx) %*% as.matrix(object$Covariance) %*% t(cbind(1, newx)) ) )
    prediction = cbind(prediction, se)
    prediction = as.data.frame(prediction)
    colnames(prediction) = c("Predicted Value", "StdErr of Predicted Value")
  }
  
  return(list(prediction = prediction, newx = newx))
}


##### Logistic #####

### Basic ###
WEE_binary <-function(formula, D, data, pd_pop, iter=5, boot=0, ...) {
  mf = model.frame(formula, data=data)
  y = model.response(mf, "numeric")
  namesx = all.vars(formula)[-1]
  xx = model.matrix(attr(mf, "terms"), data=mf)
  x = matrix(xx[,-1], nrow=nrow(data), byrow=F); colnames(x) = namesx
  temp.data = data.frame(cbind(D, y, x))
  f = update(formula, y  ~ . )
  n1 = sum(D == 1); n0 = sum(D == 0)      
  
  # compute the weight p(D|X)
  gamma = coef(glm(D~x, family="binomial"))
  PO = function(gamma0){
    gamma[1] = gamma0
    (mean( exp(xx%*%gamma)/(1+exp(xx%*%gamma)))  -pd_pop)^2
  }; gamma[1] = suppressWarnings(optim(gamma[1],PO)$par)
  temp.data$estpx = exp(xx%*%gamma)/(1+exp(xx%*%gamma))
  
  # estimate P(Y=1) in cases and controls separately
  pyD1 = glm(f, family="binomial", data=temp.data[which(temp.data$D==1),])
  pyD0 = glm(f, family="binomial", data=temp.data[which(temp.data$D==0),])
  pred1 = predict(pyD0, newdata=data.frame(temp.data[which(temp.data$D==1),])) # generate pseudo control
  pred0 = predict(pyD1, newdata=data.frame(temp.data[which(temp.data$D==0),])) # generate pseudo case
  py1 = exp(pred1)/(1+exp(pred1)) 
  py0 = exp(pred0)/(1+exp(pred0)) 
  
  # generate pseudo observations for $iter$ times and get the averaged $iter$ estimates as coefficient
  pseudo = NULL
  for (iiter in 1:iter) {
    pseudo1 = rbinom(n1, 1, py1)
    pseudo0 = rbinom(n0, 1, py0)
    data1 = cbind( rep(0, n1), pseudo1, temp.data[which(D==1), c(namesx, "estpx")] ) 
    data0 = cbind( rep(1, n0), pseudo0, temp.data[which(D==0), c(namesx, "estpx")] )
    colnames(data1)[1:2] = c("D",  "y")
    colnames(data0)[1:2] = c("D", "y")
    alldat = rbind(temp.data, data1, data0)
    alldat$estpx[which(alldat$D == 0)] = 1- alldat$estpx[which(alldat$D == 0)]
    ## add optional arguments in glm
    dots <- list(...)
    arg_cur = modifyList(list(f = f,family = "binomial", data = alldat, weights = alldat$estpx),dots)
    # pseudo = rbind(pseudo, suppressWarnings(glm(f, family="binomial", data=alldat, weights=estpx)$coef))
    pseudo = rbind(pseudo,suppressWarnings(do.call(glm,arg_cur)$coef))
  }
  # the point estmate
  coef=apply(pseudo, 2, mean)
  
  # bootstrap SE	
  if(boot==0){
    TAB = list(coefficients =coef, Oddsratio = exp(coef))
    TAB
  } else
  {
    sample_cases = temp.data[which(temp.data$D == 1),]
    sample_controls = temp.data[which(temp.data$D == 0),]
    
    bootcoef = NULL	
    
    nCores=detectCores()
    registerDoParallel(nCores - 1)
    
    bootcoef = foreach ( iboot=1:boot, .combine=rbind ) %dopar% {
      
      boot_cases_sample = sample_cases[sample(n1, n1, replace=T),]
      boot_controls_sample = sample_controls[sample(n0, n0, replace=T),]
      bootsample = rbind(boot_cases_sample, boot_controls_sample)
      bootmf = model.frame(f, data=bootsample)
      bootxx = model.matrix(attr(bootmf, "terms"), data=bootmf)
     
      # compute the weight p(D|X)
      gamma = coef(glm(update(formula, D  ~ . ), family="binomial", data=bootsample))
      PO = function(gamma0){
        gamma[1] = gamma0
        (mean( exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma)))  -pd_pop)^2
      }; gamma[1] = suppressWarnings(optim(gamma[1],PO)$par)
      bootsample$estpx = exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma))
      
      pyD1 = glm(f, family="binomial", data= boot_cases_sample)
      pyD0 = glm(f, family="binomial", data= boot_controls_sample)
      pred1 = predict(pyD0, newdata=data.frame(boot_cases_sample)) # pseudo control
      pred0 = predict(pyD1, newdata=data.frame(boot_controls_sample))
      py1 = exp(pred1)/(1+exp(pred1)) 
      py0 = exp(pred0)/(1+exp(pred0)) 
      
      # generate pseudo observations for T(=10) times and get the averaged T estimates as coefficient
      pseudo = NULL
      for (iiter in 1:iter) {
        pseudo1 = rbinom(n1, 1, py1)  
        pseudo0 = rbinom(n0, 1, py0)
        data1 = cbind( rep(0, n1), pseudo1, boot_cases_sample[, c(namesx,  "estpx")] ) 
        data0 = cbind( rep(1, n0), pseudo0, boot_controls_sample[, c(namesx, "estpx")] )
        colnames(data1)[1:2] = c("D",  "y")
        colnames(data0)[1:2] = c("D",  "y")
        
        alldat = rbind(bootsample, data1, data0)
        alldat$estpx[which(alldat$D == 0)] = 1- alldat$estpx[which(alldat$D == 0)]
        
        pseudo = rbind(pseudo, suppressWarnings(glm(f, family="binomial", data= alldat, weights = estpx)$coef))
      }
      
      bootcoef = apply(pseudo, 2, mean)
      
    }
    
    var = apply(bootcoef, 2, var)
    cov = cov(bootcoef)
    wald = coef^2/var
    pvalue = pchisq(wald, df=1, lower.tail=F)
  
    TAB = list(coefficients=coef, Oddsratio = exp(coef), StdErr=sqrt(var), Wald=wald, p.value=pvalue, Covariance=cov)
    TAB    
  }
}

### Make it generic ###
WEE.logistic<-function(x, ...) UseMethod("WEE.logistic")

WEE.logistic<-function(formula, D, data, pd_pop, iter=5, boot=0, ...) {
  D<-as.numeric(D)
  data<-as.data.frame(data)
  pd_pop<-as.numeric(pd_pop)
  
  est <- WEE_binary(formula, D, data, pd_pop, iter=iter, boot=boot)
  est$call <- match.call()
  
  class(est) <- "WEE.logistic"
  est
}

print.WEE.logistic <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

### Basic summary ###
summary.WEE.logistic <- function(object, ...) {
  TAB = cbind(Estimate=object$coefficients, Oddsratio = exp(object$coefficients), StdErr=object$StdErr, Wald=object$Wald, p.value=object$p.value)
  if (length(object) == 6) {
    res = list(call=object$call,
               coefficients=TAB,
               var.cov=object$Covariance) 
  }
  else {
    res = list(call=object$call,
               coefficients=TAB)
  }
  class(res) = "summary.WEE.logistic"
  res
}

print.summary.WEE.logistic <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat("Coefficients:\n")
  if (ncol(x$coefficients) > 1) {
    printCoefmat(x$coefficients,P.values = TRUE, has.Pvalue = TRUE)
  } 
  else {
    print(x$coefficients)
  }
  cat("\n")
  
  if (length(x) == 3){
    cat("Covariance:\n")
    print(as.data.frame(x$var.cov))
  }
}

### Prediction based on the model fitted ###
predict.WEE.logistic <- function(object, newx, ...){
  len = length(object$coefficients) - 1
  newx = matrix(as.numeric(as.matrix(newx)),ncol = len)
  colnames (newx) = names(object$coefficients)[-1]
  
  linear = cbind(1, newx) %*% as.matrix(object$coefficients)
  response = exp(linear)/(1+exp(linear)) #prediction of response
  prediction = cbind(linear, response)
  colnames(prediction) = c("Linear Predictor", "Predicted Response")
  
  if (length(object) == 6){
    se_linear = sqrt( diag( cbind(1,newx) %*% object$Covariance %*% t(cbind(1, newx)) ) )
    se_response = exp(2 * linear) / (1 + exp(linear))^4 * se_linear
    prediction = cbind(prediction, se_linear, se_response)
    colnames(prediction)[3:4] = c("StdErr of Linear Predictor", "StdErr of Predicted Response")
  }
  
  return(list(prediction = prediction, newx = newx))
}

##### Quantile #####

### Basic ###
WEE_rq <- function(formula, D, data, pd_pop, tau, iter=5, boot=0,  ...) {
  mf = model.frame(formula, data=data)
  y = model.response(mf, "numeric")
  namesx = all.vars(formula)[-1]
  xx = model.matrix(attr(mf, "terms"), data=mf)
  x = matrix(xx[,-1], nrow=nrow(data), byrow=F); colnames(x)=namesx; p=dim(x)[2]
  temp.data = data.frame(cbind(D, y, x))
  f = update(formula, y  ~ . )
  n1 = sum(D == 1); n0 = sum(D == 0)      
  ltau = length(tau); flag = seq(1, ltau)  
  
  # compute the weight p(D|X)
  gamma = coef(glm(D~x, family="binomial"))
  PO = function(gamma0){
    gamma[1] = gamma0
    (mean( exp(xx%*%gamma)/(1+exp(xx%*%gamma)))  -pd_pop)^2
  }; gamma[1] = suppressWarnings(optim(gamma[1],PO)$par)
  temp.data$estpx = exp(xx%*%gamma)/(1+exp(xx%*%gamma))
  
  # estimate py in cases and controls separately
  pyD1 = quantreg::rq(f,  tau=-1, data=temp.data[which(temp.data$D == 1),])
  pyD0 = quantreg::rq(f,  tau=-1, data=temp.data[which(temp.data$D == 0),])
  xcase = xx[which(D == 1),]; xcontrol = xx[which(D == 0),]
  pred1 =	as.matrix(xcase) %*% pyD0$sol[4:length(pyD0$sol[,1]), ] # pseudo control
  pred0 = as.matrix(xcontrol) %*% pyD1$sol[4:length(pyD0$sol[,1]), ]
  
  step1=apply(cbind(pred1[,1], pred1), 1, function(x) rearrange(stepfun(pyD0$sol[1,], x)))
  step0=apply(cbind(pred0[,1], pred0), 1, function(x) rearrange(stepfun(pyD1$sol[1,], x)))

  # generate pseudo observations for $iter$ times and get the averaged $iter$ estimates as coefficient
  pseudo = NULL
  for (iiter in 1:iter) {
    pseudo1 = unlist(lapply(step1, function(x) x(runif(1))))
    pseudo0 = unlist(lapply(step0, function(x) x(runif(1))))
    data1 = cbind( rep(0, n1), pseudo1, temp.data[which(D == 1), c(namesx, "estpx")] ) 
    data0 = cbind( rep(1, n0), pseudo0, temp.data[which(D == 0), c(namesx, "estpx")] )
    colnames(data1)[1:2] = c("D", "y")
    colnames(data0)[1:2] = c("D", "y")
    alldat = rbind(temp.data, data1, data0)
    alldat$estpx[which(alldat$D==0)] = 1 - alldat$estpx[which(alldat$D == 0)]
    ## Optional arguments
    # pseudo = rbind(pseudo, quantreg::rq(f, tau=tau, data=alldat, weights=estpx, method=method)$coef)
    dots <- list(...)
    arg_cur = modifyList(list(f = f, tau = tau, data = alldat, weights = alldat$estpx),dots)
    pseudo = rbind(pseudo, do.call(quantreg::rq,arg_cur)$coef)
  }
  # the point estmate
  if (ltau == 1) {
    coef = list(apply(pseudo, 2, mean))
  } else  
  {  
    coef_est = do.call("rbind", lapply(  (seq(0:p)%%(p+1)), function(x) apply(pseudo[which(seq(nrow(pseudo))%%(p+1) == x),], 2, mean)))
    coef = split(coef_est, col(coef_est))   
  } 
  
  # bootstrap SE	
  if(boot == 0){ 
    names(coef) = paste("tau=", tau)
    for (iflag in flag) {
      names(coef[[iflag]]) = c("intercept", namesx)
    }
    TAB = vector("list", ltau)
    for (iflag in flag) {
      TAB[[iflag]] = list(coefficients = coef[[iflag]])					
    }
    names(TAB)=paste('tau=', tau)
  } else
  {
    sample_cases = temp.data[which(temp.data$D == 1),]
    sample_controls = temp.data[which(temp.data$D == 0),]
    
    bootcoef=NULL
    
    nCores=detectCores()
    registerDoParallel(nCores - 1)
    
    bootcoef = foreach ( iboot=1:boot, .combine=rbind, .packages='quantreg' ) %dopar% {
  
      boot_cases_sample = sample_cases[sample(n1, n1, replace=T),]
      boot_controls_sample = sample_controls[sample(n0, n0, replace=T),]
      
      bootsample = rbind(boot_cases_sample, boot_controls_sample)
      bootmf = model.frame(f, data=bootsample)
      bootxx = model.matrix(attr(bootmf, "terms"), data=bootmf)
      xcase = bootxx[which(bootsample$D == 1), ]
      xcontrol = bootxx[which(bootsample$D == 0), ]
      
      # compute p(D|X)
      gamma = coef(glm(update(formula, D  ~ . ), family="binomial", data=bootsample))
      PO = function(gamma0){
        gamma[1] = gamma0
        (mean( exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma)))  -pd_pop)^2
      }; gamma[1] = suppressWarnings(optim(gamma[1],PO)$par)
      bootsample$estpx = exp(bootxx%*%gamma)/(1+exp(bootxx%*%gamma))
      
      pyD1 = quantreg::rq(f, tau=-1, data=boot_cases_sample)
      pyD0 = quantreg::rq(f, tau=-1, data=boot_controls_sample)
      pred1 = as.matrix(xcase) %*% pyD0$sol[4:length(pyD0$sol[,1]), ]
      pred0 = as.matrix(xcontrol) %*% pyD1$sol[4:length(pyD0$sol[,1]), ]
      
      step1 = apply(cbind(pred1[,1], pred1), 1, function(x) rearrange(stepfun(pyD0$sol[1,], x)))
      step0 = apply(cbind(pred0[,1], pred0), 1, function(x) rearrange(stepfun(pyD1$sol[1,], x)))
      
      # generate pseudo observations for T(=10) times and get the averaged T estimates as coefficient
      pseudo = NULL
      for (iiter in 1:iter) {
        
        pseudo1 = unlist(lapply(step1, function(x) x(runif(1))))
        pseudo0 = unlist(lapply(step0, function(x) x(runif(1))))
        
        data1 = cbind( rep(0, n1), pseudo1, boot_cases_sample[, c(namesx,  "estpx")] ) 
        data0 = cbind( rep(1, n0), pseudo0, boot_controls_sample[, c(namesx, "estpx")] )
        colnames(data1)[1:2] = c("D", "y")
        colnames(data0)[1:2] = c("D", "y")
        
        alldat = rbind(bootsample, data1, data0)
        alldat$estpx[which(alldat$D == 0)] = 1- alldat$estpx[which(alldat$D == 0)]
        
        dots <- list(...)
        arg_cur = modifyList(list(f = f, tau = tau, data = alldat, weights = alldat$estpx),dots)
        pseudo = rbind(pseudo, do.call(quantreg::rq,arg_cur)$coef)
        #pseudo = rbind(pseudo, quantreg::rq(f, tau=tau, data=alldat, weights=estpx)$coef)
        
      }
      
      
      if (ltau == 1 ) {
        bootcoef = apply(pseudo, 2, mean)
      } else {  
        bootcoef = t(do.call("rbind", lapply(  (seq(0:p)%%(p+1)), function(x) apply(pseudo[which(seq(nrow(pseudo))%%(p+1)== x),], 2, mean)))  )  
        
      }
    }
    
    
    if (ltau == 1) {
      bootcoef = list(bootcoef)
    } else{
      temp = lapply(   (1:ltau)%%ltau, function(x) bootcoef[which(seq(nrow(bootcoef))%%ltau== x), ] )
      bootcoef = temp
    }
    
    var = lapply(seq(1:ltau), function(x) apply(bootcoef[[x]], 2, var))
    cov = lapply(seq(1:ltau), function(x) cov(bootcoef[[x]]))
    wald = lapply(seq(1:ltau), function(x) coef[[x]]^2/var[[x]])
    pvalue = lapply(seq(1:ltau), function(x) pchisq(wald[[x]], df=1, lower.tail=F))
    
    TAB = vector("list", ltau)
    for (iflag in flag) {
        TAB[[iflag]] = list(coefficients=coef[[iflag]], StdErr=sqrt(var[[iflag]]), 
                            Wald=wald[[iflag]], p.value=pvalue[[iflag]], Covariance=cov[[iflag]])	
        names(TAB[[iflag]]$coefficients) = c("intercept", namesx)
        names(TAB[[iflag]]$StdErr) = c("intercept", namesx)
        names(TAB[[iflag]]$Wald) = c("intercept", namesx)
        names(TAB[[iflag]]$p.value) = c("intercept", namesx)
        rownames(TAB[[iflag]]$Covariance) = colnames(TAB[[iflag]]$Covariance) = c("intercept", namesx)
    }
    names(TAB) = paste('tau=', tau) 
  }	
  LIS = list(coefficients=TAB, tau=tau)
  LIS
}

### Make it generic ###
WEE.quantile <- function(formula, D, data, pd_pop, tau, iter=5, boot=0, ...) UseMethod("WEE.quantile")  

WEE.quantile <- function(formula, D, data, pd_pop, tau, iter=5, boot=0, ...){
  D = as.numeric(D)
  data = as.data.frame(data)
  pd_pop = as.numeric(pd_pop)
  
  est = WEE_rq(formula, D, data, pd_pop, tau, iter=iter, boot=boot)
  est$call = match.call()
  class(est) = "WEE.quantile"
  est
}

### ??? how to use print(object) or just object ###
print.WEE.quantile <- function(x, ...){
  x$TAB = x$coefficients
  tau_len = length(x$tau)
  boot_len = length(x$TAB[[1]])-1 #=0 if boot = 0
  object = x$TAB
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  
  if (tau_len == 1){
    if (boot_len == 0)
    { out = as.data.frame(object[[1]])
    names(out) = "Estimates"
    print(out)}
    else {
      if (length(object[[1]]) == 4)
      {out = as.data.frame(cbind(Estimates = object[[1]]$coefficients,StdErr = object[[1]]$StdErr,
                                 Wald = object[[1]]$Wald,p.value = object[[1]]$p.value))
      print(out)}
      else
      {out = as.data.frame(cbind(Estimates = object[[1]]$coefficients,StdErr = object[[1]]$StdErr,
                                 Wald = object[[1]]$Wald,p.value = object[[1]]$p.value))
      covariance = as.data.frame(object[[1]]$Covariance)
      out = list(Coefficients = out, Covariance = covariance)
      print(out)}
    }
  }
  
  else{ # tau_len > 1
    if (boot_len != 0)
    { 
      if (length(object[[1]]) == 5) {
        for (i in 1:tau_len){
          out = as.data.frame(cbind(Estimates = object[[i]]$coefficients,StdErr = object[[i]]$StdErr,Wald = object[[i]]$Wald,p.value = object[[i]]$p.value))
          cat(names(object)[i],"\n",sep = "")
          out = list(out, var.cov = object[[i]]$Covariance)
          names(out) = c("Coefficients","Covariance")
          print(out) }}
      else {
        for (i in 1:tau_len)
        {cat(names(object)[i],"\n",sep = "")
          
          out = as.data.frame(cbind(Estimates = object[[i]]$coefficients,StdErr = object[[i]]$StdErr,Wald = object[[i]]$Wald,p.value = object[[i]]$p.value))
          print(out)
        }
      }
    }
    else {
      if (length(object[[1]]) == 5) {
        for (i in 1:tau_len){
          out = as.data.frame(object[[i]]$coefficients)
          cat(names(object)[i],"\n",sep = "")
          out = list(out, var.cov = object[[1]]$Covariance)
          names(out) = c("Coefficients","Covariance")
          print(out) }}
      else {
        for (i in 1:tau_len)
        {
          cat(names(object)[i],"\n",sep = "")
          out = as.data.frame(object[[i]]$coefficients)
          names(out) = "Estimates"
          print(out)
        } }
    }
  }
}

### Basic summary ###
summary.WEE.quantile <- function(object, ...){
  tau_cur = object$tau
  call = object$call
  object = object$coefficients
  # dim, number of taus
  taul = length(tau_cur)
  
  # boot = 0
  if (length(object[[1]]) == 1)
  {
    if (taul < 2){
      TAB =  list(Coefficients = as.data.frame(object[[1]]))
      names(TAB) = "Coefficients"
    } else
    {
      TAB = vector(mode = "list",length = taul)
      for (dimi in 1:taul){
        TAB[[dimi]] = list(Coefficients = as.data.frame(object[[dimi]]))
      }
      names(TAB) = paste('tau=',tau_cur)
    }
  }
  # boot > 0 
  else
  {
    if (taul < 2)
    {
      TAB = list(cbind(Estimate=object[[1]]$coefficients, 
                       StdErr=object[[1]]$StdErr, 
                       Wald=object[[1]]$Wald,
                       p.value=object[[1]]$p.value))
      names(TAB) = "Coefficients"
      if (!is.null(object[[1]]$Covariance)){
        TAB = list(as.data.frame(TAB$Coefficients), object[[1]]$Covariance)
        names(TAB) = c("Coefficients","Covariance")
      }
    }
    else{
      TAB = vector(mode = "list", length = taul)
      for (dimi in 1:taul){
        TAB[[dimi]] = list(cbind(Estimate=object[[dimi]]$coefficients, StdErr=object[[dimi]]$StdErr, 
                                 Wald=object[[dimi]]$Wald, p.value=object[[dimi]]$p.value))
        names(TAB[[dimi]]) = "Coefficients"
        if (!is.null(object[[dimi]]$Covariance)){
          TAB[[dimi]] = list(as.data.frame(TAB[[dimi]]$Coefficients), object[[dimi]]$Covariance)
          names(TAB[[dimi]]) = c("Coefficients","Covariance") 
        } 
      }
      names(TAB) = paste('tau =',tau_cur)
    }
  }
  res = list(call=call,
             coefficients = TAB,
             tau=tau_cur)
  class(res) = "summary.WEE.quantile"
  res
}

print.summary.WEE.quantile <- function (x, ...){
  
  # ntau: number of taus
  ntau = length(x$tau) 
  cat("Call:\n")
  print(x$call)
  cat("\n") 
  
  # boot: whether employed boots
  boot = ncol(as.data.frame(x$coefficients[[1]]))
  
  if (boot == 1) # no boots
  { 
    if (ntau == 1)
      print(as.data.frame(x$coefficients))
    
    else {
      for (i in 1:ntau){
        cat(names(x$coefficients)[i])
        cat("\n")
        print(as.data.frame(x$coefficients[[i]]))
        cat("\n")}
    }
  }
  
  else {
    
    if (ntau == 1)
    {cat ("Coefficients:\n")
      printCoefmat(as.data.frame(x$coefficients$Coefficients), 
                   P.values = TRUE, has.Pvalue = TRUE)
      cat ("\n")
      if (length(x$coefficients)==2) {
        cat ("Covariance:\n")
        print(as.data.frame(x$coefficients$Covariance))
        cat ("\n")}
    }
    
    else {
      for (i in 1:ntau)
      { cat(names(x$coefficients)[i])
        cat("\n")
        cat("Coefficients:\n")
        printCoefmat(data.frame(x$coefficients[[i]]$Coefficients),
                     P.values = TRUE, has.Pvalue = TRUE)
        cat("\n")
        if (length(x$coefficients[[1]])==2) {
          cat("Covariance:\n")
          print(as.data.frame(x$coefficients[[i]]$Covariance))
          cat("\n")}}
    }
  }
}

### Predict based on the fitted model ###
predict.WEE.quantile <- function(object, newx, ...){
  
  object$TAB = object$coefficients
  name_tau = paste("tau=", object$tau, sep="")
  len = length(name_tau)
  
  len_coef = length(object$TAB[[1]]$coefficients)
  newx = matrix(as.numeric(as.matrix(newx)),ncol = (len_coef-1))
  colnames (newx) = names(object$TAB[[1]]$coefficients)[-1]
  
  prediction = matrix(0, ncol=len, nrow=nrow(newx))
  for (ii in 1:len){
    prediction[, ii] = cbind(1, newx) %*% object$TAB[[ii]]$coefficients 
  }
  colnames(prediction) = name_tau
  p = list(predicted_value = prediction, newx = newx)
  
  if (!is.null(object$TAB[[1]]$Covariance)) {
    StdErr = matrix(0, ncol=len, nrow=nrow(newx))
    for (ii in 1:len){
      StdErr[, ii] = sqrt( diag( cbind(1, newx) %*% object$TAB[[ii]]$Covariance %*% t(cbind(1, newx)) ) )
    }
    colnames(StdErr) = name_tau
    p = list(predicted_value = prediction, StdErr_predicted_value = StdErr, newx = newx)
  } 
  
  class(p) = "predict.WEE.quantile"
  p
}

### Plot estimated coefficients against taus ###
plot.WEE.quantile <- function(x, CI = FALSE, level=.95, index=1, ...){
  x$TAB = x$coefficients
  
  len = length(x$tau)
  
  len_coef = length(x$TAB[[1]]$coefficients)
  item = names(x$TAB[[1]]$coefficients)
  
  coef =  matrix(0, ncol=len, nrow=len_coef)
  for (ii in 1:len){
    coef[, ii] = x$TAB[[ii]]$coefficients
  }
  
  if (CI){
    if (is.null(x$TAB[[1]]$StdErr)) {print("No confidence interval is plotted when zero boot in fitting step")}
    else{
      record = matrix(0, ncol=len, nrow=2*len_coef)
      for (ii in 1:len){
        record[(1:len_coef), ii] = x$TAB[[ii]]$coefficients + qnorm(1-(1-level)/2)*x$TAB[[ii]]$StdErr
        record[((len_coef+1):(2*len_coef)), ii] = x$TAB[[ii]]$coefficients - qnorm(1-(1-level)/2)*x$TAB[[ii]]$StdErr
      }
    }
  }
  
  if (!exists("record")){
    for (jj in 1:length(index)){
      dots <- list(...)
      arg_cur = modifyList(list(x$tau, coef[index[jj],], main=paste0("Estimated coefficient of ", item[index[jj]]), 
           xlab="taus", ylab= paste0("Estimated coefficient of ",item[index[jj]])), dots)
      do.call(plot,arg_cur)
    }
  } else{
    for (jj in 1:length(index)){
      lim=c(  min(record[(index[jj]+len_coef), ]), max(record[index[jj], ]) )
       plot(x$tau, coef[index[jj], ], main=paste0("Estimated coefficient of ", item[index[jj]]), 
            type="n", xlab="taus", ylab= paste0("Estimated coefficient of ",item[index[jj]]), ylim=lim) 
       polygon(c(x$tau,rev(x$tau)),c(record[index[jj], ],rev(record[(index[jj]+len_coef), ])),col = "grey75", border = FALSE)
      dots <- list(...)
      arg_cur = modifyList(list(x$tau, coef[index[jj],]),dots)
      do.call(lines,arg_cur)
    }
  }
  
}

### Plot predicted values against taus ###
plot.predict.WEE.quantile <- function(x, CI= FALSE, level=.95, index=1, ...){
  
  name_tau = colnames(x$predicted_value)
  len = length(name_tau)
  
  used.tau = NULL
  for (ii in 1:len){
    used.tau[ii] = as.numeric(substring(name_tau[ii], 6))
  }
  
  x$newx = as.matrix(x$newx)
  newx = x$newx
  num_x = nrow(x$newx)
  
  if (CI){
    if (length(x) == 2) print("No confidence interval is plotted when zero-boot in the fitting step")
    else
    {
      U = x$predicted_value + qnorm(1-(1-level)/2)*x$StdErr_predicted_value
      L = x$predicted_value - qnorm(1-(1-level)/2)*x$StdErr_predicted_value
    }
  }
  
  if ((CI)&(length(x) != 2)){
    for (ii in 1:length(index)){
      lim=c(  min(L[index[ii], ]), max(U[index[ii], ]) )
      plot(used.tau, x$predicted_value[index[ii], ], main = paste0("Prediction based on covariates: ",paste0(round(newx[index[ii],], digits=2),collapse = " ")), type="n", xlab="tau", ylab="", ylim=lim)
      polygon(c(used.tau,rev(used.tau)),c(as.numeric(U[index[ii],]),as.numeric(rev(L[index[ii],]))),col = "grey75", border = FALSE)
      dots = list(...)
      arg_cur = modifyList(list(used.tau,x$predicted_value[index[ii], ]),dots)
      do.call(lines,arg_cur)
    }
  }
  else {
    for (ii in 1:length(index)) {
      dots = list(...)
      arg_cur = modifyList(list(used.tau, x$predicted_value[index[ii], ], main = paste0("Prediction based on covariates: ",
                          paste0(round(newx[index[ii],], digits=2),collapse = " ")), xlab="tau", ylab=""),dots)
      do.call(plot,arg_cur)
    }}
}
