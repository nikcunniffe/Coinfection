###########################################################
# Functions to return equilibrium densities of the models #
###########################################################

#
# calcEqmNiSP()
#     Calculates equilibrium densities as predicted by the full NiSP model
#       (i.e. Eqns 74,75 and 77 in the Supplementary Information)
#     Solves the two-parameter model with specific clearance (gamma>0)
#     One parameter NiSP model without specific clearance as used in main text is nested within this model (gamma=0)
#
# Inputs
#     beta = (scaled) infection rate
#     gamma = (scaled) clearance rate
#     n = total number of pathogens
#
# Return value 
#     n+1 dimensional vector (M_0,M_1,...,M_n) at equilibrium
#
calcEqmNiSP <- function(beta,gamma,n)
{
  # calculate derived parameter F (Eqn 72)
  F <- beta - (gamma + 1)
  
  # create space for the coefficient matrix (H) and the rhs vector (b) (Eqn 80)
  H <- matrix(0, nrow = n + 1, ncol = n + 1)
  b <- matrix(0, nrow = n + 1, ncol = 1)
  
  # set element of b (Eqn 79)
  b[1] <- -1
  
  # set elements of H (Eqn 81)
  H[1,1] <- -(n * F + 1)
  H[1,2] <- gamma
  
  for(k in 1:(n-1))
  {
    H[k+1,k] <- (n - k + 1) * F
    H[k+1,k+1] <- -((n - k) * F + k * gamma + 1)
    H[k+1,k+2] <- (k + 1) * gamma
  }
  
  H[n+1,n] <- F
  H[n+1,n+1] <- -(n * gamma + 1)

  # solve to find v (Eqn 80)
  eqmDen <- solve(H, b)
  return(eqmDen)
}

#
# calcEqmNiDP()
#     Calculates equilibrium densities as predicted by the full NiDP model
#         (i.e. Eqns 61 and 62 in the Supplementary Information)
#     Solves the 2n parameter model with specific clearance for each pathogen (gamma_i>0)
#     n parameter model without specific clearance as used in main text is nested within this model (gamma_i=0)
#
# Inputs
#     vBeta = vector of n (scaled) infection rates
#     vGamma = vector of n (scaled) clearance rates
#     n = total number of pathogens
#
# Return value
#   A list with two elements
#     Element 1 = 2^n dimensional vector of equilibrium densities
#     Element 2 = 2^n by n+1 "infection type" matrix
#       - first n columns show presence or absence of each pathogen carried by individuals in class corresponding to this row
#       - column n+1 is the number of distinct pathogens carried by individuals in class corresponding to this row
#
calcEqmNiDP <- function(vBeta,vGamma,n)
{
  # calculate vector of derived parameters, F_i (Eqn 60)
  vF <- vBeta - (vGamma + 1)
  
  # create space for infection type matrix
  vInfType <- 0:((2^n)-1)
  nInfType <- length(vInfType) 
  vInfMatrix <- matrix(0, nrow = nInfType, ncol = n+1)

  # populate the infection type matrix
  for(i in 1:nInfType)
  {
    for(j in 1:n)
    {
      if(bitwAnd(vInfType[i],2^(j-1))>0)
      {
        vInfMatrix[i,j] <- 1
      }
    }
    vInfMatrix[i,n+1] <- sum(vInfMatrix[i,1:n])
  }

  # create space for the coefficient matrix (H) and the rhs vector (b) (Eqn 65)
  H <- matrix(0, nrow = nInfType, ncol = nInfType)
  b <- matrix(0, nrow = nInfType, ncol = 1)
  
  # set element of b (Eqn 64)
  b[1] <- -1
  
  # set elements of H (matrix directly after Eqn 65)
  H[1,1] <- -(sum(vF)+1)

  singleInfIDX <- which(vInfMatrix[,n+1]==1)  
  for(idx in singleInfIDX)
  {
    for(k in 1:n)
    {
      if(vInfMatrix[idx,k])
      {
        H[1,idx] <- vGamma[k]
      }
    }
  }

  for(i in 1:n)
  {
    vThisManyIDX <- which(vInfMatrix[,n+1]==i)
    for(j in vThisManyIDX)
    {
      thisCoeff <- 0
      for(k in 1:n)
      {
        if(vInfMatrix[j,k])
        {
          # Pathogen in set Omega
          newIDX <- j - 2^(k-1)
          H[j,newIDX] <- vF[k]
          thisCoeff <- thisCoeff + vGamma[k]
          
        }else{
          # Pathogen in set Lambda  
          newIDX <- j + 2^(k-1)
          H[j,newIDX] <- vGamma[k]
          thisCoeff <- thisCoeff + vF[k]
        }
      }
      thisCoeff <- -thisCoeff-1
      H[j,j] <- thisCoeff
    }
  }

  # solve to find v (Eqn 65)
  eqmDen <- solve(H, b)

  # return both the equilibrium densities and the infection type matrix
  return(list(eqmDen,vInfMatrix))
}

#################################################
# Functions to fit the two parameter NiSP model #
#################################################

#
# trueParams_NiSP_TwoParam()
#   The two parameter NiSP model is fitted in transformed form with
#   parameters constrained to: i) be positive; ii) correspond to values R0 > 1
#   This function inverts the parameter transformation
#
# Inputs
#   p1 = log(R0-1) [R0 = beta/(gamma+1) for the two parameter NiSP model]
#   p2 = log(gamma)
#   n = total number of pathogens
#
# Return value 
#     2 dimensional vector of meaningful parameters (beta,gamma)
#
trueParams_NiSP_TwoParam <- function(p1,p2,n)
{
  gamma <- exp(p2)
  r0 <- 1+exp(p1)
  beta <- (1+gamma)*r0
  
  return(c(beta,gamma))
}

#
# calcEqmNiSP_Transformed_TwoParam()
#   The two parameter NiSP model is fitted in transformed form with
#   parameters constrained to: i) be positive; ii) correspond to values R0 > 1
#
# Inputs
#   p1 = log(R0-1) [R0 = beta/(gamma+1) for the two parameter NiSP model]
#   p2 = log(gamma)
#   n = total number of pathogens
#
# Return value 
#     n+1 dimensional vector (M_0,M_1,...,M_n) at equilibrium
#
calcEqmNiSP_Transformed_TwoParam <- function(p1,p2,n)
{
  v <- trueParams_NiSP_TwoParam(p1,p2,n)
  calcEqmNiSP(v[1],v[2],n)
}

#
# NiSP_LL_TwoParam()
#   Calculates the log likelihood according to the two parameter NiSP model
#
# Inputs
#   vP = vector of two transformed parameters
#   n = total number of pathogens
#   vPathogens = (data) vector of numbers of distinct pathogens
#   vCounts = (data) vector of counts carrying each number of distinct pathogens
#   
# Return value 
#   Log-likelihood of the data in vPathogens, vCounts given the parameters in vP
#
NiSP_LL_TwoParam <- function(vP,n,vPathogens,vCounts)
{
  allP <- calcEqmNiSP_Transformed_TwoParam(vP[1],vP[2],n)
  thisLL <- 0
  for(i in 1:length(vPathogens))
  {
    thisLL <- thisLL + vCounts[i] * log(allP[vPathogens[i]+1])
  }
  return(thisLL)  
}

#
# fit_NiSP_TwoParam()
#   Fits the two parameter NiSP model in transformed form
#   The model is fitted a number of times with random starting values and best-fitting model returned
#
# Inputs
#   n = total number of pathogens
#   vPathogens = (data) vector of numbers of distinct pathogens
#   vCounts = (data) vector of counts carrying each number of distinct pathogens
#   nFits = number of times to fit the model
#   minP1...maxP2 = (optional) bounds on starting parameter values
#
# Return value 
#   A fitted model (note the parameters in $par are in transformed form)
#
fit_NiSP_TwoParam <- function(n,vPathogens,vCounts,nFits,minP1=-5,maxP1=5,minP2=-5,maxP2=5)
{
  p1 <- runif(1,min=minP1,max=maxP1)
  p2 <- runif(1,min=minP2,max=maxP2)
  currentBestFit <- optim(c(p1,p2), NiSP_LL_TwoParam, control = list(fnscale = -1,maxit = 20000), n = n, vPathogens = vPathogens, vCounts = vCounts)
  nLeft <- nFits - 1
  while(nLeft > 0)
  {
    p1 <- runif(1,min=minP1,max=maxP1)
    p2 <- runif(1,min=minP2,max=maxP2)
    thisFit <- optim(c(p1,p2), NiSP_LL_TwoParam, control = list(fnscale = -1,maxit = 20000), n = n, vPathogens = vPathogens, vCounts = vCounts)
    if(thisFit$value > currentBestFit$value)
    {
      currentBestFit <- thisFit
    }
    nLeft <- nLeft - 1
  }
  return(currentBestFit)
}

#################################################
# Functions to fit the one parameter NiSP model #
#################################################

#
# trueParams_NiSP_OneParam()
#   The one parameter NiSP model is fitted in transformed form to ensure R0>1
#   This function inverts the parameter transformation
#
# Inputs
#   p1 = log(R0-1) [R0 = beta for the one parameter NiSP model]
#   n = total number of pathogens
#
# Return value 
#     beta
#
trueParams_NiSP_OneParam <- function(p1,n)
{
  beta <- 1 + exp(p1)
  return(beta)
}

#
# calcEqmNiSP_Transformed_OneParam()
#   The one parameter NiSP model is fitted in transformed form to ensure R0>1
#
# Inputs
#   p1 = log(R0-1) = log(beta - 1) [since for this model R0 = beta]
#   n = total number of pathogens
#
# Return value 
#     n+1 dimensional vector (M_0,M_1,...,M_n) at equilibrium
#
calcEqmNiSP_Transformed_OneParam <- function(p1,n)
{
  calcEqmNiSP(trueParams_NiSP_OneParam(p1,n),0,n)
}

#
# NiSP_LL_OneParam()
#   Calculates the log likelihood according to the one parameter NiSP model
#
# Inputs
#   p = single transformed parameter
#   n = total number of pathogens
#   vPathogens = (data) vector of numbers of distinct pathogens
#   vCounts = (data) vector of counts carrying each number of distinct pathogens
#   
# Return value 
#   Log-likelihood of the data in vPathogens, vCounts given the parameter p
#
NiSP_LL_OneParam <- function(p,n,vPathogens,vCounts)
{
  allP <- calcEqmNiSP_Transformed_OneParam(p,n)
  thisLL <- 0
  for(i in 1:length(vPathogens))
  {
    thisLL <- thisLL + vCounts[i] * log(allP[vPathogens[i]+1])
  }
  return(thisLL)  
}

#
# fit_NiSP_OneParam()
#   Fits the one parameter NiSP model in transformed form
#   Since this is a one dimensional optimisation problem, and so is well-posed, only need to fit the model once
#
# Inputs
#   n = total number of pathogens
#   vPathogens = (data) vector of numbers of distinct pathogens
#   vCounts = (data) vector of counts carrying each number of distinct pathogens
#   minP1...maxP1 = (optional) bounds on starting parameter value
#
# Return value 
#   A fitted model (note the parameter in $par is in transformed form)
#
fit_NiSP_OneParam <- function(n,vPathogens,vCounts,minP1=-5,maxP1=5)
{
  p1 <- runif(1,min=minP1,max=maxP1)
  currentBestFit <- optim(p1, NiSP_LL_OneParam, method="Brent", lower = minP1, upper = maxP1, control = list(fnscale = -1,maxit = 20000), n = n, vPathogens = vPathogens, vCounts = vCounts)
  return(currentBestFit)
}

#######################################
# Functions to fit the binomial model #
#######################################

#
# binomial_LL()
#   Calculates the log likelihood according to the binomial model
#
# Inputs
#   p = probability of infection
#   n = total number of pathogens
#   vPathogens = (data) vector of numbers of distinct pathogens
#   vCounts = (data) vector of counts carrying each number of distinct pathogens
#   
# Return value 
#   Log-likelihood of the data in vPathogens, vCounts given the parameter p
#
binomial_LL <- function(p,n,vPathogens,vCounts)
{
  thisLL <- 0
  for(i in 1:length(vPathogens))
  {
    thisLL <- thisLL + vCounts[i] * (lchoose(n, vPathogens[i]) + vPathogens[i] * log(p) + (n - vPathogens[i]) * log(1 - p))
  }
  return(thisLL)  
}

#
# fitBinomialModel()
#   Fits the binomial model
#
# Inputs
#   n = total number of pathogens
#   vPathogens = (data) vector of numbers of distinct pathogens
#   vCounts = (data) vector of counts carrying each number of distinct pathogens
#
# Return value 
#   A fitted binomial model
#
fitBinomialModel <- function(n,vPathogens,vCounts)
{
  return(optim(runif(1), binomial_LL, method = "Brent", lower = 0, upper = 1, control = list(fnscale = -1), n = n, vPathogens = vPathogens, vCounts = vCounts))
}

############################################################################
# Functions to fit the n parameter NiDP model (without specific clearance) #
############################################################################

#
# trueParams_NiDP_nParam()
#   The n parameter NiDP model is fitted in transformed form to ensure R0>1
#   This function inverts the parameter transformation
#
# Inputs
#   vp1 = vector of values of log(R0_i-1) [R0_i = beta_i for the n parameter NiDP model]
#   n = total number of pathogens
#
# Return value 
#   vector of n values of beta
#
trueParams_NiDP_nParam <- function(vp1,n)
{
  vBeta <- 1 + exp(vp1)
  return(vBeta)
}

#
# calcEqmNiDP_nParam_Transformed()
#   The n parameter NiDP model is fitted in transformed form to ensure R0>1
#
# Inputs
#   vp1 = vector of values of log(R0_i-1) [R0_i = beta_i for the n parameter NiDP model]
#   n = total number of pathogens
#
# Return value 
#   Equilibrium values of the model in the same list format as calcEqmNiDP()
#
calcEqmNiDP_nParam_Transformed <- function(vp1,n)
{
  calcEqmNiDP(trueParams_NiDP_nParam(vp1,n),rep(0,n),n)
}

#
# NiDP_nParam_LL()
#   Calculates the log likelihood according to the n parameter NiDP model
#
# Inputs
#   vP = vector of values of log(R0_i-1)
#   n = total number of pathogens
#   mPathCodes = (data) matrix of pathogen codes
#   vCounts = (data) vector of counts of each pathogen code
#   nCounts = number of distinct pathogen codes
#   
# Return value 
#   Log-likelihood of the data in vPathogens, vCounts given the parameters in vP
#
NiDP_nParam_LL <- function(vP,n,mPathCodes,vCounts,nCounts)
{
  modList <- calcEqmNiDP_nParam_Transformed(vP,n)
  allP <- modList[[1]]
  returnedCodes <- modList[[2]][,1:n]
  thisLL <- 0
  for(i in 1:nCounts)
  {
    thisRow <- row.match(mPathCodes[i,],returnedCodes)
    if(!is.na(thisRow))
    {
      thisLL <- thisLL + vCounts[i] * log(allP[thisRow])
    }else{
      print("Error!")
    }
  }
  return(thisLL)  
}

#
# fit_NiDP_nParam()
#   Fits the n parameter NiDP model in transformed form 
#   This is the model without pathogen specific clearance
#   The parameters are therefore n values of beta_i (in transformed form)
#
# Inputs
#   n = number of pathogens
#   mPathCodes = (data) matrix of pathogen codes
#   vCounts = (data) vector of counts of each pathogen code
#   nCounts = number of distinct pathogen codes
#   minP1...maxP1 = (optional) bounds on starting parameter value
#
# Return value 
#   A fitted model (note the parameter in $par is in transformed form)
#
fit_NiDP_nParam <- function(n,mPathCodes,vCounts,nCounts,nFits,minP1=-10,maxP1=10)
{
  vp1 <- runif(n,min=minP1,max=maxP1)
  currentBestFit <- optim(vp1, NiDP_nParam_LL, control = list(fnscale = -1,maxit = 20000), n = n, mPathCodes = mPathCodes, vCounts = vCounts, nCounts = nCounts)
  nLeft <- nFits - 1
  while(nLeft > 0)
  {
    vp1 <- runif(n,min=minP1,max=maxP1)
    thisFit <- optim(vp1, NiDP_nParam_LL, control = list(fnscale = -1,maxit = 20000), n = n, mPathCodes = mPathCodes, vCounts = vCounts, nCounts = nCounts)
    if(thisFit$value > currentBestFit$value)
    {
      currentBestFit <- thisFit
    }
    nLeft <- nLeft - 1
  }
  return(currentBestFit)
}

#
# disp_NiDP_ZeroGamma()
#   Displays a fitted n parameter NiDP model
#
disp_NiDP_ZeroGamma <- function(n, vMod)
{
  cat("\tLL:",vMod$value,"\n")

  vBeta <- trueParams_NiDP_nParam(vMod$par,n)

  cat("\tR0:  ",vBeta,"\n")
  cat("\tEqm: ",calcEqmNiDP(vBeta,rep(0,n),n)[[1]],"\n")
  
  retList <- calcEqmNiDP(vBeta,rep(0,n),n)
  retList[[3]] <- vMod$value
  retList[[4]] <- vBeta
  return(retList)
}

##########################################
# Functions to fit the multinomial model #
##########################################


#
# trueParams_Multinomial()
#   The multinomial model is fitted in transformed form to ensure 0<p_i<1
#   This function inverts the parameter transformation
#
# Inputs
#   vp1 = vector of values of the logit of the probabilities
#
# Return value 
#   vector of n probabilities
#
trueParams_Multinomial <- function(vp1)
{
  vPInfect <- 1/(1 + exp(-vp1))
  return(vPInfect)
}

#
# fit_Multinomial()
#   Fits the n parameter multinomial in transformed form 
#   The parameters are n probabilities of infection (in transformed form)
#
# Inputs
#   n = number of pathogens
#   mPathCodes = (data) matrix of pathogen codes
#   vCounts = (data) vector of counts of each pathogen code
#   nCounts = number of distinct pathogen codes
#   minLogit...maxLogit = (optional) bounds on starting parameter value
#
# Return value 
#   A fitted model (note the parameter in $par is in transformed form)
#
fit_Multinomial <- function(n,mPathCodes,vCounts,nCounts,nFits,minLogit=-10,maxLogit=10)
{
  vp <- runif(n,min=minLogit,max=maxLogit)
  currentBestFit <- optim(vp, Multinomial_LL_Full, control = list(fnscale = -1,maxit = 20000), n = n, mPathCodes = mPathCodes, vCounts = vCounts, nCounts = nCounts)
  nLeft <- nFits - 1
  while(nLeft > 0)
  {
    vp <- runif(n,min=minLogit,max=maxLogit)
    thisFit <- optim(vp, Multinomial_LL_Full, control = list(fnscale = -1,maxit = 20000), n = n, mPathCodes = mPathCodes, vCounts = vCounts, nCounts = nCounts)
    if(thisFit$value > currentBestFit$value)
    {
      currentBestFit <- thisFit
    }
    nLeft <- nLeft - 1
  }
  return(currentBestFit)
}

#
# Multinomial_LL_Full()
#   Calculates the log likelihood according to the multinomial model
#
# Inputs
#   vP = vector of values of logits of infection probabilities
#   n = total number of pathogens
#   mPathCodes = (data) matrix of pathogen codes
#   vCounts = (data) vector of counts of each pathogen code
#   nCounts = number of distinct pathogen codes
#   
# Return value 
#   Log-likelihood of the data in vPathogens, vCounts given the parameters in vP
#
Multinomial_LL_Full <- function(vP,n,mPathCodes,vCounts,nCounts)
{
  modList <- calcEqm_MultinomialModel(vP,n)
  allP <- modList[[1]]
  returnedCodes <- modList[[2]][,1:n]
  thisLL <- 0
  for(i in 1:nCounts)
  {
    thisRow <- row.match(mPathCodes[i,],returnedCodes)
    if(!is.na(thisRow))
    {
      thisLL <- thisLL + vCounts[i] * log(allP[thisRow])
    }else{
      print("Error!")
    }
  }
  return(thisLL)  
}

#
# calcEqm_MultinomialModel()
#   The multinomial model is fitted in transformed form to ensure 0<p_i<1
#
# Inputs
#   vParams = vector of values of logits of proababilities
#   n = total number of pathogens
#
# Return value 
#   Equilibrium values of the model in the same list format as calcEqmMultinomial()
#
calcEqm_MultinomialModel <- function(vParams,n)
{
  calcEqmMultinomial(trueParams_Multinomial(vParams),n)
}

#
# calcEqmMultinomial()
#     Calculates equilibrium densities as predicted by the multinomial model
#
# Inputs
#     vPInfect = vector of n infection probabilities
#     n = total number of pathogens
#
# Return value
#   A list with two elements
#     Element 1 = 2^n dimensional vector of equilibrium densities
#     Element 2 = 2^n by n+1 "infection type" matrix
#       - first n columns show presence or absence of each pathogen carried by individuals in class corresponding to this row
#       - column n+1 is the number of distinct pathogens carried by individuals in class corresponding to this row
#
calcEqmMultinomial <- function(vPInfect,n)
{
  # set aside space for matrix specifying what pathogens correspond to which position
  vInfType <- 0:((2^n)-1)
  nInfType <- length(vInfType) 
  vInfMatrix <- matrix(0, nrow = nInfType, ncol = n + 1)
  
  # fill in entries of the what pathogen in which position matrix
  for(i in 1:nInfType)
  {
    for(j in 1:n)
    {
      if(bitwAnd(vInfType[i],2^(j-1))>0)
      {
        vInfMatrix[i,j] <- 1
      }
    }
    vInfMatrix[i,n+1] <- sum(vInfMatrix[i,1:n])
  }
  eqmDen <- numeric(2^n)
  for(i in 1:(2^n))
  {
    thisP <- 1
    for(j in 1:n)
    {
      if(vInfMatrix[i,j]==1)
      {
        thisP <- thisP * vPInfect[j]
      }else{
        thisP <- thisP * (1-vPInfect[j])
      }
    }
    eqmDen[i] <- thisP
  }
  return(list(eqmDen,vInfMatrix))
}

#
# disp_Multinomial()
#   Displays a fitted multinomial model
#
disp_Multinomial <- function(n, vMod)
{
  cat("\tLL:",vMod$value,"\n")
  
  vProb <- trueParams_Multinomial(vMod$par)
  
  cat("\tp:", vProb,"\n")
  cat("\tEqm:",calcEqm_MultinomialModel(vMod$par,n)[[1]],"\n")
  
  retList <- calcEqm_MultinomialModel(vMod$par,n)
  retList[[3]] <- vMod$value
  retList[[4]] <- vProb
  
  return(retList)
}