#
# Delete all data
#
rm(list=ls())

#
# Load the utility routines for model fitting
#
source("modelFitting.R")

#
# Load the package for monte-carlo goodness of fit testing
#
library(XNomial)

#
# Set resType to select the data set to use for model fitting
#
resType <- 2

if(resType == 1)
{ 
  #
  # Human papillomavirus 
  # As reported in Chaturvedi et al. (2011)
  #
  n <- 25
  vPathogens <- 0:8
  vCounts <- c(2933,1409,646,267,102,39,12,2,2)
}

if(resType == 2)
{
  #
  # Barley and cereal yellow dwarf virus 
  # As reported in Seabloom et al. (2009)
  #
  n <- 5
  vPathogens <- 0:5
  vCounts <- c(1570,224,69,17,6,4)
}

if(resType == 3)
{
  #
  # Plasmodium vivax 
  # As reported in Koepfli et al. (2011)
  #
  n <- 57
  vPathogens <- 0:9
  vCounts <- c(1023,404,291,208,118,50,16,5,1,1)
}

if(resType == 4)
{
  #
  # Anther smut
  # As reported by Lopez-Villavicencio et al. (2007)
  #
  n <- 102
  vPathogens <- 0:9
  vCounts <- c(285,74,60,32,14,3,3,2,1,1)
}

if(resType == 5)
{
  #
  # Respiratory virus 
  # As reported in Nickbakhsh et al. (2016)
  #
  n <- 11
  vPathogens <- 0:5
  vCounts <- c(17630,8568,964,105,15,2)
}

if(resType == 6)
{
  # 
  # Borrelia afzelii on voles
  # As reported by Andersson et al. (2013)
  #
  n <- 7
  vPathogens <- 0:6
  vCounts <- c(807,33,26,13,10,11,6)
}

if(resType == 7)
{
  #
  # Pathogens of Ixodes ricinus ticks
  # As reported by Moutailler et al. (2016)
  #
  n <- 37
  vPathogens <- 0:5
  vCounts <- c(147,66,24,18,5,2)
}


#
# Fit the one and two parameter NiSP models, as well as the binomial model
#
oneParamModel <- fit_NiSP_OneParam(n, vPathogens, vCounts)
twoParamModel <- fit_NiSP_TwoParam(n, vPathogens, vCounts, nFits = 25)
binomialModel <- fitBinomialModel(n,vPathogens,vCounts)

#
# Report the three fits
#
sOneParam <- sprintf("Fit.  One parameter NiSP. beta=%.3f, L=%.3f", 
                            trueParams_NiSP_OneParam(oneParamModel$par,n),
                            oneParamModel$value)
print(sOneParam)

sTwoParam <- sprintf("Fit.  Two parameter NiSP. beta=%.3f, gamma=%.3f, L=%.3f", 
                            trueParams_NiSP_TwoParam(twoParamModel$par[1],twoParamModel$par[2],n)[1],
                            trueParams_NiSP_TwoParam(twoParamModel$par[1],twoParamModel$par[2],n)[2],
                            twoParamModel$value)
print(sTwoParam)

sBinomial <- sprintf("Fit.  Binomial. p=%.3f, L=%.3f", 
                     binomialModel$par,
                     binomialModel$value)
print(sBinomial)

#
# Report comparison between the one and two parameter NiSP models
# (via a chi-squared test on the likelihood ratio)
#
diffLL <- twoParamModel$value - oneParamModel$value
modelCompPVal <- 1-pchisq(diffLL,df=1)
if(modelCompPVal < 0.05)
{
  bestNiSP <- twoParamModel
  bestNiSPPVals <- calcEqmNiSP_Transformed_TwoParam(bestNiSP$par[1],bestNiSP$par[2],n)
  sNiSPComp <- sprintf("Comp. Of the NiSP models, the two parameter model is the better fit (X2(1)=%.3f, p=%.3f)", diffLL, modelCompPVal)
}else{
  bestNiSP <- oneParamModel
  bestNiSPPVals <- calcEqmNiSP_Transformed_OneParam(bestNiSP$par[1],n)
  sNiSPComp <- sprintf("Comp. Of the NiSP models, the one parameter model is the better fit (X2(1)=%.3f, p=%.3f)", diffLL, modelCompPVal)
}
print(sNiSPComp)

#
# Report comparison between the best-fitting NiSP model and the binomial model
# (via the AIC)
#

bestNiSPAIC <- 2*length(bestNiSP$par) - 2*bestNiSP$value
binomAIC <- 2*1 - 2*binomialModel$value
deltaAIC <- binomAIC - bestNiSPAIC 
if(deltaAIC>0) # NiSP is better fit
{
  sNiSPBinomialComp = sprintf("Comp. The (best) NiSP model is a better fit than the binomial model (deltaAIC=%.3f)", deltaAIC)
  gofPred <- sum(vCounts)*bestNiSPPVals[vPathogens+1]
}else{ # Binomial is better fit
  sNiSPBinomialComp = sprintf("Comp. The binomial model is a better fit than the (best) NiSP model (deltaAIC=%.3f)", deltaAIC)
  gofPred <- sum(vCounts)*dbinom(vPathogens,n,binomialModel$par)[vPathogens+1]
}
print(sNiSPBinomialComp)

#
# Do the Monte Carlo goodness-of-fit test for the best-fitting model
#
gofCounts <- vCounts
if(length(gofCounts) < (n+1))
{
  # If some possible pathogen counts are not represented in the data then 
  # amalgamate all of these into a single class with a zero count
  gofCounts <- c(vCounts,0)
  gofPred <- c(gofPred,sum(vCounts)-sum(gofPred))
}
v <- xmonte(gofCounts, gofPred, histobins=F, statName="LLR", ntrials=1000000, detail=0)

sMC <- sprintf("GoF.  Monte-Carlo p=%.3f", v$pLLR)
print(sMC)

