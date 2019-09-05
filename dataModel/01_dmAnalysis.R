# File: 01_dmAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: differential methylation analysis using binomial model for probe list
# Date: 5/9/2019


source('dataModel/header.R')
setwd('dataModel/')

###################################################### 
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegression4RandomEffects.stan')
dfData = read.csv('dataExternal/NTD methylation_99probes and covariates for Umar.csv', header=T)

###################################### function to fit model
modelFunction = function(lData){
  fs = tryCatch(sampling(stanDso, data=lData, iter=2000, chains=4, pars=c(#'intercept',
                                                                          'nCoefFactor1',
                                                                          #'nCoefFactor2', 
                                                                          #'nCoefFactor3',
                                                                          #'nCoefFactor4',
                                                                          'sigmaFactor1',
                                                                          'sigmaFactor2',
                                                                          'sigmaFactor3',
                                                                          'sigmaFactor4'),
                         cores=4), error=function(e) NULL)
  return(fs)
}

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

extractResult = function(fit.stan, base='diet1', deflection='diet2', lev=levels(dfData$diet)){
  if (is.null(fit.stan)) return(NULL)
  mCoef = extract(fit.stan)$nCoefFactor1
  colnames(mCoef) = lev
  dif = getDifference(mCoef[,deflection], mCoef[,base])
  r = data.frame(coef.base=mean(mCoef[,base]), 
                 coef.deflection=mean(mCoef[,deflection]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
}

dfData.sub = dfData[,-c(1:5)]
cn = colnames(dfData.sub)
cn = gsub('_methylated|_unmethylated', '', cn)
cn = factor(cn)
lIndex = split(1:ncol(dfData.sub), cn)

## fit model for each corresponding positions in the data frame
lModel = lapply(lIndex, function(x){
  lStanData = list(Ntotal=nrow(dfData),
                   Nlevels1 = nlevels(dfData$diet),
                   NfactorMap1 = as.numeric(dfData$diet),
                   Nlevels2 = nlevels(dfData$sex),
                   NfactorMap2 = as.numeric(dfData$sex),
                   Nlevels3 = nlevels(dfData$genotype),
                   NfactorMap3 = as.numeric(dfData$genotype),
                   Nlevels4 = nlevels(dfData$batch),
                   NfactorMap4 = as.numeric(dfData$batch),
                   y=dfData.sub[,x[1]],
                   Ntrials=dfData.sub[,x[1]] + dfData.sub[,x[2]])
  return(modelFunction(lStanData))
})

##################################### extract results

lResults = lapply(lModel, extractResult)
dfResults = do.call(rbind, lResults)
write.csv(round(dfResults,3), file='temp/diet2vsdiet1.csv')

lResults = lapply(lModel, extractResult, deflection='diet 3')
dfResults = do.call(rbind, lResults)
write.csv(round(dfResults,3), file='temp/diet3vsdiet1.csv')

