# Name: 03_concordanceAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: various clustering and correlation plots for samples based on methylation profiles
# Date: 23/11/2017


source('header.R')

load('results/GRangesMethylatedPositionsWithMetaData.rds')

oGRbis
names(metadata(oGRbis))
dfSamples = metadata(oGRbis)$samples
mMethylated = metadata(oGRbis)$methylated
mTotal = metadata(oGRbis)$total

head(dfSamples)
head(mMethylated)
head(mTotal)
dim(mMethylated)
length(oGRbis)
oGRbis$mapper = 1:length(oGRbis)

## select a smaller sample to look at the data
i = sample(1:nrow(mMethylated), 1000, replace = F)

mMethylated.sam = mMethylated[i,]
## add a one to avoid zero divisions
mTotal.sam = mTotal[i,] + 1

mProp.sam = matrix(NA, nrow = nrow(mMethylated.sam), ncol = ncol(mMethylated.sam))
colnames(mProp.sam) = colnames(mMethylated.sam)

for (i in 1:ncol(mProp.sam)){
  mProp.sam[,i] = mMethylated.sam[,i] / mTotal.sam[,i]
}

## calculate model fitted covariance
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='twoLevelBinomialHierarchicalModel.stan')

## set up the data 
mProp.mod = matrix(NA, nrow = nrow(mMethylated.sam), ncol = ncol(mMethylated.sam))
colnames(mProp.mod) = colnames(mMethylated.sam)

dfData = data.frame(methylated=mMethylated.sam[1,],
                    total=mTotal.sam[1,],
                    grouping=factor(factor(dfSamples$group1):factor(dfSamples$group2)))
nrow(dfData)

for (i in 1:nrow(mProp.mod)){
  ## set up stan data
  lStanData = list(Ntotal=22, NgroupsLvl1=6, 
                   NgroupsLvl2Map=as.numeric(dfData$grouping),
                   y=mMethylated.sam[i,],
                   N=mTotal.sam[i,])
  
  fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2,
                      cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
  #print(fit.stan, digits=3)
  mProp.mod[i,] = colMeans(extract(fit.stan)$theta)
  cat(i, 'done\n')
}

### load the mouse genome BSGenome object
library(BSgenome.Mmusculus.UCSC.mm10)
seqlengths(BSgenome.Mmusculus.UCSC.mm10)
oGRgenome = GRangesForBSGenome('mm10')
iBins = 20

## for each chromosome of interest create an overlapping binned range
###### chr1
b = f_bin_vector(start(oGRgenome[seqnames(oGRgenome) == 'chr1']), end(oGRgenome[seqnames(oGRgenome) == 'chr1']), bins = iBins)

gr = GRanges('chr1', ranges = IRanges(b$start, b$end), strand='*')
## subset the bis data 
oGRbis.sub = subsetByOverlaps(oGRbis, gr)

## for each of the ranged bins count the number of methylation signals
mMethylated.bin = matrix(NA, nrow = iBins, ncol = nrow(dfSamples), dimnames = list(NULL, dfSamples$id))
mTotal.bin = matrix(NA, nrow = iBins, ncol = nrow(dfSamples), dimnames = list(NULL, dfSamples$id))

for (i in 1:length(gr)){
  c = subsetByOverlaps(oGRbis.sub, gr[i])
  ## use this to search into the matrix row
  mMethylated.bin[i,] = colSums(mMethylated[c$mapper,])
  mTotal.bin[i,] = colSums(mTotal[c$mapper,]) + 1
}

## raw proportions
mProp.bin = matrix(NA, nrow = nrow(mMethylated.bin), ncol = ncol(mMethylated.bin))
colnames(mProp.bin) = colnames(mMethylated.bin)

for (i in 1:ncol(mProp.bin)){
  mProp.bin[,i] = mMethylated.bin[,i] / mTotal.bin[,i]
}

############ fit model to calculate probabilities under binomial sampling model
mProp.binMod = matrix(NA, nrow = nrow(mMethylated.bin), ncol = ncol(mMethylated.bin))
colnames(mProp.binMod) = colnames(mMethylated.bin)

dfData = data.frame(methylated=mMethylated.bin[1,],
                    total=mTotal.bin[1,],
                    grouping=factor(factor(dfSamples$group1):factor(dfSamples$group2)))
nrow(dfData)

for (i in 1:nrow(mProp.binMod)){
  ## set up stan data
  lStanData = list(Ntotal=22, NgroupsLvl1=6, 
                   NgroupsLvl2Map=as.numeric(dfData$grouping),
                   y=mMethylated.bin[i,],
                   N=mTotal.bin[i,])
  
  fit.stan = sampling(stanDso, data=lStanData, iter=500, chains=2,
                      cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
  #print(fit.stan, digits=3)
  mProp.binMod[i,] = colMeans(extract(fit.stan)$theta)
  cat(i, 'done\n')
}



