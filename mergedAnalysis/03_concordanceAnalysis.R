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
metadata(oGRbis) = list()
gc(reset = T)
# ## select a smaller sample to look at the data
# i = sample(1:nrow(mMethylated), 1000, replace = F)
# 
# mMethylated.sam = mMethylated[i,]
# ## add a one to avoid zero divisions
# mTotal.sam = mTotal[i,] + 1
# 
# mProp.sam = matrix(NA, nrow = nrow(mMethylated.sam), ncol = ncol(mMethylated.sam))
# colnames(mProp.sam) = colnames(mMethylated.sam)
# 
# for (i in 1:ncol(mProp.sam)){
#   mProp.sam[,i] = mMethylated.sam[,i] / mTotal.sam[,i]
# }

## calculate model fitted covariance
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='twoLevelBinomialHierarchicalModel.stan')

# ## set up the data 
# mProp.mod = matrix(NA, nrow = nrow(mMethylated.sam), ncol = ncol(mMethylated.sam))
# colnames(mProp.mod) = colnames(mMethylated.sam)
# 
# dfData = data.frame(methylated=mMethylated.sam[1,],
#                     total=mTotal.sam[1,],
#                     grouping=factor(factor(dfSamples$group1):factor(dfSamples$group2)))
# nrow(dfData)
# 
# for (i in 1:nrow(mProp.mod)){
#   ## set up stan data
#   lStanData = list(Ntotal=22, NgroupsLvl1=6, 
#                    NgroupsLvl2Map=as.numeric(dfData$grouping),
#                    y=mMethylated.sam[i,],
#                    N=mTotal.sam[i,])
#   
#   fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2,
#                       cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
#   #print(fit.stan, digits=3)
#   mProp.mod[i,] = colMeans(extract(fit.stan)$theta)
#   cat(i, 'done\n')
# }

### load the mouse genome BSGenome object
library(BSgenome.Mmusculus.UCSC.mm10)
seqlengths(BSgenome.Mmusculus.UCSC.mm10)
oGRgenome = GRangesForBSGenome('mm10')
iBins = 1000
## choose the sequence names over which to calculate coverage
cvSeqnames = c(paste0('chr', 1:19))#, 'chrX', 'chrY')

# lChromosomes = vector(mode='list', length=length(cvSeqnames))
# names(lChromosomes) = cvSeqnames
# 
# lChromosomes.raw = vector(mode='list', length=length(cvSeqnames))
# names(lChromosomes.raw) = cvSeqnames
gr = GRanges()

## for each chromosome of interest create an overlapping binned range
for (outer in 1:length(cvSeqnames)){
  ###### choose chr
  b = f_bin_vector(start(oGRgenome[seqnames(oGRgenome) == cvSeqnames[outer]]), 
                   end(oGRgenome[seqnames(oGRgenome) == cvSeqnames[outer]]), bins = iBins)
  
  gr = append(gr, GRanges(cvSeqnames[outer], ranges = IRanges(b$start, b$end), strand='*'))
}

## subset these bins baesd on how many overlaps it has
f = countOverlaps(gr, oGRbis)
gr = gr[f >= 350]
table(f >= 350)
## randomly select 1000 regions
# i = sample(1:length(gr), 1000, replace = F)
# gr = gr[i]

## for each of the ranged bins count the number of methylation signals
mMethylated.bin = matrix(NA, nrow = length(gr), ncol = nrow(dfSamples), dimnames = list(NULL, dfSamples$id))
mTotal.bin = matrix(NA, nrow = length(gr), ncol = nrow(dfSamples), dimnames = list(NULL, dfSamples$id))
c = vector('logical', length=length(oGRbis))

# ## fill this matrix of counts by going over the bisulphite data
# l = lapply(seq_along(1:length(gr)), function(i){
#   c = overlapsAny(oGRbis, gr[i])
#   mMethylated.bin[i,] = colSums(mMethylated[oGRbis[c]$mapper,])
#   mTotal.bin[i,] = colSums(mTotal[oGRbis[c]$mapper,]) + 1
# })

## get the methylated signal over each bin
for (i in 1:length(gr)){
  c = overlapsAny(oGRbis, gr[i])
  mMethylated.bin[i,] = colSums(mMethylated[oGRbis[c]$mapper,])
  mTotal.bin[i,] = colSums(mTotal[oGRbis[c]$mapper,]) + 1
  cat(i, ' out of ', length(gr), '\n')
}

## calculate the raw values
## raw proportions
mProp.raw = matrix(NA, nrow = nrow(mMethylated.bin), ncol = ncol(mMethylated.bin))
colnames(mProp.raw) = colnames(mMethylated.bin)

for (i in 1:ncol(mProp.raw)){
  mProp.raw[,i] = mMethylated.bin[,i] / mTotal.bin[,i]
}

############ fit model to calculate probabilities under binomial sampling model
mProp.mod = matrix(NA, nrow = nrow(mMethylated.bin), ncol = ncol(mMethylated.bin))
colnames(mProp.mod) = colnames(mMethylated.bin)


fGrouping = factor(factor(dfSamples$group1):factor(dfSamples$group2))
nlevels(fGrouping)

## save the first level parameters as well
mProp.modLevel1 = matrix(NA, nrow = nrow(mMethylated.bin), ncol = nlevels(fGrouping))
colnames(mProp.modLevel1) = as.character(levels(fGrouping))

for (i in 1:nrow(mProp.mod)){
  ## set up stan data
  lStanData = list(Ntotal=22, NgroupsLvl1=6, 
                   NgroupsLvl2Map=as.numeric(fGrouping),
                   y=mMethylated.bin[i,],
                   N=mTotal.bin[i,])
  
  fit.stan = sampling(stanDso, data=lStanData, iter=500, chains=2,
                      cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
  #print(fit.stan, digits=3)
  mProp.mod[i,] = colMeans(extract(fit.stan)$theta)
  mProp.modLevel1[i,] = colMeans(extract(fit.stan)$omega1)
  
  cat(i, ' out of ', length(gr), '\n')
}

####### perform a conjugate analysis at the group level
## functions to use
getalphabeta = function(m, v){
  al.be = (m * (1-m) / v) - 1
  al = al.be * m
  be = al.be * (1-m)
  return(c(alpha=al, beta=be))
}

getFittedTheta = function(param){
  if (param[1] <= 0) param[1] = 0.5
  if (param[2] <= 0) param[2] = 0.5
  iAlpha = param[1];
  iBeta = param[2];
  iSuc = param[3];
  iFail = param[4];
  return(mean(rbeta(1000, iSuc+iAlpha, iFail+iBeta)))
}

## the first level parameters using conjugate model
mProp.modConj = matrix(NA, nrow = nrow(mMethylated.bin), ncol = nlevels(fGrouping))
colnames(mProp.modConj) = as.character(levels(fGrouping))

for (i in 1:nrow(mProp.modConj)){
  ## get the hyperparameters using all the data
  suc = tapply(mMethylated.bin[i,], fGrouping, sum)
  tot = tapply(mTotal.bin[i,], fGrouping, sum)
  fail = tot - suc
  ivThetas.data = suc/tot
  ## get hyperparameters using population data
  l = getalphabeta(mean(ivThetas.data), var(ivThetas.data))
  # put data in matrix form to cycle the function over each row
  m = cbind(l['alpha'], l['beta'], suc, fail)
  mProp.modConj[i,] = apply(m, 1, getFittedTheta)
  
  cat(i, ' out of ', length(gr), '\n')
}

lProportions = list(model=mProp.mod, model.level1=mProp.modLevel1, raw=mProp.raw, model.conj = mProp.modConj)

save(lProportions, file='temp/lProportions.binsize.1000.cov350.rds')

### some diagnostic plots
library(downloader)
library(car)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

colnames(mProp.mod) = dfSamples$title
colnames(mProp.raw) = dfSamples$title

oDiag.mod = CDiagnosticPlots(mProp.mod, 'Model - bin 1000')
oDiag.raw = CDiagnosticPlots(mProp.raw, 'Raw - bin 1000')
l = CDiagnosticPlotsGetParameters(oDiag.raw)

fBatch = fGrouping

l$PCA.jitter = F; l$HC.jitter = F;
#l$PCA.scaleSubjects = F;
l$PCA.scaleVariables = F;
#l$HC.scaleSubjects = F;
l$HC.scaleVaribles = F;
oDiag.mod = CDiagnosticPlotsSetParameters(oDiag.mod, l)
oDiag.raw = CDiagnosticPlotsSetParameters(oDiag.raw, l)

par(mfrow=c(2,1))
plot.mean.summary(oDiag.mod, fBatch)
plot.mean.summary(oDiag.raw, fBatch)

plot.sigma.summary(oDiag.mod, fBatch)
plot.sigma.summary(oDiag.raw, fBatch)

boxplot.median.summary(oDiag.mod, fBatch)
boxplot.median.summary(oDiag.raw, fBatch)

m = oDiag.mod@lData$mean[,'m']
s = oDiag.mod@lData$sigma[,'m']

par(mfrow=c(1,2))
col.p = rainbow(nlevels(fBatch))
col = col.p[as.numeric(fBatch)]
plot(m, s, pch=20, col=col, main='Model - bin 1000', xlab='mean', ylab='sigma')
text(m, s, names(m), pos=1, cex=0.6)


m = oDiag.raw@lData$mean[,'m']
s = oDiag.raw@lData$sigma[,'m']

col.p = rainbow(nlevels(fBatch))
col = col.p[as.numeric(fBatch)]
plot(m, s, pch=20, col=col, main='Raw - bin 1000', xlab='mean', ylab='sigma')
text(m, s, names(m), pos=1, cex=0.6)

par(mfrow=c(1,1))
plot.PCA(oDiag.mod, fBatch)
plot.PCA(oDiag.raw, fBatch)

plot.dendogram(oDiag.mod, fBatch, labels_cex = 0.8)
plot.dendogram(oDiag.raw, fBatch, labels_cex = 0.8)

oDiag.lvl1 = CDiagnosticPlots(mProp.modLevel1, 'Group Level bin 1000')
oDiag.lvl1 = CDiagnosticPlotsSetParameters(oDiag.lvl1, l)
fBatch.2 = factor(levels(fGrouping))
plot.mean.summary(oDiag.lvl1, fBatch.2)
plot.sigma.summary(oDiag.lvl1, fBatch.2)
plot.PCA(oDiag.lvl1, fBatch.2)
plot.dendogram(oDiag.lvl1, fBatch.2, labels_cex = 0.8)

oDiag.lvl1 = CDiagnosticPlots(mProp.modConj, 'Conjugate bin 1000 Cov 600')
oDiag.lvl1 = CDiagnosticPlotsSetParameters(oDiag.lvl1, l)
fBatch.2 = factor(levels(fGrouping))
plot.mean.summary(oDiag.lvl1, fBatch.2)
plot.sigma.summary(oDiag.lvl1, fBatch.2)
plot.PCA(oDiag.lvl1, fBatch.2)
plot.dendogram(oDiag.lvl1, fBatch.2, labels_cex = 0.8)




