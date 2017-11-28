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
iBins = 500
## choose the sequence names over which to calculate coverage
cvSeqnames = c(paste0('chr', 1:19))#, 'chrX', 'chrY')

lChromosomes = vector(mode='list', length=length(cvSeqnames))
names(lChromosomes) = cvSeqnames

lChromosomes.raw = vector(mode='list', length=length(cvSeqnames))
names(lChromosomes.raw) = cvSeqnames

## for each chromosome of interest create an overlapping binned range
for (outer in 1:length(lChromosomes)){
  ###### choose chr
  b = f_bin_vector(start(oGRgenome[seqnames(oGRgenome) == cvSeqnames[outer]]), 
                   end(oGRgenome[seqnames(oGRgenome) == cvSeqnames[outer]]), bins = iBins)
  
  gr = GRanges(cvSeqnames[outer], ranges = IRanges(b$start, b$end), strand='*')
  # ## take a sample from this ranges
  # sam = sample(1:length(gr), 100, replace = F)
  # gr = gr[sam]
  ## subset the bis data 
  oGRbis.sub = subsetByOverlaps(oGRbis, gr)
  
  ## for each of the ranged bins count the number of methylation signals
  mMethylated.bin = matrix(NA, nrow = iBins, ncol = nrow(dfSamples), dimnames = list(NULL, dfSamples$id))
  mTotal.bin = matrix(NA, nrow = iBins, ncol = nrow(dfSamples), dimnames = list(NULL, dfSamples$id))
  
  for (i in 1:length(gr)){
    c = subsetByOverlaps(oGRbis.sub, gr[i])
    if (length(c) < 2) next;
    ## use this selected set of ranges overlapping with the bin
    ## to search into the corresponding matrix row to find number of overlaps
    mMethylated.bin[i,] = colSums(mMethylated[c$mapper,])
    mTotal.bin[i,] = colSums(mTotal[c$mapper,]) + 1
  }
  
  mMethylated.bin = na.omit(mMethylated.bin)
  mTotal.bin = na.omit(mTotal.bin)
  
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
  lChromosomes[[cvSeqnames[outer]]] = mProp.binMod
  lChromosomes.raw[[cvSeqnames[outer]]] = mProp.bin
  cat(outer, 'done ', cvSeqnames[outer], '\n')
}

save(lChromosomes, file='temp/lChromosome.rds')
save(lChromosomes.raw, file='temp/lChromosome.raw.rds')

### some diagnostic plots
library(downloader)
library(car)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.mod = CDiagnosticPlots(logit(do.call(rbind, lChromosomes)), 'Model Estimated')
oDiag.raw = CDiagnosticPlots(logit(do.call(rbind, lChromosomes.raw)), 'Raw Proportions')
l = CDiagnosticPlotsGetParameters(oDiag.mod)

fBatch = factor(dfSamples$group3)

l$PCA.jitter = F; l$HC.jitter = F;
l$PCA.scaleSubjects = F;
l$PCA.scaleVariables = F;
l$HC.scaleSubjects = F;
l$HC.scaleVaribles = F;
oDiag.mod = CDiagnosticPlotsSetParameters(oDiag.mod, l)
oDiag.raw = CDiagnosticPlotsSetParameters(oDiag.raw, l)

par(mfrow=c(1,2))
plot.mean.summary(oDiag.mod, fBatch)
plot.mean.summary(oDiag.raw, fBatch)

plot.sigma.summary(oDiag.mod, fBatch)
plot.sigma.summary(oDiag.raw, fBatch)

boxplot.median.summary(oDiag.mod, fBatch)
boxplot.median.summary(oDiag.raw, fBatch)

plot.PCA(oDiag.mod, fBatch)
plot.PCA(oDiag.raw, fBatch)

plot.dendogram(oDiag.mod, fBatch, labels_cex = 0.8)
plot.dendogram(oDiag.raw, fBatch, labels_cex = 0.8)

par(mfrow=c(1,1))
plot.dendogram(oDiag.mod, dfData$grouping, labels_cex = 0.8)
plot.dendogram(oDiag.raw, dfData$grouping, labels_cex = 0.8)






