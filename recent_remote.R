# File: recent_remote.R
# Desc: Data set from tb recent and remote
# Auth: u.niazi@imperial.ac.uk
# Date: 01/06/2016

library(annotate)
library(limma)
library(org.Hs.eg.db)
library(lumi)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CGraphClust/master/CGraphClust.R'
download(url, 'CGraphClust.R')

# load the required packages
source('CGraphClust.R')
# delete the file after source
unlink('CGraphClust.R')

## internal function
# Function: f_lGetPCAClusterCount
# Desc: Takes first two components of the PCA and counts possible clusters. The function does this by 
#       binning the vector of data into X bins, assigning a class label to each bin, counting how many
#       observations in each bin, any total number of bins with at least one observations. This is calculated
#       for both the components of the pca matrix, and the max number of bins with at least one observation, in
#       first or second dimension is reported along with a data.frame with cluster labels.
# Args: pr.out = principal component object returned by prcomp function
# Rets: returns list with 2 elements: 
#       1 - cluster.count = possible number of clusters in the data
#       2 - cluster.label = data.frame with cluster labels
f_lGetPCAClusterCount = function(pr.out){
  # how many clusters in data, using first 2 components
  x1 = pr.out$x[,1]
  x2 = pr.out$x[,2]
  # bin the data from the 2 components
  h1 = hist(x1, plot=F)
  # give a class label to each bin
  c1 = cut(x1, h1$breaks, labels = 1:(length(h1$mids)))
  h2 = hist(x2, plot=F)
  c2 = cut(x2, h2$breaks, labels = 1:(length(h2$mids)))
  # labels for vectors and the class labels
  dfClust = data.frame(lab=names(x1), c1, c2)
  # get contingency table
  mClust = as.matrix(table(c1 = dfClust$c1, c2 = dfClust$c2))
  # count the max of row and col sums that are not zero
  ir = length(which(rowSums(mClust) != 0))
  ic = length(which(colSums(mClust) != 0))
  iClust.count = ifelse(ir > ic, ir, ic)
  lRet = list(cluster.count=iClust.count, cluster.label=dfClust)
  return(lRet)
}

f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}
### end internal functions

# global variables
p.old = par()


## load the data
dfSamples = read.csv('Data_external/Recent_remote/samples.csv', header=T)
dfData = read.csv('Data_external/Recent_remote/rr_matrix_norm.txt', sep = ' ', header=T)

# data sorting
cn = colnames(dfData)#[6:(ncol(dfData))]
cn = gsub('X(.+)', replacement = '\\1', cn, perl = T)

colnames(dfData) = cn

mData = as.matrix(dfData)
mData = t(mData)

## put the two in the same order i.e. sample annotations and data
table(rownames(mData) %in% dfSamples$sample)
i = match(rownames(mData), dfSamples$sample)
dfSamples = dfSamples[i,]

# create grouping factor
fSamples = as.character(dfSamples$group)
fSamples = factor(fSamples, levels = c('remote', 'recent'))
table(fSamples)
head(fSamples)
## sanity check
table(rownames(mData) %in% dfSamples$sample)

## create expression set object
dim(mData)
# make sure samples are in columns and probes/genes in rows
mData = t(mData)
oLumi = ExpressionSet(mData)

# pheno data
dfSamples$fSamples = fSamples
rownames(dfSamples) = dfSamples$sample
all(rownames(pData)==colnames(mData))

pData(oLumi) = dfSamples

# make lumi object from this
oLumi = addNuID2lumi(oLumi, lib.mapping='lumiHumanIDMapping')

### perform DE analysis
# remove unannotated genes
mDat = exprs(oLumi)
cvSym = getSYMBOL(rownames(mDat), 'lumiHumanAll.db')
# number of genes not annotated
table(is.na(cvSym))
mDat = mDat[!is.na(cvSym),]
# sanity check
dim(mDat)
fSamples = oLumi$fSamples

design = model.matrix(~fSamples)
colnames(design) = levels(fSamples)
head(design)

fit = lmFit(mDat, design)
fit = eBayes(fit)

# get annotation
df = select(lumiHumanAll.db, keys = rownames(mDat), columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# sanity check
nrow(mDat) == nrow(df)
# add annotation to limma object
fit$genes = df
topTable(fit, adjust='BH')

# look at top tables for each comparison
for (i in 2:length(levels(fSamples))){
  print(paste(levels(fSamples)[i], i))
  print(topTable(fit, coef=i, adjust='BH'))
}

# get the list of genes for each comparison i.e. each coefficient compared to base line
lSigGenes.adj = vector('list', length = length(levels(fSamples))-1)
names(lSigGenes.adj) = levels(fSamples)[2:length(levels(fSamples))]

for (i in 2:length(levels(fSamples))){
  p.adj = p.adjust(fit$p.value[,i], method = 'BH')
  lSigGenes.adj[[i-1]] = names(p.adj)[p.adj < 0.1]
}

#cvSigGenes.adj = unique(cvSigGenes.adj)
sapply(lSigGenes.adj, length)


######### Volcano plots
# plot volcano plots
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  f_plotVolcano(dfGenes, paste(names(n[i])), fc.lim = c(-2, 2), p.adj.cut = 0.1)
}

###### names of genes
dfRes = topTable(fit, adjust='BH', number=Inf, p.value=0.1)
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
cvCommonGenes = NULL
for (i in seq_along(n)) {
  cvCommonGenes = append(cvCommonGenes, lSigGenes.adj[[names(n[i])]])
}
cvCommonGenes = unique(cvCommonGenes)
mCommonGenes = sapply(seq_along(n), function(x) cvCommonGenes %in% lSigGenes.adj[[names(n[x])]])
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(n)

## all overexpressed genes if interested in
i = 1:nrow(mCommonGenes)

m1 = as.matrix(mCommonGenes[i,])
m1 = mDat[rownames(m1),]

fGroups = fSamples
colnames(m1) = fGroups
m1 = m1[,order(fGroups)]
fGroups = fGroups[order(fGroups)]

# ignore step if stabalization not required
m1 = t(apply(m1, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(m1) = fGroups

# scale across rows
rownames(m1) = dfRes[rownames(m1), 'SYMBOL']
mCounts = t(m1)
mCounts = scale(mCounts)
mCounts = t(mCounts)
# threshhold the values
mCounts[mCounts < -3] = -3
mCounts[mCounts > 3] = 3

library(NMF)
# draw the heatmap  color='-RdBu:50'
aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)


# write results csv files
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  dfGenes.2 = dfGenes[lSigGenes.adj[[names(n[i])]],]
  rownames(dfGenes.2) = NULL
  f = paste('Results/', 'Significant_genes_at_10pcFDR_Combined', names(n[i]), '.csv', sep='')
  dfGenes.2 = dfGenes.2[,c(2, 3, 4, 5, 6, 8, 9)]
  write.csv(dfGenes.2, file=f)
}







