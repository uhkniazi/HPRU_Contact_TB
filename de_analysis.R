# File: de_analysis.R
# Desc: Data set from tb contacts
# Auth: u.niazi@imperial.ac.uk
# Date: 15/06/2016

library(annotate)
library(limma)
library(org.Hs.eg.db)

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
dfSamples = read.csv('Data_external/Contact_DWC_data/Sample_annotation.csv', header=T)
dfData = read.csv('Data_external/Contact_DWC_data/Contact_normalised_filtered_expr.txt', sep = '\t', header=T)

# data sorting
cn = colnames(dfData)[6:(ncol(dfData))]
cn = gsub('X(.+)', replacement = '\\1', cn, perl = T)

dfData.sub = dfData[,6:(ncol(dfData))]
colnames(dfData.sub) = cn

mData = as.matrix(dfData.sub)
rownames(mData) = dfData$EntrezID
mData = t(mData)

## put the two in the same order i.e. sample annotations and data
table(rownames(mData) %in% dfSamples$ArrayID)
i = match(rownames(mData), dfSamples$ArrayID)
dfSamples = dfSamples[i,]

# create grouping factor
fSamples = as.character(dfSamples$Group)
i = which(fSamples %in% c('TB_INFECTION', 'TB_RESISTANCE'))
fSamples = factor(fSamples[i])
mData = mData[i,]

## sanity check
table(rownames(mData) %in% dfSamples$ArrayID[i])

### perform DE analysis
mDat = t(mData)
cvSym = select(org.Hs.eg.db, keys = rownames(mDat), columns = c('SYMBOL'), keytype = 'ENTREZID')$SYMBOL

# number of genes not annotated
table(is.na(cvSym))
# remove unannotated genes
mDat = mDat[!is.na(cvSym),]

design = model.matrix(~fSamples)
colnames(design) = levels(fSamples)
head(design)

fit = lmFit(mDat, design)
fit = eBayes(fit)

# get annotation
df = select(org.Hs.eg.db, keys = rownames(mDat), columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
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
  p.adj = p.adjust(fit$p.value[,i], method = 'none')
  lSigGenes.adj[[i-1]] = names(p.adj)[p.adj < 0.01]
}

#cvSigGenes.adj = unique(cvSigGenes.adj)
sapply(lSigGenes.adj, length)

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

## make heatmap
i = 1:nrow(mCommonGenes)

m1 = as.matrix(mCommonGenes[i,])
m1 = mDat[rownames(m1),]

source('../CGraphClust/CGraphClust.R')
library(NMF)
fGroups = fSamples
colnames(m1) = fGroups
m1 = m1[,order(fGroups)]
fGroups = fGroups[order(fGroups)]

# ignore step if stabalization not required
m1 = t(apply(m1, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(m1) = fGroups

# scale across rows
rn = f_dfGetGeneAnnotation(rownames(m1))
rownames(m1) = rn$SYMBOL
mCounts = t(m1)
mCounts = scale(mCounts)
mCounts = t(mCounts)
# threshhold the values
mCounts[mCounts < -3] = -3
mCounts[mCounts > 3] = 3

# draw the heatmap  color='-RdBu:50'
aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)


# write results csv files
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  dfGenes.2 = dfGenes[lSigGenes.adj[[names(n[i])]] %in% dfGenes$ID,]
  rownames(dfGenes.2) = NULL
  f = paste('Temp/', 'Significant_genes_at_p_0.01', names(n[i]), '.csv', sep='')
  dfGenes.2 = dfGenes.2[,c(2, 3, 4, 5, 6, 8, 9)]
  write.csv(dfGenes.2, file=f)
}

############################ variable selection

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')


## format data for variable selection
mData = mDat[rownames(mCommonGenes),]
fGroups = fSamples
mData = t(mData)
table(fGroups)

cn = f_dfGetGeneAnnotation(colnames(mData))
colnames(mData) = cn$SYMBOL

dfData = data.frame(mData)
## random forest step
oVar.r = CVariableSelection.RandomForest(dfData, groups = fGroups, boot.num=100, big.warn = F)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
cvTopGenes = rownames(dfRF)[1:20]
f_dfGetGeneAnnotation(gsub('X(.+)', '\\1', cvTopGenes, perl=T))

dfData = dfData[,colnames(dfData) %in% cvTopGenes]
oVar.s = CVariableSelection.ReduceModel(dfData, fGroups, boot.num=100)
plot.var.selection(oVar.s)

## select test set
test = sample(1:nrow(dfData), nrow(dfData) * 0.2, replace = F)
table(fGroups[test]); table(fGroups[-test])

dfPrint = NULL
par(mfrow=c(1,2))
## 10 fold cv
## 10 fold nested cross validation with various variable combinations
# try models of various sizes with CV
for (i in 1:6){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.s, i)
  dfData.train = as.data.frame(dfData[-test, cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  dfData.test = data.frame(dfData[test, cvTopGenes.sub])
  colnames(dfData.test) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = (dfData.test), train.dat = (dfData.train), test.groups = fGroups[test],
                             train.groups = fGroups[-test], level.predict = 'TB_RESISTANCE', boot.num = 500, k.fold = 5)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  x = getAUCVector(oCV)
  vc = (paste('Variable Count', i))
  gn = paste(cvTopGenes.sub, collapse = ' ')
  sig = (signif(quantile(x, probs = c(0.025, 0.975)), 2))
  dfPrint = rbind(dfPrint, c(vc, gn, sig))
}

dfPrint





