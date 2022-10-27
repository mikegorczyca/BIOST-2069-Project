# From https://github.com/braunr/TimeSignatR/blob/master/example/TSexample.R

# set working directory corresponding to location in github repo in link above
directory = "C:/Users/14843/Desktop/TimeSignatR-master"
setwd(directory)

RNGversion("3.5.1")
set.seed(194) 

library(limma)
library(glmnet)

source("R/TimeStampFns.R")

# This file contains expression data for all three datasets,
# time, subject, and condition metadata, as well as previously
# obtained ZeitZeiger and PLSR predictions for comparison
load("example/DATA/TSexampleData.Rdata")









DE_limmaF = function(ctx, full_model, label, out.dir){
  fit.ctx = eBayes(contrasts.fit(full_model, ctx))%>%topTable(n = nrow(full_model$coefficient))
  fit.ctx$ID <- rownames(fit.ctx)
  fit.ctx$GeneSymbol = names(symbol[match(fit.ctx$ID, symbol)])
  # symbol[match(rownames(fit.ctx), names(symbol))]   # genesBM$mgi_symbol[match(rownames(fit.ctx), genesBM$ensembl_gene_id)]
  write.csv(fit.ctx, paste0(out.dir, "/DE_limmaF_", label, ".csv"))
  return(fit.ctx)
}

df$Region <- df$BrainRegion
regions = unique(df$Region)
library(edgeR)

all.expr <- all.expr[,order(colnames(all.expr))]
all.meta <- all.meta[order(all.meta$samp),]

all.meta <- all.meta[which(!is.na(all.meta$CPhrs)),]
all.expr <- all.expr[,which(colnames(all.expr) %in% all.meta$samp)]




library('annotate')
library('data.table')
library('foreach')
library('ggplot2')
library('limma')
library('limorhyde')
# library('org.Mm.eg.db')
library('qs')
library(splines)
recalibrateExprs <- function(exprMat,subjectIDs){
  # exprMat is a genes * subjects matrix
  # subjectIDs a vector of subject IDs
  exprs <- as.data.frame(t(exprMat))
  exprs <- split(exprs,subjectIDs)
  FCs <- lapply(exprs,scale,center=T,scale=F)
  FCs <- t(do.call(rbind,FCs))[,colnames(exprMat)]
  return(FCs)
}

for (a in unique(all.meta$study)) {
  a.meta = all.meta[which(all.meta$study==a),]
  cond = factor(a.meta$cond) # , levels = c(0, 1))
  a.meta$time_sin = sin(a.meta$LocalTime*pi/12)
  a.meta$time_cos = cos(a.meta$LocalTime*pi/12)
  a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
  a.expr <- recalibrateExprs(a.expr, a.meta$sID)
  # times = limorhyde(all.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
  rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
    design = model.matrix(~ time_cos + time_sin, data = a.meta[cond == condNow,])
    fit = lmFit(a.expr[, a.meta$cond == condNow], design)
    fit = eBayes(fit, trend = TRUE)
    rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
    setnames(rhyNow, 'rn', 'gene_id')
    rhyNow[, cond := condNow]}
  
  rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
  rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(rhyLimmaSummary, 'adj.P.Val')
  # write.csv(rhyLimmaSummary, paste0(, "_rhy_limma.csv"))
  a.meta$time_sin = sin(a.meta$CPhrs*pi/12)
  a.meta$time_cos = cos(a.meta$CPhrs*pi/12)
  # a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
  # times = limorhyde(all.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
  rhyLimma2 = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
    design = model.matrix(~ time_cos + time_sin, data = a.meta[cond == condNow,])
    fit = lmFit(a.expr[, a.meta$cond == condNow], design)
    fit = eBayes(fit, trend = TRUE)
    rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
    setnames(rhyNow, 'rn', 'gene_id')
    rhyNow[, cond := condNow]}
  
  rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
  rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(rhyLimmaSummary2, 'adj.P.Val')
  a1 = length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a2 = length(union(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a3 = length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$P.Value <0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$P.Value<0.01)]))
  a4 = length(union(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$P.Value <0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$P.Value<0.01)]))
  print(c(a1/a2, a3/a4))  # print(c(a3/a1, a3/a2, a6/a4, a6/a5))
  
  
}

print(c(a3/a1, a3/a2))
0.9743071 0.9088507 0.9749271 0.9181768
0.9752672 0.8777589 0.9712575 0.8892544
0.9792686 0.8178988 0.9718554 0.8426882
0.9809200 0.7518029 0.9744439 0.7936019


for (a in unique(all.meta$study)) {
  a.meta = all.meta[which(all.meta$study==a),]
  cond = factor(a.meta$cond) # , levels = c(0, 1))
  # times <- periodicSpline(a.meta$LocalTime*pi/12, 7, period = 2*pi, ord = 4L)
  times = limorhyde(a.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 7, intercept = FALSE)
  a.meta$t1 = times[,1]
  a.meta$t2 = times[,2]
  a.meta$t3 = times[,3]
  a.meta$t4 = times[,4]
  a.meta$t5 = times[,5]
  a.meta$t6 = times[,6]
  a.meta$t7 = times[,7]
  a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
  a.expr <- recalibrateExprs(a.expr, a.meta$sID)
  rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
    design = model.matrix(~ t1+t2+t3+t4+t5+t6+t7, data = a.meta[cond == condNow,])
    fit = lmFit(a.expr[, a.meta$cond == condNow], design)
    fit = eBayes(fit, trend = TRUE)
    rhyNow = data.table(topTable(fit, coef = 2:8, number = Inf), keep.rownames = TRUE)
    setnames(rhyNow, 'rn', 'gene_id')
    rhyNow[, cond := condNow]}
  
  rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
  rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(rhyLimmaSummary, 'adj.P.Val')
  # write.csv(rhyLimmaSummary, paste0(, "_rhy_limma.csv"))
  times = limorhyde(a.meta$CPhrs, period = 24, sinusoid = FALSE, nKnots = 7, intercept = FALSE)
  a.meta$t1 = times[,1]
  a.meta$t2 = times[,2]
  a.meta$t3 = times[,3]
  a.meta$t4 = times[,4]
  a.meta$t5 = times[,5]
  a.meta$t6 = times[,6]
  a.meta$t7 = times[,7]
  #a.meta$time_sin = sin(times*pi/12)
  #a.meta$time_cos = cos(times*pi/12)
  # a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
  # times = limorhyde(all.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
  rhyLimma2 = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
    design = model.matrix(~ t1+t2+t3+t4+t5+t6+t7, data = a.meta[cond == condNow,])
    fit = lmFit(a.expr[, a.meta$cond == condNow], design)
    fit = eBayes(fit, trend = TRUE)
    rhyNow = data.table(topTable(fit, coef = 2:8, number = Inf), keep.rownames = TRUE)
    setnames(rhyNow, 'rn', 'gene_id')
    rhyNow[, cond := condNow]}
  
  rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
  rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(rhyLimmaSummary2, 'adj.P.Val')
  
  a1 = length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a2 = length(union(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a3 = length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$P.Value <0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$P.Value<0.01)]))
  a4 = length(union(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$P.Value <0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$P.Value<0.01)]))
  print(c(a1/a2, a3/a4))  # print(c(a3/a1, a3/a2, a6/a4, a6/a5))
  
}

0: 0.8875783 0.8969957
3: 0.8586455 0.8664530
5: 0.8039778 0.8226132
7: 0.7409674 0.7774212


for (a in unique(all.meta$study)) {
  a.meta = all.meta[which(all.meta$study==a),]
  cond = factor(a.meta$cond) # , levels = c(0, 1))
  times = limorhyde(a.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 5, intercept = FALSE)
  a.meta$t1 = times[,1]
  a.meta$t2 = times[,2]
  a.meta$t3 = times[,3]
  a.meta$t4 = times[,4]
  a.meta$t5 = times[,5]
  a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
  a.expr <- recalibrateExprs(a.expr, a.meta$sID)
  rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
    design = model.matrix(~ t1+t2+t3+t4+t5, data = a.meta[cond == condNow,])
    fit = lmFit(a.expr[, a.meta$cond == condNow], design)
    fit = eBayes(fit, trend = TRUE)
    rhyNow = data.table(topTable(fit, coef = 2:6, number = Inf), keep.rownames = TRUE)
    setnames(rhyNow, 'rn', 'gene_id')
    rhyNow[, cond := condNow]}
  
  rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
  rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(rhyLimmaSummary, 'adj.P.Val')
  # write.csv(rhyLimmaSummary, paste0(, "_rhy_limma.csv"))
  times = limorhyde(a.meta$CPhrs, period = 24, sinusoid = FALSE, nKnots = 5, intercept = FALSE)
  a.meta$t1 = times[,1]
  a.meta$t2 = times[,2]
  a.meta$t3 = times[,3]
  a.meta$t4 = times[,4]
  a.meta$t5 = times[,5]
  #a.meta$time_sin = sin(times*pi/12)
  #a.meta$time_cos = cos(times*pi/12)
  # a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
  # times = limorhyde(all.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
  rhyLimma2 = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
    design = model.matrix(~ t1+t2+t3+t4+t5, data = a.meta[cond == condNow,])
    fit = lmFit(a.expr[, a.meta$cond == condNow], design)
    fit = eBayes(fit, trend = TRUE)
    rhyNow = data.table(topTable(fit, coef = 2:6, number = Inf), keep.rownames = TRUE)
    setnames(rhyNow, 'rn', 'gene_id')
    rhyNow[, cond := condNow]}
  
  rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
  rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(rhyLimmaSummary2, 'adj.P.Val')
  setorderv(rhyLimmaSummary2, 'adj.P.Val')
  a1 = length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a2 = length(union(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a3 = length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$P.Value <0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$P.Value<0.01)]))
  a4 = length(union(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$P.Value <0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$P.Value<0.01)]))
  print(c(a1/a2, a3/a4)) 
}


for (a in unique(all.meta$study)) {
  a.meta = all.meta[which(all.meta$study==a),]
  cond = factor(a.meta$cond) # , levels = c(0, 1))
  times = limorhyde(a.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
  a.meta$t1 = times[,1]
  a.meta$t2 = times[,2]
  a.meta$t3 = times[,3]
  a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
  a.expr <- recalibrateExprs(a.expr, a.meta$sID)
  rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
    design = model.matrix(~ t1+t2+t3, data = a.meta[cond == condNow,])
    fit = lmFit(a.expr[, a.meta$cond == condNow], design)
    fit = eBayes(fit, trend = TRUE)
    rhyNow = data.table(topTable(fit, coef = 2:4, number = Inf), keep.rownames = TRUE)
    setnames(rhyNow, 'rn', 'gene_id')
    rhyNow[, cond := condNow]}
  
  rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
  rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(rhyLimmaSummary, 'adj.P.Val')
  # write.csv(rhyLimmaSummary, paste0(, "_rhy_limma.csv"))
  times = limorhyde(a.meta$CPhrs, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
  a.meta$t1 = times[,1]
  a.meta$t2 = times[,2]
  a.meta$t3 = times[,3]
  #a.meta$time_sin = sin(times*pi/12)
  #a.meta$time_cos = cos(times*pi/12)
  # a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
  # times = limorhyde(all.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
  rhyLimma2 = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
    design = model.matrix(~ t1+t2+t3, data = a.meta[cond == condNow,])
    fit = lmFit(a.expr[, a.meta$cond == condNow], design)
    fit = eBayes(fit, trend = TRUE)
    rhyNow = data.table(topTable(fit, coef = 2:4, number = Inf), keep.rownames = TRUE)
    setnames(rhyNow, 'rn', 'gene_id')
    rhyNow[, cond := condNow]}
  
  rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
  rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(rhyLimmaSummary2, 'adj.P.Val')
  setorderv(rhyLimmaSummary2, 'adj.P.Val')
  a1 = length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a2 = length(union(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a3 = length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$P.Value <0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$P.Value<0.01)]))
  a4 = length(union(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$P.Value <0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$P.Value<0.01)]))
  print(c(a1/a2, a3/a4)) 
}



cosinor: 0.9737425 0.9353095 0.9725640 0.9382419
3-spline: 0.9675269 0.9099823 0.9645578 0.9223947
5-spline: 0.9676054 0.8680647 0.9617689 0.8826192
7-spline: 0.9698977 0.8094122 0.9638509 0.8385204

for (a in unique(all.meta$study)) {
  a.meta = all.meta[which(all.meta$study==a),]
  cond = factor(a.meta$cond) # , levels = c(0, 1))
  times = limorhyde(a.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 2, intercept = FALSE)
  a.meta$t1 = times[,1]
  a.meta$t2 = times[,2]
  a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
  a.expr <- recalibrateExprs(a.expr, a.meta$sID)
  rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
    design = model.matrix(~ t1+t2, data = a.meta[cond == condNow,])
    fit = lmFit(a.expr[, a.meta$cond == condNow], design)
    fit = eBayes(fit, trend = TRUE)
    rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
    setnames(rhyNow, 'rn', 'gene_id')
    rhyNow[, cond := condNow]}
  
  rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
  rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(rhyLimmaSummary, 'adj.P.Val')
  # write.csv(rhyLimmaSummary, paste0(, "_rhy_limma.csv"))
  times = limorhyde(a.meta$CPhrs, period = 24, sinusoid = FALSE, nKnots = 2, intercept = FALSE)
  a.meta$t1 = times[,1]
  a.meta$t2 = times[,2]
  #a.meta$time_sin = sin(times*pi/12)
  #a.meta$time_cos = cos(times*pi/12)
  # a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
  # times = limorhyde(all.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
  rhyLimma2 = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
    design = model.matrix(~ t1+t2+t3, data = a.meta[cond == condNow,])
    fit = lmFit(a.expr[, a.meta$cond == condNow], design)
    fit = eBayes(fit, trend = TRUE)
    rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
    setnames(rhyNow, 'rn', 'gene_id')
    rhyNow[, cond := condNow]}
  
  rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
  rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(rhyLimmaSummary2, 'adj.P.Val')
  setorderv(rhyLimmaSummary2, 'adj.P.Val')
  length(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)])
  length(rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)])
  length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a1 = length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a2 = length(union(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$adj.P.Val<0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$adj.P.Val<0.01)]))
  a3 = length(intersect(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$P.Value <0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$P.Value<0.01)]))
  a4 = length(union(rhyLimmaSummary$gene_id[which(rhyLimmaSummary$P.Value <0.01)], rhyLimmaSummary2$gene_id[which(rhyLimmaSummary2$P.Value<0.01)]))
  print(c(a1/a2, a3/a4)) 
}



