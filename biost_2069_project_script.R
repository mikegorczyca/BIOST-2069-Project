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

#df$Region <- df$BrainRegion
#regions = unique(df$Region)
library(edgeR)

all.expr <- all.expr[,order(colnames(all.expr))]
all.meta <- all.meta[order(all.meta$samp),]

all.meta <- all.meta[which(!is.na(all.meta$CPhrs)),]
all.expr <- all.expr[,which(colnames(all.expr) %in% all.meta$samp)]


study <- all.expr[,which(all.meta$study=="TrTe")]
info <- all.meta[which(all.meta$study=="TrTe"),]
info$time_cos <- cos(info$CPhrs*pi/12) # LocalTime*pi/12) # CPhrs*pi/12)
info$time_sin <- sin(info$CPhrs*pi/12) # LocalTime*pi/12) # CPhrs*pi/12)
design = model.matrix(~ time_cos + time_sin, data = info)
fit = lmFit(study, design)
fit = eBayes(fit) #, trend = TRUE)
rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
setnames(rhyNow, 'rn', 'gene_id')

phase <- atan2(rhyNow$time_cos, rhyNow$time_sin)
library(overlap)
w <- densityFit(phase+pi, seq(0, 2*pi, 0.01), 100)
plot(seq(0, 2*pi, 0.01), w, type = "l", col = "red", 
     xlab = "Phase Estimate", ylab = "Density", lwd = 2,
     main = "Study 3 ICT")
# w <- vm.kde(phase+pi,  thumb = "rot")


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

for (a in c("TrTe", "V1")) {
  for (p in c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3)) {
    
    a.meta = all.meta[which(all.meta$study==a),]
    cond = factor(a.meta$cond) # , levels = c(0, 1))
    a.meta$cond <- as.factor(a.meta$cond)
    a.meta$time_sin = sin(a.meta$LocalTime*pi/12)
    a.meta$time_cos = cos(a.meta$LocalTime*pi/12)
    a.expr <- all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
    a.expr <- recalibrateExprs(a.expr, a.meta$sID)
    rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
      design = model.matrix(~ time_cos + time_sin, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimma$name <- paste0(rhyLimma$gene_id, rhyLimma$cond)
    rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
    rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(rhyLimmaSummary, 'adj.P.Val')
    
    
    #############################################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    a = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma = data.table(topTable(fit, coef = 5:6, number = Inf), keep.rownames = TRUE)
    setnames(drLimma, 'rn', 'gene_id')
    drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= p]$gene_id]
    drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma, 'adj.P.Val')
    #############################################################################

    
    # write.csv(rhyLimmaSummary, paste0(, "_rhy_limma.csv"))
    a.meta$time_sin = sin(a.meta$CPhrs*pi/12)
    a.meta$time_cos = cos(a.meta$CPhrs*pi/12)
    # a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
    # times = limorhyde(all.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
    rhyLimma2 = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
      design = model.matrix(~ time_cos + time_sin, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
    rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    rhyLimma2$name <- paste0(rhyLimma2$gene_id, rhyLimma2$cond)
    
    #############################################
    design = model.matrix(~ cond * (time_cos + time_sin), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma2 = data.table(topTable(fit, coef = 5:6, number = Inf), keep.rownames = TRUE)
    setnames(drLimma2, 'rn', 'gene_id')
    drLimma2 = drLimma2[gene_id %in% rhyLimmaSummary2[adj.P.Val <= p]$gene_id]
    drLimma2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma2, 'adj.P.Val')
    #############################################################################

    a1 = length(intersect(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a2 = length(union(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a3 = length(intersect(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    a4 = length(union(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    #a1 = length(intersect(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a2 = length(union(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a3 = length(intersect(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    #a4 = length(union(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$adj.P.Val<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$adj.P.Val<p)] <- 1
    # a5 <- adj.rand.index(cl1, cl2)
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$P.Value<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$P.Value<p)] <- 1
    #a6 <- adj.rand.index(cl1, cl2)
    or <- c(a1/a2, a3/a4, a1, a2, a3, a4, 
                  length(which(rhyLimma$adj.P.Val<p)),
                  length(which(rhyLimma$P.Value<p)),
                  length(which(rhyLimma2$adj.P.Val<p)),
                  length(which(rhyLimma2$P.Value<p)))
    ###############################################################################################
    
    #a.meta = all.meta[which(all.meta$study==a),]
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
    #a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
    #a.expr <- recalibrateExprs(a.expr, a.meta$sID)
    rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
      design = model.matrix(~ t1+t2+t3+t4+t5+t6+t7, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:8, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3+t4+t5+t6+t7), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma = data.table(topTable(fit, coef = 10:16, number = Inf), keep.rownames = TRUE)
    setnames(drLimma, 'rn', 'gene_id')
    # drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= p]$gene_id]
    drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma, 'adj.P.Val')
    #############################################################################
    
    rhyLimma$name <- paste0(rhyLimma$gene_id, rhyLimma$cond)
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
      design = model.matrix(~ t1+t2+t3+t4+t5+t6+t7, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:8, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimma2$name <- paste0(rhyLimma2$gene_id, rhyLimma2$cond)
    rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
    rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3+t4+t5+t6+t7), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma2 = data.table(topTable(fit, coef = 10:16, number = Inf), keep.rownames = TRUE)
    setnames(drLimma2, 'rn', 'gene_id')
    # drLimma2 = drLimma2[gene_id %in% rhyLimmaSummary2[adj.P.Val <= p]$gene_id]
    drLimma2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma2, 'adj.P.Val')
    #############################################################################
    
    
    
    
    a1 = length(intersect(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a2 = length(union(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a3 = length(intersect(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    a4 = length(union(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    #a1 = length(intersect(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a2 = length(union(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a3 = length(intersect(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    #a4 = length(union(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$adj.P.Val<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$adj.P.Val<p)] <- 1
    # a5 <- adj.rand.index(cl1, cl2)
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$P.Value<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$P.Value<p)] <- 1
    #a6 <- adj.rand.index(cl1, cl2)
    or <- rbind(or,
                c(a1/a2, a3/a4, a1, a2, a3, a4, 
                  length(which(rhyLimma$adj.P.Val<p)),
                  length(which(rhyLimma$P.Value<p)),
                  length(which(rhyLimma2$adj.P.Val<p)),
                  length(which(rhyLimma2$P.Value<p)))) 
    ###############################################################################################
    
    #a.meta = all.meta[which(all.meta$study==a),]
    cond = factor(a.meta$cond) # , levels = c(0, 1))
    times = limorhyde(a.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 5, intercept = FALSE)
    a.meta$t1 = times[,1]
    a.meta$t2 = times[,2]
    a.meta$t3 = times[,3]
    a.meta$t4 = times[,4]
    a.meta$t5 = times[,5]
    #a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
    #a.expr <- recalibrateExprs(a.expr, a.meta$sID)
    rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
      design = model.matrix(~ t1+t2+t3+t4+t5, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:6, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
    # rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
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
      design = model.matrix(~ t1+t2+t3+t4+t5, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:6, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
    rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    a1 = length(intersect(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a2 = length(union(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a3 = length(intersect(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    a4 = length(union(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    #a1 = length(intersect(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a2 = length(union(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a3 = length(intersect(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    #a4 = length(union(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$adj.P.Val<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$adj.P.Val<p)] <- 1
    # a5 <- adj.rand.index(cl1, cl2)
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$P.Value<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$P.Value<p)] <- 1
    #a6 <- adj.rand.index(cl1, cl2)
    
    
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3+t4+t5), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma = data.table(topTable(fit, coef = 8:12, number = Inf), keep.rownames = TRUE)
    setnames(drLimma, 'rn', 'gene_id')
    # drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= p]$gene_id]
    drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma, 'adj.P.Val')
    #############################################################################
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3+t4+t5), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma2 = data.table(topTable(fit, coef = 8:12, number = Inf), keep.rownames = TRUE)
    setnames(drLimma2, 'rn', 'gene_id')
    # drLimma2 = drLimma2[gene_id %in% rhyLimmaSummary2[adj.P.Val <= p]$gene_id]
    drLimma2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma2, 'adj.P.Val')
    #############################################################################
    a1 = length(intersect(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a2 = length(union(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a3 = length(intersect(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    a4 = length(union(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    rhyLimma$name <- paste0(rhyLimma$gene_id, rhyLimma$cond)
    rhyLimma2$name <- paste0(rhyLimma2$gene_id, rhyLimma2$cond)
    #a1 = length(intersect(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a2 = length(union(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a3 = length(intersect(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    #a4 = length(union(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    or <- rbind(or,
                c(a1/a2, a3/a4, a1, a2, a3, a4, 
                  length(which(rhyLimma$adj.P.Val<p)),
                  length(which(rhyLimma$P.Value<p)),
                  length(which(rhyLimma2$adj.P.Val<p)),
                  length(which(rhyLimma2$P.Value<p)))) 
    ###############################################################################################
    
    #a.meta = all.meta[which(all.meta$study==a),]
    cond = factor(a.meta$cond) # , levels = c(0, 1))
    times = limorhyde(a.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
    a.meta$t1 = times[,1]
    a.meta$t2 = times[,2]
    a.meta$t3 = times[,3]
    #a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
    #a.expr <- recalibrateExprs(a.expr, a.meta$sID)
    rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
      design = model.matrix(~ t1+t2+t3, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
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
      design = model.matrix(~ t1+t2+t3, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:4, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
    rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma = data.table(topTable(fit, coef = 6:8, number = Inf), keep.rownames = TRUE)
    setnames(drLimma, 'rn', 'gene_id')
    # drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= p]$gene_id]
    drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma, 'adj.P.Val')
    #############################################################################
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma2 = data.table(topTable(fit, coef = 6:8, number = Inf), keep.rownames = TRUE)
    setnames(drLimma2, 'rn', 'gene_id')
    # drLimma2 = drLimma2[gene_id %in% rhyLimmaSummary2[adj.P.Val <= p]$gene_id]
    drLimma2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma2, 'adj.P.Val')
    #############################################################################
    rhyLimma$name <- paste0(rhyLimma$gene_id, rhyLimma$cond)
    rhyLimma2$name <- paste0(rhyLimma2$gene_id, rhyLimma2$cond)
    a1 = length(intersect(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a2 = length(union(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a3 = length(intersect(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    a4 = length(union(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    rhyLimma$name <- paste0(rhyLimma$gene_id, rhyLimma$cond)
    rhyLimma2$name <- paste0(rhyLimma2$gene_id, rhyLimma2$cond)
    #a1 = length(intersect(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a2 = length(union(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a3 = length(intersect(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    #a4 = length(union(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$adj.P.Val<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$adj.P.Val<p)] <- 1
    # a5 <- adj.rand.index(cl1, cl2)
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$P.Value<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$P.Value<p)] <- 1
    #a6 <- adj.rand.index(cl1, cl2)
    or <- rbind(or,
                c(a1/a2, a3/a4, a1, a2, a3, a4, 
                  length(which(rhyLimma$adj.P.Val<p)),
                  length(which(rhyLimma$P.Value<p)),
                  length(which(rhyLimma2$adj.P.Val<p)),
                  length(which(rhyLimma2$P.Value<p)))) #a5, a6, a1, a2, a3, a4)) #
    print(c("study", a, "thresh", p))
    print(or[c(1,4,3,2), c(1, 3, 4, 7, 9,
                           2, 5, 6, 8, 10)])
  }
}



















for (a in c("TrTe", "V1")) {
  for (p in c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3)) {
    
    a.meta = all.meta[which(all.meta$study==a),]
    a.meta$c1 <- 0
    a.meta$c1[which(a.meta$cond==1)] <- 1
    a.meta$c2 <- 0
    a.meta$c2[which(a.meta$cond==2)] <- 1
    cond = factor(a.meta$cond) # , levels = c(0, 1))
    a.meta$cond <- as.factor(a.meta$cond)
    a.meta$time_sin = sin(a.meta$LocalTime*pi/12)
    a.meta$time_cos = cos(a.meta$LocalTime*pi/12)
    a.expr <- all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
    a.expr <- recalibrateExprs(a.expr, a.meta$sID)
    rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
      design = model.matrix(~ time_cos + time_sin, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimma$name <- paste0(rhyLimma$gene_id, rhyLimma$cond)
    rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
    rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(rhyLimmaSummary, 'adj.P.Val')
    
    
    #############################################
    # design = model.matrix(~ cond * (time_cos + time_sin), data = a.meta)
    design = model.matrix(~ c1:(time_cos + time_sin) + 
                            c2:(time_cos + time_sin)+c1+c2+0, data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma = data.table(topTable(fit, coef = 5:6, number = Inf), keep.rownames = TRUE)
    setnames(drLimma, 'rn', 'gene_id')
    drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= p]$gene_id]
    drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma, 'adj.P.Val')
    #############################################################################
    
    
    # write.csv(rhyLimmaSummary, paste0(, "_rhy_limma.csv"))
    a.meta$time_sin = sin(a.meta$CPhrs*pi/12)
    a.meta$time_cos = cos(a.meta$CPhrs*pi/12)
    # a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
    # times = limorhyde(all.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
    rhyLimma2 = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
      design = model.matrix(~ time_cos + time_sin, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
    rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    rhyLimma2$name <- paste0(rhyLimma2$gene_id, rhyLimma2$cond)
    
    #############################################
    design = model.matrix(~ cond * (time_cos + time_sin), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma2 = data.table(topTable(fit, coef = 5:6, number = Inf), keep.rownames = TRUE)
    setnames(drLimma2, 'rn', 'gene_id')
    drLimma2 = drLimma2[gene_id %in% rhyLimmaSummary2[adj.P.Val <= p]$gene_id]
    drLimma2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma2, 'adj.P.Val')
    #############################################################################
    
    a1 = length(intersect(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a2 = length(union(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a3 = length(intersect(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    a4 = length(union(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    #a1 = length(intersect(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a2 = length(union(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a3 = length(intersect(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    #a4 = length(union(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$adj.P.Val<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$adj.P.Val<p)] <- 1
    # a5 <- adj.rand.index(cl1, cl2)
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$P.Value<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$P.Value<p)] <- 1
    #a6 <- adj.rand.index(cl1, cl2)
    or <- c(a1/a2, a3/a4, a1, a2, a3, a4, 
            length(which(rhyLimma$adj.P.Val<p)),
            length(which(rhyLimma$P.Value<p)),
            length(which(rhyLimma2$adj.P.Val<p)),
            length(which(rhyLimma2$P.Value<p)))
    ###############################################################################################
    
    #a.meta = all.meta[which(all.meta$study==a),]
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
    #a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
    #a.expr <- recalibrateExprs(a.expr, a.meta$sID)
    rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
      design = model.matrix(~ t1+t2+t3+t4+t5+t6+t7, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:8, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3+t4+t5+t6+t7), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma = data.table(topTable(fit, coef = 10:16, number = Inf), keep.rownames = TRUE)
    setnames(drLimma, 'rn', 'gene_id')
    # drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= p]$gene_id]
    drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma, 'adj.P.Val')
    #############################################################################
    
    rhyLimma$name <- paste0(rhyLimma$gene_id, rhyLimma$cond)
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
      design = model.matrix(~ t1+t2+t3+t4+t5+t6+t7, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:8, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimma2$name <- paste0(rhyLimma2$gene_id, rhyLimma2$cond)
    rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
    rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3+t4+t5+t6+t7), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma2 = data.table(topTable(fit, coef = 10:16, number = Inf), keep.rownames = TRUE)
    setnames(drLimma2, 'rn', 'gene_id')
    # drLimma2 = drLimma2[gene_id %in% rhyLimmaSummary2[adj.P.Val <= p]$gene_id]
    drLimma2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma2, 'adj.P.Val')
    #############################################################################
    
    
    
    
    a1 = length(intersect(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a2 = length(union(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a3 = length(intersect(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    a4 = length(union(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    #a1 = length(intersect(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a2 = length(union(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a3 = length(intersect(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    #a4 = length(union(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$adj.P.Val<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$adj.P.Val<p)] <- 1
    # a5 <- adj.rand.index(cl1, cl2)
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$P.Value<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$P.Value<p)] <- 1
    #a6 <- adj.rand.index(cl1, cl2)
    or <- rbind(or,
                c(a1/a2, a3/a4, a1, a2, a3, a4, 
                  length(which(rhyLimma$adj.P.Val<p)),
                  length(which(rhyLimma$P.Value<p)),
                  length(which(rhyLimma2$adj.P.Val<p)),
                  length(which(rhyLimma2$P.Value<p)))) 
    ###############################################################################################
    
    #a.meta = all.meta[which(all.meta$study==a),]
    cond = factor(a.meta$cond) # , levels = c(0, 1))
    times = limorhyde(a.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 5, intercept = FALSE)
    a.meta$t1 = times[,1]
    a.meta$t2 = times[,2]
    a.meta$t3 = times[,3]
    a.meta$t4 = times[,4]
    a.meta$t5 = times[,5]
    #a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
    #a.expr <- recalibrateExprs(a.expr, a.meta$sID)
    rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
      design = model.matrix(~ t1+t2+t3+t4+t5, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:6, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimmaSummary = rhyLimma[, .(P.Value = min(P.Value)), by = gene_id]
    # rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
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
      design = model.matrix(~ t1+t2+t3+t4+t5, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:6, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
    rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    a1 = length(intersect(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a2 = length(union(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a3 = length(intersect(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    a4 = length(union(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    #a1 = length(intersect(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a2 = length(union(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a3 = length(intersect(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    #a4 = length(union(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$adj.P.Val<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$adj.P.Val<p)] <- 1
    # a5 <- adj.rand.index(cl1, cl2)
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$P.Value<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$P.Value<p)] <- 1
    #a6 <- adj.rand.index(cl1, cl2)
    
    
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3+t4+t5), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma = data.table(topTable(fit, coef = 8:12, number = Inf), keep.rownames = TRUE)
    setnames(drLimma, 'rn', 'gene_id')
    # drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= p]$gene_id]
    drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma, 'adj.P.Val')
    #############################################################################
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3+t4+t5), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma2 = data.table(topTable(fit, coef = 8:12, number = Inf), keep.rownames = TRUE)
    setnames(drLimma2, 'rn', 'gene_id')
    # drLimma2 = drLimma2[gene_id %in% rhyLimmaSummary2[adj.P.Val <= p]$gene_id]
    drLimma2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma2, 'adj.P.Val')
    #############################################################################
    a1 = length(intersect(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a2 = length(union(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a3 = length(intersect(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    a4 = length(union(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    rhyLimma$name <- paste0(rhyLimma$gene_id, rhyLimma$cond)
    rhyLimma2$name <- paste0(rhyLimma2$gene_id, rhyLimma2$cond)
    #a1 = length(intersect(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a2 = length(union(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a3 = length(intersect(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    #a4 = length(union(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    or <- rbind(or,
                c(a1/a2, a3/a4, a1, a2, a3, a4, 
                  length(which(rhyLimma$adj.P.Val<p)),
                  length(which(rhyLimma$P.Value<p)),
                  length(which(rhyLimma2$adj.P.Val<p)),
                  length(which(rhyLimma2$P.Value<p)))) 
    ###############################################################################################
    
    #a.meta = all.meta[which(all.meta$study==a),]
    cond = factor(a.meta$cond) # , levels = c(0, 1))
    times = limorhyde(a.meta$LocalTime, period = 24, sinusoid = FALSE, nKnots = 3, intercept = FALSE)
    a.meta$t1 = times[,1]
    a.meta$t2 = times[,2]
    a.meta$t3 = times[,3]
    #a.expr = all.expr[, which(colnames(all.expr) %in% a.meta$samp)]
    #a.expr <- recalibrateExprs(a.expr, a.meta$sID)
    rhyLimma = foreach(condNow = unique(a.meta$cond), .combine = rbind) %do% {
      design = model.matrix(~ t1+t2+t3, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
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
      design = model.matrix(~ t1+t2+t3, data = a.meta[which(a.meta$cond==condNow),])
      fit = lmFit(a.expr[, which(a.meta$cond==condNow)], design)
      fit = eBayes(fit) #, trend = TRUE)
      rhyNow = data.table(topTable(fit, coef = 2:4, number = Inf), keep.rownames = TRUE)
      setnames(rhyNow, 'rn', 'gene_id')
      rhyNow[, cond := condNow]}
    
    rhyLimmaSummary2 = rhyLimma2[, .(P.Value = min(P.Value)), by = gene_id]
    rhyLimmaSummary2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    setorderv(rhyLimmaSummary2, 'adj.P.Val')
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma = data.table(topTable(fit, coef = 6:8, number = Inf), keep.rownames = TRUE)
    setnames(drLimma, 'rn', 'gene_id')
    # drLimma = drLimma[gene_id %in% rhyLimmaSummary[adj.P.Val <= p]$gene_id]
    drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma, 'adj.P.Val')
    #############################################################################
    #############################################
    design = model.matrix(~ cond * (t1+t2+t3), data = a.meta)
    fit = lmFit(a.expr, design)
    fit = eBayes(fit) #, trend = TRUE)
    drLimma2 = data.table(topTable(fit, coef = 6:8, number = Inf), keep.rownames = TRUE)
    setnames(drLimma2, 'rn', 'gene_id')
    # drLimma2 = drLimma2[gene_id %in% rhyLimmaSummary2[adj.P.Val <= p]$gene_id]
    drLimma2[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
    setorderv(drLimma2, 'adj.P.Val')
    #############################################################################
    rhyLimma$name <- paste0(rhyLimma$gene_id, rhyLimma$cond)
    rhyLimma2$name <- paste0(rhyLimma2$gene_id, rhyLimma2$cond)
    a1 = length(intersect(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a2 = length(union(drLimma$gene_id[which(drLimma$adj.P.Val<p)], drLimma2$gene_id[which(drLimma2$adj.P.Val<p)]))
    a3 = length(intersect(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    a4 = length(union(drLimma$gene_id[which(drLimma$P.Value <p)], drLimma2$gene_id[which(drLimma2$P.Value<p)]))
    rhyLimma$name <- paste0(rhyLimma$gene_id, rhyLimma$cond)
    rhyLimma2$name <- paste0(rhyLimma2$gene_id, rhyLimma2$cond)
    #a1 = length(intersect(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a2 = length(union(rhyLimma$name[which(rhyLimma$adj.P.Val<p)], rhyLimma2$name[which(rhyLimma2$adj.P.Val<p)]))
    #a3 = length(intersect(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    #a4 = length(union(rhyLimma$name[which(rhyLimma$P.Value <p)], rhyLimma2$name[which(rhyLimma2$P.Value<p)]))
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$adj.P.Val<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$adj.P.Val<p)] <- 1
    # a5 <- adj.rand.index(cl1, cl2)
    cl1 <- rep(0, nrow(rhyLimmaSummary))
    cl1[which(rhyLimmaSummary$P.Value<p)] <- 1
    cl2 <- rep(0, nrow(rhyLimmaSummary2))
    cl2[which(rhyLimmaSummary2$P.Value<p)] <- 1
    #a6 <- adj.rand.index(cl1, cl2)
    or <- rbind(or,
                c(a1/a2, a3/a4, a1, a2, a3, a4, 
                  length(which(rhyLimma$adj.P.Val<p)),
                  length(which(rhyLimma$P.Value<p)),
                  length(which(rhyLimma2$adj.P.Val<p)),
                  length(which(rhyLimma2$P.Value<p)))) #a5, a6, a1, a2, a3, a4)) #
    print(c("study", a, "thresh", p))
    print(or[c(1,4,3,2), c(1, 3, 4, 7, 9,
                           2, 5, 6, 8, 10)])
  }
}






fit = getModelFit(a.expr, a.meta)
fit = getPosteriorFit(fit)
rhyStats = getRhythmStats(fit, features = c('13170', '13869'))






library(GEOquery)
# dfl = getGSEDataTables("GSE3494")
#gse <- getGEO("gse39445") # "gse48113")
gse <- getGEO("gse39445")
test = gse[[1]] # get just the first element in the list
head(fData(test))
symbols = fData(test)[,'GENE_SYMBOL']
expr <- gse$GSE39445_series_matrix.txt.gz@assayData$exprs # gse$GSE48113_series_matrix.txt.gz@assayData$exprs
rownames(expr) <- symbols
mn <- c()
for (i in 1:nrow(expr)) {
  mn <- c(mn, mean(na.omit(expr[i,])))
}
# mn <-apply(expr, 1, mean, na.omit=TRUE)
for (i in 1:length(mn)) {
  expr[i,which(is.na(expr[i,]))] <- mn[i]
}
# expr[is.na(expr)] <- mean(na.omit(expr))
expr <- expr[,order(colnames(expr))]
# expr <- expr[,-which(colnames(expr)=="GSM1168837")]
# expr < -
od <-apply(expr, 1, sd)
expr <- expr[od>0.5, ] # 22570 genes

# archer - 12895
# try this
# this has been changed
# Filtered_Data = raw.counts[which(countCPM>0.5),] 
# two consecutive missing time points removed
# 
gene_variance = apply(expr, 1, function(a){mean(a)/sd(a)})

#df2 = as.data.frame(list(exp = unlist(log2(expr+1)),
#                         Sample_name = rep(1:287, each = nrow(expr))))
#library(ggplot2)
#ggplot(aes(y = exp, x = Sample_name), data = df2) +
#  geom_boxplot(aes(fill = Region, color = Region))+
#  coord_flip()+
#  ggtitle(paste("Raw Counts of Filtered Genes (log2 scale)"))

library(limma)
d = DGEList(counts = Filtered_Data)
design = model.matrix(~treatment) #  + df$PMI) # + df$pH)
# https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/variancePartition/inst/doc/dream.html
y = limma::voom(d, design = design) # ,

















a1 = length(intersect(rhyLimma$gene_id[which(rhyLimma$adj.P.Val<p)], rhyLimma2$gene_id[which(rhyLimma2$adj.P.Val<p)]))
a2 = length(union(rhyLimma$gene_id[which(rhyLimma$adj.P.Val<p)], rhyLimma2$gene_id[which(rhyLimma2$adj.P.Val<p)]))
a3 = length(intersect(rhyLimma$gene_id[which(rhyLimma$P.Value <p)], rhyLimma2$gene_id[which(rhyLimma2$P.Value<p)]))
a4 = length(union(rhyLimma$gene_id[which(rhyLimma$P.Value <p)], rhyLimma2$gene_id[which(rhyLimma2$P.Value<p)]))
cl1 <- rep(0, nrow(rhyLimma))
cl1[which(rhyLimma$adj.P.Val<p)] <- 1
cl2 <- rep(0, nrow(rhyLimma2))
cl2[which(rhyLimma2$adj.P.Val<p)] <- 1
# a5 <- adj.rand.index(cl1, cl2)
cl1 <- rep(0, nrow(rhyLimma))
cl1[which(rhyLimma$P.Value<p)] <- 1
cl2 <- rep(0, nrow(rhyLimma2))
cl2[which(rhyLimma2$P.Value<p)] <- 1
#a6 <- adj.rand.index(cl1, cl2)
or <- c(a1/a2, a3/a4, a1, a2, a3, a4, 
        length(which(rhyLimma$adj.P.Val<p)),
        length(which(rhyLimma$P.Value<p)),
        length(which(rhyLimma2$adj.P.Val<p)),
        length(which(rhyLimma2$P.Value<p)))



# pathway enrichment analysis -- what biological function. Use functional class scoring
# gene ontology

# https://www.bioconductor.org/packages/devel/bioc/vignettes/annotate/inst/doc/GOusage.pdf

# 2-table comparison (Fisher's exact test) # usually a few hundredg genes selection
# -- encrichment analysis

# GSEA weighted version of KS-test -- high # indicates symbol (e.g., correlation with outcome, p-value, t-statistics)


