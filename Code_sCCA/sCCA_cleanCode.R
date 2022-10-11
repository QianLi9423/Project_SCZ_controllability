# clean code for sCCA analysis 2022.10.11
# 0. load packages ####
library(lattice)
library(ggplot2)
library(caret)
library(readxl)
library(cowplot)
library(parallel)
library(imputeR)
library(PMA)
library(emdbook)
library(R.matlab)
library(MASS)
library(permute)
library(matrixStats)
library(scales)
library(ggrepel)
library(raster)
library(sp)
library(rasterVis)
library(plyr)
library(psych)
library(reshape2)
library(CORElearn)
require(mgcv)
require(visreg)
library(circlize)
source("F:\\Projects_SCZ\\ccafunction.R")
set.seed(3456)


# 1. load data ####
# load controllability
netpath <- paste("F:\\Projects\\ave_PT.mat")
#netpath <- paste("F:\\Projects\\modal_PT.mat")
mat <- readMat(netpath)
SCH_fea<- mat$ave.PT.125
#SCH_fea<- mat$modal.PT.125

# load clinical data
infopath <- paste("F:\\Projects\\info_PT_gp_ns_ps.mat")
infomat <- readMat(infopath)
behavior.data <- infomat$info.gp.ns.ps
data <- list(brain = SCH_fea, behavior = behavior.data) # Note: check your data

# 2. obtain best regularization parameters ####
SCH_id <- 1:dim(behavior.data)[1]
sampleid <- createResample(SCH_id, list = T, times = 300)
brain_sample <- mclapply(sampleid, function(id) data$brain[id,])
behavior_sample <- mclapply(sampleid, function(id) data$behavior[id,])

x_pen <- seq(0.1,1,length.out=10)
y_pen <- seq(0.1,1,length.out=10)

scca.gs <- ccaDWfoldgs(brain_sample, behavior_sample, x_pen, y_pen)
gs.mat <- matrix(scca.gs$GS[, 'COR_MEAN'], nrow = 10, ncol = 10)
rownames(gs.mat) <- x_pen
colnames(gs.mat) <- y_pen

sample_corPlot <- corPlot(gs.mat, n = 100, zlim = c(min(gs.mat),max(gs.mat)), numbers = F, gr = colorRampPalette(c("#B52127", "white", "#2171B5")))

px <- scca.gs$PENX
py <- scca.gs$PENY


# 3. sCCA analysis ####
px<-0.5 #according to step2
py<-0.8 #according to step2
dim_num<-3
scca.entirety <- ccaDW(data$brain, data$behavior, px, py, dim_num) 
#sCCA overall correlation based on bootstrap
cor_SCH <- scca_corSE(sampleid, px, py, dim_num, scca.entirety, data, "Correlations")#
cor_SCH$plot
cor_SCH$cor.df

#permutation test
SCH_perm <- perm_scca(data, px, py, dim_num, scca.entirety, cor_SCH$cor.df)
data.pValue<-SCH_perm$pval
data.pValue_adjusted = p.adjust(data.pValue, method='bonferroni')
data.pValue_adjusted_fdr = p.adjust(data.pValue, method='fdr')
pass<-c(1,2) #according to data.pValue_adjusted

#obtain the features whose confidence intervals did not cross zero
SCH_id <- 1:dim(behavior.data)[1]
SCH_boot <- boot_scca(SCH_id, data, px, py, dim_num, scca.entirety, pass)

# 4. GAM model ####
SCH<-read_excel("F:\\Projects\\data_match_PT_gam.xlsx") #load age sex
SCH_net<-data$brain
fea.std <- apply(SCH_net,2,scale)
brain.df <- data.frame(barin_score=fea.std%*%scca.entirety$u[, SCH_perm$pass[1]],age = SCH$Age, sex = as.ordered(as.factor(SCH$Gender)))

brain_gam <- gam(barin_score~s(age)+sex, data = brain.df, method="REML")
gam_plot(brain_gam)

pval <- data.frame(age = summary(brain_gam)$s.table[,'p-value'] ,sex = summary(brain_gam)$p.table['sex.L','Pr(>|t|)'])
adj.pval<-p.adjust(pval,method = 'fdr')
gam_barin<-list(brain_gam=brain_gam, pval=pval, adj.pval=adj.pval)

gam_barin[[2]]
gam_barin[[3]]

# 5. Reliability analysis ####
# 5.1 reliability of overall sCCA correlation
# 5.1.1 mean and sd of overall sCCA correlation 
bl_dataset_impute<- cbind(SCH_fea,behavior.data)
df1<-bl_dataset_impute
temppsych = list()
temppsych$X1 = df1[,c(1:227)]
temppsych$Y1 = df1[,c(228:230)]
dim(temppsych$Y1)
set.seed(123)
perm=1:1000
rr_func=function(perm) {
  training<- sample(nrow(df1), nrow(df1), replace=T )
  dftraining<-df1[training,]
  temppsych = list()
  temppsych$X1 = as.matrix(dftraining[,c(1:227)])
  temppsych$Y1 = as.matrix(dftraining[,c(228:230)])
  result.rgcca = ccaDW(temppsych$X1, temppsych$Y1, px, py, 3)
  return(list(sort(result.rgcca$cors,decreasing = TRUE)))
}
results <- mclapply(perm, rr_func)
results_resample_cor<-array(0,c(1000,3))
for (i in 1:1000) {
  results_resample_cor[i,]<- results[[i]][[1]]
}
mean_results_resample_cor<- colMeans(abs(results_resample_cor));mean_results_resample_cor
Sds_results_resample_cor<- colSds(abs(results_resample_cor));Sds_results_resample_cor

# 5.1.2 sample size effect
data <- list(brain = SCH_fea, behavior = behavior.data) # Note: check your data
df1<- cbind(SCH_fea,behavior.data) #added by LQ
#df1<-scale(df1)
temppsych = list()
temppsych$X1 = df1[,c(1:227)]
temppsych$Y1 = df1[,c(228:230)]
px<-0.5
py<-0.8
set.seed(123)
perm=1:1000
boot_func<-function(perm) {
  cor_boot_thickness_baseline<-matrix(0,15,1)
  sampleseq<- c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5)
  for (j in sampleseq) {
    training<- sample(nrow(df1), ceiling(j*nrow(df1)), replace=T )
    temppsych = list()
    temppsych$X1 = df1[training,c(1:227)]
    temppsych$Y1 = df1[training,c(228:230)]
    rownames(temppsych$X1)<-seq(1,ceiling(j*nrow(df1)))
    rownames(temppsych$Y1)<-seq(1,ceiling(j*nrow(df1)))
    result.rgcca <- ccaDW(temppsych$X1, temppsych$Y1, px, py, 1)
    
    cor_boot_thickness_baseline[which(sampleseq==j),]<-result.rgcca$cors
  }
  return(cor_boot_thickness_baseline)
}

results <- mclapply(perm, boot_func)
absresults<-lapply(results,abs)
thickness_baseline_resample<-t(matrix(unlist(absresults),15,1000))

bb<-apply(thickness_baseline_resample,2,mean);bb
plotdf<- melt(t(thickness_baseline_resample),id='thickness_baseline_resample') 
sampleseq<- rep(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5),1000)
plotdf_sample <- data.frame(sampleseq = sampleseq,cor_sample_Mean = plotdf$value)

# 5.2 for u v reliatbility 
# 5.2.1 mean and sd for u and v
SCH_boot <- boot_scca(SCH_id, data, px, py, dim_num, scca.entirety, 1)
SCH_boot$u;SCH_boot$v;
bottv_mean<-apply(abs(matrix(unlist(SCH_boot$scca.boot.v),nrow=3)),1,mean)
bottu_mean<-apply(abs(matrix(unlist(SCH_boot$scca.boot.u),nrow=227)),1,mean)
bottu_mean[c(7,11,13,14,18,20,23,26,38,64,65,75,92,145,164,171,188,199,218,219,221,222)] # obtain means and sd for significant features
bottv_sd<-apply(abs(matrix(unlist(SCH_boot$scca.boot.v),nrow=3)),1,sd)
bottu_sd<-apply(abs(matrix(unlist(SCH_boot$scca.boot.u),nrow=227)),1,sd)
bottu_sd[c(7,11,13,14,18,20,23,26,38,64,65,75,92,145,164,171,188,199,218,219,221,222)]

org<-scca.entirety$v[,1]; boot<-SCH_boot$scca.boot.v[[1]]
btst_v <- bootstats2(org,boot,0.975,0.025);btst_v
org<-scca.entirety$u[,1]; boot<-SCH_boot$scca.boot.u[[1]]
btst_u <- bootstats2(org,boot,0.995,0.005);btst_u
btst_u[c(7,11,13,14,18,20,23,26,38,64,65,75,92,145,164,171,188,199,218,219,221,222),]

# 5.2.2 leave one-out analysis for u and v
df1<- cbind(SCH_fea,behavior.data) #added by LQ
temppsych = list()
temppsych$X1 = df1[,c(1:227)]
temppsych$Y1 = df1[,c(228:230)]
u.cor<-matrix(0,dim(df1)[1],3)
v.cor<- matrix(0,dim(df1)[1],3)

for (i in 1:dim(df1)[1]){
  trainsamp<- i
  testlist = list()
  testlist$X1 = temppsych$X1[trainsamp, ]
  testlist$Y1 = temppsych$Y1[trainsamp, ]

  trainlist = list()
  trainlist$X1 = temppsych$X1[-trainsamp, ]
  trainlist$Y1 = temppsych$Y1[-trainsamp, ]
  result.rgcca = ccaDW(trainlist$X1, trainlist$Y1, px, py, 3)
  u.cor[i,]<-diag(cor(result.rgcca$u,u))
  v.cor[i,]<-diag(cor(result.rgcca$v,v))
}

mean.u.cor<-colMeans(u.cor);mean.u.cor
mean.v.cor<-colMeans(v.cor);mean.v.cor