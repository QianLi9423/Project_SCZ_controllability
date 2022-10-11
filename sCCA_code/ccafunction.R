#source function for sCCA analysis
bootstats2 <- function(org,bootdata, uplim,lwlim){
  boot.stats <- sapply(seq_along(1:length(org)), function(i) org[i] - quantile(bootdata[i,] - org[i],c(uplim,lwlim),na.rm =T))
  boot.stats <- as.data.frame(t(boot.stats))
  colnames(boot.stats)<-c("low","high")
  boot.stats$ci <- boot.stats$high - boot.stats$low
  boot.stats$load <- org
  boot.stats$fea <- 1:length(org)
  boot.stats
  }

ccaDWpermorder <- function(X,Y,pen_x,pen_y,rank,cca_org){
  perm.mode<-PMA::CCA(x=X, z=Y, typex=c("standard"),typez=c("standard"), penaltyx=pen_x, penaltyz=pen_y, K=rank, niter=20, trace=FALSE)
  perm.mode.reorder<-reorderCCA(perm.mode,cca_org,rank)
}

ccaDW <- function(X,Y,pen_x,pen_y,rank){
  mode<-PMA::CCA(x=X, z=Y, typex=c("standard"),typez=c("standard"), penaltyx=pen_x, penaltyz=pen_y, K=rank, niter=150, trace=FALSE)
}

ccaDWfoldgs <-function(Xlist,Ylist,pen_xseq,pen_yseq){
  #loop through x penalty sequence
  gs_x<-function(X,Y,pen_xseq,pen_yseq){
    gs_x <- mclapply(pen_xseq,function(pen_x) ccaDWfold(X,Y,pen_x,pen_yseq))
  }
  #loop through y penalty sequence
  gs_y <- mclapply(pen_yseq,function(pen_y) gs_x(Xlist,Ylist,pen_xseq,pen_y))
  
  #summarize the result
  gs2 <- unlist(gs_y,recursive = F)
  gs.out <- t(sapply(gs2, function(x) unlist(x)))
  
  #find the parameters with the largest correlations
  best.para <- gs.out[which.max(gs.out[,'COR_MEAN']),]
  best.cor <- as.numeric(best.para[1])
  best.penx <- as.numeric(best.para[2])
  best.peny <- as.numeric(best.para[3])
  out <-list(GS = gs.out,COR= best.cor, PENX = best.penx, PENY = best.peny)
}

ccaDWfold <- function(Xlist,Ylist,pen_x,pen_y){
  result.cca<-mclapply(seq_along(Xlist),function(i) ccaDW(Xlist[[i]],Ylist[[i]],pen_x,pen_y,1))
  result.cor <- sapply(seq_along(Xlist), function(i) result.cca[[i]]$cors)
  result.cor.mean <- mean(result.cor)
  out <-list(COR_MEAN = result.cor.mean, PEN_X= pen_x, PEN_Y = pen_y)
}

perm.plot <-function(perm_file,cor_file,pval,dim){
  perm_file <- data.frame(cor_perm = perm_file[,dim])
  p <- ggplot(perm_file,aes(cor_perm))+
    geom_histogram(binwidth = 0.005, fill = "blue", alpha = 0.5)  +
    geom_vline(xintercept = cor_file$cor[dim], colour = "red", linetype = "longdash") +
    labs(x = "Correlations") +
    annotate("text", x = median(perm_file$cor_perm,na.rm = T), y = c(10,5),label = c("Permuted Data","(1000 times)"),size =10,colour = "black" ) +
    annotate("text",x = -Inf, y = Inf, hjust = -0.1,vjust = 1,label = paste("Cor=",round(cor_file$cor[dim],2),", p<",round(pval[dim],3)), size = 10, colour = "red" ) +
    theme_classic(base_size = 35) + 
    scale_y_continuous(expand = c(0, 0))
  p
}


reorderCCA <- function(res,org,k){
  u <- org$u #base selected 
  v <- org$u
  cor.match<-abs(cor(v,res$u[,1:k]))

  order.match<-order(-abs(org$cors))

  res.match <- 0
  res.cor<-0
  for(i in 1:k){
      res.match[order.match[i]]<-which.max(cor.match[order.match[i],])
      res.cor[order.match[i]]<-max(cor.match[order.match[i],])
      cor.match[,res.match[order.match[i]]]<- -1
  }

  u.od <- res$u[,res.match]
  v.od <- res$v[,res.match]
  cors.final <- res$cors[res.match]
  
  u.org<-org$u[,res.match]
  v.org<-org$v[,res.match]

  for(i in 1:length(res.match)){
      a<-0
      for(j in 1:dim(org$u)[1]){
        a<-a+u.od[j,i]*u.org[j,i]
      }
      if(a<0){
        u.od[,i]<-(-1)*u.od[,i]
        v.od[,i]<-(-1)*v.od[,i]
      }
  }
  
  u.final <- u.od
  v.final <- v.od
  
  res.one.reorder <- list(u= u.final ,v= v.final, cors = cors.final, pos = res.match, dimcor = res.cor)

  out <- res.one.reorder
}


scca_corSE<-function(sampleid,px,py,dim_num,scca.entirety,data,title){
 brain_sample <- mclapply(sampleid, function(id) data$brain[id,])
 behavior_sample <- mclapply(sampleid, function(id) data$behavior[id,])
 scca.single  <- mclapply(seq_along(sampleid),function(i) ccaDW(brain_sample[[i]],behavior_sample[[i]],px,py,dim_num))
 scca.cca.ro <- sapply(scca.single,function(x) reorderCCA(x,scca.entirety,dim_num))
 scca.cca.cor <- (rowMeans(simplify2array(scca.cca.ro['cors',]),na.rm =T))
 scca.cca.cor.se <- (rowSds(simplify2array(scca.cca.ro['cors',]),na.rm =T)/sqrt(dim(scca.cca.ro)[2]))

 cor.df <- data.frame(modenum = as.factor(1:dim_num), cor = scca.cca.cor, se = scca.cca.cor.se)
 cor.lim <- aes(ymax = cor + se, ymin = cor - se)

 p.cor <- ggplot(cor.df,aes(1:length(modenum), cor, label = round(cor,2))) +
  geom_bar(width = 0.75, stat = 'identity',  fill = c(rep("red",dim_num)),alpha=0.6) +
  geom_errorbar(cor.lim,  width=0.25) +
  geom_text(size = 4, position = position_dodge(width = 0.9), vjust=-1,color='gray')+
  scale_x_discrete(name ="Mode", limits=as.character(c(1:dim_num)) ) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,1),name = title, breaks=seq(0,1,length=5)) +
  theme_classic(base_size = 20) +
  coord_cartesian(ylim=c(0.2,1)) +
  theme(legend.position = 'none')

 out<-list(cor.df=cor.df,plot=p.cor)
}

perm_scca<-function(data,px,py,dimnum,scca.entirety,cor.df){

 num.perm <- 1000
 behavior.perm <- rlply(num.perm,data$behavior[sample(nrow(data$behavior)),])
 scca.perm.cca<-sapply(behavior.perm, function(y_perm){ out<-ccaDWpermorder(data$brain,y_perm,px,py,dimnum,scca.entirety)} )  
 perm.cor <- simplify2array(scca.perm.cca['cors',])
 perm.pval <- sapply(seq_along(cor.df$cor),function(x) (length(which(perm.cor[x,] >= cor.df$cor[x])) ) / length(which(is.na(perm.cor[x,]) == FALSE)))
 perm.cor.df <- as.data.frame(t(perm.cor))
 perm.pass <- which(perm.pval < 0.05)
 permplots <-lapply(perm.pass,function(x) perm.plot(perm.cor.df,cor.df,perm.pval,x))

 out<-list(pass=perm.pass,pval=perm.pval,plot=permplots)
}


boot_u <- function(org, boot){
  btst <- bootstats2(org,boot,0.995,0.005)
  
  fea <- which(btst$low * btst$high >0)
  load <- org[fea]
  
  p <- ggplot(btst,aes(fea,load))+
    geom_point(aes(colour = low * high > 0)) +
    geom_errorbar(aes(ymax = high, ymin = low),  width=0.25) +
    theme_classic()

  out <- list(fea = fea,load = load,plot = p)
}

boot_v <- function(org, boot){
  btst <- bootstats2(org,boot,0.975,0.025)
  
  p <- ggplot(btst,aes(fea,load))+
    geom_point(aes(colour = low * high > 0)) +
    geom_errorbar(aes(ymax = high, ymin = low),  width=0.25) +
    theme_classic()
  
  fea <- which(btst$low * btst$high >0)
  load <- org[fea]
  ##fea <- fea[order(-abs(load))]
  ##load <- load[order(-abs(load))]
  out <- list(plot = p, fea = fea, load = load)
}

#SCH_boot <- boot_scca(SCH_id, data, px, py, dim_num, scca.entirety, SCH_perm$pass)
boot_scca<-function(tscore,data,px,py,dimnum,scca.entirety,perm.pass){
 bootnum = 1000
 bootid<-createResample(tscore, list = T, times = bootnum) #create lists of subjs for bootstrap

 brain_boot <- lapply(bootid, function(id) data$brain[id,]) #creat bootstrap samples for connectivity features
 behavior_boot <- lapply(bootid, function(id) data$behavior[id,]) #creat bootstrap samples for clinical features 

 scca.boot<- mclapply(seq_along(bootid),function(i) ccaDW(brain_boot[[i]],behavior_boot[[i]],px,py,dimnum),mc.cores = 1) #run scca on these bootstrap sample

 scca.boot.ro<- lapply(1:bootnum,function(i) reorderCCA(scca.boot[[i]],scca.entirety,dimnum)) #reorder to match components across samples

 scca.boot.u <- lapply(perm.pass, function(x) sapply(1:bootnum, function(i) scca.boot.ro[[i]]$u[,x])) #extract loadings on connectivity features
 scca.boot.v <- lapply(perm.pass, function(x) sapply(1:bootnum, function(i) scca.boot.ro[[i]]$v[,x])) #extract loadings on clinical features

 u.boot <- lapply(seq_along(perm.pass), function(x) boot_u(scca.entirety$u[,perm.pass[x]], scca.boot.u[[x]] ))
 v.boot <- lapply(seq_along(perm.pass), function(x) boot_v(scca.entirety$v[,perm.pass[x]], scca.boot.v[[x]] ))

 out<-list(u=u.boot,v=v.boot,scca.boot.v=scca.boot.v,scca.boot.u=scca.boot.u)

}

gam_scca <- function(fea,loadings,SCH){
 fea.std <- apply(fea,2,scale)
 brain.df <- data.frame(barin_score=fea.std%*%loadings,age = SCH$age, sex = as.ordered(as.factor(SCH$gender)))
 
 brain_gam <- gam(barin_score~s(age)+sex, data = brain.df, method="REML")

 pval <- data.frame(age = summary(brain_gam)$s.table[,'p-value'] ,sex = summary(brain_gam)$p.table['sex.L','Pr(>|t|)'])
 adj.pval<-p.adjust(pval,method = 'bonferroni')

 out<-list(brain_gam=brain_gam, pval=pval, adj.pval=adj.pval) 
}

gam_plot <- function(gammodel){
 op<-par(mfrow = c(1, 2))
 visreg(gammodel,ylab="brain_score")
 par(op)
}

scca_fea.con <- function(loadings,idx,boot_fea){
 u_loading.relief <- loadings
 u_loading.relief[-boot_fea] <- 0
 uloading_mat <- array(0,c(264,264))
 uloading_mat[upper.tri(uloading_mat,diag=F)][idx] <- u_loading.relief
 uloading_mat <- sna::symmetrize(uloading_mat,rule = "upper")
 uloading_vevtor <- uloading_mat[upper.tri(uloading_mat, diag = F)]
 fea <- which(uloading_vevtor!=0)
 sccaloadings <- uloading_vevtor[fea]

 out<-list(uloading_vevtor=uloading_vevtor,uloading_mat=uloading_mat,fea=fea,sccaloadings=sccaloadings)
}

scca_fea <- function(loadings,idx){
 uloading_mat <- array(0,c(264,264))
 uloading_mat[upper.tri(uloading_mat,diag=F)][idx] <- loadings
 uloading_mat <- sna::symmetrize(uloading_mat,rule = "upper")
 uloading_vevtor <- uloading_mat[upper.tri(uloading_mat, diag = F)]
 fea <- which(uloading_vevtor!=0)
 sccaloadings <- uloading_vevtor[fea]

 out<-list(uloading_vevtor=uloading_vevtor,uloading_mat=uloading_mat,fea=fea,sccaloadings=sccaloadings)
}