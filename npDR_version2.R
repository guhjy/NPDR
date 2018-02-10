# load libraries/functions
args=(commandArgs(TRUE))
#setwd(".")

packages <- c("foreach","doParallel","doRNG","boot",
              "rmutil","mvtnorm","gam","sandwich",
              "devtools","glmnet","data.table","rpart",
              "ranger","nnet","arm","earth","e1071","xgboost",
              "foreach")
#userLib <-  "~/R/R_LIBS_USER"
#.libPaths(userLib)

#ifelse(!dir.exists(userLib), dir.create(userLib), FALSE)

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN') #lib=userLib, 
  }
}

for (package in packages) {
  library(package, character.only=T)#, lib.loc=.libPaths()
}

devtools::install_github("ecpolley/SuperLearner") #,lib=userLib
library(SuperLearner)#, lib.loc=.libPaths()

# https://stackoverflow.com/questions/12872699/error-unable-to-load-installed-packages-just-now/25932828#25932828
install.packages("rJava",repos='http://lib.stat.cmu.edu/R/CRAN')#lib=userLib, 
library(rJava)#, lib.loc=.libPaths()

devtools::install_github("ck37/ck37r")#,lib=userLib
library(ck37r)#, lib.loc=.libPaths()

install.packages("tmle",repos='http://lib.stat.cmu.edu/R/CRAN')#lib=userLib,
library(tmle)

time<-proc.time()
print(args)
set.seed(as.numeric(args[[1]]))
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }
nsim <- as.numeric(args[[2]])

cols <- c("ipwPMT","ipwNPT","ipwPMF","ipwNPF",
          "regPMT","regNPT","regPMF","regNPF",
          "aipwPMT","aipwNPT","aipwPMF","aipwNPF",
          "tmlePMT","tmleNPT","tmlePMF","tmleNPF")
res.est <- data.frame(matrix(nrow=1,ncol=length(cols))); #nsim*5
colnames(res.est) <- cols
res.se <- res.est

##  true value
true<-6;true

source("~/Dropbox/SLwrappers/create.SL.glmnet.R")
source("~/Dropbox/SLwrappers/create.SL.gam.Wrapper.R")
source("~/Dropbox/SLwrappers/rangerWrapper.R")
source("~/Dropbox/SLwrappers/rangerWrapper1.R")
source("~/Dropbox/SLwrappers/rangerWrapper2.R")
source("~/Dropbox/SLwrappers/rangerWrapper3.R")
source("~/Dropbox/SLwrappers/rangerWrapper4.R")
source("~/Dropbox/SLwrappers/rangerWrapper5.R")
source("~/Dropbox/SLwrappers/create.SL.xgboost.R")

tune = list(ntrees = c(100, 500), max_depth = c(1, 2), minobspernode = 10,
            shrinkage = c(0.1, 0.01, 0.001))
xgb_grid = create.SL.xgboost(tune = tune) #, env = sl_env
create.SL.glmnet(alpha = c(0,.5)) #setting alpha = 0 gives LASSO, 1 gives elastic net.
create.SL.gam(deg.gam = c(3,4,5)) 
sl.lib<-c("SL.glmnet.0","SL.glmnet.0.5","SL.glmnet","SL.rpartPrune",
          "SL.gam.3","SL.gam.4","SL.gam.5","SL.glm.interaction","SL.earth",
          "SL.ranger","SL.ranger1","SL.ranger2","SL.ranger3","SL.ranger4","SL.ranger5",
          "SL.bayesglm","SL.xgboost","SL.mean",xgb_grid$names)

sl.lib<-c("SL.glmnet.0","SL.glmnet.0.5","SL.glmnet","SL.rpartPrune",
          "SL.gam.3","SL.gam.4","SL.gam.5","SL.earth",
          "SL.ranger","SL.ranger1","SL.ranger2",
          "SL.bayesglm","SL.mean")

npDR<-function(counter,bs=F,bootNum=100){
  # data management
  i<-counter
  samp<-rep(1:5,nsim*5)[counter]
  ss<-c(100,200,600,1200,2000)
  n<-ss[samp]
  cat("Now running iteration",i,"with a sample size of",n,'\n');flush.console()
  
  # confounders
  sigma<-matrix(0,nrow=4,ncol=4);diag(sigma)<-4
  x <- rmvnorm(n, mean=rep(0,4), sigma=sigma)
  
  #GGally::ggpairs(data.frame(x))
  
  z<-x
  z[,1]<-exp(x[,1]/2)
  z[,2]<-x[,2]/(1+exp(x[,1]))+10
  z[,3]<-(x[,1]*x[,3]/25+.6)^3
  z[,4]<-(x[,2]*x[,4]+20)^2
  
  #GGally::ggpairs(data.frame(z))
  
  # design matrix for exposure and outcome model
  dMatT<-model.matrix(as.formula(paste("~(",paste("x[,",1:ncol(x),"]",collapse="+"),")"))) ## removed all interactions
  
  beta<-c(120,3.5,2.5,-1,3.5) # ,2,2.5,1.5,1.5,1.5,1
  theta<-c(-2,log(2),log(2.5),log(.5),log(1.5)) # ,log(1.75),log(1.5),log(1.25),log(1.25),log(1.25),log(1.25)
  
  # design matrix for propensity score model
  mu <- dMatT%*%beta
  # propensity score model
  pi <- expit(dMatT%*%theta);r<-rbinom(n,1,pi)
  # outcome model: true expsoure effect = 6
  y <- r*6 + mu + rnorm(n,0,60)

  # induce misspecification
  dMat<-model.matrix(as.formula(paste("~(",paste("z[,",1:ncol(x),"]",collapse="+"),")")))
  
  # correct specification
  dat <- data.frame(x,r,y); colnames(dat)[1:ncol(x)] <- paste("x",1:ncol(x),sep="")
  W <- dMatT[,-1];colnames(W) <- paste("W",1:ncol(dMatT[,-1]), sep="");A <- r;Y <- y
  tQForm<-as.formula(paste0("Y~A+", paste(paste0("W",1:ncol(dMatT[,-1])), collapse="+")))
  tgForm<-as.formula(paste0("A~", paste(paste0("W",1:ncol(dMatT[,-1])), collapse="+")))
  
  tmlePMT <- tmle(Y,A,W,family="gaussian",Qform=tQForm,gform=tgForm)
  
  folds<-c(2,2,3,5,5)[samp]
  cat("Number of cross-validation folds is",folds,'\n');flush.console()
  
  tmleNPT <- tmle(Y,A,data.frame(W),family="gaussian",Q.SL.library=sl.lib,g.SL.library=sl.lib)

  # setup_parallel_tmle(parallel="multicore",max_cores=20)
  # set.seed(1, "L'Ecuyer-CMRG")
  # tmleNPT <- tmle_parallel(Y=Y,A=A,W=data.frame(W),
  #                          family="gaussian",V=folds,
  #                          Q.SL.library=sl.lib,g.SL.library=sl.lib)

  # SL_out.monT<-data.table(N=n,model="outcome",t(tmleNPT$Qinit$coef))
  # SL_exp.monT<-data.table(N=n,model="exposure",t(tmleNPT$g$coef))
  # 
  # if(i==1&samp==1){
  #   write.table(SL_out.monT,"SL_outDat_True.txt",sep="\t",row.names=F)
  #   write.table(SL_exp.monT,"SL_expDat_True.txt",sep="\t",row.names=F)
  # } else{
  #   write.table(SL_out.monT,"SL_outDat_True.txt",sep="\t",row.names=F,col.names=F,append=T)
  #   write.table(SL_exp.monT,"SL_expDat_True.txt",sep="\t",row.names=F,col.names=F,append=T)
  # }
  
  pihatPT <- tmlePMT$g$g1W
  #swPT<-dat$r*(mean(dat$r)/tmlePMT$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmlePMT$g$g1W))
  pihatNPT <- tmleNPT$g$g1W
  #swNPT<-dat$r*(mean(dat$r)/tmleNPT$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmleNPT$g$g1W))
  muhatPT  <- tmlePMT$Qinit$Q[,2]*A+tmlePMT$Qinit$Q[,1]*(1-A)
  muhatPT1 <- tmlePMT$Qinit$Q[,2];muhatPT0 <- tmlePMT$Qinit$Q[,1]
  muhatNPT  <- tmleNPT$Qinit$Q[,2]*A+tmleNPT$Qinit$Q[,1]*(1-A)
  muhatNPT1 <- tmleNPT$Qinit$Q[,2];muhatNPT0 <- tmleNPT$Qinit$Q[,1]
  
  # misspecified
  dat <- data.frame(z,r,y); colnames(dat)[1:ncol(z)] <- paste("z",1:ncol(x),sep="")
  W <- dMat[,-1];colnames(W) <- paste("W",1:ncol(dMat[,-1]), sep="");A <- r;Y <- y
  tQForm<-as.formula(paste0("Y~A+", paste(paste0("W",1:ncol(dMat[,-1])), collapse="+")))
  tgForm<-as.formula(paste0("A~", paste(paste0("W",1:ncol(dMat[,-1])), collapse="+")))
  tmlePMF <- tmle(Y,A,W,family="gaussian",Qform=tQForm,gform=tgForm)
  
  folds<-c(2,2,3,5,5)[samp]
  tmleNPF <- tmle(Y,A,data.frame(W),family="gaussian",Q.SL.library=sl.lib,g.SL.library=sl.lib)
  
  # SL_out.monF<-data.table(N=n,model="outcome",t(tmleNPF$Qinit$coef))
  # SL_exp.monF<-data.table(N=n,model="exposure",t(tmleNPF$g$coef))
  # 
  # if(i==1&samp==1){
  #   write.table(SL_out.monF,"SL_outDat_MisSpec.txt",sep="\t",row.names=F)
  #   write.table(SL_exp.monF,"SL_expDat_MisSpec.txt",sep="\t",row.names=F)
  # } else{
  #   write.table(SL_out.monF,"SL_outDat_MisSpec.txt",sep="\t",row.names=F,col.names=F,append=T)
  #   write.table(SL_exp.monF,"SL_expDat_MisSpec.txt",sep="\t",row.names=F,col.names=F,append=T)
  # }
  
  pihatPF <- tmlePMF$g$g1W
  #swPF<-dat$r*(mean(dat$r)/tmlePMF$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmlePMF$g$g1W))
  pihatNPF <- tmleNPF$g$g1W
  #swNPF<-dat$r*(mean(dat$r)/tmleNPF$g$g1W) + (1-dat$r)*((1-mean(dat$r))/(1-tmleNPF$g$g1W))
  muhatPF  <- tmlePMF$Qinit$Q[,2]*A+tmlePMF$Qinit$Q[,1]*(1-A)
  muhatPF1 <- tmlePMF$Qinit$Q[,2];muhatPF0 <- tmlePMF$Qinit$Q[,1]
  muhatNPF  <- tmleNPF$Qinit$Q[,2]*A+tmleNPF$Qinit$Q[,1]*(1-A)
  muhatNPF1 <- tmleNPF$Qinit$Q[,2];muhatNPF0 <- tmleNPF$Qinit$Q[,1]
  
  # compute estimators
  res.est$ipwPMT <- sum((dat$r*dat$y/pihatPT)/sum(dat$r/pihatPT)) -  sum((1-dat$r)*dat$y/(1-pihatPT))/sum((1-dat$r)/(1-pihatPT))
  res.est$ipwNPT <- sum((dat$r*dat$y/pihatNPT)/sum(dat$r/pihatNPT)) -  sum((1-dat$r)*dat$y/(1-pihatNPT))/sum((1-dat$r)/(1-pihatNPT))
  res.est$regPMT <- mean(muhatPT1 - muhatPT0)
  res.est$regNPT <- mean(muhatNPT1 - muhatNPT0)
  
  res.est$aipwPMT  <- mean((((2*dat$r-1)*(dat$y - muhatPT))/((2*dat$r-1)*pihatPT + (1-dat$r)) + muhatPT1 - muhatPT0))
  res.est$aipwNPT  <- mean((((2*dat$r-1)*(dat$y - muhatNPT))/((2*dat$r-1)*pihatNPT + (1-dat$r)) + muhatNPT1 - muhatNPT0))
  
  res.est$tmlePMT <- tmlePMT$estimates$ATE$psi
  res.est$tmleNPT <- tmleNPT$estimates$ATE$psi
  
  res.est$ipwPMF <- sum((dat$r*dat$y/pihatPF)/sum(dat$r/pihatPF)) -  sum((1-dat$r)*dat$y/(1-pihatPF))/sum((1-dat$r)/(1-pihatPF))
  res.est$ipwNPF <- sum((dat$r*dat$y/pihatNPF)/sum(dat$r/pihatNPF)) -  sum((1-dat$r)*dat$y/(1-pihatNPF))/sum((1-dat$r)/(1-pihatNPF))
  res.est$regPMF <- mean(muhatPF1 - muhatPF0)
  res.est$regNPF <- mean(muhatNPF1 - muhatNPF0)
  
  res.est$aipwPMF  <- mean((((2*dat$r-1)*(dat$y - muhatPF))/((2*dat$r-1)*pihatPF + (1-dat$r)) + muhatPF1 - muhatPF0))
  res.est$aipwNPF  <- mean((((2*dat$r-1)*(dat$y - muhatNPF))/((2*dat$r-1)*pihatNPF + (1-dat$r)) + muhatNPF1 - muhatNPF0))
  
  res.est$tmlePMF<-tmlePMF$estimates$ATE$psi
  res.est$tmleNPF<-tmleNPF$estimates$ATE$psi
  
  # compute closed-form SEs
  #mod1<-lm(y~r,data=dat,weights=swPT)
  #mod2<-lm(y~r,data=dat,weights=swNPT)
  res.se$ipwPMT <- sd((dat$r*dat$y/pihatPT)/sum(dat$r/pihatPT) -  ((1-dat$r)*dat$y/(1-pihatPT))/sum((1-dat$r)/(1-pihatPT)))
  res.se$ipwNPT <- sd((dat$r*dat$y/pihatNPT)/sum(dat$r/pihatNPT) -  ((1-dat$r)*dat$y/(1-pihatNPT))/sum((1-dat$r)/(1-pihatNPT)))
  res.se$aipwPMT <- sd((((2*dat$r-1)*(dat$y - muhatPT))/((2*dat$r-1)*pihatPT + (1-dat$r)) + muhatPT1 - muhatPT0))
  res.se$aipwNPT <- sd((((2*dat$r-1)*(dat$y - muhatNPT))/((2*dat$r-1)*pihatNPT + (1-dat$r)) + muhatNPT1 - muhatNPT0))
  res.se$tmlePMT <- sqrt(tmlePMT$estimates$ATE$var.psi)*sqrt(n)
  res.se$tmleNPT <- sqrt(tmleNPT$estimates$ATE$var.psi)*sqrt(n)
  
  #mod1<-lm(y~r,data=dat,weights=swPF)
  #mod2<-lm(y~r,data=dat,weights=swNPF)
  res.se$ipwPMF <- sd((dat$r*dat$y/pihatPF)/sum(dat$r/pihatPF) -  ((1-dat$r)*dat$y/(1-pihatPF))/sum((1-dat$r)/(1-pihatPF)))
  res.se$ipwNPF <- sd((dat$r*dat$y/pihatNPF)/sum(dat$r/pihatNPF) -  ((1-dat$r)*dat$y/(1-pihatNPF))/sum((1-dat$r)/(1-pihatNPF)))
  res.se$aipwPMF <- sd((((2*dat$r-1)*(dat$y - muhatPF))/((2*dat$r-1)*pihatPF + (1-dat$r)) + muhatPF1 - muhatPF0))
  res.se$aipwNPF <- sd((((2*dat$r-1)*(dat$y - muhatNPF))/((2*dat$r-1)*pihatNPF + (1-dat$r)) + muhatNPF1 - muhatNPF0))
  res.se$tmlePMF <- sqrt(tmlePMF$estimates$ATE$var.psi)*sqrt(n)
  res.se$tmleNPF <- sqrt(tmleNPF$estimates$ATE$var.psi)*sqrt(n)
  
  # compute bootstrap CIs
  datB<- data.frame(z,x,r,y); colnames(datB)[1:(2*ncol(z))] <- c(paste("z",1:ncol(x),sep=""),paste("x",1:ncol(x),sep=""))
  cat("Bootstrapping Parametric",'\n');flush.console()
  plugin1 <- function(d,j,dim){
    x<-d[j,c(paste0("x",1:4))]
    mMat<-model.matrix(as.formula(paste("~(",paste("x[,",1:dim,"]",collapse="+"),")")))[,-1]
    dat<-data.table(y=d[j,]$y,r=d[j,]$r,mMat);mumod<-glm(y~.,data=dat)
    dat$r<-1;muhat1 <- predict(mumod,newdata=dat)
    dat$r<-0;muhat0 <- predict(mumod,newdata=dat)
    meanPT<-mean(muhat1-muhat0)
    x<-d[j,c(paste0("z",1:4))]
    mMat<-model.matrix(as.formula(paste("~(",paste("x[,",1:dim,"]",collapse="+"),")")))[,-1]
    dat<-data.table(y=d[j,]$y,r=d[j,]$r,mMat);mumod<-glm(y~.,data=dat)
    dat$r<-1;muhat1 <- predict(mumod,newdata=dat)
    dat$r<-0;muhat0 <- predict(mumod,newdata=dat)
    meanPF<-mean(muhat1-muhat0)
    res<-c(meanPT,meanPF)
  }
  bs1 <- boot(datB,plugin1,R=bootNum,dim=4)
  res.se$regPMT <- sd(bs1$t[,1])*sqrt(n)
  res.se$regPMF <- sd(bs1$t[,2])*sqrt(n)
  if (bs==T){
    cat("Bootstrapping Nonparametric",'\n');flush.console()
    plugin2<-function(d,dim,f){
      j<-sample(1:nrow(d),nrow(d),replace=T);dat<-d[j,]
      x<-dat[,c(paste0("x",1:p))]
      W<-model.matrix(as.formula(paste("~(",paste("x[,",1:dim,"]",collapse="+"),")")))[,-1]
      Y<-dat$y;A<-as.matrix(dat$r) 
      X=data.frame(cbind(A,W))
      names(X)<-c("A",paste0("x",1:ncol(W)))
      mumod<-SuperLearner(Y,X,family=gaussian,SL.library=sl.lib,cvControl=list(V=f))
      X$A<-1;muhat1 <- predict(mumod,newdata=X,onlySL=T)$pred
      X$A<-0;muhat0 <- predict(mumod,newdata=X,onlySL=T)$pred
      gT<-mean(muhat1-muhat0)
      x<-dat[,c(paste0("z",1:p))]
      W<-model.matrix(as.formula(paste("~(",paste("x[,",1:dim,"]",collapse="+"),")")))[,-1]
      Y<-dat$y;A<-as.matrix(dat$r) 
      X=data.frame(cbind(A,W))
      names(X)<-c("A",paste0("x",1:ncol(W)))
      mumod<-SuperLearner(Y,X,family=gaussian,SL.library=sl.lib,cvControl=list(V=f))
      X$A<-1;muhat1 <- predict(mumod,newdata=X,onlySL=T)$pred
      X$A<-0;muhat0 <- predict(mumod,newdata=X,onlySL=T)$pred
      gF<-mean(muhat1-muhat0)
      g<-c(gT,gF)
      return(g)
    }
    coreNum<-detectCores()
    cl <- makeCluster(coreNum)
    registerDoParallel(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    clusterExport(cl, sl.lib)
    pak<-c("tmle","SuperLearner","glmnet","rpart","ranger","nnet","arm","earth")
    results<-foreach(i=1:bootNum,.packages=pak) %dorng% {
      b<-plugin2(dat,dim=p,f=folds);b
    }
    stopCluster(cl)
    results<-do.call(rbind,results)
    res.se$regNPT <- sd(results[,1])*sqrt(n)
    res.se$regNPF <- sd(results[,2])*sqrt(n)
  } else{
    res.se$regNPT <- res.se$regNPF <- NA
  }
  
  res.cov <- res.est-1.96*res.se/sqrt(n) < true & true < res.est+1.96*res.se/sqrt(n)
  res.width <- (res.est+1.96*res.se/sqrt(n)) - (res.est-1.96*res.se/sqrt(n))
  
  tmp<-data.frame(rbind(res.est-true,(res.est-true)^2,res.cov,res.width))
  tmp<-cbind(n,tmp)
  rownames(tmp)<-c("bias","rmse","cov","width");colnames(tmp)[1]<-"N";
  cat('\n')
  cat("Printing Results",'\n')
  print(round(tmp,3));cat('\n');flush.console()
  setDT(tmp, keep.rownames = TRUE)[];colnames(tmp)[1] <- "type"
  
  if(i==1&samp==1){
    write.table(tmp,"results.txt",sep="\t",row.names=F)
  } else{
    write.table(tmp,"results.txt",sep="\t",row.names=F,col.names=F,append=T)
  }
  return(tmp)
}


cores<-20 #detectCores()
#print(cores)
#results<-lapply(1:(nsim*5), function(x) npDR(x,bs=F))
results<-mclapply(1:(nsim*5), function(x) npDR(x,bs=F),mc.cores=cores)
# results<-do.call(rbind,results)
#write.table(results,"results.txt",sep="\t",row.names=F)

# cl <- makeCluster(cores)
# registerDoParallel(cl)
# clusterCall(cl, function(x) .libPaths(x), .libPaths())
# a <- 5; b <- nsim
# rng <- RNGseq(a*b, as.numeric(args[[1]]))
# pak<-c("SuperLearner","tmle","ranger","rpart",
#        "nnet","arm","earth","e1071","xgboost",
#        "glmnet","gam","earth","rmutil")
# obj<-c('res.est','res.se','true')
# results<-foreach(i=1:a,
#                  .packages=pak,
#                  .export=obj) %:% 
#   foreach(j=1:b,
#           .packages=pak,
#           .export=obj,
#           r=rng[(i-1)*b + 1:b],
#           .combine='rbind') %dopar% {
#             rngtools::setRNG(r)
#             r <- nparDR(i,j,bs=T)
#             r}
# stopCluster(cl)
# results<-do.call(rbind,results)
# tt<-(proc.time()-time)/60;tt
# write.table(results,"results.txt",sep="\t",row.names=F)
# write.table(tt,"time.txt",sep="\t",row.names=F)

proc.time()-time