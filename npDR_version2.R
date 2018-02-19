# load libraries/functions
args=(commandArgs(TRUE))

packages <- c("foreach","doParallel","doRNG","boot","ranger","GGally",
              "mvtnorm","devtools","glmnet","data.table","broom","sandwich",
              "xgboost","foreach","ggplot2","tabplot","gridExtra")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN') 
  }
}

for (package in packages) {
  library(package, character.only=T)
}

devtools::install_github("ecpolley/SuperLearner") 
library(SuperLearner)

thm <- theme_classic() + 
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

# devtools::install_github("imbs-hl/ranger")
# library(ranger)

# don't need rJava, as ck37r is no longer being used
# https://stackoverflow.com/questions/12872699/error-unable-to-load-installed-packages-just-now/25932828#25932828
#install.packages("rJava",repos='http://lib.stat.cmu.edu/R/CRAN') 
#library(rJava)

# install.packages("bartMachine",repos='http://lib.stat.cmu.edu/R/CRAN')
# options(java.parameters = "-Xmx5g")
# library(bartMachine)

#devtools::install_github("ck37/ck37r")
#library(ck37r)

# install.packages("tmle",repos='http://lib.stat.cmu.edu/R/CRAN')
# library(tmle)

# previous version
packageurl <- "https://cran.r-project.org/src/contrib/Archive/tmle/tmle_1.2.0-4.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
library(tmle)

time<-proc.time()
print(args)
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }
nsim <- as.numeric(args[[2]])

cols <- c("ipwNPT","ipwPMT","ipwPMF","ipwNPF",
          "regNPT","regPMT","regPMF","regNPF")
res.est <- data.frame(matrix(nrow=1,ncol=length(cols)));
colnames(res.est) <- cols
res.se <- res.est

xgboost_learner = create.Learner("SL.xgboost", tune = list(ntrees = c(500),
                                                           max_depth = c(10,20),
                                                           shrinkage = c(0.01,0.001)))

ranger_learner = create.Learner("SL.ranger",param=list(num.trees=500),
                                            tune=list(min.node.size=c(5,25,50),mtry=c(2,3)))

glmnet_learner = create.Learner("SL.glmnet", tune=list(alpha=seq(0,1,.1))) 

sl.lib <- c(ranger_learner$names,xgboost_learner$names)

init<-data.frame(type=NA,run=NA,N=NA,res.est,max_swPT=NA,max_swNPT=NA,mean_swNPT=NA,mean_swPT=NA)
write.table(init,"./results_ipw.txt",sep="\t",row.names=F)

##  true value
true<-3
npDR<-function(counter,bs=F,corr=as.numeric(args[[3]]),
               varX=as.numeric(args[[4]]),bootNum=100){
  # data management
  set.seed(counter)
  i<-counter
  samp<-rep(1:4,nsim*4)[counter]
  ss<-c(100,200,600,1200)
  n<-ss[samp]
  cat("Now running iteration",i,"with a sample size of",n,'\n');flush.console()
  
  # confounders
  sigma<-matrix(corr,nrow=4,ncol=4);diag(sigma)<-1
  x <- rmvnorm(n, mean=rep(0,4), sigma=sigma)
  
  # design matrix for exposure and outcome model
  dMatT<-model.matrix(as.formula("~x[,1]+x[,2]+x[,3]+x[,4]"))
  
  beta<-c(120,4,4,4,4) 
  theta<-c(-1,log(1.75),log(1.75),log(1.75),log(1.75))
  
  # design matrix for propensity score model
  mu <- dMatT%*%beta
  # propensity score model
  pi <- expit(dMatT%*%theta)
  r<-rbinom(n,1,pi)
  # outcome model: true expsoure effect = 3
  y <- r*true + mu + rnorm(n,0,1)
  
  W <- x;colnames(W) <- paste("W",1:4,sep="");A <- r;Y <- y
  W1 <- x[,1]
  W2 <- x[,2]
  W3 <- x[,3]
  W4 <- x[,4]
  tQForm<-as.formula("Y~A+W1+W2+W3+W4")
  tgForm<-as.formula("A~W1+W2+W3+W4")
  
  cat("Running correct parametric TMLE",'\n');flush.console()
  g.predPMT <- glm(A~W1+W2+W3+W4,family=binomial)$fitted.values
  Q.PMT <- glm(Y~A+W1+W2+W3+W4,family=gaussian)

  pihatPT <- g.predPMT
  swPT <- A*(mean(A)/pihatPT) + (1-A)*(mean(1-A)/(1-pihatPT))
  max_swPT=max(swPT)
  mean_swPT=mean(swPT)
  
  folds=5
  index<-split(1:n,1:folds)
  
  Q.mod<-SuperLearner(Y,data.frame(A,W),family=gaussian,
                      SL.library = sl.lib,
                      cvControl=list(V=folds,validRows=index))
  newD1<-data.frame(A=1,W);newD0<-data.frame(A=0,W)
  Q1.NPT<-predict(Q.mod,newdata=newD1,onlySL=T)$pred
  Q0.NPT<-predict(Q.mod,newdata=newD0,onlySL=T)$pred
  
  Q1.PT<-predict(Q.PMT,newdata=newD1)
  Q0.PT<-predict(Q.PMT,newdata=newD0)
  
  g.predNPT<-SuperLearner(A,data.frame(W),family=binomial,
                          cvControl=list(V=folds,validRows=index),
                          SL.library = sl.lib)$SL.predict
  
  g.predNPT <- pmax(.01,g.predNPT)
  print(summary(g.predNPT))
  
  swNPT <- A*(mean(A)/g.predNPT) + (1-A)*(mean(1-A)/(1-g.predNPT))
  max_swNPT=max(swNPT)
  mean_swNPT=mean(swNPT)
  
  # induce misspecification
  z<-x
  # z[,1]<-(1-r)*(exp(x[,1]/2))            + (r)*(x[,1]+2*x[,2])
  # z[,2]<-(1-r)*exp(x[,2])/(x[,3]+5)      + (r)*x[,3]*exp(x[,2])/5
  # z[,3]<-(1-r)*((x[,1]*x[,3]/5))         + (r)*(x[,3] - x[,1])
  # z[,4]<-(1-r)*exp(x[,4])                + (r)*x[,4]
  z[,1]<-(1-r)*(x[,1])                   + (r)*(exp(x[,1]))
  z[,2]<-(1-r)*(x[,2]*x[,1])                   + (r)*(x[,2]*x[,2])
  z[,3]<-(1-r)*(2*exp(x[,3]))                  + (r)*(x[,3]*x[,1])
  z[,4]<-(1-r)*(exp(x[,4])/(x[,3]+5))                   + (r)*(x[,4]*x[,2])
  
  W <- z;colnames(W) <- paste("W",1:4,sep="");A <- r;Y <- y
  W1 <- z[,1]
  W2 <- z[,2]
  W3 <- z[,3]
  W4 <- z[,4]
  
  cat("Running correct parametric TMLE",'\n');flush.console()
  g.predPMF <- glm(A~W1+W2+W3+W4,family=binomial)$fitted.values
  Q.PMF <- glm(Y~A+W1+W2+W3+W4,family=gaussian)
  
  pihatPF <- g.predPMF
  swPF <- A*(mean(A)/pihatPF) + (1-A)*(mean(1-A)/(1-pihatPF))
  max_swPF=max(swPF)
  mean_swPF=mean(swPF)
  
  Q.mod<-SuperLearner(Y,data.frame(A,W),family=gaussian,
                      SL.library = sl.lib,
                      cvControl=list(V=folds,validRows=index))
  newD1<-data.frame(A=1,W);newD0<-data.frame(A=0,W)
  Q1.NPF<-predict(Q.mod,newdata=newD1,onlySL=T)$pred
  Q0.NPF<-predict(Q.mod,newdata=newD0,onlySL=T)$pred
  
  Q1.PF<-predict(Q.PMF,newdata=newD1)
  Q0.PF<-predict(Q.PMF,newdata=newD0)
  
  g.predNPF<-SuperLearner(A,data.frame(W),family=binomial,
                          cvControl=list(V=folds,validRows=index),
                          SL.library = sl.lib)$SL.predict
  
  g.predNPF <- pmax(.01,g.predNPF)
  print(summary(g.predNPF))
  
  swNPF <- A*(mean(A)/g.predNPF) + (1-A)*(mean(1-A)/(1-g.predNPF))
  max_swNPF=max(swNPF)
  mean_swNPF=mean(swNPF)
  
  # compute estimators
  res.est$ipwPMT <- coef(lm(Y~A,weights=swPT))[2] 
  res.est$ipwNPT <- coef(lm(Y~A,weights=swNPT))[2]
  
  res.est$regPMT <- mean(Q1.PT-Q0.PT)
  res.est$regNPT <- mean(Q1.NPT-Q0.NPT)
  
  res.se$ipwPMT <- (vcovHC(lm(Y~A,weights=swPT))[2,2])*sqrt(n)
  res.se$ipwNPT <- (vcovHC(lm(Y~A,weights=swNPT))[2,2])*sqrt(n)
  
  res.est$ipwPMF <- coef(lm(Y~A,weights=swPF))[2] 
  res.est$ipwNPF <- coef(lm(Y~A,weights=swNPF))[2]
  
  res.est$regPMF <- mean(Q1.PF-Q0.PF)
  res.est$regNPF <- mean(Q1.NPF-Q0.NPF)
  
  res.se$ipwPMF <- (vcovHC(lm(Y~A,weights=swPF))[2,2])*sqrt(n)
  res.se$ipwNPF <- (vcovHC(lm(Y~A,weights=swNPF))[2,2])*sqrt(n)
  
  # compute bootstrap CIs
  datB<- data.frame(z,x,r,y)
  names(datB)<-c(paste0("z",1:4),paste0("x",1:4),"r","y")
  cat("Bootstrapping Parametric",'\n');flush.console()
  plugin1 <- function(d,j){
    # correct
    dat<-d[j,c("y","r",paste0("x",1:4))]
    mumod<-glm(y~.,data=dat,family=gaussian)
    dat$r<-1;muhat1 <- predict(mumod,newdata=dat)
    dat$r<-0;muhat0 <- predict(mumod,newdata=dat)
    meanPT<-mean(muhat1-muhat0)
    # misspecified
    dat<-d[j,c("y","r",paste0("z",1:4))]
    mumod<-glm(y~.,data=dat,family=gaussian)
    dat$r<-1;muhat1 <- predict(mumod,newdata=dat)
    dat$r<-0;muhat0 <- predict(mumod,newdata=dat)
    meanPF<-mean(muhat1-muhat0)
    res<-c(meanPT,meanPF)
  }
  bs1 <- boot(datB,plugin1,R=bootNum)
  res.se$regPMT <- sd(bs1$t[,1])*sqrt(n)
  res.se$regPMF <- sd(bs1$t[,2])*sqrt(n)
  if (bs==T){
    cat("Bootstrapping Nonparametric",'\n');flush.console()
    plugin2<-function(d){
      j<-sample(1:nrow(d),nrow(d),replace=T);dat<-d[j,]
      W<-dat[,c(paste0("x",1:4))]
      Y<-dat$y;A<-as.matrix(dat$r) 
      X=data.frame(cbind(A,W))
      names(X)<-c("A",paste0("x",1:ncol(W)))
      mumod<-SuperLearner(Y,X,family=gaussian,SL.library=sl.lib)
      X$A<-1;muhat1 <- predict(mumod,newdata=X,onlySL=T)$pred
      X$A<-0;muhat0 <- predict(mumod,newdata=X,onlySL=T)$pred
      gT<-mean(muhat1-muhat0)
      W<-dat[,c(paste0("z",1:4))]
      Y<-dat$y;A<-as.matrix(dat$r) 
      X=data.frame(cbind(A,W))
      names(X)<-c("A",paste0("x",1:ncol(W)))
      mumod<-SuperLearner(Y,X,family=gaussian,SL.library=sl.lib)
      X$A<-1;muhat1 <- predict(mumod,newdata=X,onlySL=T)$pred
      X$A<-0;muhat0 <- predict(mumod,newdata=X,onlySL=T)$pred
      gF<-mean(muhat1-muhat0)
      g<-c(gT,gF)
      return(g)
    }
    b<-lapply(1:bootNum,function(x) plugin2(datB))
    results<-do.call(rbind,b)
    res.se$regNPT <- sd(results[,1])*sqrt(n)
    res.se$regNPF <- sd(results[,2])*sqrt(n)
  } else{
    res.se$regNPT <- res.se$regNPF <- NA
  }
  
  res.cov <- res.est-1.96*res.se/sqrt(n) < true & true < res.est+1.96*res.se/sqrt(n)
  res.width <- (res.est+1.96*res.se/sqrt(n)) - (res.est-1.96*res.se/sqrt(n))
  
  tmp<-data.frame(rbind(res.est-true,(res.est-true)^2,res.cov,res.width),max_swPT,max_swNPT,mean_swNPT,mean_swPT)
  tmp<-cbind(i,n,tmp)
  rownames(tmp)<-c("bias","rmse","cov","width");colnames(tmp)[1:2]<-c("run","N");
  cat('\n')
  cat("Printing Results",'\n')
  print(round(tmp,3));cat('\n');flush.console()
  setDT(tmp, keep.rownames = TRUE)[];colnames(tmp)[1] <- "type"
  
  write.table(tmp,"./results_ipw.txt",sep="\t",row.names=F,col.names=F,append=T)
  
  return(tmp)
}

cores<-20
mclapply(1:(nsim*4), function(x) npDR(x,bs=F),mc.cores=cores)
# d<-do.call(rbind,r)
# d$N<-factor(d$N,levels=c("100","200","600","1200"))
# 
# head(d)
# 
# ggplot(subset(d,type=="bias"),aes(x=ipwPMT,y=ipwO,group=N,color=N)) + 
#   geom_point()
# 
# ggplot(subset(d,type=="bias"),aes(ipwPMT,group=N,color=N)) + 
#   geom_density() + xlab("Bias")
# 
# ggplot(subset(d,type=="bias"),aes(ipwO,group=N,color=N)) + 
#   geom_density() + xlab("Bias")

