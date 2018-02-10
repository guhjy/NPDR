library(here)
library(tidyverse)
library(reshape2)

## BIAS
a <- read_delim("./results.txt",delim="\t");a
a <- a %>% filter(type=="bias")
nrow(a)
A<-c("PMT","NPT","PMF","NPF")
vv<-c(paste("ipw",A,sep=""),paste("reg",A,sep=""),
      paste("aipw",A,sep=""),paste("tmle",A,sep=""))
res1 <- melt(a,measure.vars=vv)[,c("N","variable","value")]
head(res1);tail(res1)
res1$N<-as.factor(res1$N);res1$length<-str_length(res1$variable)
table(res1$length)
res1$estimator<-str_sub(res1$variable,1,res1$length-3)
res1$estimator<-factor(res1$estimator)
res1$estimator<-factor(res1$estimator,levels(res1$estimator)[c(3,2,1,4)])
res1$model_type<-str_sub(res1$variable,-3,-1)
res1$variable<-res1$length<-NULL
head(res1);tail(res1)

lngth<-nrow(res1[res1$N==100&res1$estimator=="reg"&res1$model_type=="PMT",])
summary(res1[res1$N==100&res1$estimator=="aipw"&res1$model_type=="PMT",]$value)

names(res1)
res1 %>% group_by(N,estimator,model_type) %>% summarize(n=n())

data.frame(res1 %>% group_by(N,estimator,model_type) %>% summarize(q=quantile(value,.9)))

res1 <- res1 %>% mutate(value=ifelse(value>60,0,value))

D <- res1 %>% group_by(N,estimator,model_type) %>% summarize(mv = mean(value),sd.mv=sd(value)/sqrt(lngth))
D <- D %>% ungroup(N,estimator,model_type)
D

facet_labs<-list(
  "NPF"="Nonpar Misspec",
  "NPT"="Nonpar Correct",
  "PMF"="Par Misspec",
  "PMT"="Par Correct"
)
type_labeller <- function(variable,value){
  return(facet_labs[value])
}

ggplot(D,aes(as.factor(N),mv,group=estimator)) + 
  theme_light() + facet_grid(model_type ~ . ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x = "Sample Size",y = "Bias",fill = "Estimator\n") +
  geom_bar(aes(fill=estimator), width = 0.5,position = position_dodge(),stat="identity") +
  #geom_errorbar(aes(ymin=mv-sd.mv, ymax=mv+sd.mv),width=.1,position = position_dodge(.5)) +
  scale_fill_grey(labels = c("g Comp", "IPW", "AIPW","TMLE")) +
  #coord_cartesian(ylim=c(-15,15)) +
  geom_hline(yintercept = 0,color="black",
             linetype="dashed",size=.25)

## MSE
a <- read_delim("./results.txt",delim="\t");a
a <- a %>% filter(type=="rmse")
nrow(a)

a %>% group_by(N) %>% summarize(n=n())

A<-c("PMT","NPT","PMF","NPF")
vv<-c(paste("ipw",A,sep=""),paste("reg",A,sep=""),
      paste("aipw",A,sep=""),paste("tmle",A,sep=""))
res1 <- melt(a,measure.vars=vv)[,c("N","variable","value")] ############# vv
head(res1);tail(res1)
res1$N<-as.factor(res1$N);res1$length<-str_length(res1$variable)
res1$estimator<-str_sub(res1$variable,1,res1$length-3)
res1$estimator<-factor(res1$estimator)
res1$estimator<-factor(res1$estimator,levels(res1$estimator)[c(3,2,1,4)])
res1$model_type<-str_sub(res1$variable,-3,-1)
res1$variable<-res1$length<-NULL

data.frame(res1 %>% group_by(N,estimator,model_type) %>% summarize(q=quantile(value,.9)))
res1 <- res1 %>% mutate(value=ifelse(value>1000,1000,value))

D <- res1 %>% group_by(N,estimator,model_type) %>% summarize(mv = median(value))
D <- D %>% ungroup(N,estimator,model_type)

ggplot(D,aes(as.factor(N),mv,group=estimator)) +
  facet_grid(model_type ~ .,labeller = type_labeller) + theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x = "Sample Size",y = "rMSE",shape = "Estimator\n") +
  scale_shape(labels = c("g Comp", "IPW", "AIPW","TMLE")) +
  scale_y_log10() + 
  geom_point(aes(shape=estimator),position=position_dodge(.5),size=3) + 
  theme(text = element_text(size=10)) 
