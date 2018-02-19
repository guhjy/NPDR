library(here)
library(tidyverse)
library(reshape2)
library(gtable)
library(grid)
library(gridExtra)

## BIAS
a <- read_delim("./results_test.txt",delim="\t");a
a <- a %>% filter(type=="bias"&!is.na(type));a
a$N<-a$n
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

data.frame(res1 %>% group_by(N,estimator,model_type) %>% summarize(q=quantile(value,.9)))

res1 <- res1 %>% mutate(value=ifelse(abs(value)>60,0,value))

D <- res1 %>% group_by(N,estimator,model_type) %>% summarize(mv = mean(value),sd.mv=sd(value)/sqrt(lngth),M=n())
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

ggplot(res1,aes(N,value,fill = estimator)) + 
  facet_grid(model_type ~ .,labeller = type_labeller) + theme_light() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x = "Sample Size",y = "Bias",fill="Estimator\n") + 
  scale_fill_grey(labels = c("g Comp", "IPW", "AIPW","TMLE")) +
  geom_boxplot(width=.5,position = position_dodge(.65)) + 
  geom_hline(yintercept = 0,color="red",linetype="dashed",size=.25)
  

ggplot(D,aes(as.factor(N),mv,group=estimator)) + 
  theme_light() +  facet_grid(model_type ~ .,labeller = type_labeller) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x = "Sample Size",y = "Bias",fill = "Estimator\n") +
  geom_boxplot(width=.25,position = position_dodge(),outlier.colour = NA)
  #geom_boxplot(aes(fill=estimator), width = 0.5,position = position_dodge(),stat="identity") +
  #geom_bar(aes(fill=estimator), width = 0.5,position = position_dodge(),stat="identity") +
  #geom_text(size = 5, position = position_dodge(width=.5)) +
  #geom_errorbar(aes(ymin=mv-sd.mv, ymax=mv+sd.mv),width=.1,position = position_dodge(.5)) +
  scale_fill_grey(labels = c("g Comp", "IPW", "AIPW","TMLE")) +
  #coord_cartesian(ylim=c(-2,2)) +
  geom_hline(yintercept = 0,color="black",
             linetype="dashed",size=.25)

# ind_seq<-c("PMT","PMF","NPT","NPF")
# plotFunc<-function(ind){
#   ggplot(D[D$model_type==ind,],aes(as.factor(N),mv,group=estimator,label=M)) + 
#     theme_light() +  facet_grid(model_type ~ .,labeller = type_labeller) + 
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank()) + 
#     labs(x = "Sample Size",y = "Bias",fill = "Estimator\n") +
#     geom_bar(aes(fill=estimator), width = 0.5,position = position_dodge(),stat="identity") +
#     geom_text(size = 5, position = position_dodge(width=.5)) +
#     #geom_errorbar(aes(ymin=mv-sd.mv, ymax=mv+sd.mv),width=.1,position = position_dodge(.5)) +
#     scale_fill_grey(labels = c("g Comp", "IPW", "AIPW","TMLE")) +
#     #coord_cartesian(ylim=c(-2,2)) +
#     geom_hline(yintercept = 0,color="black",
#                linetype="dashed",size=.25)
# }
# plot_object<-lapply(ind_seq,function(x) plotFunc(x))
# 
# legend = gtable_filter(ggplotGrob(plot_object[[1]]), "guide-box")
# thm<-theme(legend.position="none",
#            axis.title.x=element_blank(),
#            axis.title.y=element_blank(),
#            axis.ticks.x=element_blank(),
#            axis.text.x=element_blank())
# 
# grid.arrange(arrangeGrob(plot_object[[4]] + thm, 
#                          plot_object[[3]] + thm,
#                          plot_object[[2]] + thm,
#                          plot_object[[1]] + theme(legend.position="none",axis.title.y=element_blank()), 
#                          nrow = 4,
#                          left = textGrob("Bias", rot = 90, vjust = 1)),
#              legend, 
#              widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
#              ncol=2)
# 
# D
# res1 %>% filter(N==100,estimator=="ipw",model_type=="NPF")

## MSE
a <- read_delim("./results_test.txt",delim="\t");a
a$N<-a$n
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

D <- res1 %>% group_by(N,estimator,model_type) %>% summarize(mv = mean(value))
D <- D %>% ungroup(N,estimator,model_type)

plotFunc<-function(ind){
  ggplot(D[D$model_type==ind,],aes(as.factor(N),mv,group=estimator)) +
    facet_grid(model_type ~ .,labeller = type_labeller) + theme_light() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    labs(x = "Sample Size",y = "rMSE",shape = "Estimator\n") +
    scale_y_log10() +
    scale_shape(labels = c("g Comp", "IPW", "AIPW","TMLE")) +
    geom_point(aes(shape=estimator),position=position_dodge(.5),size=2.5)
}
plot_object<-lapply(ind_seq,function(x) plotFunc(x))

legend = gtable_filter(ggplotGrob(plot_object[[1]]), "guide-box")
thm<-theme(legend.position="none",
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           axis.ticks.x=element_blank(),
           axis.text.x=element_blank())

grid.arrange(arrangeGrob(plot_object[[4]] + thm, 
                         plot_object[[3]] + thm,
                         plot_object[[2]] + thm,
                         plot_object[[1]] + theme(legend.position="none",axis.title.y=element_blank()), 
                         nrow = 4,
                         left = textGrob("rMSE", rot = 90, vjust = 1)),
             legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
             ncol=2)




###### TABLE
a <- read_delim("./results_ipw.txt",delim="\t");a
a <- a %>% filter(!is.na(type));a

a$N<-factor(a$N,levels=c("100","200","600","1200"))

ggplot(subset(a,type=="bias"),aes(ipwPMT,group=N,color=N)) +
  geom_density() + xlab("Bias")

ggplot(subset(a,type=="bias"),aes(ipwNPT,group=N,color=N)) +
  geom_density() + xlab("Bias")

ggplot(subset(a,type=="bias"),aes(regNPT,group=N,color=N)) +
  geom_density() + xlab("Bias")

a %>% filter(type=="bias")  %>% group_by(N) %>% summarise(M=n(),ipwPT=mean(ipwPMT),ipwNPT=mean(ipwNPT),
                                                          ipwPF=mean(ipwPMF),ipwNPF=mean(ipwNPF),
                                                          regPT=mean(regPMT),regNPT=mean(regNPT),
                                                          regPF=mean(regPMF),regNPF=mean(regNPF))
a %>% filter(type=="rmse")  %>% group_by(N) %>% summarise(M=n(),ipwPT=mean(ipwPMT),ipwNPT=mean(ipwNPT),
                                                          ipwPF=mean(ipwPMF),ipwNPF=mean(ipwNPF),
                                                          regPT=mean(regPMT),regNPT=mean(regNPT),
                                                          regPF=mean(regPMF),regNPF=mean(regNPF))
a %>% filter(type=="cov")   %>% group_by(N) %>% summarise(M=n(),ipwPT=mean(ipwPMT),ipwNPT=mean(ipwNPT),
                                                          ipwPF=mean(ipwPMF),ipwNPF=mean(ipwNPF),
                                                          regPT=mean(regPMT),regNPT=mean(regNPT),
                                                          regPF=mean(regPMF),regNPF=mean(regNPF))
a %>% filter(type=="width") %>% group_by(N) %>% summarise(M=n(),ipwPT=mean(ipwPMT),ipwNPT=mean(ipwNPT),
                                                          ipwPF=mean(ipwPMF),ipwNPF=mean(ipwNPF),
                                                          regPT=mean(regPMT),regNPT=mean(regNPT),
                                                          regPF=mean(regPMF),regNPF=mean(regNPF))



