# In this script we will re-analyse the data behind the
# recent GARD study with a certain question I have been
# curious about for a long time - dose RSI capture enough of
# the variance in outcome such that a dose-response becomes
# visible empirically. A secondary question is whether
# an interaction between RSI and Dose is seen.

require(survival)
require(prodlim)
require(ggplot2)


load("data_df.rda")
summary(rsi.all)
# Notice one patient, see later, has a negative survival time
# We shall discuss data quality and other statsitical
# analysis issues in another analysis
summary(factor(rsi.all$Site))

sum(rsi.all$Time_OS<0,na.rm=T)

# Our interest is in OS since Local Progression is an interval-censored
# event and the intervals of assessment will be different for different
# cohorts which make a pooled analysis problematic.

# We only want all positive survival times and only those that had RT
# for our analysis
rsi.all<-rsi.all[is.na(rsi.all$Time_OS)==F &
                   rsi.all$Time_OS>0 & rsi.all$Received_RT==1,]

# Include GARD 
alpha=log(rsi.all$RSI)/(-2) - 0.05*2
rsi.all$gard=rsi.all$n*rsi.all$d*(alpha + 0.05*rsi.all$d)

# Explore each tumour site
summary(factor(rsi.all$Site))
# What is the event rate like in each study
summary(factor(rsi.all$Event_OS[rsi.all$Site=="breast_TN"]))
# too few events - think about this how comfortable do you 
# feel about 9 events! We will though include it
t1<-rsi.all[rsi.all$Site=="breast_TN",c("Time_OS","Event_OS","RSI")]
ggplot()+geom_point(data=t1,aes(Time_OS,RSI,col=factor(Event_OS)),alpha=0.5,size=2)+
  theme_bw(base_size=16)+xlab("Time (Years)")+ylab("RSI")
# no late events!
summary(factor(rsi.all$Event_OS[rsi.all$Site=="endometrial"]))
summary(factor(rsi.all$Event_OS[rsi.all$Site=="glioma"]))
summary(factor(rsi.all$Event_OS[rsi.all$Site=="lung"]))
summary(factor(rsi.all$Event_OS[rsi.all$Site=="pancreas"]))

summary(factor(rsi.all$Event_OS[rsi.all$Site=="melanoma"]))# 7 events and only 10 patients!
t1<-rsi.all[rsi.all$Site=="melanoma",c("Time_OS","Event_OS","RSI")]
ggplot()+geom_point(data=t1,aes(Time_OS,RSI,col=factor(Event_OS)),alpha=0.5,size=2)+
  theme_bw(base_size=16)+xlab("Time (Years)")+ylab("RSI")
# no late events!
# In order to do the analysis we need dose-variance


ggplot()+geom_jitter(data=rsi.all,aes(Site,TD,col=Site),height=0,width=0.2,alpha=0.3)+
  theme_bw(base_size=16)+xlab("Tumour Site")+ylab("Total RT Dose (Gy)")

# We can see that there is variance in dose in the studies but this variance is
# clearly not random - however what its informative of cannot be assessed as no
# details are provided

# lets explore the varinace in "RSI"
ggplot()+geom_jitter(data=rsi.all,aes(Site,RSI,col=Site),height=0,width=0.2,alpha=0.3)+
  theme_bw(base_size=16)+xlab("Tumour Site")+ylab("RSI")


# now we perform the assessments of interest
m0<-(coxph(Surv(Time_OS,Event_OS)~RSI,data=rsi.all[rsi.all$Site=="endometrial",]))
summary(m0)#0.55 (0.05)
# lets check for non-linearity
m1<-(coxph(Surv(Time_OS,Event_OS)~pspline(RSI),data=rsi.all[rsi.all$Site=="endometrial",]))
summary(m1)#0.55 (0.05)
# does dose become a predictor after adjusting for RSI
m2<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD,data=rsi.all[rsi.all$Site=="endometrial",]))
summary(m2)# 
anova(m0,m2)#0.671 - no
# does an interaction between RSI and dose exist
m3<-(coxph(Surv(Time_OS,Event_OS)~RSI*TD,data=rsi.all[rsi.all$Site=="endometrial",]))
summary(m3)# 
anova(m2,m3)#0.382

# move on to breast

m0<-(coxph(Surv(Time_OS,Event_OS)~RSI,data=rsi.all[rsi.all$Site=="breast_TN",]))
summary(m0)#0.72 (0.07) p - 0.0245 HR = 3261 (2.83- >10000)
# HR makes no sense lets check non-linearity

m1b<-(coxph(Surv(Time_OS,Event_OS)~pspline(RSI),data=rsi.all[rsi.all$Site=="breast_TN",]))
summary(m1b)#0.72 (0.07) p - 0.0245 HR = 3261 (2.83- >10000)
# hint even with such low event rate
summary(rsi.all$RSI[rsi.all$Site=="breast_TN"])
p1<-predict(m1b,newdata = data.frame(RSI=seq(0.3,0.5,by=0.01)),
            se=T)
plot(seq(0.3,0.5,by=0.01),p1$fit,type="l",
     ylab="log(HR)",xlab="RSI",ylim=c(-5,10))
lines(seq(0.3,0.5,by=0.01),p1$fit+p1$se.fit,lty=2)
lines(seq(0.3,0.5,by=0.01),p1$fit-p1$se.fit,lty=2)
abline(h=0,col=2)
mtext("Triple Negative Breast")


m2<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD,data=rsi.all[rsi.all$Site=="breast_TN",]))
summary(m2)# 
anova(m0,m2)#0.260
m3<-(coxph(Surv(Time_OS,Event_OS)~RSI*TD,data=rsi.all[rsi.all$Site=="breast_TN",]))
summary(m3)# 
anova(m2,m3)#0.827

m0<-(coxph(Surv(Time_OS,Event_OS)~RSI,data=rsi.all[rsi.all$Site=="melanoma",]))
summary(m0)#0.88 (0.04) HR = 41490 (15.66 - >10000)
m1<-(coxph(Surv(Time_OS,Event_OS)~pspline(RSI),data=rsi.all[rsi.all$Site=="melanoma",]))
summary(m1)#mom-linearity
summary(rsi.all$RSI[rsi.all$Site=="melanoma"])
p1<-predict(m1,newdata = data.frame(RSI=seq(0.25,0.6,by=0.01)),
            se=T)
plot(seq(0.25,0.6,by=0.01),p1$fit,type="l",
     ylab="log(HR)",xlab="RSI",ylim=c(-5,10))
lines(seq(0.25,0.6,by=0.01),p1$fit+p1$se.fit,lty=2)
lines(seq(0.25,0.6,by=0.01),p1$fit-p1$se.fit,lty=2)
abline(h=0,col=2)
mtext("Melanoma")

m1<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD,data=rsi.all[rsi.all$Site=="melanoma",]))
summary(m1)# 
anova(m0,m1)#0.208
m2<-(coxph(Surv(Time_OS,Event_OS)~RSI*TD,data=rsi.all[rsi.all$Site=="melanoma",]))
summary(m2)# 
anova(m1,m2)# there is no interaction because there is little variance in dose NA

# glioma where they also adjusted for MGMT expression - we shall first look
# without it - its important to note that for some reason this was the only
# tumour type where any adjustment was done
m0<-(coxph(Surv(Time_OS,Event_OS)~RSI,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m0)#0.51 (0.03)
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m1)# 
anova(m0,m1)#0.664
m2<-(coxph(Surv(Time_OS,Event_OS)~RSI*TD,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m2)# 
anova(m1,m2)#0.681

m3<-(coxph(Surv(Time_OS,Event_OS)~MGMT_Expression,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m3)#0.59 (0.03)
m4<-(coxph(Surv(Time_OS,Event_OS)~MGMT_Expression+RSI,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m4)#0.60 (0.03)
anova(m4,m3)#0.050
m5<-(coxph(Surv(Time_OS,Event_OS)~MGMT_Expression+RSI+TD,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m5)# 
anova(m5,m4)#0.815
m6<-(coxph(Surv(Time_OS,Event_OS)~MGMT_Expression+RSI*TD,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m6)# 
anova(m6,m5)#0.605
# After adjusting for MGMT we do see RSI explain a little of the variance but
# still not elucidating a dose-response 


m0<-(coxph(Surv(Time_OS,Event_OS)~RSI,data=rsi.all[rsi.all$Site=="lung",]))
summary(m0)#0.55 (0.04)
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD,data=rsi.all[rsi.all$Site=="lung",]))
summary(m1)# 
anova(m1,m0)#0.538
m2<-(coxph(Surv(Time_OS,Event_OS)~RSI*TD,data=rsi.all[rsi.all$Site=="lung",]))
summary(m2)# 
anova(m2,m1)#0.538

m0<-(coxph(Surv(Time_OS,Event_OS)~RSI,data=rsi.all[rsi.all$Site=="pancreas",]))
summary(m0)#0.55 (0.06)
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD,data=rsi.all[rsi.all$Site=="pancreas",]))
summary(m1)# 
anova(m1,m0)#0.585
m2<-(coxph(Surv(Time_OS,Event_OS)~RSI*TD,data=rsi.all[rsi.all$Site=="pancreas",]))
summary(m2)# 
anova(m2,m1)#0.496

rsi.all$Site<-factor(rsi.all$Site,levels=unique(rsi.all$Site))

m0<-(coxph(Surv(Time_OS,Event_OS)~RSI+strata(Site),data=rsi.all[rsi.all$Site!="melanoma"&
                                                                  rsi.all$Site!="breast_TN",]))
summary(m0)#0.52 (0.02)
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI*Site+strata(Site),data=rsi.all[rsi.all$Site!="melanoma"&
                                                                       rsi.all$Site!="breast_TN",]))
anova(m0,m1)# no interactions
m2<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD+strata(Site),data=rsi.all[rsi.all$Site!="melanoma"&
                                                                          rsi.all$Site!="breast_TN",]))
summary(m2)
anova(m0,m2)#0.584
m3<-(coxph(Surv(Time_OS,Event_OS)~RSI*TD*Site+strata(Site),data=rsi.all[rsi.all$Site!="melanoma"&
                                                                          rsi.all$Site!="breast_TN",]))
summary(m3)# 
anova(m2,m3)#0.8282

# re-do with the two problematic cohorts

m0<-(coxph(Surv(Time_OS,Event_OS)~RSI+strata(Site),data=rsi.all))
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI*Site+strata(Site),data=rsi.all))
anova(m1,m0)
m2<-(coxph(Surv(Time_OS,Event_OS)~RSI*Site+TD+strata(Site),data=rsi.all))
summary(m2)
anova(m2,m1)
m3<-(coxph(Surv(Time_OS,Event_OS)~RSI*Site*TD+strata(Site),data=rsi.all))
anova(m3,m2)

rsi.all$Site<-as.character(rsi.all$Site)
rsi.all$Site[rsi.all$Site=="glioma"]<-"Glioma"
rsi.all$Site[rsi.all$Site=="pancreas"]<-"Pancreas"
rsi.all$Site[rsi.all$Site=="lung"]<-"Lung"
rsi.all$Site[rsi.all$Site=="endometrial"]<-"Endometrial"
rsi.all$Site[rsi.all$Site=="breast_TN"]<-"TN Breast"
rsi.all$Site[rsi.all$Site=="melanoma"]<-"Melanoma"
ggplot()+geom_jitter(data=rsi.all,aes(Site,TD,col=Site),height=0,width=0.2,alpha=0.3)+
  theme_bw(base_size=16)+xlab("Tumour Site")+ylab("Dose (Gy)")
