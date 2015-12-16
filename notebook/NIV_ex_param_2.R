# Script from Marm Fitzpatrick generating the original fit for the NiV-based
# SIR Model

rm(list=ls()) # clears workspace
library(gplots)
#Count data (not used directly for model fitting)------------
m<-read.csv("c:/marm/research/other/nipah/niVFarLong.csv");m$Date=as.Date(m$Date,"%m/%d/%Y")
m1=data.frame(month=NA, Age = NA, NiV = -1);m1=m1[-1,];levels(m1$Age)=c("A","J","P")
ages=c("Adult","Juv");NiVs=c("pos","neg")
for (i in 1:nrow(m)) {
  for (j in 1:length(ages)) {
    for (k in 1:length(NiVs)) {
      ch=paste(ages[j],".",NiVs[k],sep="")
      if (m[i,ch]>0) {mT=data.frame(week=m[rep(i,m[i,ch]),"Rweek"],Age=ages[j],NiVc=NiVs[k]); m1=rbind(m1,mT) }
    }}}
m1$NiV=-9; m1$NiV[m1$NiVc=="pos"]=1;m1$NiV[m1$NiVc=="neg"]=0
plot(N~I(104+Rweek),data=m,xlab="Week",ylab="Total Count");abline(v=(1:8)*52,col = "gray", lty = "dotted") #pop counts over time
f1=lm(N~I(104+Rweek),data=m);summary(f1)
abline(a=coef(f1)[1],b=coef(f1)[2])
plot(Njuv~Rweek,data=m,type="b",ylim=c(0,250))
lines(Nadult~Rweek,data=m,col="red")
abline(v=(1:8)*52,col = "gray", lty = "dotted")
plot(Fjuv~Rweek,data=m,type="b",ylab=c("Fraction juvenile"))
abline(h=.5)
summary(lm(Fjuv~Rweek,data=m))

#----------------------------------------------------

#seroprevalence data
md<-read.csv("c:/marm/research/other/nipah/niVFarLong.csv")
md$Juv.prev=md$Juv.pos/md$Juv.tot;md$Adult.prev=md$Adult.pos/md$Adult.tot;
md$Pjse=(md$Juv.prev*(1-md$Juv.prev)/md$Juv.tot)^.5;md$Pase=(md$Adult.prev*(1-md$Adult.prev)/md$Adult.tot)^.5
plotCI(x=md$Rweek,y=md$Juv.prev,uiw=md$Pjse,type="b",col="red");abline(v=(1:8)*52,col = "gray", lty = "dotted")
plotCI(x=md$Rweek,y=md$Adult.prev,uiw=md$Pase,type="b");abline(v=(1:8)*52,col = "gray", lty = "dotted")

NiV= function(Bjj,Bja,Baj,Baa,SA,Ias){ #Model which simulated dynamics given parameter estimates
  st=104;tend=st+max(md$Rweek);#starting and ending weeks
  A=rep(NA,tend);Rw=A;Rw[]=0;Nj=A;wc=A;wc[1:(st)]=1:52;nJ=Nj;#initializes variables
  #A - #adults, Nj-# juvs; nJ-# new juveniles born that week; wc-week
#  A[1:(st+1)]=195;Nj[1:(st+1)]=100;nJ[]=0;nJ[c(27:35,79:87)]=100/9 #initial conditions 2 month reproduction
  A[1:(st+1)]=195;Nj[1:(st+1)]=100;nJ[]=0;nJ[c(22:35,74:87)]=100/14 #initial conditions 3 month reproduction
  Sa=(1-Ias)*A;Ia=0.005*A;Ra=(Ias-0.005)*A; # # of adult (a) susceptibles, infecteds, recovereds
  Sj=0.8*Nj;Ij=0.005*Nj;Rj=0.195*Nj ## # of juveniles (j) susceptibles, infecteds, recovereds
  SJ=((1.05-SA^52)/.5)^(1/52); #adult survival is a fitted parameter, juv survival is scaled to give lambda = 0.92
  g=1;w=1#gamma = recovery rate = 1/ infectious period; w = week counter
  for (t in (st+1):tend) {#loop over weeks of data
#    if(w>26&w<36) {nJ[t]= A[t-13]*0.5*(1/9)} # 2 month reproduction
    if(w>21&w<36) {nJ[t]= A[t-13]*0.5*(1/14)} # 3 month reproduction
    #    w/in age trans                     b/w age trans    recovery mortality    births  transition to adults
    dSj=-exp(Bjj)*Sj[t]*Ij[t]/Nj[t]-exp(Baj)*Sj[t]*Ia[t]/Nj[t]        -(1-SJ)*Sj[t]+nJ[t]-(nJ[t-st]*SJ^st)*Sj[t]/Nj[t]
    dIj=+exp(Bjj)*Sj[t]*Ij[t]/Nj[t]+exp(Baj)*Sj[t]*Ia[t]/Nj[t]-g*Ij[t]-(1-SJ)*Ij[t]      -(nJ[t-st]*SJ^st)*Ij[t]/Nj[t]
    dRj=+                                            g*Ij[t]-(1-SJ)*Rj[t]                -(nJ[t-st]*SJ^st)*Rj[t]/Nj[t]
    dSa=-exp(Baa)*Sa[t]*Ia[t]/A[t]-exp(Bja)*Sa[t]*Ij[t]/A[t]          -(1-SA)*Sa[t]      +(nJ[t-st]*SJ^st)*Sj[t]/Nj[t]
    dIa=+exp(Baa)*Sa[t]*Ia[t]/A[t]+exp(Bja)*Sa[t]*Ij[t]/A[t]  -g*Ia[t]-(1-SA)*Ia[t]      +(nJ[t-st]*SJ^st)*Ij[t]/Nj[t]
    dRa=+                                            g*Ia[t]-(1-SA)*Ra[t]                +(nJ[t-st]*SJ^st)*Rj[t]/Nj[t]
    Sj[t+1]=Sj[t]+dSj;Ij[t+1]=Ij[t]+dIj;Rj[t+1]=Rj[t]+dRj
    Sa[t+1]=Sa[t]+dSa;Ia[t+1]=Ia[t]+dIa;Ra[t+1]=Ra[t]+dRa
    Nj[t+1]=Sj[t+1]+Ij[t+1]+Rj[t+1];A[t+1]=Sa[t+1]+Ia[t+1]+Ra[t+1]
    wc[t]=w;w=w+1;Rw[t]=t-st;if(w==53) {w=1}}
  prev=data.frame(w=c(Rw,tend-st+1),Aprev=Ra/A,Jprev=Rj/Nj)
  return(prev)
}

ff=function(B) {
  Bjj=B[1];Bja=B[2];Baj=B[3];Baa=B[4];SA=B[5];Ias=B[6]
  prev=NiV(Bjj,Bja,Baj,Baa,SA,Ias)
  lke=lkf(prev)
  print(c(B,lke))
  return(lke)  }

lkf=function(prev) {   #This function estimates the neg log likelihood using data in data.frame m1
  ud=unique(m1$week)#unique time points
  lkA=matrix(0,(length(ud)-1),1);#likelihoods for Adult data
  lkJ=matrix(0,(length(ud)-1),1);#likelihoods for juvenile data
  for (i in 2:length(ud)) {
    rowA=sum(m1$week==ud[i]&m1$Age=="Adult")#total adults on week
#    m1[m1$week==ud[i]&m1$Age=="Adult",]
    ppA=sum(m1$NiV[m1$week==ud[i]&m1$Age=="Adult"])#observed # pos adult
    npA=rowA-ppA#observed # neg adult
    rowJ=sum(m1$week==ud[i]&m1$Age=="Juv")#total juvs on week
    ppJ=sum(m1$NiV[m1$week==ud[i]&m1$Age=="Juv"])#observed # pos Juv
    npJ=rowJ-ppJ#observed # neg Juv
    lkA[i-1]=(prev$Aprev[prev$w==ud[i]]^ppA)*(1-prev$Aprev[prev$w==ud[i]])^npA #likelihood of adult data
    lkJ[i-1]=(prev$Jprev[prev$w==ud[i]]^ppJ)*(1-prev$Jprev[prev$w==ud[i]])^npJ #likelihood of adult data
  }
  tkl=-sum(log(lkA),log(lkJ)) #sum likelihood from both
  return(tkl)			}

#-----------Fit model to data
#B=trans coeffjuv-juv  adu->juv   juv->ad   adu-adu   ad surv       initial adult seroprevalence
sp=data.frame(Bjj=.328,Baj=-13.3,Bja=-16.1,Baa=-22.6,SA=.81^(1/52),Ias=.808)
fp=optim(par=sp,fn=ff,control=list());fp#fit model with starting estimates
#note: betas are all log_e transformed so real value is exp(par)


#Simulate fitted model---------------------------------
st=104;tend=st+max(md$Rweek);tend=1000
#next line inserts fitted parameter estimates to simulate fitted model
Bjj=fp$par[1];Bja=fp$par[2];Baj=fp$par[3];Baa=fp$par[4];SA=fp$par[5];Ias=fp$par[6] #uncomment this line to simulate fitted model
#A=rep(NA,tend);Rw=A;Rw[]=0;Nj=A;wc=A;wc[1:(st)]=1:52;nJ=Nj;A[1:(st+1)]=195;Nj[1:(st+1)]=100;nJ[]=0;nJ[]=0;nJ[c(27:35,79:87)]=100/9 #initial conditions 2 month reproduction
A=rep(NA,tend);Rw=A;Rw[]=0;Nj=A;wc=A;wc[1:(st)]=1:52;nJ=Nj;A[1:(st+1)]=195;Nj[1:(st+1)]=100;nJ[]=0;nJ[]=0;nJ[c(22:35,74:87)]=100/14 #initial conditions 3 month reproduction
Sj=0.8*Nj;Ij=0.005*Nj;Rj=0.195*Nj
Sa=(1-Ias)*A;Ia=0.005*A;Ra=(Ias-0.005)*A;
SJ=((1.05-SA^52)/.5)^(1/52);g=1;w=1; #1.138 = stable; 8.8% decline/yr estimated lambda
#SA=.9^(1/52);SJ=.46^(1/52);g=1;w=1
#SA=.8^(1/52);
#SA=.9^(1/52);SJ=((1.05-SA^52)/.5)^(1/52);g=1;w=1; #1.138 = stable; 8.8% decline/yr estimated lambda
#Bjj=1.40617773;Bja=0;Baj=0;Baa=2.36570967 #uncomment this line to simulate fitted model
for (t in (st+1):tend) {
  if(w>21&w<36) {nJ[t]= A[t-13]*0.5*(1/14)} # 3 month reproduction
#  if(w>26&w<36) {nJ[t]= A[t-13]*0.5*(1/9)}  #2 month reproduction
  #    w/in age trans   b/w age trans  recovery             mortality  births  transition to adults
  dSj=-exp(Bjj)*Sj[t]*Ij[t]/Nj[t]-exp(Baj)*Sj[t]*Ia[t]/Nj[t]        -(1-SJ)*Sj[t]+nJ[t]-(nJ[t-st]*SJ^st)*Sj[t]/Nj[t]
  dIj=+exp(Bjj)*Sj[t]*Ij[t]/Nj[t]+exp(Baj)*Sj[t]*Ia[t]/Nj[t]-g*Ij[t]-(1-SJ)*Ij[t]      -(nJ[t-st]*SJ^st)*Ij[t]/Nj[t]
  dRj=+                                            g*Ij[t]-(1-SJ)*Rj[t]      -(nJ[t-st]*SJ^st)*Rj[t]/Nj[t]
  dSa=-exp(Baa)*Sa[t]*Ia[t]/A[t]-exp(Bja)*Sa[t]*Ij[t]/A[t]          -(1-SA)*Sa[t]      +(nJ[t-st]*SJ^st)*Sj[t]/Nj[t]
  dIa=+exp(Baa)*Sa[t]*Ia[t]/A[t]+exp(Bja)*Sa[t]*Ij[t]/A[t]  -g*Ia[t]-(1-SA)*Ia[t]      +(nJ[t-st]*SJ^st)*Ij[t]/Nj[t]
  dRa=+                                            g*Ia[t]-(1-SA)*Ra[t]      +(nJ[t-st]*SJ^st)*Rj[t]/Nj[t]
  Sj[t+1]=Sj[t]+dSj;Ij[t+1]=Ij[t]+dIj;Rj[t+1]=Rj[t]+dRj
  Sa[t+1]=Sa[t]+dSa;Ia[t+1]=Ia[t]+dIa;Ra[t+1]=Ra[t]+dRa
  Nj[t+1]=Sj[t+1]+Ij[t+1]+Rj[t+1];A[t+1]=Sa[t+1]+Ia[t+1]+Ra[t+1]
  wc[t]=w;w=w+1;Rw[t]=t-st;if(w==53) {w=1}}
prev=data.frame(w=c(Rw,tend-st+1),Aprev=Ra/A,Jprev=Rj/Nj)
dm=data.frame(t=1:(tend+1),wc=c(wc,0),A,Nj,nJ=c(nJ,0),Sa,Ia,Ra,Sj,Ij,Rj)

#Plot juvenile fitted model and data
par(mar=c(2.5,3,0,3),mgp = c(1.5, .5, 0),cex=1.)
plot(Rj/Nj~t,data=dm,type="l",xlab="Week",ylab="Seroprevalence",xlim=c(st,tend),ylim=c(0,.6));abline(v=(1:8)*52,col = "gray", lty = "dotted")
plotCI(x=I(st+md$Rweek),y=md$Juv.prev,uiw=md$Pjse,type="b",add=T,col="red",pch=19)
yrs=2006:2012;dts=118+52*0:6
text(dts,-0.01,labels=yrs)

plot(A~t, data=dm)
par(new=T); plot(Nj~t,data=dm,col="green")

#Plot fitted model and data
pdf(file="c:/marm/research/other/nipah/NiVProj.pdf",width=11,height=8.5)
library(gplots)
po=2; #po=1 - Juveniles only; po=2
par(mar=c(2.5,3,0,3),mgp = c(1.5, .5, 0),cex=1.)
#if(po==1) plot(Rj/Nj~t,data=dm,type="l",xlab="Week",ylab="Seroprevalence",xlim=c(st,tend),ylim=c(0,1));abline(v=(1:8)*52,col = "gray", lty = "dotted")
if(po>1) plot(Rj/Nj~t,data=dm,type="l",xlab="Week",ylab="Seroprevalence",xlim=c(st,tend),ylim=c(-0.1,1));abline(v=(1:(max(dm$t)/min(dm$t))*52),col = "gray", lty = "dotted")
plotCI(x=I(st+md$Rweek),y=md$Juv.prev,uiw=md$Pjse,add=T,pch=19)
lines(Ra/A~t,data=dm,col="red")
plotCI(x=I(st+md$Rweek),y=md$Adult.prev,uiw=md$Pase,add=T,col="red")
if(po==1) legend(117,1.04, c("Juveniles", "Adults"),col = c("black","red"),lty=1,bty="n")
#if(po>1) legend(117,1.04, c("Infected Juveniles","Juveniles", "Adults"),col = c("green","black","red"),lty=1,bty="n")
if(po>1) legend(117,1.04, c("New Juveniles","Infected Juveniles","Juvenile seroprev.", "Adult  seroprev."),col = c("blue", "green","black","red"),lty=1,bty="n") # add new juveniles
yrs=2006:2012;dts=118+52*0:6
abline(h=0)
if(po==1) text(dts,-0.02,labels=yrs)
if(po>1)  text(dts,-.12,labels=yrs)

#Add Just Ij
par(new=T);plot(3*Ij~t,data=dm,type="l",xlim=c(st,tend),ylim=c(0,35),xlab="",ylab="",xaxt="n",yaxt="n",col="green")
axis(4);mtext("Density of infected juveniles",side=4,line=1.5)
par(new=T);plot(nJ~t,data=dm,type="l",xlim=c(st,tend),ylim=c(0,35),xlab="",ylab="",xaxt="n",yaxt="n",col="blue")

#Add Ij + Adult counts
par(new=T);plot(9*Ij~t,data=dm,type="l",xlim=c(st,tend),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n",col="green")
par(new=T); plot(N~I(104+Rweek),data=m,xaxt="n",yaxt="n",xlab="",ylab="",col="blue",pch=15,xlim=c(st,tend),ylim=c(0,350));
#par(new=T); plot(N~I(104+Rweek),data=m,xaxt="n",yaxt="n",xlab="",ylab="",col="blue",pch=15,ylim=c(0,350));
par(new=T);plot(nJ~t,data=dm,type="l",xlim=c(st,tend),ylim=c(0,35),xlab="",ylab="",xaxt="n",yaxt="n",col="blue")
axis(4);mtext("Total Count and density of infected juveniles",side=4,line=1.5)
dev.off()


#plot of adult data + counts
plotCI(x=I(st+md$Rweek),y=md$Adult.prev,uiw=md$Pase,col="red")
lines(x=I(st+md$Rweek),y=md$Adult.prev)
abline(v=(1:8)*52,col = "gray", lty = "dotted")
par(new=T);
plot(N~I(104+Rweek),data=m,xaxt="n",yaxt="n",xlab="",ylab="",col="blue",pch=15,ylim=c(0,350));

#par(new=T)
#plot(7*Ia~t,data=dm,type="l",xlab="Week",ylab="Density of infected adult bats",xlim=c(st,tend));abline(v=(1:8)*52,col = "gray", lty = "dotted")
plot(3*Ij~t,data=dm,type="l",xlim=c(st,tend),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n",col="green")
par(new=T)
plot(N~I(104+Rweek),data=m,xaxt="n",yaxt="n",xlab="",ylab="",col="blue",pch=15,ylim=c(0,350));
axis(4)
mtext("Total Count",side=4,line=1.5)


write.csv(dm,file="c:/marm/research/other/Nipah/out.csv")
#Juvenile dynamics
plot(Sj~t,data=dm,type="l",xlim=c(st,tend),xlab="Week",ylab="Density of infected juvenile bats",ylim=c(0,100));abline(v=(1:8)*52,col = "gray", lty = "dotted")
lines(Ij~t,data=dm,type="l",xlim=c(st,tend),xlab="Week",ylab="Density of infected juvenile bats",col="red");abline(v=(1:8)*52,col = "gray", lty = "dotted")
lines(Rj~t,data=dm,type="l",xlim=c(st,tend),xlab="Week",ylab="Density of infected juvenile bats",col="blue");abline(v=(1:8)*52,col = "gray", lty = "dotted")
lines(Nj~t,data=dm,type="l",xlim=c(st,tend),xlab="Week",ylab="Density of infected juvenile bats",col="green");abline(v=(1:8)*52,col = "gray", lty = "dotted")
legend(350, 105, c("Susceptible", "Infected","Recovered","Total"),col = c("black","red","blue","green"),bty="n",lty=1)

#Adult dynamics
plot(Sa~t,data=dm,type="l",xlim=c(st,tend),xlab="Week",ylab="Density of infected juvenile bats",ylim=c(0,200));abline(v=(1:8)*52,col = "gray", lty = "dotted")
lines(Ia~t,data=dm,type="l",xlim=c(st,tend),xlab="Week",ylab="Density of infected juvenile bats",col="red");abline(v=(1:8)*52,col = "gray", lty = "dotted")
lines(Ra~t,data=dm,type="l",xlim=c(st,tend),xlab="Week",ylab="Density of infected juvenile bats",col="blue");abline(v=(1:8)*52,col = "gray", lty = "dotted")
lines(A~t,data=dm,type="l",xlim=c(st,tend),xlab="Week",ylab="Density of infected juvenile bats",col="green");abline(v=(1:8)*52,col = "gray", lty = "dotted")
legend(350, 190, c("Susceptible", "Infected","Recovered","Total"),col = c("black","red","blue","green"),bty="n",lty=1)
#--------------------------------------
