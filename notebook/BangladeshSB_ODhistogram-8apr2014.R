# Script provided by J. Olival as example of how to calculate cut-off values
# to determine seropositivity from Luminex data.


require(MASS)

sb <-read.csv('data/Olival_smallbats/SmallBat_cleaned-8April2014.csv')
names(sb)

NiVsG <- sb$NiVsG
NiVsF <- sb$NiVsF
NiVN <- sb$NiVN
HeVsG <- sb$HeVsG
CedVsG <- sb$CedVsG
EboVGP <- sb$EboVGP
MarVGP <- sb$MarVGP
MenVN <- sb$MenVN
SARSN <- sb$SARSN
MERSN <- sb$MERSN

plot(sb$NiVsG, sb$NiVN)
###############HISTOGRAMS WITH GAMMA DISTRIBUTION Max Likelihood Estimator of 95% AND 99% CUTOFFS####

###NiVsG
pdf("NiVsG_cutoff.pdf",width=11,height=8.5)
h = hist(NiVsG, breaks=100, prob=T, main="ALL NiVsG", xlab="Luminex, MFI values",
        ylab="Frequency", col="lightgrey") #histogram of all data minus positive controls
fitdistr(NiVsG/1000, "gamma") # use to estimate parmaters, scale data by dividing by 1000
f = function(theta,dta){
  l = -sum(dgamma(dta,rate=theta[1],shape=theta[2],log=TRUE))
}
o = optim( c(1.32,5.73),  f, dta = NiVsG)
print(o)
p = qgamma(c(0.95, 0.99), rate=o$par[1],shape=o$par[2]); print(p)
abline(v=p, col=c("blue", "green"), lwd=2 )
lines( h$breaks, dgamma(h$breaks,rate=o$par[1],shape=o$par[2]), col="red",lwd=3,lty="44")
legend("topright", c("Gamma distribution", "95% cutoff (626)", "99% cutoff (925)", "n=7 bats >99%"),
       fill=c("red", "blue", "green", "white"))
x <- NiVsG [NiVsG >925] # figure out how many individual bats above cutoff
length (x) #n=10, added to xlab
dev.off()

###NiVsF
pdf("NiVsF_cutoff.pdf",width=11,height=8.5)
h = hist(NiVsF, breaks=100, prob=T, main="ALL NiVsF", xlab="Luminex, MFI values",
         ylab="Frequency", col="lightgrey") #histogram of all data minus positive controls
fitdistr(NiVsF/1000, "gamma") # use to estimate parmaters, scale data by dividing by 1000
f = function(theta,dta){
  l = -sum(dgamma(dta,rate=theta[1],shape=theta[2],log=TRUE))
}
o = optim( c(1.87, 10.08),  f, dta = NiVsF)
print(o)
p = qgamma(c(0.95, 0.99), rate=o$par[1],shape=o$par[2]); print(p)
abline(v=p, col=c("blue", "green"), lwd=2 )
lines( h$breaks, dgamma(h$breaks,rate=o$par[1],shape=o$par[2]), col="red",lwd=3,lty="44")
legend("topright", c("Gamma distribution", "95% cutoff (449)", "99% cutoff (634)", "n=8 bats >99%"),
       fill=c("red", "blue", "green", "white"))
x <- NiVsF [NiVsF >634] # figure out how many individual bats above cutoff
length (x) #of bats >99% cutoff, added to legend
dev.off()

###NiVN
pdf("NiVN_cutoff.pdf",width=11,height=8.5)
h = hist(NiVN, breaks=100, prob=T, main="ALL NiVN", xlab="Luminex, MFI values",
         ylab="Frequency", col="lightgrey") #histogram of all data minus positive controls
fitdistr(NiVN/1000, "gamma") # use to estimate parmaters, scale data by dividing by 1000
f = function(theta,dta){
  l = -sum(dgamma(dta,rate=theta[1],shape=theta[2],log=TRUE))
}
o = optim( c(1.67, 7.64),  f, dta = NiVN)
print(o)
p = qgamma(c(0.95, 0.99), rate=o$par[1],shape=o$par[2]); print(p)
abline(v=p, col=c("blue", "green"), lwd=2 )
lines( h$breaks, dgamma(h$breaks,rate=o$par[1],shape=o$par[2]), col="red",lwd=3,lty="44")
legend("topright", c("Gamma distribution", "95% cutoff (549)", "99% cutoff (786)", "n=8 bats >99%"),
       fill=c("red", "blue", "green", "white"))
x <- NiVN [NiVN >786] # figure out how many individual bats above cutoff
length (x) #of bats >99% cutoff, added to legend
dev.off()

###HeVsG
pdf("HeVsG_cutoff.pdf",width=11,height=8.5)
h = hist(HeVsG, breaks=100, prob=T, main="ALL HeVsG", xlab="Luminex, MFI values",
         ylab="Frequency", col="lightgrey") #histogram of all data minus positive controls
fitdistr(HeVsG/1000, "gamma") # use to estimate parmaters, scale data by dividing by 1000
f = function(theta,dta){
  l = -sum(dgamma(dta,rate=theta[1],shape=theta[2],log=TRUE))
}
o = optim( c(0.77, 2.92),  f, dta = HeVsG)
print(o)
p = qgamma(c(0.95, 0.99), rate=o$par[1],shape=o$par[2]); print(p)
abline(v=p, col=c("blue", "green"), lwd=2 )
lines( h$breaks, dgamma(h$breaks,rate=o$par[1],shape=o$par[2]), col="red",lwd=3,lty="44")
legend("topright", c("Gamma distribution", "95% cutoff (865)", "99% cutoff (1385)", "n=5 bats >99%"),
       fill=c("red", "blue", "green", "white"))
x <- HeVsG [HeVsG >1385] # figure out how many individual bats above cutoff
length (x) #of bats >99% cutoff, added to legend
dev.off()

###CedVsG
pdf("CedVsG_cutoff.pdf",width=11,height=8.5)
h = hist(CedVsG, breaks=100, prob=T, main="ALL CedVsG", xlab="Luminex, MFI values",
         ylab="Frequency", col="lightgrey") #histogram of all data minus positive controls
fitdistr(CedVsG/1000, "gamma") # use to estimate parmaters, scale data by dividing by 1000
f = function(theta,dta){
  l = -sum(dgamma(dta,rate=theta[1],shape=theta[2],log=TRUE))
}
o = optim( c(0.96, 4.14),  f, dta = CedVsG)
print(o)
p = qgamma(c(0.95, 0.99), rate=o$par[1],shape=o$par[2]); print(p)
abline(v=p, col=c("blue", "green"), lwd=2 )
lines( h$breaks, dgamma(h$breaks,rate=o$par[1],shape=o$par[2]), col="red",lwd=3,lty="44")
legend("topright", c("Gamma distribution", "95% cutoff (704)", "99% cutoff (1088)", "n=8 bats >99%"),
       fill=c("red", "blue", "green", "white"))
x <- CedVsG [CedVsG >1088] # figure out how many individual bats above cutoff
length (x) #of bats >99% cutoff, added to legend
dev.off()

###EboVGP
pdf("EboVGP_cutoff.pdf",width=11,height=8.5)
h = hist(EboVGP, breaks=100, prob=T, main="ALL EboVGP", xlab="Luminex, MFI values",
         ylab="Frequency", col="lightgrey") #histogram of all data minus positive controls
fitdistr(EboVGP/1000, "gamma") # use to estimate parmaters, scale data by dividing by 1000
f = function(theta,dta){
  l = -sum(dgamma(dta,rate=theta[1],shape=theta[2],log=TRUE))
}
o = optim( c(0.57, 0.75),  f, dta = EboVGP)
print(o)
p = qgamma(c(0.95, 0.99), rate=o$par[1],shape=o$par[2]); print(p)
abline(v=p, col=c("blue", "green"), lwd=2 )
lines( h$breaks, dgamma(h$breaks,rate=o$par[1],shape=o$par[2]), col="red",lwd=3,lty="44")
legend("topright", c("Gamma distribution", "95% cutoff (2152)", "99% cutoff (3223)", "n=19 bats >99%"),
       fill=c("red", "blue", "green", "white"))
x <- EboVGP [EboVGP >3223] # figure out how many individual bats above cutoff
length (x) #of bats >99% cutoff, added to legend
dev.off()

###MarVGP
pdf("MarVGP_cutoff.pdf",width=11,height=8.5)
h = hist(MarVGP, breaks=100, prob=T, main="ALL MarVGP", xlab="Luminex, MFI values",
         ylab="Frequency", col="lightgrey") #histogram of all data minus positive controls
fitdistr(MarVGP/1000, "gamma") # use to estimate parmaters, scale data by dividing by 1000

f = function(theta,dta){
  l = -sum(dgamma(dta,rate=theta[1],shape=theta[2],log=TRUE))
}
o = optim( c(2.81, 23.58),  f, dta = MarVGP)
print(o)
p = qgamma(c(0.95, 0.99), rate=o$par[1],shape=o$par[2]); print(p)
abline(v=p, col=c("blue", "green"), lwd=2 )
lines( h$breaks, dgamma(h$breaks,rate=o$par[1],shape=o$par[2]), col="red",lwd=3,lty="44")
legend("topright", c("Gamma distribution", "95% cutoff (254)", "99% cutoff (342)", "n=16 bats >99%"),
       fill=c("red", "blue", "green", "white"))
x <- MarVGP [MarVGP >342] # figure out how many individual bats above cutoff
length (x) #of bats >99% cutoff, added to legend
dev.off()

###MenVN
pdf("MenVN_cutoff.pdf",width=11,height=8.5)
h = hist(MenVN, breaks=100, prob=T, main="ALL MenVN", xlab="Luminex, MFI values",
         ylab="Frequency", col="lightgrey") #histogram of all data minus positive controls
fitdistr(MenVN/1000, "gamma") # use to estimate parmaters, scale data by dividing by 1000

f = function(theta,dta){
  l = -sum(dgamma(dta,rate=theta[1],shape=theta[2],log=TRUE))
}
o = optim( c(1.16, 0.51),  f, dta = MenVN)
print(o)
p = qgamma(c(0.95, 0.99), rate=o$par[1],shape=o$par[2]); print(p)
abline(v=p, col=c("blue", "green"), lwd=2 )
lines( h$breaks, dgamma(h$breaks,rate=o$par[1],shape=o$par[2]), col="red",lwd=3,lty="44")
legend("topright", c("Gamma distribution", "95% cutoff (6445)", "99% cutoff (9691)", "n=11 bats >99%"),
       fill=c("red", "blue", "green", "white"))
x <- MenVN [MenVN >9691] # figure out how many individual bats above cutoff
length (x) #of bats >99% cutoff, added to legend
dev.off()

###SARSN
pdf("SARSN_cutoff.pdf",width=11,height=8.5)
h = hist(SARSN, breaks=100, prob=T, main="ALL SARSN", xlab="Luminex, MFI values",
         ylab="Frequency", col="lightgrey") #histogram of all data minus positive controls
fitdistr(SARSN/1000, "gamma") # use to estimate parmaters, scale data by dividing by 1000

f = function(theta,dta){
  l = -sum(dgamma(dta,rate=theta[1],shape=theta[2],log=TRUE))
}
o = optim( c(1.94, 9.22),  f, dta = SARSN)
print(o)
p = qgamma(c(0.95, 0.99), rate=o$par[1],shape=o$par[2]); print(p)
abline(v=p, col=c("blue", "green"), lwd=2 )
lines( h$breaks, dgamma(h$breaks,rate=o$par[1],shape=o$par[2]), col="red",lwd=3,lty="44")
legend("topright", c("Gamma distribution", "95% cutoff (502)", "99% cutoff (706)", "n=12 bats >99%"),
       fill=c("red", "blue", "green", "white"))
x <- SARSN [SARSN >706] # figure out how many individual bats above cutoff
length (x) #of bats >99% cutoff, added to legend
dev.off()

###MERSN
pdf("MERSN_cutoff.pdf",width=11,height=8.5)
h = hist(MERSN, breaks=100, prob=T, main="ALL MERSN", xlab="Luminex, MFI values",
         ylab="Frequency", col="lightgrey") #histogram of all data minus positive controls
fitdistr(MERSN/1000, "gamma") # use to estimate parmaters, scale data by dividing by 1000

f = function(theta,dta){
  l = -sum(dgamma(dta,rate=theta[1],shape=theta[2],log=TRUE))
}
o = optim( c(0.93, 1.63),  f, dta = MERSN)
print(o)
p = qgamma(c(0.95, 0.99), rate=o$par[1],shape=o$par[2]); print(p)
abline(v=p, col=c("blue", "green"), lwd=2 )
lines( h$breaks, dgamma(h$breaks,rate=o$par[1],shape=o$par[2]), col="red",lwd=3,lty="44")
legend("topright", c("Gamma distribution", "95% cutoff (1757)", "99% cutoff (2729)", "n=13 bats >99%"),
       fill=c("red", "blue", "green", "white"))
x <- MERSN [MERSN >2729] # figure out how many individual bats above cutoff
length (x) #of bats >99% cutoff, added to legend
dev.off()


#compare MERS vs SARS seropositive individuals
pdf("MERSvsSARS_individualbats.pdf",width=11,height=8.5)
plot(sb$MERSN, sb$SARSN, ylab="SARS-N, MFI", xlab="MERS-N, MFI", main="Individual bat values (all data), MERS vs SARS")
abline(a=0, b=1, col="red")
cor.test(sb$MERSN, sb$SARSN)
text(1000,6000, label="r=0.529662", cex=1.3)
dev.off()

################### Comparing MFI values from first vs. second test runs ##################

pdf("NiV_HeV_1stvs2ndTest.pdf",width=11,height=8.5)
par(mfrow=c(2,2), oma=c(0,0,2,0))

plot(sb$NiVN, sb$R.NiVN, ylab="Secondary test, MFI", xlab="Primary test, MFI", main="NiV N, All Data")
abline(a=0, b=1, col="red")
cor.test(sb$NiVN, sb$R.NiVN)
text(1000,8000, label="r=0.9839", cex=1.3)

plot(sb$NiVsF, sb$R.NiVsF, ylab="Secondary test, MFI", xlab="Primary test, MFI", main="NiV sF, All Data")
abline(a=0, b=1, col="red")
cor.test(sb$NiVsF, sb$R.NiVsF)
text(800,7500, label="r=0.8504", cex=1.3)

plot(sb$NiVsG, sb$R.NiVsG, ylab="Secondary test, MFI", xlab="Primary test, MFI", main="NiV sG, All Data")
abline(a=0, b=1, col="red")
cor.test(sb$NiVsG, sb$R.NiVsG)
text(1000,5000, label="r=0.810", cex=1.3)

plot(sb$HeVsG, sb$R.HeVsG, ylab="Secondary test, MFI", xlab="Primary test, MFI", main="HeV sG, All Data")
abline(a=0, b=1, col="red")
cor.test(sb$HeVsG, sb$R.HeVsG)
text(2000,20000, label="r=0.994", cex=1.3)

title("Compare initial vs. retested MFI values for same individuals", outer=TRUE)
dev.off()

## MenV and CedV
pdf("MenVvsCedV_1stvs2ndTest.pdf",width=11,height=8.5)
par(mfrow=c(1,2))

plot(sb$CedVsG, sb$R.CedVsG, ylab="Secondary test, MFI", xlab="Primary test, MFI", main="CedV sG, All Data")
abline(a=0, b=1, col="red")
cor.test(sb$CedVsG, sb$R.CedVsG)
text(1000,8000, label="r=0.924", cex=1.3)

plot(sb$MenVN, sb$R.MenVN, ylab="Secondary test, MFI", xlab="Primary test, MFI", main="MenV N, All Data")
abline(a=0, b=1, col="red")
cor.test(sb$MenVN, sb$R.MenVN)
text(1000,15000, label="r=0.964", cex=1.3)
dev.off()
#EBOV and MARV
pdf("EBOVvsMARV_1stvs2ndTest.pdf",width=11,height=8.5)
par(mfrow=c(1,2))
names(sb)
summary(sb$EboVGP)
names(sb)
plot(sb$EboVGP, sb$R.EbovGP, ylab="Secondary test, MFI", xlab="Primary test, MFI",
     main="EboV GP, All Data", xlim=c(0,25000))
abline(a=0, b=1, col="red")

plot(sb$MarVGP, sb$R.MarVGP, ylab="Secondary test, MFI", xlab="Primary test, MFI",
     main="MarV GP, All Data", xlim=c(0,1500))
abline(a=0, b=1, col="red")
dev.off()

#########SPECIES SPECIFIC ###########
unique(sb$Species)  #generate data sets for each species seperately
R = sb[sb$Species=="Rousettus leschenaultii", ]
C = sb[sb$Species=="Cynopterus sphinx", ]
M = sb[sb$Species=="Megaderma lyra", ]
