##### SCRIPT INFORMATION

# SMB_Simulator.R

# A script to simulate SMB data under a direct transfer kinetic model

##########################

##### DEVELOPER NOTES

# none

##########################

##### SCRIPT USER INPUT

### Environment

rm(list = ls())
setwd('~/path/to/file/')

### Key variables

# Experiment parameters
time=60*5 # total reaction time, s
dt=150e-3 # exposure time-step, s
n=2558 # number of particle traces to simulate (2558)
prots=2 # number of binding sites in functional protein
m=100 # number of data sets to simulate

# Model constants
Kd=8.9e-9 # dissociation constant, M
kn1=1.1e-2 # dissociation rate constant, 1/s
ktheta=1500 # direct transfer rate constant, 1/M/s

# Initial values
Et=5e-12 # total predator, M
L1t=10e-9 # total ligand 1, M
L2t=10e-9 # total ligand 2, M

# Plotting parameters
err=0.1 # relative error to be introduced into particle traces
plot.hist=T # would you like a histogram of ASIs, T/F
picks='dt' # which particle traces would you like to plot, all/1/2/co/dt/none
save.pdf=F # would you like to export graphs to a pdf file, T/F
save.data=T # would you like to export data to an RData file, T/F
obs.N=36 # how many transfer events did you see in a comparable experiment

##########################

##### BEGIN SCRIPT

### Semi-autonomous script functions

# Remaining model constants
k1=kn1/Kd # association rate constant, 1/M/s

# Initial values: advanced
E.0=Et-(Et+(L1t+L2t)+Kd-sqrt((Et+(L1t+L2t)+Kd)^2-4*Et*(L1t+L2t)))/2 # initial free predator, M
L1.0=L1t-L1t/(L1t+L2t)*(Et+(L1t+L2t)+Kd-sqrt((Et+(L1t+L2t)+Kd)^2-4*Et*(L1t+L2t)))/2 # initial free decoy, M
L2.0=L2t-L2t/(L1t+L2t)*(Et+(L1t+L2t)+Kd-sqrt((Et+(L1t+L2t)+Kd)^2-4*Et*(L1t+L2t)))/2 # initial free prey, M
EL1.0=L1t/(L1t+L2t)*(Et+(L1t+L2t)+Kd-sqrt((Et+(L1t+L2t)+Kd)^2-4*Et*(L1t+L2t)))/2 # initial predator-decoy complex, M
EL2.0=L2t/(L1t+L2t)*(Et+(L1t+L2t)+Kd-sqrt((Et+(L1t+L2t)+Kd)^2-4*Et*(L1t+L2t)))/2 # initial predator-prey complex, M

# Load accessory packages
#none#

# Load custom functions
bind.change <- function(state,probs){
  if(state==-1){
    new.state=sample(c(-1,0,1),1,T,c(1-probs[['D']]-probs[['T.21']],probs[['D']],probs[['T.21']]))
  }
  if(state==0){
    new.state=sample(c(-1,0,1),1,T,c(probs[['B.L2']],1-probs[['B.L1']]-probs[['B.L2']],probs[['B.L1']]))
  }
  if(state==1){
    new.state=sample(c(-1,0,1),1,T,c(probs[['T.12']],probs[['D']],1-probs[['D']]-probs[['T.21']]))
  }
  return(new.state)
}
iter.diff <- function(a,b){
  Ma=matrix(a,nrow = length(a),ncol = length(b),byrow = F)
  Mb=matrix(b,nrow = length(a),ncol = length(b),byrow = T)
  diffs=abs(c(Ma-Mb))
  return(diffs)
}
change.dist <- function(jump,t){
  up1=t[which(jump[,1]>0)]
  down1=t[which(jump[,1]<0)]
  up2=t[which(jump[,2]>0)]
  down2=t[which(jump[,2]<0)]
  dists=NA
  if(length(up1)>0 & length(down2)>0){
    dists=append(dists,iter.diff(up1,down2))
  }
  if(length(up2)>0 & length(down1)>0){
    dists=append(dists,iter.diff(up2,down1))
  }
  if(length(dists)>1){
    dists=dists[-1]
  }
  return(dists)
}
c.pd <- function(input,range){
  apply(X = as.array(range[1]:range[2]),MARGIN = 1,FUN = function(x,input){sum(input<=x)/length(input)},input=input)
}

### Autonomous script functions

# Save simulation parameters
constants=list('k1'=k1,'kn1'=kn1,'ktheta'=ktheta)
initials=list('E'=E.0,'L1'=L1.0,'L2'=L2.0,'EL1'=EL1.0,'EL2'=EL2.0)
probs=list('B.L1'=pexp(dt,k1*L1.0),'B.L2'=pexp(dt,k1*L2.0),'D'=pexp(dt,kn1),'T.12'=pexp(dt,ktheta*L2.0),'T.21'=pexp(dt,ktheta*L1.0))
parameters=mget(ls())
if (m>1){STATS=list(NULL)}

for (m.i in 1:m){

# Generate initial data variables
t=seq(0,time,dt)
data=array(NA,dim=c(length(t),n,prots))
for (i in 1:prots){
  data[1,,i]=sample(c(-1,0,1),n,T,c(EL2.0,E.0,EL1.0))
}

# Simulate binding states
for (i in 2:length(t)){
  data[i,,]=apply(data[i-1,,],1:2,FUN = bind.change,probs=probs)
}

# Transform to particle traces
trace.1=apply((data>0),MARGIN = c(1,2),FUN = sum)+array(err*rnorm(n*length(t)),dim=c(length(t),n))
trace.2=apply((data<0),MARGIN = c(1,2),FUN = sum)+array(err*rnorm(n*length(t)),dim=c(length(t),n))

# Calculate times between state changes
changes.1=rbind(rep(NA,times=ncol(apply((data>0),MARGIN = c(1,2),FUN = sum))),((apply((data>0),MARGIN = c(1,2),FUN = sum))[-1,]-(apply((data>0),MARGIN = c(1,2),FUN = sum))[-length(t),]))
changes.2=rbind(rep(NA,times=ncol(apply((data<0),MARGIN = c(1,2),FUN = sum))),((apply((data<0),MARGIN = c(1,2),FUN = sum))[-1,]-(apply((data<0),MARGIN = c(1,2),FUN = sum))[-length(t),]))
interstate.times=NA
for (i in 1:ncol(changes.1)){
  interstate.times=append(interstate.times,change.dist(cbind(changes.1[,i],changes.2[,i]),t))
}
ITs=na.omit(interstate.times)
for (i in 1:5){
  if(m==1){show(paste0('From ',n,' total traces, among ',sum(((colSums(apply((data>0),MARGIN = c(1,2),FUN = sum))>=1) & (colSums(apply((data<0),MARGIN = c(1,2),FUN = sum))>=1))),' co-binding traces, there are ',sum(ITs<=((i-1)*dt)),' anti-correlated state changes within ',i*dt,'s'))}
}
if (m>1){
  STATS[[m.i]]=ITs
}

# Prepare graph file output
if (save.pdf==T & (m==1 | (m>1 & m.i==m))){
  pdf(file = 'GRAPHS.pdf',paper = 'letter')
}

# Plot distribution of times between anti-correlated state changes
if (plot.hist==T & m==1){
  par(mfrow=c(2,1))
  hist(ITs,xlab='Time Between Anti-Correlated State Changes (s)',main=paste0('Distribution: kθ = ',ktheta))
}

# Plot particle traces
if (picks=='all'){
  plot.me=rep(T,times=ncol(trace.1))
}
if (picks=='1'){
  plot.me=(colSums(apply((data>0),MARGIN = c(1,2),FUN = sum))>=1)
}
if (picks=='2'){
  plot.me=(colSums(apply((data<0),MARGIN = c(1,2),FUN = sum))>=1)
}
if (picks=='co'){
  plot.me=((colSums(apply((data>0),MARGIN = c(1,2),FUN = sum))>=1) & (colSums(apply((data<0),MARGIN = c(1,2),FUN = sum))>=1))
}
if (picks=='dt'){
  plot.me=(colSums(((changes.1>0 & changes.2<0) | (changes.1<0 & changes.2>0)),na.rm = TRUE)>0)
}
par(mfrow=c(6,3),mar=c(4.5,5,2,1))
for (i in 1:ncol(trace.1)){
  if(plot.me[i]==T & picks!='none' & m==1){
    plot(NULL,NULL,main=paste0(i,' (DTs = ',(colSums(((changes.1>0 & changes.2<0) | (changes.1<0 & changes.2>0)),na.rm = TRUE))[i],')'),xlab='time (s)',ylab='Signal (AU)',xlim=range(t),ylim=c(0,prots))
    lines(t,trace.1[,i],col='green')
    lines(t,trace.2[,i],col='red')
    temp=(((changes.1>0 & changes.2<0) | (changes.1<0 & changes.2>0))[,i]);temp[is.na(temp)]=F;points(t[temp],rep(prots,times=(colSums(((changes.1>0 & changes.2<0) | (changes.1<0 & changes.2>0)),na.rm = TRUE))[i]),pch='v',cex=3)
  }
}

# Plot probability density plots
if (m>1 & m.i==m){
  STATS.summ=matrix(NA,nrow=m,ncol = 6)
  for (i in 1:ncol(STATS.summ)){
    STATS.summ[,i]=apply(X = as.array(1:m),MARGIN = 1,FUN = function(x,dt,i,STATS){sum(STATS[[x]]<=((i-1)*dt))},dt=dt,i=i,STATS=STATS)
  }
  STATS.probs=(apply(X = STATS.summ,MARGIN = 2,FUN = c.pd,range=c(0,max(STATS.summ)+1)))
  STATS.probs.x=seq(0,max(STATS.summ)+1,1)
  p.05=apply(STATS.probs,2,function(x){((1:length(x))[x>=0.975])[1]})
  p.05b=apply(STATS.probs,2,function(x){(rev((1:length(x))[x<=0.025]))[1]})
  #
  par(mfrow=c(2,1),mar=c(4.5,5,2,1))
  plot(NULL,NULL,main='Probability of ≤N DTs in Experiment',xlab='N',ylab='Probability',xlim=c(0,max(STATS.summ)+1),ylim=0:1,cex.main=2,cex.lab=2,cex.axis=2)
  abline(v=obs.N,col='orange',lwd=4,lty='dashed')
  for (i in 1:ncol(STATS.summ)){
    lines(STATS.probs.x,STATS.probs[,i],col=i,lwd=2,type='s')
  }
  legend('bottomright',legend = paste0('DT∆ ≤ ',dt*(1:ncol(STATS.summ)),' s'),fill = 1:ncol(STATS.summ),col = 1:ncol(STATS.summ),cex=1.5)
  text(x = 0,y = 1,labels = paste0('N = ',obs.N),col='orange',cex=3,adj = c(0,1))
  #
  plot(NULL,NULL,main='Expected N (p > 0.05) for DT∆ Thresholds',xlab='DT∆ (s)',ylab='N',xlim=range(dt*(1:length(p.05))),ylim=c(0,max(p.05)),cex.main=2,cex.lab=2,cex.axis=2)
  arrows(dt*(1:length(p.05)),p.05b,y1 = p.05,angle = 90,code = 3,lwd=5)
  abline(h=obs.N,col='orange',lwd=4,lty='dashed')
  text(x = min(dt*(1:length(p.05))),y = max(p.05),labels = paste0('N = ',obs.N),col='orange',cex=3,adj = c(0,1))
}

# Initiate file output
if (save.pdf==T){
  dev.off()
}
if (save.data==T & m==1){
  save(list=ls(),file = 'DATA.RData')
}
if (save.data==T & m>1 & m.i==m){
  save(list=ls(),file = 'bigDATA.RData')
}

}



##### END SCRIPT




