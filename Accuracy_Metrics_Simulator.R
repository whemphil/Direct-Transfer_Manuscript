##### SCRIPT INFORMATION

# Accuracy-Metrics_Simulator.R

# A script to simulate FP data under a displacement-transfer kinetic model

# Associated with manuscript "Direct Ligand Transfer in TREX1 and Other Proteins: Implications for Nucleic Acid Binding Proteins" 

##########################

##### DEVELOPER NOTES

# The manuscript this code is associated with only utilizes variation in the 'Environment' and 'Key variables' sections of SCRIPT USER INPUT

##########################

##### SCRIPT USER INPUT

### Environment

rm(list = ls())
setwd('~/path/to/output/')

# Accessory packages
library(deSolve)
library(FPalyze) # github.com/whemphil/FPalyze

### Key variables

# Reaction parameters
time=60*60*2 # total reaction time, s

# Model constants
Kd.P=10e-9 # prey dissociation constant, M
Kd.D=10e-9 # decoy dissociation constant, M
kn1.P=1e-3 # prey dissociation rate constant, 1/s
kn1.D=1e-3 # decoy dissociation rate constant, 1/s
ktheta.P=0 # prey displacement-transfer rate constant, 1/M/s
ktheta.D=0 # decoy displacement-transfer rate constant, 1/M/s

# Initial values
Et=c(3*Kd.P/1e-9)*1e-9 # total predator, M
Dt=c(1e4/2^c(0:10),0)*1e-9 # total decoy, M
Pt=c(5)*1e-9 # total prey, M

# Data error and sampling
feaux.start=60*1 # start-read delay of fake export data, s
read.freq=30 # read frequency of fake export data, s
experimental.error=0.05 # proportion of variance to simulate data variance, 0-0.3 suggested
count=20 # number of fake data sets to generate

# Analysis parameters
manual.baseline=F # should true baseline signal be provided to analysis software

### Extraneous variables

# Simulation parameters
dt=1e0 # integration time-step, s

# Reaction parameters
pre.equilibrated=T # start from predator-prey equilibrium binding state, T/F

# Data polarization value adjustment parameters
kT=5e-4 # decay constant for relative change in mP values versus temperature/time, 1/s
Mx.H=100 # highest mP value for saturated prey binding, mP
Mx.L=100 # lowest mP value for saturated prey binding, mP
Mn.H=0 # highest mP value for unbound prey, mP
Mn.L=0 # lowest mP value for unbound prey, mP

# Data error and sampling
fake.data.type='COMP' # what type of fake data are you simulating, Kd/COMP/none
data.size='full' # is the experiment size half (2 replicates) or full (4 replicates), half/full

# Saved data features
fake.id='NULLtest' # id tag for the fake data replicates generated
save.id='Sim-NULL' # identifier for save files
save.pdf='no' # decision to save pdf of graphs, yes/no
plotting=T # should plots be produced, T/F

##########################

##### BEGIN SCRIPT

### Semi-autonomous script functions

# Remaining model constant calculations
k1.P=kn1.P/Kd.P # prey association rate constant, 1/M/s
k1.D=kn1.D/Kd.D # decoy association rate constant, 1/M/s

# Miscellaneous additions
if (fake.data.type!='none' & data.size=='full'){
  Feaux.Data=matrix(0,nrow=length(seq(0,time+feaux.start,dt)),ncol = 49)
  colnames(Feaux.Data)=c('Time',paste(rep(c('A','B'),each=24),rep(1:24,times=2),sep = ''))
}
if (fake.data.type!='none' & data.size=='half'){
  Feaux.Data=matrix(0,nrow=length(seq(0,time+feaux.start,dt)),ncol = 25)
  colnames(Feaux.Data)=c('Time',paste(rep(c('A'),each=24),rep(1:24,times=1),sep = ''))
}

# Initial values: advanced
if (pre.equilibrated==T){
  E.0=Et-(Et+Pt+Kd.P-sqrt((Et+Pt+Kd.P)^2-4*Et*Pt))/2 # initial free predator, M
  D.0=Dt # initial free decoy, M
  P.0=Pt-(Et+Pt+Kd.P-sqrt((Et+Pt+Kd.P)^2-4*Et*Pt))/2 # initial free prey, M
  ED.0=0e-9 # initial predator-decoy complex, M
  EP.0=(Et+Pt+Kd.P-sqrt((Et+Pt+Kd.P)^2-4*Et*Pt))/2 # initial predator-prey complex, M
}
if (pre.equilibrated==F){
  E.0=Et # initial free predator, M
  D.0=Dt # initial free decoy, M
  P.0=rep(Pt,times=length(Et)) # initial free prey, M
  ED.0=0e-9 # initial predator-decoy complex, M
  EP.0=rep(0e-9,times=length(Et)) # initial predator-prey complex, M
}

### Autonomous script functions

# Save simulation parameters
parameters=mget(ls())

# Plot bins
if (save.pdf=='yes'){
  pdf(file = paste(save.id,'.pdf',sep=''))
}
if (length(D.0)>1 & plotting==T){
  par(mfrow=c(round(length(D.0)/2+0.1),2),mar=c(2,2,2,2))
}
if (length(D.0)==1 & plotting==T){
  par(fig = c(0,1,0,1))
}

# Simulation results storage
SIM.Results=list()

# RNA range
for (j in 1:length(D.0)){

  # Plot frame
  if (plotting==T){
    plot(NULL,NULL,main = paste('Predator-Prey Binding Curve: [Decoy] = ',D.0[j]*1e9,' nM',sep = ''),xlim = c(0,time/60),ylim = c(0,1),xlab = 'Time (min)',ylab = 'Fraction Prey Bound')
  }
  
  # Enzyme range
  for (k in 1:length(E.0)){
  
    # Record times
    t=seq(0,time+feaux.start,dt)

    # Record constants
    constants=c(k1.P=k1.P,kn1.P=kn1.P,k1.D=k1.D,kn1.D=kn1.D,ktheta.P=ktheta.P,ktheta.D=ktheta.D)

    # Record initial values
    initials=c(E=E.0[k], D=D.0[j], P=P.0[k], ED=ED.0, EP=EP.0[k])

    # Record Equations
    Equations=function(t,initials,constants){
        with(as.list(c(initials,constants)),{
            
        dE=kn1.P*EP+kn1.D*ED-k1.P*E*P-k1.D*E*D
        dD=kn1.D*ED+ktheta.P*ED*P-k1.D*E*D-ktheta.D*EP*D
        dP=kn1.P*EP+ktheta.D*EP*D-k1.P*E*P-ktheta.P*ED*P
        dED=k1.D*E*D+ktheta.D*EP*D-kn1.D*ED-ktheta.P*ED*P
        dEP=k1.P*E*P+ktheta.P*ED*P-kn1.P*EP-ktheta.D*EP*D

        list(c(dE,dD,dP,dED,dEP))
        })
    }

    # Numerical evaluation
    sim.sim=as.data.frame(ode(y=initials,times=t,func=Equations,parms=constants))
    
    # Clean up results
    RXN.Results=list('t'=sim.sim$time,'E'=sim.sim$E,'D'=sim.sim$D,'P'=sim.sim$P,'ED'=sim.sim$ED,'EP'=sim.sim$EP,'mP'=(1-sim.sim$P/Pt)*(exp(-sim.sim$time*kT)*(Mx.H+Mn.L-Mx.L-Mn.H)+Mx.L-Mn.L)+exp(-sim.sim$time*kT)*(Mn.H-Mn.L)+Mn.L)
    SIM.Results[[paste(E.0[k],D.0[j],sep = '_')]]=RXN.Results
    
    # Plot results
    if (plotting==T){
      lines(RXN.Results[['t']][seq(1,time/dt+1,0.1/dt)]/60,1-(RXN.Results[['P']][seq(1,time/dt+1,0.1/dt)]/Pt),col=k)
    }
    
    if (fake.data.type=='Kd'){
      Feaux.Data[,c(1)]=t
      if (data.size=='full'){
        Feaux.Data[,c(2*(k-1)+2,2*(k-1)+3,2*(k-1)+26,2*(k-1)+27)]=rep(RXN.Results[['mP']],times=4)
      }
      if (data.size=='half'){
        Feaux.Data[,c(2*(k-1)+2,2*(k-1)+3)]=rep(RXN.Results[['mP']],times=2)
      }
    }
    
  }

  # Storing fake data
  if (fake.data.type=='COMP'){
    Feaux.Data[,c(1)]=t
    if (data.size=='full'){
      Feaux.Data[,c(2*(j-1)+2,2*(j-1)+3,2*(j-1)+26,2*(j-1)+27)]=rep(RXN.Results[['mP']],times=4)
    }
    if (data.size=='half'){
      Feaux.Data[,c(2*(j-1)+2,2*(j-1)+3)]=rep(RXN.Results[['mP']],times=2)
    }
  }
  
}

# Save simulation data
SIM.Results[['Parameters']]=parameters
if (fake.data.type!='none'){
  for (i in 1:count){
    SIM.Results[[paste0('Feaux Data ',i)]]=Feaux.Data[seq(1+feaux.start/dt,nrow(Feaux.Data),read.freq/dt),]+rnorm(length(Feaux.Data[seq(1+feaux.start/dt,nrow(Feaux.Data),read.freq/dt),]),0,experimental.error*(Mx.H-Mn.L))
    write.table(SIM.Results[[paste0('Feaux Data ',i)]],file = paste0('FakeData_',fake.id,'-',i,'.txt'),sep = '\t',col.names = TRUE,row.names = FALSE,quote = FALSE)
  }
}
save(SIM.Results,file=paste0(save.id,'.RData'))
if (save.pdf=='yes'){
  dev.off()
}
#rm(list=setdiff(ls(),'SIM.Results'))

# Run example analysis
temp=matrix(NA,nrow = count,ncol = 3)
for (i in 1:count){
  sink('NULLout.txt');suppressMessages(try(test.run<-FPalyze(experiment.type = fake.data.type,path.to.file = '',file.name = paste0('FakeData_',fake.id,'-',i,'.txt'),save.data = F,plot.pdf = T,save.console = F,time.step = read.freq,t.zero = feaux.start/60,incubation.time = time/60,default.mP.values = T,data.size = data.size,manual.fEP.adjust = manual.baseline,FP.baseline = Mn.L)));sink()
  if(exists('test.run')==T & fake.data.type=='COMP'){
    temp[i,1:2]=coefficients(test.run[['fun.model.opt2']])[c(1,4)]
    temp[i,3]=test.run[['delta.BIC']]
    rm(test.run)
  }
  if(exists('test.run')==T & fake.data.type=='Kd'){
    temp[i,1]=signif(summary(test.run[['model.std']])[['coefficients']][1,1],2)
    temp[i,2]=signif(summary(test.run[['model.hill']])[['coefficients']][1,1],2)
    temp[i,3]=signif(summary(test.run[['model.hill']])[['coefficients']][4,1],2)
    rm(test.run)
  }
}
temp=na.omit(temp)
if(ktheta.D==0 & fake.data.type=='COMP'){
  summaries=data.frame('Kd.P'=Kd.P,'Kd.D'=Kd.D,'kn1.D'=kn1.D,'ktheta.P'=ktheta.P,'kn1.P'=kn1.P,'ktheta.D'=ktheta.D,'kn1.avg'=mean(temp[,1]),'kn1.sd'=sd(temp[,1]),'ktheta.avg'=mean(temp[,2]),'ktheta.sd'=sd(temp[,2]),'BIC.acc'=round(mean(temp[,3]>=0),2),'kn1.close'=(abs(kn1.P-mean(temp[,1]))<=3*sd(temp[,1])),'ktheta.close'=(abs(ktheta.D-mean(temp[,2]))<=3*sd(temp[,2])))
}
if(ktheta.D!=0 & fake.data.type=='COMP'){
  summaries=data.frame('Kd.P'=Kd.P,'Kd.D'=Kd.D,'kn1.D'=kn1.D,'ktheta.P'=ktheta.P,'kn1.P'=kn1.P,'ktheta.D'=ktheta.D,'kn1.avg'=mean(temp[,1]),'kn1.sd'=sd(temp[,1]),'ktheta.avg'=mean(temp[,2]),'ktheta.sd'=sd(temp[,2]),'BIC.acc'=round(mean(temp[,3]<0),2),'kn1.close'=(abs(kn1.P-mean(temp[,1]))<=3*sd(temp[,1])),'ktheta.close'=(abs(ktheta.D-mean(temp[,2]))<=3*sd(temp[,2])))
}
if(fake.data.type=='Kd'){
  summaries=data.frame('Kd.P'=Kd.P/1e-9,'Kd.std.avg'=mean(temp[,1]),'Kd.std.sd'=sd(temp[,1]),'kd.hill.avg'=mean(temp[,2]),'kd.hill.sd'=sd(temp[,2]),'n.hill.avg'=mean(temp[,3]))
}
show(summaries)
show(binom.test(if(ktheta.D==0) sum(temp[,3]>=0) else sum(temp[,3]<0),length(temp[,3])))
View(summaries)

##### END SCRIPT









