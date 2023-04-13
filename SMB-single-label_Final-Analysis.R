rm(list = ls())

setwd('/path/to/files/0uM/200ms/')

load(file = 'r1/Refined-Particle_Data.RData')
e1.0uM.200=residence.data$Residence
load(file = 'r2/Refined-Particle_Data.RData')
e1.0uM.200=c(e1.0uM.200,residence.data$Residence)
load(file = 'r3/Refined-Particle_Data.RData')
e1.0uM.200=c(e1.0uM.200,residence.data$Residence)
load(file = 'r4/Refined-Particle_Data.RData')
e1.0uM.200=c(e1.0uM.200,residence.data$Residence)
example.trace=refined.particle.trace.rolls[,1]
ex.time=seq(0.2,1487*0.2,0.2)/60

setwd('~/path/to/files/0uM/500ms/')

load(file = 'r1/Refined-Particle_Data.RData')
e1.0uM.500=residence.data$Residence
load(file = 'r2/Refined-Particle_Data.RData')
e1.0uM.500=c(e1.0uM.500,residence.data$Residence)
load(file = 'r3/Refined-Particle_Data.RData')
e1.0uM.500=c(e1.0uM.500,residence.data$Residence)
load(file = 'r4/Refined-Particle_Data.RData')
e1.0uM.500=c(e1.0uM.500,residence.data$Residence)

setwd('/path/to/files/10uM/200ms/')

load(file = 'r1/Refined-Particle_Data.RData')
e1.10uM.200=residence.data$Residence
load(file = 'r2/Refined-Particle_Data.RData')
e1.10uM.200=c(e1.10uM.200,residence.data$Residence)
load(file = 'r3/Refined-Particle_Data.RData')
e1.10uM.200=c(e1.10uM.200,residence.data$Residence)
load(file = 'r4/Refined-Particle_Data.RData')
e1.10uM.200=c(e1.10uM.200,residence.data$Residence)

setwd('/path/to/files/10uM/500ms/')

load(file = 'r1/Refined-Particle_Data.RData')
e1.10uM.500=residence.data$Residence
load(file = 'r2/Refined-Particle_Data.RData')
e1.10uM.500=c(e1.10uM.500,residence.data$Residence)
load(file = 'r3/Refined-Particle_Data.RData')
e1.10uM.500=c(e1.10uM.500,residence.data$Residence)
load(file = 'r4/Refined-Particle_Data.RData')
e1.10uM.500=c(e1.10uM.500,residence.data$Residence)

setwd('/path/to/files/0uM/200ms/')

load(file = 'r1/Refined-Particle_Data.RData')
e2.0uM.200=residence.data$Residence
load(file = 'r2/Refined-Particle_Data.RData')
e2.0uM.200=c(e2.0uM.200,residence.data$Residence)
load(file = 'r3/Refined-Particle_Data.RData')
e2.0uM.200=c(e2.0uM.200,residence.data$Residence)
load(file = 'r4/Refined-Particle_Data.RData')
e2.0uM.200=c(e2.0uM.200,residence.data$Residence)

setwd('/path/to/files/0uM/500ms/')

load(file = 'r1/Refined-Particle_Data.RData')
e2.0uM.500=residence.data$Residence
load(file = 'r2/Refined-Particle_Data.RData')
e2.0uM.500=c(e2.0uM.500,residence.data$Residence)
load(file = 'r3/Refined-Particle_Data.RData')
e2.0uM.500=c(e2.0uM.500,residence.data$Residence)
load(file = 'r4/Refined-Particle_Data.RData')
e2.0uM.500=c(e2.0uM.500,residence.data$Residence)

setwd('/path/to/files/10uM/200ms/')

load(file = 'r1/Refined-Particle_Data.RData')
e2.10uM.200=residence.data$Residence
load(file = 'r2/Refined-Particle_Data.RData')
e2.10uM.200=c(e2.10uM.200,residence.data$Residence)
load(file = 'r3/Refined-Particle_Data.RData')
e2.10uM.200=c(e2.10uM.200,residence.data$Residence)
load(file = 'r4/Refined-Particle_Data.RData')
e2.10uM.200=c(e2.10uM.200,residence.data$Residence)

setwd('/path/to/files/10uM/500ms/')

load(file = 'r1/Refined-Particle_Data.RData')
e2.10uM.500=residence.data$Residence
load(file = 'r2/Refined-Particle_Data.RData')
e2.10uM.500=c(e2.10uM.500,residence.data$Residence)
load(file = 'r3/Refined-Particle_Data.RData')
e2.10uM.500=c(e2.10uM.500,residence.data$Residence)
load(file = 'r4/Refined-Particle_Data.RData')
e2.10uM.500=c(e2.10uM.500,residence.data$Residence)


par(fig = c(0.01,0.4,0.6,1))
usable=(ex.time>=2 & ex.time<=4)
event=c(808/300,876/300)
plot(ex.time[usable],example.trace[usable],type = 'l',col='grey',main = 'Representative Particle Trace',xlab = 'Time (min)',ylab = 'Signal (A.U.)',cex.lab=1.4,cex.main=1.4,cex.axis=1.4)
lines(c(0,event[1]),rep(0,times=2),col='black',lwd=3)
lines(event,rep(7,times=2),col='green',lwd=3)
lines(c(event[2],5),rep(0,times=2),col='black',lwd=3)
legend('topright',legend = c('Bound','Unbound'),col = c('green','black'),fill = c('green','black'),cex=1)
arrows(x0 = 3,x1 = 5.6,y0 = 12e3,y1 = 12e3,code = 3,lwd = 1,angle = 90,length = 0.15)
mtext(expression(tau*' = 13.6 s'),side = 3,line = -10.5,at = 2.4,cex = 1.3)


par(fig = c(0.01,1,0.205,0.6),new=TRUE)
fit.0=MASS::fitdistr(c(e1.0uM.200,e1.0uM.500,e2.0uM.200,e2.0uM.500),densfun = 'exponential')
fit.0sd=sd(c((MASS::fitdistr(c(e1.0uM.200),densfun = 'exponential'))$estimate,(MASS::fitdistr(c(e1.0uM.500),densfun = 'exponential'))$estimate,(MASS::fitdistr(c(e2.0uM.200),densfun = 'exponential'))$estimate,(MASS::fitdistr(c(e2.0uM.500),densfun = 'exponential'))$estimate))
hist.0=hist(c(e1.0uM.200,e1.0uM.500,e2.0uM.200,e2.0uM.500),plot = FALSE,breaks = 20)
fit.10=MASS::fitdistr(c(e1.10uM.200,e1.10uM.500,e2.10uM.200,e2.10uM.500),densfun = 'exponential')
fit.10sd=sd(c((MASS::fitdistr(c(e1.10uM.200),densfun = 'exponential'))$estimate,(MASS::fitdistr(c(e1.10uM.500),densfun = 'exponential'))$estimate,(MASS::fitdistr(c(e2.10uM.200),densfun = 'exponential'))$estimate,(MASS::fitdistr(c(e2.10uM.500),densfun = 'exponential'))$estimate))
hist.10=hist(c(e1.10uM.200,e1.10uM.500,e2.10uM.200,e2.10uM.500),plot = FALSE,breaks = 20)
plot(hist.0$breaks,c(hist.0$counts,hist.0$counts[length(hist.0$counts)])/max(hist.0$counts),type ='s',lwd=1,col='blue',xlab = 'Residence Time (s)',ylab = 'Relative Likelihood',main = expression('Residence Time Distributions: TREX1 + [Cy5]d(N)'[5]*' + d(N)'[5]),ylim=0:1,xlim = c(0,100),cex.main=2,cex.lab=1.5,cex.axis=1.5)
lines(hist.10$breaks,c(hist.10$counts,hist.10$counts[length(hist.10$counts)])/max(hist.10$counts),type='s',lwd=1,col='red')
lines(seq(0,300,0.1),dexp(seq(0,300,0.1),rate = fit.10$estimate)/max(dexp(seq(0,300,0.1),rate = fit.10$estimate)),col='red',lwd=5)
lines(seq(0,300,0.1),dexp(seq(0,300,0.1),rate = fit.0$estimate)/max(dexp(seq(0,300,0.1),rate = fit.0$estimate)),col='blue',lwd=5)
legend('topright',legend = c('Competitor = 0 µM','Competitor = 10 µM'),fill = c('blue','red'),col= c('blue','red'),cex = 1.5)
mtext(paste0('n = ',c(length(c(e1.0uM.200,e1.0uM.500,e2.0uM.200,e2.0uM.500)),length(c(e1.10uM.200,e1.10uM.500,e2.10uM.200,e2.10uM.500)))),side = 3,line = c(-2,-4),at = 60,cex = 2,col = c('blue','red'))

par(fig = c(0.01,1,0,0.3),new=TRUE)
boxplot(e1.0uM.200,e1.0uM.500,e2.0uM.200,e2.0uM.500,e1.10uM.200,e1.10uM.500,e2.10uM.200,e2.10uM.500,
        at = c(4:-3),
        names = paste0(rep(rep(c('A-','B-'),each=2),times=2),rep(c('0.2s','0.5s'),times=4)),
        las = 2,
        col = rep(c('blue','red'),each=4),
        border = rep(c('blue','red'),each=4),
        range=0,
        varwidth = TRUE,
        horizontal = TRUE,
        notch = TRUE,
        xaxt='n',
        yaxt='n',
        ylim=c(0,100),
        axes=FALSE
)
mtext(paste0(rep(rep(c('A - ','B - '),each=2),times=2),rep(c('0.2s','0.5s'),times=4)),side = 1,line = rev(-1*seq(1.2,6,(6-1.2)/7)),at = -5,cex = 0.8,col = rep(c('blue','red'),each=4))


par(fig = c(0.405,1,0.65,1),new=TRUE)
plot(NULL,NULL,xlim = c(-1,11),ylim = c(0,max(1/c(mean(e1.10uM.200),mean(e1.10uM.500),mean(e2.10uM.200),mean(e2.10uM.500)))),xaxt='n',xlab='',ylab='',main = expression('TREX1 + [Cy5]d(N)'[5]*' + d(N)'[5]),axes=FALSE,cex.main=2)
axis(side = 2,at = seq(0,20,5)*1e-2,cex.axis=1.5)
mtext(expression('k'['off']*''^obs*' (s'^-1*')'),side = 2,line = 2.5,at = 1e-1,cex = 2)
axis(side = 1,at = c(0,10),labels = c(0,10),cex.axis=1.5)
mtext(paste0('[Competitor] (µM)'),side = 1,line = 2.5,at = 5,cex = 1.6)
points(c(0,10),c(fit.0$estimate,fit.10$estimate),col=c('blue','red'),pch = '-',cex = 5)
arrows(x0 = c(0,10),x1 = c(0,10),y0 = c(fit.0$estimate-fit.0sd,fit.10$estimate-fit.10sd), y1 = c(fit.0$estimate+fit.0sd,fit.10$estimate+fit.10sd),code = 3,col = c('blue','red'),lwd = 1,angle = 90,length = 0.5)
temp.1=signif(fit.0$estimate,2)
temp.2=signif(diff(c(fit.0$estimate,fit.10$estimate))/10e-6,2)
mtext(substitute(paste('k'[-1],' = ',temp.1,' (s'^-1,')'),list(temp.1 = temp.1)),side = 1,line = 6,at = 1.5,cex = 1.6)
mtext(substitute(paste('k'[theta],' = ',temp.2,' (M'^-1,'s'^-1,')'),list(temp.2 = temp.2)),side = 1,line = 6,at = 8.5,cex = 1.6)

