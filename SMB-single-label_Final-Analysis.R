
rm(list = ls())

setwd('~/path/to/file/')

load(file = 'no-decoy_rep-1_refine.particles.RData')
null300.RT=residence.data$Residence
load(file = 'no-decoy_rep-2_refine.particles.RData')
null150.RT=residence.data$Residence
load(file = 'high-decoy_rep-1_refine.particles.RData')
decoy300.RT=residence.data$Residence
load(file = 'high-decoy_rep-2_refine.particles.RData')
decoy150.RT=residence.data$Residence

null300.fit=density(null300.RT)
null150.fit=density(null150.RT)
decoy300.fit=density(decoy300.RT)
decoy150.fit=density(decoy150.RT)
par(fig = c(0,1,0.6,1))
plot(null300.fit$x,null300.fit$y/max(null300.fit$y),type='l',xlab = 'Residence Time (s)',ylab = 'Relative Likelihood',main = expression('Residence Time Distributions: TREX1 + [Cy5]d(N)'[5]*' + d(N)'[5]),col='blue',ylim = c(0,1),xlim = c(0,400),cex.main=2)
lines(null150.fit$x,null150.fit$y/max(null150.fit$y),col='blue')
lines(decoy300.fit$x,decoy300.fit$y/max(decoy300.fit$y),col='red')
lines(decoy150.fit$x,decoy150.fit$y/max(decoy150.fit$y),col='red')
legend('topright',legend = c('Decoy = 0 µM','Decoy = 10 µM'),fill = c('blue','red'),col= c('blue','red'),cex = 1.5)
par(fig = c(0,1,0.4,0.7),new=TRUE)
boxplot(null300.RT,null150.RT,decoy300.RT,decoy150.RT,
        at = c(4:1),
        names = c("", "", "", ""),
        las = 2,
        col = c("blue","blue","red","red"),
        border = c("blue","blue","red","red"),
        range=0,
        varwidth = TRUE,
        horizontal = TRUE,
        notch = TRUE,
        xaxt='n',
        yaxt='n',
        ylim=c(0,400),
        axes=FALSE
)
mtext(c('Rep-1','Rep-2','Rep-1','Rep-2'),side = 1,line = rev(c(-1.5,-2.9,-4.5,-5.9)),at = -25,cex = 1.5,col = rev(c('red','red','blue','blue')))
par(fig = c(0.005,1,0.1,0.5),new=TRUE)
plot(NULL,NULL,xlim = c(-1,11),ylim = c(0,max(1/c(mean(null300.RT),mean(null150.RT),mean(decoy300.RT),mean(decoy150.RT)))),xaxt='n',xlab='',ylab='',main = expression('Apparent Dissociation Rates: TREX1 + [Cy5]d(N)'[5]*' + d(N)'[5]),axes=FALSE,cex.main=2)
axis(side = 2,at = c(0:3)*1e-2,cex.axis=1.2)
mtext(expression('k'['off']*' (s'^-1*')'),side = 2,line = 2.5,at = 1.5e-2,cex = 1.5)
axis(side = 1,at = c(0,10),labels = c(0,10),cex.axis=1.25)
mtext(paste0('[Decoy] (µM)'),side = 1,line = 2.5,at = 5,cex = 1.5)
points(c(0,10),c(mean(1/c(mean(null300.RT),mean(null150.RT))),mean(1/c(mean(decoy300.RT),mean(decoy150.RT)))),col=c('blue','red'),pch = '-',cex = 5)
arrows(x0 = c(0,10),x1 = c(0,10),y0 = c(min(1/c(mean(null300.RT),mean(null150.RT))),min(1/c(mean(decoy300.RT),mean(decoy150.RT)))), y1 = c(max(1/c(mean(null300.RT),mean(null150.RT))),max(1/c(mean(decoy300.RT),mean(decoy150.RT)))),code = 3,col = c('blue','red'),lwd = 1,angle = 90,length = 0.5)
temp.1=signif(mean(1/c(mean(null300.RT),mean(null150.RT))),2)
temp.2=signif((mean(1/c(mean(decoy300.RT),mean(decoy150.RT)))-mean(1/c(mean(null300.RT),mean(null150.RT))))/10e-6,2)
mtext(substitute(paste('k'[-1],' = ',temp.1,' (s'^-1,')'),list(temp.1 = temp.1)),side = 1,line = 6,at = 2.5,cex = 2)
mtext(substitute(paste('k'[theta],' = ',temp.2,' (M'^-1,'s'^-1,')'),list(temp.2 = temp.2)),side = 1,line = 6,at = 7.5,cex = 2)

