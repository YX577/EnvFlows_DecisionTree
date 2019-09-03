# Plot null and alternative hypothesis distributions

xx <- seq(-4,8,0.01)
d_0 <- dnorm(xx,0,1)
d_thresh <- dnorm(xx,3,1)
d_alt_high <- dnorm(xx,4.5,1)
d_alt_low <- dnorm(xx,1.5,1)

# Plot distributions
par(mfrow=c(3,1))
par(mar=c(3,4,3,2)+0.1)

plot(xx,d_thresh,typ="l",col="blue",
     xlab="",
     ylab="",
     xaxt='n',
     yaxt='n',
     bty='n')
lines(xx,d_alt_high,typ="l",col="red")
text(-3,0.32,"Decision tree: \nAlteration greater \nthan threshold",cex=1.1)
text(7,0.3,pos=4,"(a)")
#abline(v=4.5,col="gray80",lty=2)
axis(side=1, at=c(-4,3,4.5,8), labels=c("","Thr %","Alt %",""), pos=0, lty=1, col="black")

plot(xx,d_thresh,typ="l",col="blue",
     xlab="",
     ylab="",
     xaxt='n',
     yaxt='n',
     bty='n')
lines(xx,d_alt_low,typ="l",col="red")
text(-3,0.32,"Decision tree: \nAlteration less \nthan threshold",cex=1.1)
text(7,0.3,pos=4,"(b)")
#abline(v=1.5,col="gray80",lty=2)
axis(side=1, at=c(-4,1.5,3,8), labels=c("","Alt %","Thr %",""), pos=0, lty=1, col="black")

par(mar=c(5,4,3,2)+0.1)  # To enable x-axis label to appear
plot(xx,d_0,typ="l",col="blue",
     xlab="Extent of Alteration",
     ylab="",
     xaxt='n',
     yaxt='n',
     bty='n')
lines(xx,d_thresh,typ="l",col="red")
text(-3,0.32,"A priori \neffect size \npower analysis",cex=1.1)
text(7,0.3,pos=4,"(c)")
#abline(v=1.5,col="gray80",lty=2)
axis(side=1, at=c(-4,0,3,8), labels=c("","0%","Thr%",""), pos=0, lty=1, col="black")
legend(4.8,0.35,c("Null Dist.","Alt. Dist."),col=c("blue","red","gray80"),lty=c(1,1,2),bty="n")

par(mfrow=c(1,1))
par(mar=c(5,4,4,2)+0.1)

