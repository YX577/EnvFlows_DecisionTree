# Figure showing MWW test without a threshold

par(mfrow=(c(2,1)))

x <- seq(0,400,1)
h_pre <- dnorm(x,100,20)
h_post <- dnorm(x,80,16)

plot(x,h_pre,typ="l",col="black", lty=1,xlim=c(0,400),ylim=c(0,0.03),
     xlab=expression(paste("Discharge (",m^3,"/s)",sep="")),
     ylab="Probability",las=1,cex.axis=0.8,cex.lab=0.8) # MAKE THESE AXIS LABELS GO AWAY
title("Hypothetical AFDC quantile distributions",cex.main=1.00)
text(5,0.025,"(a)",cex=0.7)
lines(x,h_pre)
lines(x,h_post,col="gray",lty=1)
legend("topright",c("Pre-dam","Post-dam"),lty=c(1,1),col=c("black","gray"),bty="n",cex=0.7)

# Sample from these distributions
m_size <- 10
n_size <- 10

U_max <- m_size*n_size
U_rng <- seq(0,U_max,1)  

norm_pre <- rnorm(m_size,100,20) 
norm_post <- rnorm(n_size,80,16)

# Run MWW test
U_MWWtest <- MWW_test(log(norm_post),log(norm_pre),alternative="l",alt_dist="snorm",paired=F, exact=F, correct=F)
U_mu_null <- as.numeric(U_MWWtest[7])  
U_sigma2_null <- as.numeric(U_MWWtest[8]) 
U_mu_alt <- as.numeric(U_MWWtest[9]) 
U_sigma2_alt <- as.numeric(U_MWWtest[10]) 

# Generate and plot test statistic probability distributions
U_dist_null <- dnorm(U_rng,U_mu_null,sqrt(U_sigma2_null))  
U_dist_alt <- dnorm(U_rng,U_mu_alt,sqrt(U_sigma2_alt))

plot(U_rng,U_dist_null,type="l",lty=2,xlim=c(0,max(U_rng)),ylim=c(0,1.2*max(U_dist_alt)),
     xlab="U test statistic",ylab="Probability",las=1,cex.axis=0.8,cex.lab=0.8)
title("Distribution of U test statistic",cex.main=1.00)
text(1,1.1*max(U_dist_alt),"(b)",cex=0.7)
lines(U_rng,U_dist_alt,lty=3)
legend("topright",c("Null dist","Alt dist"),lty=c(2,3),cex=0.7,bty="n")

par(mfrow=(c(1,1)))

# Now do this with a threshold ----------------------------------------------------------------------------------

par(mfrow=(c(2,1)))

# Plot scaled pre-dam flow distribution with pre- and post-dam flows
h_pre_scaled <- dnorm(x,70,14)
plot(x,h_pre,typ="l",col="black", lty=1,xlim=c(0,400),ylim=c(0,0.03),
     xlab=expression(paste("Discharge (",m^3,"/s)",sep="")),
     ylab="Probability",
     cex.axis=0.8,cex.lab=0.8,las=1) # MAKE THESE AXIS LABELS GO AWAY
title("Hypothetical AFDC quantile distributions",cex.main=1.00)
text(1.5,0.025,"(a)",cex=0.7)
lines(x,h_pre_scaled,lty=2)
lines(x,h_post,col="gray",lty=1)
legend("topright",c("Pre-dam","Post-dam","Pre-dam scaled"),lty=c(1,1,2),cex=0.7,col=c("black","gray","gray40"),bty="n")

# Scale pre-dam flows
norm_pre_scaled <- 0.8*norm_pre

# Run MWW test
#U_MWWtest <- MWW_test(log(norm_post),log(norm_pre_scaled),alternative="l",alt_dist="snorm",paired=F, exact=F, correct=F)
U_MWWtest <- MWW_test(log(norm_post),log(norm_pre_scaled),alternative="l",alt_dist="snorm",paired=F, exact=F, correct=F)
U_mu_null <- as.numeric(U_MWWtest[7])  
U_sigma2_null <- as.numeric(U_MWWtest[8]) 
U_mu_alt <- as.numeric(U_MWWtest[9]) 
U_sigma2_alt <- as.numeric(U_MWWtest[10]) 

# Generate and plot test statistic probability distributions
U_dist_null <- dnorm(U_rng,U_mu_null,sqrt(U_sigma2_null))  
U_dist_alt <- dnorm(U_rng,U_mu_alt,sqrt(U_sigma2_alt))

plot(U_rng,U_dist_null,type="l",lty=2,xlim=c(0,max(U_rng)),ylim=c(0,1.2*max(U_dist_alt)),
     xlab="U test statistic",ylab="Probability",cex.axis=0.8,cex.lab=0.8,las=1)
title("Distributions of U test statistic",cex.main=1.00)
lines(U_rng,U_dist_alt,lty=3)
text(1,1.1*max(U_dist_alt),"(b)",cex=0.7)
legend("topright",c("Null dist","Alt dist"),lty=c(2,3),cex=0.7,bty="n")

par(mfrow=(c(1,1)))

