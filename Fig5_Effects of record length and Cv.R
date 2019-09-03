# FIGURE 5: 

# Effects of record length and Cv on hypothesis testing errors and loss probabilities
# for a hypothetical threshold exceedance 


# Set parameters for record length and Cv
MM <- 50
N_length <- c(10,37,100)
Cv <- c(0.25, 0.5, 1)

# Create a vector with a range of alphas to be tested
a <- c(0.0000001,seq(0.001,0.009,0.001),seq(0.01,0.99,0.01),seq(0.991,0.999,0.001),0.9999999)
b <- array(NA,dim=c(length(a),length(N_length),length(Cv),MM))

for (m in 1:MM){
  
  for (l in 1:length(N_length)){
    
    for (v in 1:length(Cv)){
      
      # Create streamflow data
      X <- rnorm(N_length[l],1000,Cv[v]*1000)
      Y <- X * 1.1
      
      # Change negative values to zero
      X[X < 0] = 0
      Y[Y < 0] = 0
      
      for (i in 1:length(a)){
        
        # Test sensitivity of Type II error to Type I error tolerance
        
        # Step 1: Create a set of standard normal variates for application of Eq 4b in Shieh et al (2008)
        Z_0_1 <- rnorm(10000,0,1)
        mu_null <- length(X)*length(Y)/2
        sigma2_null <- length(X)*length(Y)*(length(X)+length(Y)+1)/12
        
        # Step 2: Standardized data assuming two data sets have the same distribution
        # Use mean and standard deviation from X 
        Y_stdz <- (Y-mean(X))/sd(X)
        X_stdz <- (X-mean(X))/sd(X)
        
        # Step 3: Compute Type II error and power of statistical test
        theta <- mean(Y_stdz)-mean(X_stdz) # Assuming shift in central tendency only
        if (theta > 3.719016) {theta = 3.719016}
        if (theta < -3.719016) {theta = -3.719016}
        if (a[i] < 0.0001) {a[i] = 0.0001}
        if (a[i] > 0.9999) {a[i] = 0.9999}
        z_alpha <- qnorm(1-a[i])
        if (z_alpha > 3.719016) {z_alpha = 3.719016}
        if (z_alpha < -3.719016) {z_alpha = -3.719016}
        pr1 <- pnorm(theta/sqrt(2))
        pr2 <- mean(pnorm(Z_0_1+theta)^2)
        pr3 <- pr2
        mu_alt <- length(X)*length(Y)*pr1
        sigma2_alt <- length(X)*length(Y)*pr1*(1-pr1)+length(X)*length(Y)*length(Y-1)*(pr2-pr1^2)+length(X)*length(Y)*length(X-1)*(pr3-pr1^2)
        if(sigma2_alt <= 0){sigma2_alt=0.0001}
        stat_power <- pnorm((mu_alt-mu_null-z_alpha*sqrt(sigma2_null))/sqrt(sigma2_alt))
        beta <- 1 - stat_power
        b[i,l,v,m] = beta
      }
      
    }
    
  }
  
}

# Average simulation results
b_avg <- apply(b,c(1,2,3),mean)

# Plot effects of record length for different Cvs
ltype <- c(2,1,4)
lwidth <- c(1,2,1)

par(mfrow=c(2,2))

# For Cv = 0.25
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),
     xlab=expression(paste("Type I error probability (",alpha,")")),
     ylab=expression(paste("Type II error probability ( ",beta,")")),
     cex.axis=0.85)
for (v in 1:length(N_length)){
  lines(a,b_avg[,v,1],lty=ltype[v],lwd=lwidth[v])
}
legend("topright",legend=c("10","37","100"),bty="n",lty=c(2,1,4),lwd=c(1,2,1),cex = 0.8,
       title="Record length (Yrs)")
text(0.9,0.3,"(a)",cex=0.75)
title("Effect of record length \n for Cv = 0.25",cex.main=0.9)

# For Cv = 1.00
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),   
     xlab=expression(paste("Type I error probability (",alpha,")")),
     ylab=expression(paste("Type II error probability ( ",beta,")")),
     cex.axis=0.85)
for (v in 1:length(N_length)){
  lines(a,b_avg[,v,3],lty=ltype[v],lwd=lwidth[v])
}
legend("topright",legend=c("10","37","100"),bty="n",lty=c(2,1,4),lwd=c(1,2,1),cex = 0.8,
       title="Record length (Yrs)")
text(0.9,0.3,"(b)",cex=0.75)
title("Effect of record length \n for Cv = 1.00",cex.main=0.9)

# Loss probs for Cv = 0.25
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),   
     xlab="Hydropower loss probability",
     ylab="Ecological loss probability",
     cex.axis=0.85)
for (v in 1:length(N_length)){
  lines(a/(a+1-b_avg[,v,1]),b_avg[,v,1]/(b_avg[,v,1]+1-a),lty=ltype[v],lwd=lwidth[v])
}
lines(c(0,0.5,0.5),c(0.5,0.5,0),col="gray80",lwd=0.5,lty=6)
legend("topright",legend=c("10","37","100"),bty="n",lty=c(2,1,4),lwd=c(1,2,1),cex = 0.8,
       title="Record length (Yrs)")
legend("bottomright",legend="Upper bound",col="gray80",lty=6,bty="n",cex=0.8)
text(0.9,0.3,"(c)",cex=0.75)
title("Effect of record length \n for Cv = 0.25",cex.main=0.9)

# Loss probs for Cv = 1.00
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),   
     xlab="Hydropower loss probability",
     ylab="Ecological loss probability",
     cex.axis=0.85)
for (v in 1:length(N_length)){
  lines(a/(a+1-b_avg[,v,3]),b_avg[,v,3]/(b_avg[,v,3]+1-a),lty=ltype[v],lwd=lwidth[v])
}
lines(c(0,0.5,0.5),c(0.5,0.5,0),col="gray80",lwd=0.5,lty=6)
legend("topright",legend=c("10","37","100"),bty="n",lty=c(2,1,4),lwd=c(1,2,1),cex = 0.8,
       title="Record length (Yrs)")
legend("bottomright",legend="Upper bound",col="gray80",lty=6,bty="n",cex=0.8)
text(0.9,0.3,"(d)",cex=0.75)
title("Effect of record length \n for Cv = 1.00",cex.main=0.9)

par(mfrow=c(1,1))
