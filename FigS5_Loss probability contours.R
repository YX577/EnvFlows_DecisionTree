aa <- seq(0.001,0.999,0.001)
bb <- seq(0.001,0.999,0.001)

P_NA_CNA <- array(NA,dim=c(length(aa),length(bb)))
P_NA_CA <- array(NA,dim=c(length(aa),length(bb)))
P_A_CNA <- array(NA,dim=c(length(aa),length(bb)))
P_A_CA <- array(NA,dim=c(length(aa),length(bb)))

# For a non-informative prior
k = 0.5

for (i in 1:length(aa)){
  for (j in 1:length(bb)){
    P_NA_CNA[i,j] = (1-k)*(1 - aa[i])/((1-k)*(1-aa[i])+k*bb[j])
    P_NA_CA[i,j] = (1-k)*aa[i]/((1-k)*aa[i]+k*(1-bb[j]))
    P_A_CNA[i,j] = k*bb[j]/((1-k)*(1-aa[i])+k*bb[j])
    P_A_CA[i,j] = k*(1-bb[j])/((1-k)*aa[i]+k*(1-bb[j]))
  }
}



# 2 x 2 figure of contour plots
par(mfrow=c(2,2))

contour(aa,bb,P_NA_CNA,xlab=expression(paste("Type I Error Prob.",(alpha))),ylab=expression(paste("Type II Error Prob.",(beta)))
        ,col="gray67",main="Maintain operations \n P(NA|CNA)")

contour(aa,bb,P_A_CNA,xlab=expression(paste("Type I Error Prob.",(alpha))),ylab=expression(paste("Type II Error Prob.",(beta)))
        ,main="Ecological regret \n P(A|CNA)")

contour(aa,bb,P_NA_CA,xlab=expression(paste("Type I Error Prob.",(alpha))),ylab=expression(paste("Type II Error Prob.",(beta)))
        ,main="Hydropower regret \n P(NA|CA)")

contour(aa,bb,P_A_CA,xlab=expression(paste("Type I Error Prob.",(alpha))),ylab=expression(paste("Type II Error Prob.",(beta)))
        ,main="Change operations \n P(A|CA)",col="gray67")

par(mfrow=c(1,1))


# Create masks to prevent contours > 0.5
P_NA_CA_masked <- P_NA_CA
P_NA_CA_masked[P_NA_CA_masked > 0.5] <- 0.501

P_A_CNA_masked <- P_A_CNA
P_A_CNA_masked[P_A_CNA_masked > 0.5] <- 0.501


# 1 x 3 figure of contour plots for publication
par(mfrow=c(1,2))

contour(aa,bb,P_NA_CA_masked,nlevels=5,xlab=expression(paste("Type I Error Prob."," ",(alpha),sep="")),ylab=expression(paste("Type II Error Prob."," ",(beta),sep=""))
        ,main="Hydropower loss prob. \n P(NA|CA)")

contour(aa,bb,P_A_CNA_masked,nlevels=5,xlab=expression(paste("Type I Error Prob."," ",(alpha),sep="")),ylab=expression(paste("Type II Error Prob."," ",(beta),sep=""))
        ,main="Ecological loss prob. \n P(A|CNA)")

#contour(aa,bb,P_NA_CA_masked + P_A_CNA_masked,nlevels=5,xlab=expression(paste("Type I Error Prob.","",(alpha),sep="")),ylab=expression(paste("Type II Error Prob.","",(beta),sep=""))
#        ,main="Sum of loss prob. \n P(NA|CA) + P(A|CNA)")

par(mfrow=c(1,1))
