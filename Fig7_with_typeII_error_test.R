# Type I vs. Type II error tradeoffs

# Effect size
ES <- seq(-3,3,0.1)

# Type I and II error cost percentages
loss_costfrac <- seq(0,1,0.05)

# Number of test runs
MM <- 1000
         
# Generate null distribution sample
NN <- 37
ny_pre <- NN
ny_post <- NN

# Create output matrices
aa <- vector(mode = "numeric",length=length(ES))
bb <- vector(mode = "numeric",length=length(ES))

pNA_CA <- vector(mode = "numeric",length=length(ES))
pA_CNA <- vector(mode = "numeric",length=length(ES))

mu_alt <- vector(mode = "numeric",length=length(ES))
sigma2_alt <- vector(mode = "numeric",length=length(ES))

ER_aa <- array(NA,dim=c(length(ES),length(loss_costfrac),MM))
ER_bb <- array(NA,dim=c(length(ES),length(loss_costfrac),MM))

#Create a set of standard normal variates for application of type II error distin Shieh et al (2006)
Z_0_1 <- qnorm((seq(1,10000,1)-0.5)/10000,0,1)

# Compute mu and sigma of the null distribution of W 
mu_null <- ny_pre*ny_post/2
sigma2_null <- ny_pre*ny_post*(ny_pre+ny_post+1)/12

for(m in 1:MM){

  for(i in 1:length(ES)){
    
    # Generate null distribution
    rand_null <- rnorm(NN,0,1)
    
    # Generate alternative distribution with a specified effect size
    rand_alt <- rnorm(NN,0,1) + ES[i]
    
    # Compute type I error using t distribution
    if(ES[i] >= 0){
      aa[i] = wilcox.test(rand_alt,rand_null,alternative="greater",paired=F,exact=F,correct=T,conf.int=F)$p.value
    } else {
      aa[i] = wilcox.test(-rand_alt,rand_null,alternative="greater",paired=F,exact=F,correct=T,conf.int=F)$p.value
    }
    
    # Compute type II error using standard normal assumption from Shieh et al. (2006)
    if(ES[i] >= 0){
      theta <- ES[i]
    } else {
      theta <- -ES[i]
    }
    
    pr1 <- pnorm(theta/sqrt(2))
    pr2 <- mean(pnorm(Z_0_1+theta)^2)
    pr3 <- pr2
    
    mu_alt[i] = ny_pre*ny_post*pr1
    sigma2_alt[i] = ny_pre*ny_post*pr1*(1-pr1) + ny_pre*ny_post*(ny_post-1)*(pr2-pr1^2) + ny_pre*ny_post*(ny_pre-1)*(pr3-pr1^2)
    z_alpha <- qnorm(1-aa[i])
    stat_power_up <- pnorm((mu_alt[i] - mu_null - z_alpha * sqrt(sigma2_null))/sqrt(sigma2_alt[i]))
    bb[i] = 1 - stat_power_up
  
    # Compute regret probabilities
    pNA_CA[i] = aa[i]/(aa[i]+1-bb[i])
    pA_CNA[i] = bb[i]/(bb[i]+1-aa[i])
    
    # Compute expected regrets for different error consequences
    if(ES[i]>=0){
      for(j in 1:length(loss_costfrac)){
        ER_aa[i,j,m] = pNA_CA[i] * loss_costfrac[j]
        ER_bb[i,j,m] = pA_CNA[i] * (1 - loss_costfrac[j])
      }
    } else {
      for(j in 1:length(loss_costfrac)){
        ER_aa[i,j,m] = pA_CNA[i] * loss_costfrac[j] 
        ER_bb[i,j,m] = pNA_CA[i] * (1 - loss_costfrac[j])
      }
    }
      
  }

}

# Average ER over mm
ER_AA <- apply(ER_aa,c(1,2),mean)
ER_BB <- apply(ER_bb,c(1,2),mean)

# Plot power curve
plot(aa,bb)

# Plot decision recommendation
ER_diff <- ER_AA - ER_BB
#ER_diff[ER_diff < 0] 
#ER_diff[ER_diff >= 0]

'

# Plot error cost-weighted analysis with type I and II errors
filled.contour(ES,loss_costfrac,ER_diff,
               xlab=expression(paste("Effect size (",sigma,")",sep="")),
               ylab="Hydropower Loss Cost Fraction",cex.axis=0.8,
               levels = c(-1,0,1),col=c("black","white"),
               main="Decision Tree - 37 yrs")

'

# Compute critical value of t stat
t_crit <- qt(0.95,NN)

# Create a matrix of values to perturb
ES_mat_alpha <- matrix(rep(ES,length(loss_costfrac)),nrow=length(loss_costfrac),ncol=length(ES),byrow=T)
ES_mat_alpha[ES_mat_alpha  <= t_crit] <- -0.5
ES_mat_alpha[ES_mat_alpha  > t_crit] <- 0.5

# Determine values for detemrinistic method
ES_mat_det <- matrix(rep(ES,length(loss_costfrac)),nrow=length(loss_costfrac),ncol=length(ES),byrow=T)
ES_mat_det[ES_mat_det > 0] <- 0.5
ES_mat_det[ES_mat_det <= 0] <- -0.5


# ----
'

# Make a panel plot showing the differences
par(mfrow=c(3,1))
#par(mar=c(3,4,2,2)+0.1)


# Plot error cost-weighted analysis with type I and II errors
filled.contour(ES,loss_costfrac,ER_diff,
               xlab="Effect size",
               ylab="Hydropower Loss Cost Fraction",
               levels = c(-1,0,1),col=c("white","black"),
               main="Decision Tree")

# Plot analysis using NHST
filled.contour(ES,loss_costfrac,t(ES_mat_alpha),
               xlab=expression(paste("Effect size (",sigma,")",sep="")),
               ylab="Hydropower Loss Cost Fraction",cex.axis=0.8,
               levels = c(-1,0,1),col=c("white","black"),
               main=expression(paste("NHST (", alpha,") = 0.05 - 37 Yrs")))
#text(-0.5,0.5,"t = 1.812461",cex=0.7)

# Plot analysis using deterministic threshold

filled.contour(ES,loss_costfrac,t(ES_mat_det),
               xlab=expression(paste("Effect size (",sigma,")",sep="")),
                ylab="Hydropower Loss Cost Fraction",cex.axis=0.8,
                levels = c(-1,0,1),col=c("white","black"),
                main="Deterministic - 37 Yrs")


#par(mar=c(5,4,4,2)+0.1)
par(mfrow=c(1,1))

'

# Plot null and alternative hypothesis distributions
for(k in 1:13){
  x_null <- seq(0,NN^2,0.1)
  d_null <- dnorm(x_null,mu_null,sqrt(sigma2_null))
  plot(x_null,d_null,typ="l",col="gray50",
       xlab="U test statistic",
       ylab="Density",
       ylim=c(0,0.1))
  abline(v = qnorm(1-aa[(k-1)*5+1],mu_null,sqrt(sigma2_null)),col="red")
  x_alt <- seq(0,NN^2,0.1)
  d_alt <- dnorm(x_alt,mu_alt[(k-1)*5+1],sqrt(sigma2_alt[(k-1)*5+1]))
  lines(x_alt,d_alt)
}
'

