# Mann-Whitney-Wilcoxon test function 

# This function improves upon the built-in wilcox.test function by also computing the power
# of the test using methods from Shieh et al. (2008).

MWW_test <- function(Y,X,alternative="g",alt_dist="snorm",paired=F, exact=F, correct=F){
  # Y = impacted/treated sample
  # X = baseline sample
  # alternative = greater than, less than, two-sided
  # alt_dist = assumed distribution of alternative
     # Options include:
        # snorm (Standard normal)
        # Add others from Shieh et al. (2008) later
  # paired = indicates whether you want a paired test, default = false
  # exact = indicates whether exact p-value should be computed, default = false, which
  #         permits use of the normal approximation
  # correct = indicates whether continuity correction should be applied to p-values, default = false
  #
  
  # Initialize a list for tracking output
  output <- list()
  # 1 U statistic
  # 2 Common language effect size (CLES)
  # 3 MWW test p-value
  # 4 theta 
  # 5 z_alpha (1 - p-value)
  # 6 beta
  # 7 mu_null
  # 8 sigma2_null
  # 9 mu_alt
  #10 sigma2_alt
  
  # COMPUTE TYPE I ERROR: Apply wilcox.test function 
  WT_output <- wilcox.test(Y,X,alternative="g",exact=F)
  U_stat <- as.numeric(WT_output[1])
  pval <- as.numeric(WT_output[3])
  if (pval < 0.0001) {pval = 0.0001}
  if (pval > 0.9999) {pval = 0.9999}
  output$U_stat <- U_stat
  output$CLES <- U_stat/(length(X)*length(Y)) # Common language effect size (% pairs where Y > X)
  output$pval <- pval
  
  # COMPUTE Type II error based on Type I error (or allow alpha input?)
  
  # Step 1: Create a set of standard normal variates for application of Eq 4b in Shieh et al (2008)
  Z_0_1 <- qnorm((seq(1,10000,1)-0.5)/10000,0,1)
  mu_null <- length(X)*length(Y)/2
  sigma2_null <- length(X)*length(Y)*(length(X)+length(Y)+1)/12
  
  # Step 2: Standardized data assuming two data sets have the same distribution
  Y_stdz <- (Y-mean(X))/sd(X)
  X_stdz <- (X-mean(X))/sd(X)
  
  # Step 3: Compute Type II error and power of statistical test
  theta <- mean(Y_stdz) - mean(X_stdz) # Assuming shift in central tendency only
  #if (theta > 3.719016) {theta = 3.719016}
  #if (theta < -3.719016) {theta = -3.719016}
  z_alpha <- qnorm(1-pval)
  pr1 <- pnorm(theta/sqrt(2))
  pr2 <- mean(pnorm(Z_0_1+theta)^2)
  pr3 <- pr2
  mu_alt <- length(X)*length(Y)*pr1
  sigma2_alt <- length(X)*length(Y)*pr1*(1-pr1)+length(X)*length(Y)*length(Y-1)*(pr2-pr1^2)+length(X)*length(Y)*length(X-1)*(pr3-pr1^2)
  if(sigma2_alt <= 0){sigma2_alt=0.001}
  stat_power <- pnorm((mu_alt-mu_null-z_alpha*sqrt(sigma2_null))/sqrt(sigma2_alt))
  beta <- 1 - stat_power
  
  # Record output
  output$theta <- theta
  output$z_alpha <- z_alpha
  output$beta <- beta
  output$mu_null <- mu_null
  output$sigma2_null <- sigma2_null
  output$mu_alt <- mu_alt
  output$sigma2_alt <- sigma2_alt
  
  return(output) 
}

# Test function 

Y <- rnorm(10,1200,100)
X <- rnorm(10,1000,100)
MWW_test(Y,X,alternative="g",exact="F")
  