# MANN-WHITNEY-WILCOXON TEST WITH HP_RES OUTPUT

# Source other R scripts as needed

# Source FDC_Maker function
source("C:/Users/joryh/ownCloud/documents/Tufts/Hydro_Alteration/HA_HypTest/Codes/Misc_Functions.R")

# Reference pre-dam and post-dam daily flow time series
#load("Q_out_results.Rdata")

# Define annual exceedance probabilities of interest p_key
p <- seq(1,99)/100
p_key = c(5,95)

# Define vector of percent deviation thresholds
devs <- seq(0,0.95,by = 0.05)

#Create a set of Type I error values for which to compute the Type II errors
alphas <- c(0.00000001,0.0001,seq(0.001,0.009,0.001),seq(0.01,0.99,0.01),seq(0.991,0.999,0.001),0.9999,0.99999999)

# Subset param_set to choose parameter combinations to run from hp_res_v*
params <- read.csv("params.csv")
param_set <- params[params[,1]==1,]

# Initialize vectors for storing results of overall significance
t1_ovrl_error_up <- array(NA,dim=c(length=length(devs),length=length(alphas),length=nrow(param_set)))
t2_ovrl_error_up <- array(NA,dim=c(length=length(devs),length=length(alphas),length=nrow(param_set)))
t1_ovrl_error_down <- array(NA,dim=c(length=length(devs),length=length(alphas),length=nrow(param_set)))
t2_ovrl_error_down <- array(NA,dim=c(length=length(devs),length=length(alphas),length=nrow(param_set)))
t1_ovrl_error_bload_hp <- array(NA,dim=c(length=length(devs),length=length(alphas),length=nrow(param_set)))
t2_ovrl_error_bload_hp <- array(NA,dim=c(length=length(devs),length=length(alphas),length=nrow(param_set)))

# Initialize vectors for storing results for one-tailed tests in each direction
# Note that only error_up cases matter for low-flows (Q95) and error_down cases matter for high
# flows (Q5)

t1_error_up = array(NA,dim=c(length(p_key),length(devs),nrow(param_set)))
t1_error_down = array(NA,dim=c(length(p_key),length(devs),nrow(param_set)))
t2_error_up <- array(NA,dim=c(length(p_key),length(devs),length(alphas),nrow(param_set))) 
t2_error_down <- array(NA,dim=c(length(p_key),length(devs),length(alphas),nrow(param_set)))

# For simple deterministic median analysis 
# (Consider changing to U stat to make it consistent with MWW test?)
diff_pre_post_up_res <- array(NA,dim=c(length(p_key),length(devs),nrow(param_set)))
diff_pre_post_down_res <- array(NA,dim=c(length(p_key),length(devs),nrow(param_set)))
diff_pre_post_up_all_p_key <- array(NA,dim=c(length(devs),nrow(param_set)))
diff_pre_post_down_all_p_key <- array(NA,dim=c(length(devs),nrow(param_set)))
diff_pre_post_bload_hp <- array(NA,dim=c(length(devs),nrow(param_set)))

# Compute AFDCs for pre- and post-dam periods using outflows from different alternative reservoir operation policies
for (b in 1:nrow(param_set)){

  AFDC.pre <- FDC.annual(Q_in,as.POSIXlt(dates),p,WY=T) 
  AFDC.post <- FDC.annual(Q_out_res[,b],as.POSIXlt(dates),p,WY=T)
  if(min(AFDC.pre) < 0){stop("Negative pre-dam flow value")}
  if(min(Q_out_res[,b]) < 0){stop("Negative post-dam flow value")}
  
  # Compute pre-and post-dam flow statistics ------------------------------------------------------------
  
  # Mean Q
  AFDC.pre_means <- apply(AFDC.pre,2,mean)
  AFDC.post_means <- apply(AFDC.post,2,mean)
    
  # Cv 
  AFDC.pre_cv <- apply(AFDC.pre,2,sd)/AFDC.pre_means
  AFDC.post_cv <- apply(AFDC.post,2,sd)/AFDC.post_means

  '
  # Create a plot to compare pre- and post-dam FDCs
  par(mfrow=c(1,2))
  par(mar=c(5.1,5.1,4.1,2.1))

  # Plot pre-dam AFDCs 
  plot(p*100,AFDC.pre[1,]/35.31,type="l",log="y",col="gray69",
       #xlim=c(0,100),ylim=(c(1000,250000)),
       xlim = c(0,100), ylim = c(10,5000),
       xlab="Exceedance Probability (%)",ylab = expression(Discharge ~ (~m^{3}/s)),cex.lab=0.9,
       main="Pre-dam AFDCs", cex.main = 1.05)
  
  for (i in 2:37){
  lines(p*100,AFDC.pre[i,]/35.31,col="gray69")
  } 
  
  abline(v=5,lty=2,col="black")
  abline(v=95,lty=2,col="black")
  text(15,5000,"Q5",cex=0.75)
  text(85,5000,"Q95",cex=0.75)
  
  # Plot post-dam AFDCs 
  
  # Start with plot of first year in each sample
  plot(p*100,AFDC.post[1,]/35.31,type="l",log="y",col="gray69",
       
       #xlim=c(0,100),ylim=(c(1000,250000)),
       xlim = c(0,100), ylim = c(10,5000),
       xlab="Exceedance Probability (%)",ylab = expression(Discharge ~ (~m^{3}/s)),cex.lab=0.9,
       main="Post-dam AFDCs", cex.main = 1.05)
  
  
  # Then plot remaining years
  for (i in 2:37){
    lines(p*100,AFDC.post[i,]/35.31,col="gray69")
  }
  
  abline(v=5,lty=2,col="black")
  abline(v=95,lty=2,col="black")
  text(15,5000,"Q5",cex=0.75)
  text(85,5000,"Q95",cex=0.75)

  par(mfrow=c(1,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  '
  
  # Compute AFDCs for exceedance probabilities of concern------------------------------------------------
  
  # Create arrays
  AFDC.pre_up = array(NA,dim=c(nrow(AFDC.pre),length(p_key),length(devs)))
  AFDC.pre_down = array(NA,dim=c(nrow(AFDC.pre),length(p_key),length(devs)))

  # Analyze data
  for (d in 1:length(devs)){
    # To test hypotheses that pre-dam flows are less than post-dam flows
    AFDC.pre_up[,,d]=(1+devs[d])*AFDC.pre[,p_key] #AFDC.pre_up[Year,Exc Prob,d]  
    # To test hypotheses that pre-dam flows are greater than post-dam flows
    AFDC.pre_down[,,d]=(1-devs[d])*AFDC.pre[,p_key] #AFDC.pre_down[Year,Exc Prob,d]
  }
  
  # Compute properties of these shifted AFDCs
  AFDC.pre_up_means <- apply(AFDC.pre_up,3,mean)
  AFDC.pre_up_sds <- apply(AFDC.pre_up,3,sd)
  AFDC.pre_down_means <- apply(AFDC.pre_down,3,mean)
  AFDC.pre_down_sds <- apply(AFDC.pre_down,3,sd)
  
  # Recreate AFDC.post so that it only has key exc prob
  AFDC.post = AFDC.post[,p_key]
  
  # Compute logs of AFDCs
  log_AFDC.pre_up <- log(AFDC.pre_up)
  log_AFDC.pre_down <- log(AFDC.pre_down)
  log_AFDC.post <- log(AFDC.post)

  # Apply Mann-Whitney-Wilcoxon test to each AFDC quantile for each percent deviation to 
  # compute Type I errors 
  W_results_up = list()
  W_results_down = list()
  Uval_up <- array(NA,dim=c(length(p_key),length(devs)))  # U stat values for increases in flow
  Uval_down <- array(NA,dim=c(length(p_key),length(devs))) # U stat values for decreases in flow
  pval_up <- array(NA,dim=c(length(p_key),length(devs))) #p values for increases in flow
  pval_down <- array(NA,dim=c(length(p_key),length(devs))) #p values for decreases in flow
 
  for (d in 1:length(devs)){ 
    for (f in 1:length(p_key)){
      # Test to see if post-dam flows are greater than pre-dam flows
      W_results_up = wilcox.test(AFDC.post[,f],AFDC.pre_up[,f,d],alternative="greater",paired=F,exact=F,correct=T,conf.int=F)
      Uval_up[f,d] = as.numeric(W_results_up[1])
      pval_up[f,d] = as.numeric(W_results_up[3])
      # Test to see if pre-dam flows are greater than post-dam flows
      W_results_down = wilcox.test(AFDC.pre_down[,f,d],AFDC.post[,f],alternative="greater",paired=F,exact=F,correct=T,conf.int=F)
      Uval_down[f,d] = as.numeric(W_results_down[1])
      pval_down[f,d] = as.numeric(W_results_down[3])
    }
  }
  
  '
  # Plot histograms
  par(mfrow=c(2,1))
  
  hist(AFDC.post[,2],
       xlim=c(0,10000),
       ylim=c(0,30),
       breaks=seq(0,10000,500),
       xlab=expression(paste("Annual Q95 Discharge (",m^3,"/s)",sep="")),
       ylab="Years",
       main="Distribution of Post-Dam Annual Q95")
  hist(AFDC.pre_up[,2,1],
       xlim=c(0,10000), 
       ylim=c(0,30),
       breaks=seq(0,10000,500),
       xlab=expression(paste("Annual Q95 Discharge (",m^3,"/s)",sep="")),
       ylab="Years",
       main="Distribution of Pre-Dam Annual Q95")
  
  par(mfrow=c(1,1))
  '
  
  # Store all type I errors for further analysis
  t1_error_up[,,b] = pval_up
  t1_error_down[,,b] = pval_down 
  
  # Compute distribution of alternative hypothesis ------------------------------
  
  # Estimate Type II errors assuming normal FDC exceedance prob distributions following Shieh et al. (2008)

  # Compute pre- and post-dam record length in years
  ny_pre <- nrow(AFDC.pre)
  ny_post <- nrow(AFDC.post)

  # Compute mu and sigma of the null distribution of W 
  mu_null <- ny_pre*ny_post/2
  sigma2_null <- ny_pre*ny_post*(ny_pre+ny_post+1)/12

  #Create a set of standard normal variates for application of Eq 4b in Shieh et al (2008)
  Z_0_1 <- qnorm((seq(1,10000,1)-0.5)/10000,0,1)

  # Create new arrays to store calculations
  AFDC.pre_up_std <- array(NA,dim=c(nrow(AFDC.pre),length(p_key),length(devs)))
  AFDC.post_up_std <-  array(NA,dim=c(nrow(AFDC.post),length(p_key),length(devs)))
  siegel.tukey_pval_up <-  array(NA,dim=c(length(p_key),length(devs)))  
  ab_pval_up <-  array(NA,dim=c(length(p_key),length(devs))) #Ansari-Bradley rank-based test for scale changes  
  theta_up <-  array(NA,dim=c(length(p_key),length(devs)))
  z_alpha_up <-  array(NA,dim=length(alphas))
  pr1_up <-  array(NA,dim=c(length(p_key),length(devs)))
  pr2_up <-  array(NA,dim=c(length(p_key),length(devs)))
  pr3_up <-  array(NA,dim=c(length(p_key),length(devs)))
  mu_alt_up <-  array(NA,dim=c(length(p_key),length(devs)))
  sigma2_alt_up <-  array(NA,dim=c(length(p_key),length(devs)))
  stat_power_up <-  array(NA,dim=c(length(p_key),length(devs),length(alphas)))

  AFDC.pre_down_std <-  array(NA,dim=c(nrow(AFDC.pre),length(p_key),length(devs)))
  AFDC.post_down_std <-  array(NA,dim=c(nrow(AFDC.post),length(p_key),length(devs)))
  siegel.tukey_pval_down <-  array(NA,dim=c(length(p_key),length(devs)))
  ab_pval_down <-  array(NA,dim=c(length(p_key),length(devs))) #Ansari-Bradley rank-based test for scale changes   
  theta_down <-  array(NA,dim=c(length(p_key),length(devs)))
  z_alpha_down <-  array(NA,dim=length(alphas))
  pr1_down <-  array(NA,dim=c(length(p_key),length(devs)))
  pr2_down <-  array(NA,dim=c(length(p_key),length(devs)))
  pr3_down <-  array(NA,dim=c(length(p_key),length(devs)))
  mu_alt_down <-  array(NA,dim=c(length(p_key),length(devs)))
  sigma2_alt_down <-  array(NA,dim=c(length(p_key),length(devs)))
  stat_power_down <-  array(NA,dim=c(length(p_key),length(devs),length(alphas)))

  # Compute power using equation 3 in Shieh et al. (2007)
  
  # Test whether post-dam flows are greater than pre-dam flows 
  # Standardize post-dam and pre-dam flows by pre-dam mean and sd
  for (d in 1: length(devs)){
    for (f in 1:length(p_key)){
      AFDC.pre_up_std[,f,d] = (AFDC.pre_up[,f,d]-mean(AFDC.pre_up[,f,d]))/sd(AFDC.pre_up[,f,d])
      AFDC.post_up_std[,f,d] = (AFDC.post[,f]-mean(AFDC.pre_up[,f,d]))/sd(AFDC.pre_up[,f,d])
      theta_up[f,d] = abs(mean(AFDC.post_up_std[,f,d])-mean(AFDC.pre_up_std[,f,d]))
      pr1_up[f,d] = pnorm(theta_up[f,d]/sqrt(2))
      pr2_up[f,d] = mean(pnorm(Z_0_1+theta_up[f,d])^2)
      pr3_up[f,d] = pr2_up[f,d]
      mu_alt_up[f,d] = ny_pre*ny_post*pr1_up[f,d]
      sigma2_alt_up[f,d] = ny_pre*ny_post*pr1_up[f,d]*(1-pr1_up[f,d])+ny_pre*ny_post*(ny_post-1)*(pr2_up[f,d]-pr1_up[f,d]^2)
        +ny_pre*ny_post*(ny_pre-1)*(pr3_up[f,d]-pr1_up[f,d]^2)
    
      for (a in 1:length(alphas)){
        z_alpha_up[a] = qnorm(1-alphas[a])
        stat_power_up[f,d,a] = pnorm((mu_alt_up[f,d]-mu_null-z_alpha_up[a]*sqrt(sigma2_null))/sqrt(sigma2_alt_up[f,d]))
     
      }
    }
  }

  #Store results for further analysis
  t2_error_up[,,,b] = 1-stat_power_up
  
  # Test whether post-dam flows are less than pre-dam flows
  # In other words, whether pre-dam flows are greater than post-dam flows
  # Standardize post-dam and pre-dam flows by pre-dam mean and sd
  for (d in 1:length(devs)){
    for (f in 1:length(p_key)){
      #Standardized data by pre-dam flows
      AFDC.pre_down_std[,f,d] = (AFDC.pre_down[,f,d]-mean(AFDC.pre_down[,f,d]))/sd(AFDC.pre_down[,f,d])
      AFDC.post_down_std[,f,d] = (AFDC.post[,f]-mean(AFDC.pre_down[,f,d]))/sd(AFDC.pre_down[,f,d])
      theta_down[f,d] = abs(mean(AFDC.pre_down_std[,f,d])-mean(AFDC.post_down_std[,f,d])) #Tests pre > post
      pr1_down[f,d] = pnorm(theta_down[f,d]/sqrt(2))
      pr2_down[f,d] = mean(pnorm(Z_0_1+theta_down[f,d])^2)
      pr3_down[f,d] = pr2_down[f,d]
      mu_alt_down[f,d] = ny_pre*ny_post*pr1_down[f,d]
      sigma2_alt_down[f,d] = ny_pre*ny_post*pr1_down[f,d]*(1-pr1_down[f,d])+ny_pre*ny_post*(ny_pre-1)*(pr2_down[f,d]-pr1_down[f,d]^2)+ny_pre*ny_post*(ny_post-1)*(pr3_down[f,d]-pr1_down[f,d]^2)

      for (a in 1:length(alphas)){
        z_alpha_down[a] = qnorm(1-alphas[a]) #For greater than hypothesis test (pre > post)
        stat_power_down[f,d,a] = pnorm((mu_alt_down[f,d]-mu_null-z_alpha_down[a]*sqrt(sigma2_null))/sqrt(sigma2_alt_down[f,d]))
      }
    }
  }
  
  #Store results for further analysis
  t2_error_down[,,,b] = 1-stat_power_down
  
} 


##################################### END OF LOOP ##########################################


