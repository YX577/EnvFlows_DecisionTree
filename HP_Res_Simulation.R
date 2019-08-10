# HYDROPOWER RESERVOIR SIMULATION MODEL 

# by Jory Hecht, last updated on July 13, 2019

# This program is designed to simulate alternative operation policies for conventional 
# hydropower reservoirs generating baseload electricity according to prescribed 
# operation rules. 

# CREDITS: Some of the code for fititng the 7Q10 low flows was developed from a code
# base provided by Annalise Blum

################################ IMPORT DATA ####################################

# Set working directory for laptop
setwd("C:/Users/joryh/ownCloud/documents/Tufts/Hydro_Alteration/HA_HypTest/Codes")

# Load packages
library("lmomco")  # For computing 7-day, 10-year low flow (7Q10)
library("lubridate") # For handling dates nicely

# Import inflow time series
Q_in_table <- read.csv("Q_predam_02080500.csv")

# Format dates
Q_in_table$Date <- as.Date(Q_in_table$Date,format="%m/%d/%Y")

# Remove NA's
Q_in_table <- Q_in_table[!is.na(Q_in_table$Date),]

# Limit table to observations within pre-dam record (10/1/1912 - 9/30/1949)
Q_in_table <- Q_in_table[Q_in_table$Date >= as.Date("10/1/1912",format="%m/%d/%Y") &
              Q_in_table$Date <= as.Date("9/30/1949",format="%m/%d/%Y"),]

# Create vectors from data frame
dates <- Q_in_table$Date

# Create inflow vector
Q_in <-  Q_in_table$Q_mean

# Check for missing data
Q_in_missing_bool <- is.na(Q_in)
Q_in_missing_num <- as.numeric(Q_in_missing_bool)
if (sum(Q_in_missing_bool) > 0) {
  stop ('Missing data! Check input.')   
}

# Compute mean flow
Q_mean <- mean(Q_in)

# Daily average inflows and precipitation onto reservoir surface 
plot(Q_in,type="l")

# Number of days
nd <- length(Q_in)
day_por <- seq(1,nd,1)  #Date in time series for plotting purposes


####################### Compute 7Q10 low flow ##############################

# The 7Q10 (7-day low flow with a 10-year recurrence interval) is the minimum release 
# allowed from the dam unless there is an extreme drought

# Compute 7-day mean flows
Q_in_7d <- vector("numeric",length(length(Q_in)-6))
for (i in 1:(length(Q_in)-6)){
  Q_in_7d[i] = mean(Q_in[i:(i+6)])
}

# Truncate dates to remove dates not in center (day 4) of any 7-day period in record
Q_in_7d_dates <- data.frame(year(dates[4:13511]),dates[4:13511],Q_in_7d) 
names(Q_in_7d_dates) <- c("Year","Date","Q_in_7d")


#dates_7d <- dates
#dates_7d_yr <- dates_7d$year
#Q_in_7d_dates <- cbind(dates_7d_yr[4:13511],Q_in_7d)

# Extract minimum 7-day low flows for each calendar year
tmp <- split(Q_in_7d_dates$Q_in_7d,Q_in_7d_dates$Year)
tmp_min <- lapply(tmp,min)
Q_min_7d_yrs <- matrix(unlist(tmp_min), ncol = 1, byrow = TRUE)

# Remove first and last year, which are not complete
Q_min_7d_yrs <- Q_min_7d_yrs[2:37,] #length of Q_min_7d_yrs = 38

# Sort annual seven-day low flows 
Q_min_7d_ranked <- sort(Q_min_7d_yrs)

# Compute Weibull non-exceedance probabilities
p_min_7d_ranked <- rank(Q_min_7d_ranked)/(length(Q_min_7d_ranked)+1)

# Compute 10-year, 7-day low flow using GP2, GP3, kappa, and Wakeby distributions

# Compute L-moments
Q_min_7d_lmoms <- lmoms(Q_min_7d_yrs,nmom=5)

# Make function to get pargpa for xi=0 -->GP2
parGP2<-function(x){ 
  pargpa2<-pargpa(x,xi=0)
  return(pargpa2)
}

par_lf_gp2 <- parGP2(Q_min_7d_lmoms)
par_lf_gp3 <- pargpa(Q_min_7d_lmoms)
par_lf_kap <- parkap(Q_min_7d_lmoms)
par_lf_wak <- parwak(Q_min_7d_lmoms)

# Check for parameter feasibility
parGP2valid <-are.pargpa.valid(par_lf_gp2) #True
parGP3valid <-are.pargpa.valid(par_lf_gp3) #True
parkapvalid <-are.parkap.valid(par_lf_kap) #True
parwakvalid <-are.parwak.valid(par_lf_wak) #True

# Compute quantile functions for each distribution
FGP2_est <- quagpa(p_min_7d_ranked,par_lf_gp2)
FGP3_est <- quagpa(p_min_7d_ranked,par_lf_gp3) # Note: Can use quagpa for 2- ad 3-param GP distributions
FKap_est <- quakap(p_min_7d_ranked,par_lf_kap)
FWak_est <- quawak(p_min_7d_ranked,par_lf_wak)

# Compute goodness-of-fit in terms of PPCC
FGP2_QQ_PPCC <- cor(FGP2_est,Q_min_7d_ranked) #0.897
FGP3_QQ_PPCC <- cor(FGP3_est,Q_min_7d_ranked) #0.970
FKap_QQ_PPCC <- cor(FKap_est,Q_min_7d_ranked) #0.987
FWak_QQ_PPCC <- cor(FWak_est,Q_min_7d_ranked) #0.988

plot(FGP3_est,Q_min_7d_ranked)
plot(FKap_est,Q_min_7d_ranked) 

# Compute 7Q10
FGP2_7Q10 <- quagpa(0.1,par_lf_gp2) # 650
FGP3_7Q10 <- quagpa(0.1,par_lf_gp3) # 967
FKap_7Q10 <- quakap(0.1,par_lf_kap) # 1012
FWak_7Q10 <- quawak(0.1,par_lf_wak) # 923

Q90_pre <- quantile(Q_in,0.10)

# ratio_7Q10_Qmean <- FKap_7Q10/mean(Q_in)
ratio_7Q10_Qmean <- Q90_pre/mean(Q_in)


########################## COMPUTE MONTHLY FLOW REQS ##############################

# Based on Tessmann (1980) cited from Pastor et al. (2014)

# Pastor et al. (2014) 
# https://www.hydrol-earth-syst-sci.net/18/5041/2014/hess-18-5041-2014.pdf

# Tessmann (1980): Tessmann, S.: Environmental assessment, technical appendix e in
# environmental use sector reconnaissance elements of the western
# dakotas region of south dakota study. South dakota state university, Water Resources Institute, South Dakota State University,
# Brookings, South Dakota, 1980.

# Identify low, intermediate and high flow months

# Compute mean monthly flows (discharges)
df_Q_in <- data.frame(dates,Q_in)
month(df_Q_in$dates)

MMF <- vector(mode="numeric",length=12)
for(j in 1:12){
  MMF[j] = mean(df_Q_in[month(df_Q_in$dates) == j,]$Q_in)
}

MAF <- mean(Q_in)

MMF_MAF <- MMF/MAF

EFR_frac <- vector(mode="numeric",12)

for(j in 1:12){
  if(MMF_MAF[j] <= 0.4){
    EFR_frac[j] = 0.6
  }
  if(MMF_MAF[j] > 0.4 & MMF_MAF[j] <= 0.8){
    EFR_frac[j] = 0.5
  } else {
    EFR_frac[j] = 0.4
  }
}

df_Q_in$Q_out_min <- EFR_frac[month(df_Q_in$dates)]*MMF[month(df_Q_in$dates)]

Q_out_min <- df_Q_in$Q_out_min

# Note: no low-flow months where MMF < 0.4 * MAF

###################### ON-LAKE PRECIPITATION AND EVAPORATION ###########################

# Precipitation and evaporation
P_lake <- rep(0,length=length(Q_in))
E_pan <- rep(0,length=length(Q_in))

# Evaporation parameters 
pan_evap = 0;
pan_to_lake = 0;

########################## CHOOSE RESERVOIR DESIGN PARAMETERS ########################

# Dead storage ratio        
SR_dead = 0.00; 

# Live storage ratio
SR_live = 1.00; 

# Flood storage ratio
SR_flood = 0.20; 
# FIX: Add error statement to make sure SR_flood > SR_npool

# Low-flow outlet ratio (outlet flow capacity/mean annual flow) 
lfo_ratio = 0.1; 

# Turbine ratio (total turbine flow capacity/mean annual flow)
turb_ratio = 1; 

# Spillway ratio (spillway discharge capacity/mean annual flow)
spway_ratio = 100;

# Efficiency of power generation
eff = 0.80; 

# Min turb flow fraction that prevents releases through turbines when storage available 
# for release is less than a given percentage of the total turbine flow capacity
# This parameter can be set lower if there are many turbines in a dam, i.e. have only 2 out of 10 
# turbines open

turb_min_pct = 0.5;

# DAILY RESERVOIR SIMULATION MODEL ----------------------------------------

# Initialize vectors 
Q_lfo <- vector(mode = "numeric", length=length(Q_in))
Q_turb <- vector(mode = "numeric", length=length(Q_in))
Q_spway <- vector(mode = "numeric", length=length(Q_in))
Q_out <- vector(mode = "numeric", length=length(Q_in))
SA_lake <- vector(mode = "numeric", length=length(Q_in))
E_lake <- vector(mode = "numeric", length=length(Q_in))
head <- vector(mode = "numeric", length=length(Q_in))
hpower <-vector(mode= "numeric", length=length(Q_in))
S <- vector(mode = "numeric", length=length(Q_in)+1) # 1 longer to have initial storage
S_live_pct <- vector(mode = "numeric", length=length(Q_in))
WL <- vector(mode = "numeric", length=length(Q_in))

# Variable parameters
op_rule <- c(1,2)
SR_live <- c(0.01,0.243) #storage ratio measured in years - includes both live (conservation) and dead storage
# 0.243 is based on John H. Kerr
turb_ratio <- c(0.05,1.1)
lfo_ratio<- c(0.05,ratio_7Q10_Qmean) # Replaced later with Q_out_min
turb_min_pct <- c(0.05,0.2)
head_min <- c(20,50) #Leave these after changing other units from ft to m

# Initialize parameter array
params <- matrix(nrow=prod(length(op_rule),length(SR_live),length(turb_ratio),length(turb_min_pct),length(lfo_ratio),length(head_min)),ncol=6)

rr = 0; #Change to expand.grid? 

for (hh in 1:length(op_rule)){
  for (ii in 1:length(SR_live)){
    for (jj in 1:length(turb_ratio)){
      for (kk in 1:length(turb_min_pct)){
        for (ll in 1:length(lfo_ratio)){
          for (mm in 1:length(head_min)){
            rr = rr + 1
            params[rr,1]=op_rule[hh]
            params[rr,2]=SR_live[ii]*Q_mean*86400*365.24
            params[rr,3]=turb_ratio[jj]*Q_mean
            params[rr,4]=turb_min_pct[kk]*turb_ratio[jj]*Q_mean
            params[rr,5]=lfo_ratio[ll]*Q_mean
            params[rr,6]=head_min[mm]
          }
        }
      }
    }
  }
}

# Export params as csv
write.csv(params,"params.csv",row.names)

# Create arrays for output variables
Q_out_res <- array(NA,dim=c(nd,nrow(params)))
Q_lfo_res <- array(NA,dim=c(nd,nrow(params)))
Q_turb_res <- array(NA,dim=c(nd,nrow(params)))
Q_spway_res <- array(NA,dim=c(nd,nrow(params)))
SA_lake_res <- array(NA,dim=c(nd+1,nrow(params)))
E_lake_res <- array(NA,dim=c(nd,nrow(params)))
WL_res <- array(NA,dim=c(nd+1,nrow(params)))
head_res <- array(NA,dim=c(nd,nrow(params)))
S_res <- array(NA,dim=c(nd+1,nrow(params)))
hpower_res <- array(NA,dim=c(nd,nrow(params)))

# Designate reservoir operation rule
#op_rule = 2;
# 1 = energy maximization, # 2 = flow variability preservation
foi_min = 1; # Sets minimum outflow as a fraction of inflow
foi_max = 1 # Sets maximum outflow as a fraction of inflow
# Downstream flow


############################## RUN SIMULATION ####################################

'
Turbine releases under Operating Rule 1 (energy maximization with constant min flow)   

Release maximum flow through turbines when inflow is greater than turbine capacity
Otherwise, release decisions take reservoir storage into account

When storage is greater than 10% of the live storage capacity + the minimum daily turbine 
discharge, then release the minimum of (a) the turbine discharge capacity, and (b) the sum 
of the inflow and storage available for release

If the storage available for release is less than the minimum that can be passed through
turbines, then do not release any water through the turbines

If the storage is below 10%, do not release any water through turbines 
to maintain storage in reservoir 
'

for (b in 1:nrow(params)){    #Combination of parameters

  # These lines rename the parameters to make the reserovir model more readable
  op_rule <- params[b,1]
  S_live_cap <- params[b,2]
  Q_turb_cap <- params[b,3] 
  turb_min_pct <- params[b,4]
  Q_lfo_cap <- params[b,5]
  head_min <- params[b,6]
  
  # Initialize storage - assume reservoir is full at the start of operations
  S[1] = 1*S_live_cap
  
  for(t in 1:nd){
  
  # Calculate flow passing through turbines for operating rule 1 (energy max with constant min flow)

  if(op_rule == 1){
    if(S[t] > 0.3234*S_live_cap + turb_min_pct*86400){
      Q_turb[t] = min(Q_turb_cap,(S[t]-0.3234*S_live_cap)/86400)
    } else {
      Q_turb[t] = 0
    }
  }
    
 
  # For operating rule 2 (preserving flow variability)
  if (op_rule == 2){
    if (Q_in[t] > Q_turb_cap){
      Q_turb[t] = Q_turb_cap
    } else if (S[t] > 0.3234*S_live_cap + turb_min_pct*86400){
      Q_turb[t] <- ifelse(Q_in[t]+(S[t]-0.3234*S_live_cap)/86400 > turb_min_pct,
                          min(Q_turb_cap,foi_max*Q_in[t],Q_in[t]+(S[t]-0.3234*S_live_cap)/86400),0)
    } else {Q_turb[t] = 0}
  }

    
  # Calculate flow passing through spillway  
    
    if ((S[t]+ (Q_in[t]-Q_turb[t])*86400 + P_lake[t] - E_lake[t] - S_live_cap)/86400 > 20*spway_ratio*Q_mean){
      stop ('Spillway capacity insufficient. Dam might fail.')
    }
    else if (S[t]+(Q_in[t]-Q_turb[t])*86400+P_lake[t]-E_lake[t] > S_live_cap){
      Q_spway[t] = (S[t]+(Q_in[t]-Q_turb[t])*86400+P_lake[t]-E_lake[t]-S_live_cap)/86400
    }
    else {Q_spway[t]=0} 
    
  # Compute flow passing through the low-flow outlet
  
    # Operating rule 1: Calculate flow passing through low-flow outlet
    if (op_rule == 1){
      
      Q_lfo[t] = ifelse(Q_turb[t] + Q_spway[t] < Q_out_min[t], min(max(0,(S[t]/86400-0.01)-(Q_turb[t]+Q_spway[t])),
                                                                   Q_out_min[t]-(Q_turb[t]+Q_spway[t])),
                                                                   0)
      
      # Constant of 0.01 prevents storage from reaching zero and avoids problem with water level computation
      # a surface area of zero
      
    # This rule releases the difference between the turbine release and the low-flow outlet
    # capacity when 
      
    # Operating rule 2: Calculate flow passing through low-flow outlet 
    } else if (op_rule == 2){ 
      
      Q_lfo[t] = ifelse(Q_turb[t] + Q_spway[t] < foi_min*Q_in[t],
                        ifelse(S[t] < 0.3234*S_live_cap,
                               min(Q_lfo_cap,max(0,(S[t]-0.3234*S_live_cap)/86400)),
                               min(Q_lfo_cap,foi_min*Q_in[t]-(Q_turb[t]+Q_spway[t]))),0)
                          
    } else {stop('Error: Invalid operating rule!')}    
     
    # If the difference between storage at beginning of day t and the sum of the turbine outflow 
    # and low-flow discharge capacity is greater than zero
      #Then if the turbine outflow is less than the low-flow outlet capacity (min low-flow req),
      #release the difference between through the low-flow outlet (require constant release even if inflow lower)
      #Otherwise, the turbine release is sufficient for maintaining min flows 
    # Otherwise release the remaining live storage
    
   
      # If difference between inflow and outflow is less than or equal to a given 
      # fraction of the inflow, then there is no need to release more water through 
      # the low-flow outlet
      # Otherwise the low-flow release ensures that the minimum flow reduction percentage
      # is achieved, unless the low-flow outlet capacity limits the volume that can be \
      # released to meet this target
    
    


    # Compute outflows through various outlets
    Q_out[t] = Q_lfo[t]+Q_turb[t]+Q_spway[t]
    Q_lfo_res[t,b] = Q_lfo[t]
    Q_turb_res[t,b] = Q_turb[t]
    Q_spway_res[t,b] = Q_spway[t]
    Q_out_res[t,b] = Q_out[t]

    # Compute water level for head estimate (elevation difference changes)
    # Generic water supply reservoir surface area-volume formula from Takeuchi (1997)
    #See Excel spreadsheet
    
    S_live_pct[t] = S[t]/S_live_cap*100
    
    WL[t] = 3*10^-5*S_live_pct[t]^3 - 0.0091*S_live_pct[t]^2 + 1.2248*S_live_pct[t]+237.94
    
    WL_res[t,b] = WL[t]

    # Compute head 
    
    # Assume constant tailwater elevation
    tw_elev <- 204
  
    head[t] = WL[t] - tw_elev
  
    head_res[t,b] = head[t]
  
    #Reservoir water balance 
    S[t+1] = S[t] + (Q_in[t] - Q_out[t]) * 86400 + P_lake[t] - E_lake[t];
    S_res[t,b] = S[t+1] #Initial storage

    # Compute average power generation on day t in kWh

    hpower[t] <- 0.001*24*9.81*eff*(Q_turb[t]/35.3*1000)*(0.9*head[t]/3.281)
  
        # From http://www.renewablesfirst.co.uk/hydropower/hydropower-learning-centre/how-much-power-could-i-generate-from-a-hydro-turbine/
        # 0.001 converts watts to kilowatts
        # 24 hours per day
        # 9.81 is the rate of gravitational acceleration at the Earth's surface
        # 35.3 converts the flow from cfs to m3/s, 1000 converts the flow from m3/s to L/s
        # 0.9 reduces the gross head to account for head losses
        # 3.281 converts the head from ft to meters
  
    hpower_res[t,b] = hpower[t]
  
    # Compute installed capacity
    head_rated <- 96
    HP_cap_inst <- 0.001*9.81*eff*(Q_turb_cap/35.3*1000)*(0.9*head_rated/3.281)
    
    
  }

}

# END OF LOOP

########################## COMPUTE OUTPUT STATISTICS ################################

# Min and average daily hydropower production
min_hpower <-apply(hpower_res,2,min)
avg_hpower <-apply(hpower_res,2,mean)

# Compute percentage of discharege passing through each outlet
Q_turb_pct_tot <- apply(Q_turb_res,2,sum)/apply(Q_out_res,2,sum)
Q_lfo_pct_tot <- apply(Q_lfo_res,2,sum)/apply(Q_out_res,2,sum)
Q_spway_pct_tot <- apply(Q_spway_res,2,sum)/apply(Q_out_res,2,sum)

# Compare releases from turbine, bypass (lfo), and spillway on each day
plot(Q_spway_res[,32],Q_lfo_res[,32]) # A few days with both bypass releases and spills under 32!!!
plot(Q_turb_res[,32],Q_lfo_res[,32]) # Bypass releases during low-flow periods, bypass releases when turbine capacity exceeded
plot(Q_spway_res[,32],Q_turb_res[,32]) # Spills only occur when turbine release complete

# Hydropower regret (Operating Rule 1 vs. Operating Rule 2)
regret_min_hpower <- min_hpower[1:32]-min_hpower[33:64]
regret_avg_hpower <- avg_hpower[1:32]-avg_hpower[33:64]

# Relative differences over baseline (minimum value)
ob_min_hpower <- min_hpower-min(min_hpower)
ob_avg_hpower <- avg_hpower-min(avg_hpower)

ob_min_hpower_pct <- (min_hpower-min(min_hpower))/min(min_hpower)-1
ob_avg_hpower_pct <- (avg_hpower-min(avg_hpower))/min(avg_hpower)-1

# Compare hydropower differences between operating rules

# Plot reservoir inflows, outflows
par(mar=c(5,5,4,2)+0.1)
plot(day_por/365.24,Q_in/35.31,type="l",log="y",ylim=c(10,10000),
     xlab="Years",cex.axis=0.80,
     ylab= expression(paste("Outflow (m"^"3","/s)",sep="")),las=1,
     main="Pre- and post-dam daily flows",col="gray80") 
lines(day_por/365.24,Q_out_res[,32]/35.31,col="black")
legend("bottomleft",c("Pre-dam","Post-dam"),lty=c(1,1),col=c("gray80","black"),bty="n",cex=0.8)
par(mar=c(5,4,4,2)+0.1)

# Plot low-flow releases
plot(day_por,Q_lfo_res[,32]+0.01,type="l",log="y",xlim=c(0,nd),ylim=c(0.01,300000),
     xlab="Days",
     ylab="Discharge (cfs)",
     main="Daily low-flow releases")

# Plot water levels
plot(seq(1,nd+1,1)/365.24,WL_res[,32]/3.281,type="l",
     xlab="Year of dam operations",
     ylab="Water Level (m)",
     main="Daily reservoir water level")

# Plot daily hydropower production
plot(hpower_res[,32],type="l",
     xlab="Days",
     ylab="Energy prodcution (kWh)",
     main="Daily hydropower production")

# Compute number of days with hydropower production
days_hpower <- length(hpower_res[,32][hpower_res[,32]>0])

# Compute number of days with spills
length(Q_spway_res[,32][Q_spway_res[,32]>0])

# Compute number of days with minimum low-flow releases
length(Q_lfo_res[,32][Q_lfo_res[,32]==Q_lfo_cap])

# Compute relationship between low-flow releases and hp production
plot(Q_lfo_res[,32]+0.1,hpower_res[,32],log="x")

# Compute plant factor
sum(hpower_res[,32])/24/(HP_cap_inst*nd)

# Compute average annual energy generation in GWh
sum(hpower_res[,32])/37/10^6









