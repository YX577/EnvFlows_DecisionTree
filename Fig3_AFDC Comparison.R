# Figure 3: AFDC pre- and post-dam comparison

# Define annual exceedance probabilities of interest p_key
p <- seq(1,99)/100

# Import pre- and post-dam flows
Q_pre <- load(NAME)
Q_post <- load(NAME)

# Create a plot to compare pre- and post-dam FDCs
par(mfrow=c(1,2))
par(mar=c(5,5,4,2) + 0.1)

# Plot pre-dam AFDCs 
plot(p*100,AFDC.pre[1,],type="l",log="y",col="gray69",
xlim = c(0,100), ylim = c(10,5000),
xlab="Exceedance Probability (%)",ylab = expression(Discharge ~ (m^{3}/s)),cex.lab=0.9,
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
plot(p*100,AFDC.post[1,],type="l",log="y",col="gray69",

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