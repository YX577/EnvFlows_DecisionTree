# RUN DECISION TREE

########################### ENTER REQUIRED INPUTS #######################################

# Set prior probabilities
P_A <- 0.5
P_NA <- 1 - P_A

# Type I error
aa <- 0.2
# Type II error
bb <- 0.6

# Cost of overprotection
C_NA_CA <- 10
C_A_CNA <- 10

########################## RUN BAYESIAN DECISION TREE #################################

Bayes_DT <- function(aa,bb,C_NA_CA,C_A_CNA,P_A){
  
  # Compute P_NA
  P_NA <- 1 - P_A
  
  # Compute loss probability from 2 x 2 confusion matrix
  P_CNA_NA <- 1-aa
  P_CA_NA <- aa
  P_CNA_A <- bb
  P_CA_A <- 1-bb
  
  # Compute probabilities of concluding alteration (CA) or no alteration (CNA)
  P_CA <- P_CA_NA * P_NA + P_CA_A * P_A
  P_CNA <- P_CNA_NA * P_NA + P_CNA_A * P_A
  
  # Compute loss probabilities using Bayes theorem
  P_NA_CA <- P_CA_NA * P_NA/P_CA
  P_A_CNA <- P_CNA_A * P_A/P_CNA
  
  # Compute expected losses
  EL_NA_CA <- P_NA_CA * C_NA_CA 
  EL_A_CNA <- P_A_CNA * C_A_CNA
  
  # 
  Decision <- if(EL_NA_CA < EL_A_CNA){
    "Change operating rule"
  } else {
    "Keep operating rule"
  }
  
  # Output 
  output <- list(Decision,P_NA_CA,P_A_CNA,EL_NA_CA,EL_A_CNA)
  names(output) <- c("Decision Rule",
                     "P(NA|CA)",
                     "P(A|CNA)",
                     "Expected Loss NA|CA",
                     "Expected Loss A|CNA")
    
  return(output)
  
}

