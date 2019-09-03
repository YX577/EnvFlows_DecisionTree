FDC.overall <- function(q,p=seq(1,365)/366){            #Should we remove Feb 29 obs here?
  
  # Variables:
  # i: The integer component of (n+1)*p.
  # n: The length of the record (length of q.i array)
  # p: (INPUT) Desired exceedence probabilities. Defaults to 1/366-356/366.  
  # q: (INPUT) Streamflow time series array.
  # q.i: The q array, arranged from largest to smallest. Following Vogel & Fennessey (1994), the smallest flow is a zero (added to end of array).
  # q.order: The order of the q array, from largest to smallest.
  # Qp: (OUTPUT) The array of the q.i flows falling at p exceedence probabilities.
  # theta: The decimal component of (n+1)*p.
  
  q <- q[!is.na(q)] #removes NA from daily flow time series, but do we want time series with them?
  # Add line to count NA, stop if there are any? 
  q.order <- order(q,decreasing=T) #computes ranks of the flows in decreasing order
  q.i <- q[q.order] #creates a new vector that sorts the flows in decreasing order
  q.i <- c(q.i,0) #following Vogel and Fennessey (1994), a value of 0 is added
  
  n = length(q.i) #n is the length of the daily flow series.
  
  #   p.nat <- seq(n,1,by=-1)/(n+1)
  
  # Using the Qp,2 method from Vogel & Fennessey (1994)
  i <- floor((n+1)*p)
  theta <- ((n+1)*p-i)
  Qp <- array(NA,length(p))#initializes empty array of length p
  Qp[which(i>0)] <- (1-theta[which(i>0)])*q.i[i[which(i>0)]]+
    theta[which(i>0)]*q.i[i[which(i>0)]+1]
  
  return(list(Q=Qp,p=p)) #Qp series goes under the Q header, p seriers goes under the p header
}