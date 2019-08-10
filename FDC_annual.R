# This function orginated from work done in Croteau et al. (2016):

FDC.annual <- function(q,dates,p=seq(1,365)/366,WY=T){
  
  # Variables:
  # dates: (INPUT) POSIXlt-formatted dates array corresponding to "q".  Only needed if type=c("annual", "median"). May be numeric if YYYYMMDD or character if "YYYY-MM-DD".  Defaults to NULL.
  # n: The length of the record (length of q.i array)
  # ny: The number of years of the record.
  # p: (INPUT) Desired exceedence probabilities. Defaults to 1/366-356/366.  
  # q: (INPUT) Streamflow daily time series array.
  # q.y: List with q separated into annual segments.
  # Qp: (OUTPUT) The annual FDC matrix be returned. 
  # WY: (INPUT) Whether the annual FDC should use water year (T) or calendar year (F).  Defaults to TRUE.
  
  
  # Check to make sure dates is in POSIXlt format
  if(class(dates)[1]=="numeric"){
    dates <- as.POSIXlt(as.character(dates),format="%Y%m%d")  #need to create a dates column in INPUT
  } else if (class(dates)[1]=="character"){
    dates <- as.POSIXlt(dates,format="%Y-%m-%d")
  }
  
  n <- length(q)   #the record length in days
  if(WY){
    ny <- dates$year[n] - dates$year[1]
  } else {
    ny <- dates$year[n] - dates$year[1] + 1
  }
  # Make n the length of q
  # If using water years, then the total number of years is...
  # If using calendar years, then the total number of years is...
  # WY is an argument of the FDC function
  
  
  q.y <- list() #empty list for storing q in annual segments
  Qp <- matrix(NA,nrow=ny,ncol=length(p)) #creates empty matrix 
  for(i in seq(1,ny)){
    if(WY){
      q.y[[i]] <- q[which((dates$year==(dates$year[1]+i) & dates$mon<9)|#go through q.y element-wise
                            (dates$year==(dates$year[1]+i-1) & dates$mon>=9))] #months are from 0-11 in R! 
    }
    else{
      q.y[[i]] <- q[which(dates$year==(dates$year[1]+i-1))]
    }
    Qp[i,] <- FDC.overall(q=q.y[[i]],p=p)$Q
  } 
  
  # If using water years, q.y[[i]] adjusts the calendar year for the water year
  
  colnames(Qp) <- paste("p",p,sep="")
  #names columns for the list
  
  return(Qp)
}