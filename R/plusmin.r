# Function for elect package 

# Ardo van den Hout, Cambridge 2010 - UCL 2019
# Mei Sum Chan, UCL 2018


##################################################
# Function to compute function of two LEs:
plusmin <- function(x,index=NA,func="plus",
                          probs=c(.025,0.5,.975),digits=3){
  # <LEs> = result from life.exp()
  # <index> = indicates the LEs in column with LEs
  # <func>  = required series of "plus" and "minus"
  # <probs> = probs for quantiles
  # <digits> = number of digits in output


  # Rename object:
  LEs <- x

  # Was simulation undertaken?:
  if(is.na(LEs$sim[1])){ SIM <- FALSE }else{ SIM <- TRUE }
  if(SIM==FALSE){
    stop("\nERROR. This function requires that uncertainty of LEs is
         estimated. Choose S > 0 in elect-call.\n\n")
  }

  # Start totals:
  total <- LEs$sim[,index[1]]
  pnt   <- LEs$pnt[index[1]]

  # Apply function(s):
  for(i in 1:length(func)){
    # Applying "plus":
    if(func[i]=="plus"){
      total <- total + LEs$sim[,index[i+1]]
      pnt   <- pnt + LEs$pnt[index[i+1]]
    }
    # Applying with "minus":
    if(func[i]=="minus"){
      total <- total - LEs$sim[,index[i+1]]
      pnt   <- pnt - LEs$pnt[index[i+1]]
    }
  }

  # Summarise:
  mn <- mean(total)
  se <- sd(total)
  quants <- rep(NA,length(probs))
  for(j in 1:length(probs)){
    quants[j] <- quantile(total,probs=probs[j])
  }
  out <- c(pnt,mn,se,quants)

  # Labels:
  out <- as.data.frame(matrix(out,nrow=1,ncol=length(out)))
  names <- rep(NA,3+length(probs))
  names[1:3] <- c("pnt","mn","se")
  for(j in 4:(3+length(probs))){
    names[j] <- paste(probs[j-3],"q",sep="")
  }
  row.names(out) <- "func(LEs)"
  names(out) <- names

  # Print:
  print(round(out,digits))

  }


# For compatibility with older versions:
plusmin.elect <- plusmin