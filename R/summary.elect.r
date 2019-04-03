# Functions for elect package 

# Ardo van den Hout, Cambridge 2010 - UCL 2019
# Mei Sum Chan, UCL 2018


#######################################################
# Function to summarise life.exp results:
summary.elect <- function(object,probs=c(.025,0.5,.975),digits=3,
                        StartStateTotals=FALSE, print=TRUE,
                        sd.model=FALSE,...){
   # <LEs> = result from life.exp()
   # <probs> = probs for quantiles
   # <digits> = number of digits in output
   # <StartStateTotals> = TRUE for output on start-state totals e_r (for S>0).
   # <print> = TRUE of output directly to screen

   # Rename object:
   LEs <- object

   # Use point estimate:
   pnt <- LEs$pnt

   # Was simulation undertaken?:
   if(is.na(LEs$sim[1])){ SIM <- FALSE }else{ SIM <- TRUE }

   # Use simulation if it took place:
   if(!SIM){ out <- as.data.frame(pnt)}
   if(SIM){
     mn <- apply(LEs$sim,2,mean)
     se <- apply(LEs$sim,2,sd)
     quants <- matrix(NA,ncol(LEs$sim),length(probs))
     for(i in 1:ncol(LEs$sim)){
       for(j in 1:length(probs)){
        quants[i,j] <- quantile(LEs$sim[,i],probs=probs[j])
       }
     }
     out <- as.data.frame(cbind(pnt,mn,se,quantiles=quants))
     for(j in 4:(3+length(probs))){
       names(out)[j]<-paste(probs[j-3],"q",sep="")
     }
   }

   # To print or not to print:
    if(print){
     cat("\n-----------------------------\n")
     cat("elect summary\n")
     cat("-----------------------------\n")
     cat("Covariate values in the multi-state model:\n")
     print(unlist(LEs$covars))
     # Output the sd.model if asked for (and if fitted):
     sd.model.fitted <-  !(is.na(LEs$sd.model)[1])
     if(sd.model.fitted & sd.model==FALSE){
       # Output covariates used in sd.model:
       sd.covars <- attr(LEs$sd.model$terms,"term.labels")
       cat("Covariates in the state-distribution model:\n  ",sd.covars,"\n\n")
     }
     if(sd.model.fitted & sd.model==TRUE){
       # Output covariates used in sd.model:
       sd.covars <- attr(LEs$sd.model$terms,"term.labels")
       cat("\nFitted model for state-distribution model:\n")
       print(summary(LEs$sd.model),digits=digits)
       cat("\n")
       cat("Predicted distribution for start states\n")
       cat("given specified value(s) for ", sd.covars, ":\n")
       cat(round(LEs$sd.model.inits, digits),"\n")
       cat("\n")
       }
     if(!sd.model.fitted & sd.model==TRUE){
       # Output info on no sd.model:
       cat("\nNo state-distribution model was fitted.\n")
     }
     cat("Life expectancies:")
     if(LEs$S>0){
      cat("\nUsing simulation with ",LEs$S,"replications\n")
      cat("\nPoint estimates, and mean, SEs, and quantiles from simulation:\n")
    }else{cat("\nPoint estimates:\n")}
    print(round(out,digits))
    cat("-----------------------------\n")

    # Produce start-state totals if asked for:
    # Number of states:
    model <- LEs$model
    nstates <- nrow(model$Qmatrices$baseline)
    if(SIM && StartStateTotals){
      cat("Start-state totals:\n")
      index <- 1
      for (k in 1:(nstates-1)){
        cols <- index:(index+(nstates-2))    #; print(cols)
        summ <- apply(LEs$sim[,cols],1,sum)  #; print(summ)
        mn <- mean(summ)
        se <- sd(summ)
        quants <- matrix(NA,1,length(probs))
        for(j in 1:length(probs)){
          quants[1,j] <- quantile(summ,probs=probs[j])
        }
        pnt.start <- sum(pnt[cols])
        out2 <- as.data.frame(cbind(pnt=pnt.start,mn,se,quantiles=quants))
        # Label rows of output:
        names.start <- paste("e",k,".",sep="") #; print(names.start)
        row.names(out2) <- names.start
        # Label columns of output:
        for(j in 4:(3+length(probs))){
          names(out2)[j] <- paste(probs[j-3],"q",sep="")
        }
        # Print output:
        print(round(out2,digits))
        # Update index:
        index <- index + (nstates-1)
      }
      cat("-----------------------------\n")
    }

    }
}

