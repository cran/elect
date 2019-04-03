# Function for elect package 

# Ardo van den Hout, Cambridge 2010 - UCL 2019
# Mei Sum Chan, UCL 2018

#######################################################
# Function to plot life.exp results:
plot.elect <- function(x,which=NULL,kernel="gaussian",col="red",
                       lwd=2,cex.lab=1,...){

   ##############################################
   # <LEs> = result from elect()
   # <which> = select which of the LEs to plot (order as in summary)
   # <col> = colour for lines
   # <lwd> = lwd for lines

   # Rename object:
   LEs <- x

   ###############################################
   # Two functions for when Which is not NULL:
   # Check whether a whole number (code from R help):
   is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
     abs(x - round(x)) < tol
   }
   # Compute dimensions for two-dim graph display given A graphs:
   gdim <- function(A){
     sqA    <- sqrt(A) 
     fl.sqA <- floor(sqA)
     if(is.wholenumber(sqA)){
       # Easy when A is a square:
       dims <- c(sqA,sqA)
     }else{
       # Minimal dims when A is not a square:
       ind <- 1
       dims <- c(fl.sqA, fl.sqA + ind)
       while( (fl.sqA * (fl.sqA+ind)) < A){
         ind <- ind + 1
         dims <- c(fl.sqA, fl.sqA + ind)
       }
     }
     return(dims)
   }
   
   ###############################################
   # No plot without simulation (S=0):
   if(LEs$S==0){
     cat("\nNo simulated LEs so no plot. \n\n")
   }
   
   ###############################################
   # Plot if S > 0:
   if(LEs$S > 0){

   # Number of states:
   model   <- LEs$model
   nstates <- nrow(model$Qmatrices$baseline)
 
   # Plot for any number of states:
   if(length(LEs$pnt)==((nstates-1)^2+nstates)){
     # No selection using <which>:
     if(is.null(which)){
     if(nstates <= 4){
       mfrow.dim <- c(nstates-1,nstates+1)
     }else{
       mfrow.dim <- c(nstates-2,nstates+2)
     }
     # Specify graphical parameters: 
     opar <- par(mfrow=mfrow.dim, mex=0.8,mar=c(5,5,2,1)+.1)
   for(i in 1:((nstates-1)^2+nstates)){
     plot(density(LEs$sim[,i],kernel=kernel),main=names(LEs$pnt[i]),
          ylab="density",xlab="years",col=col,lwd=lwd)
     }
     # Back to default setting of graphical parameters:
     opar <- par(mfrow=c(1,1), mex= 1, mar=c(5, 4, 4, 2) + 0.1)
     }
     
     # Selection using <which>:
     if(!is.null(which)){
       mfrow.dim <- gdim(length(which))
       # Specify graphical parameters: 
       opar <- par(mfrow=mfrow.dim, mex=0.8,mar=c(5,5,2,1)+.1)
       for(i in which){
         plot(density(LEs$sim[,i],kernel=kernel),main=names(LEs$pnt[i]),
              ylab="density",xlab="years",col=col,lwd=lwd)
       }
       # Back to default setting of graphical parameters:
       opar <- par(mfrow=c(1,1), mex= 1, mar=c(5, 4, 4, 2) + 0.1)
     }
   }
   }

  ###############################################
  # Return nothing:
  invisible(x)
}


