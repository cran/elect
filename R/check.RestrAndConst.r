# Functions for elect package

# Ardo van den Hout, Cambridge 2010 - UCL 2019
# Mei Sum Chan, UCL 2018


#######################################################
# Function to check definition of <RestrAndConst>:
check.RestrAndConst <- function(x, RestrAndConst, PRINT=FALSE){

  # Rename object:
  model <- x


  # Model parameters (for transitions):
  nbeta <- max(which(names(model$estimates)=="qcov"))
  p <- model$estimates[1:nbeta]

  # Check:
  p.RestrAndConst <- rep(NA,length(p))
  for(i in 1:length(RestrAndConst)){
     if(RestrAndConst[i]==0){
       p.RestrAndConst[i] <- 0
     }else{
       p.RestrAndConst[i] <- model$opt$par[RestrAndConst[i]]
     }
   }
   if(PRINT){
     cat("\nCheck of definition <RestrAndConst>:\n\n")
     print(cbind(Model.parameters=p,RestrAndConst.parameters=p.RestrAndConst))
   }

   if(any(is.na(p== p.RestrAndConst))){
      OK <- FALSE
   }else{
      OK <- all(p== p.RestrAndConst)
   }
   if(OK){
      cat("\n<RestrAndConst> correctly defined.\n")
   }else{
      cat("\n<RestrAndConst> not correctly defined.\n")
      if(!PRINT){cat("Use PRINT=TRUE to see the discrepancies.\n")}
      cat("\n")
   }
   return(OK)
}

# For compatibility with older versions:
check.RestrAndConst.elect <- check.RestrAndConst
