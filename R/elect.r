# Function for elect package. 

# Ardo van den Hout, Cambridge 2010 - UCL 2019
# Mei Sum Chan, UCL 2018


################################################
# Main elect function:
elect <- function(x, b.covariates, statedistdata,
                  time.scale.msm="years", h, age.max, S=0, setseed=NULL,
                  RestrAndConst=NULL, statedist.covariates="age",
                  method="step"){

  # Explanation:
  # <model> = R object defined by msm() after fitting a model
  # <b.covariates> = list with specified covarts vals (ignore intercept)
  # <time.scale.msm> = time scale in msm() model either in set {"years",
  #     "months", "weeks") or a  value  in (0,1]. Example: if you divide
  #      variable age by 10 before fitting the model, then scale = 1/10
  # <h> = grid parameter for integration and piecewise-constant model.
  #       In scale given by <time.scale.msm>
  # <statedistdata> = the data to derive distribution of living states
  # <max.age> = assumed maximum age in years. Use scale as in <model>
  # <S> = number of replications in the estimation of the variance
  #       Choose S=0 for no variance estimation
  # <setseed> = seed for the random number generation in the simulations
  # <RestrAndConst> = vector which indexes the independent parameters in
  #        model$opt$par w.r.t. to the model parameters
  # <statedist.covariates> = names of covarts for model state prevalance
  # <method> = approx. of integral: "step" for simple step function;
  #         "MiddleRiemann" or "Simpson"

  # Covariate age has to be the first in the list <b.covariates>.
  # Use the same order of covariates as in <model>
  # In the msm() model: center = FALSE, death = TRUE

  # Rename x:
  model <- x

  # Number of states and corresponding Q matrix:
  nstates <- nrow(model$Qmatrices$baseline)


##########################
# Input checks:
if(model$center==TRUE){
  stop("\nIn msm() model use argument <center = FALSE>.\n\n")
}
if(sum(model$Qmatrices$logbaseline[nstates,1:(nstates-1)])!=0){
  stop("\nModel in msm() should be illness-death model with one death state as final state.\n\n")
}
if(model$call$formula[3]!="age()"){
    stop("\nFirst covariate in msm() model should be <age>.\n\n")
}
if(is.character(time.scale.msm)){
  if(!time.scale.msm%in%c("years","months","weeks")){
    stop("\nChoose time scale to be <years>, <months>, <weeks>, or a numeric value.\n\n")
  }
}
if(is.numeric(time.scale.msm)){
  if(time.scale.msm<0 | 1<time.scale.msm){
    stop("\nNumeric for time scale has to be between 0 and 1.\n\n")
  }
}
if(h<=0){
    stop("\nChoose <h> larger than 0.\n\n")
}
if(round(S)!=S | S<0){
    stop("\n<S> should be non-negative integer. Choose <S=0> for no simulations.\n\n")
}
if(age.max<=b.covariates$age){
    stop("\nLEs<max.age> should be larger than starting age in <b.covariates>.\n\n")
}
if(!method%in%c("step","Simpson","MiddleRiemann")){
    stop("\n<method> should be <step>, <MiddleRiemann> or <Simpson>.\n\n")
}
if(!"state"%in%names(statedistdata)){
    stop("\n<state> should be variable in <statedistdata>.\n\n")
}
if(!"age"%in%names(statedistdata)){
    stop("\n<age> should be variable in <statedistdata>.\n\n")
}
if(!"age"%in%statedist.covariates){
    stop("\n <age> should be covariate for model for state prevalence.\n\n")
}
for(i in 1:length(statedist.covariates)){
  if(!statedist.covariates[i]%in%names(statedistdata)){
   stop("\n<",statedist.covariates[i],"> should be variable in <statedistdata>.\n\n")
  }
}
if(length(statedist.covariates)>length(b.covariates)){
  stop("\nNumber of covariates for model for state prevalence should not exceed number of covariates for msm model.\n\n")
}
if(is.null(model$covmat)){
  stop("\nFitted model in msm has no covariance-variance matrix. Using ELECT is not recommended.\n\n")
}

#######################
# Scale parameter to take into account the time scale:
if(time.scale.msm=="years"){    scale <- 1}
if(time.scale.msm=="months"){   scale <- 12}
if(time.scale.msm=="weeks"){    scale <- 52}
if(is.numeric(time.scale.msm)){ scale <- 1/time.scale.msm}

# Q matrix:
Q.null  <- matrix(as.numeric(model$Qmatrices$baseline!=0),nstates,nstates)
diag(Q.null) <- 0

# Number of modelled transitions:
ntrans <- sum(Q.null)

# Number of covariates (including intercept):
ncovs <- 1+length(b.covariates)

# Intensities parameters (ignoring paramtrs for misclassification):
# (Note that #parameters for transitions; may be different from #betas)
nbeta <- max(which(names(model$estimates)=="qcov"))
# Check:
if(nbeta!=(ntrans*ncovs)){
  stop("\nAll covariates in msm model should be specified in the elect call.\n\n")
}
# If check = OK:
beta  <- matrix(model$estimates[1:nbeta],ntrans,ncovs,byrow=FALSE)

# Baseline age:
age0 <- b.covariates$age
# Rest of the covariates if there are any:
if(length(b.covariates)>1){
   rest.covs <- as.vector( unlist(b.covariates[2:(ncovs-1)]))
}

# Grid for pwc hazards model and integral:
grid <- seq(0,age.max-age0,by=h)
# Shift grid for Middle Riemann rule:
if(method=="MiddleRiemann"){grid <- grid+h/2}

# Internal function for point estimates of LE:
msm.le.internal <- function(beta=beta){
  # Compute integrands:
  e <- matrix(NA,length(grid)-1,(nstates-1)^2)
  P <- diag(nstates)
  for(t in 2:length(grid)){
    # Q and P matrix:
    if(length(b.covariates)>1){
       cova <- c(1,age0+grid[t-1],rest.covs)
    }else{
       cova <- c(1,age0+grid[t-1])
    }
   # Compute Q matrix:
   Q <- Q.null
   index <- 1
   for(i in 1:(nstates-1)){
    for(j in 1:nstates){
      if(Q.null[i,j]){
        Q[i,j] <- exp(cova%*%beta[index,])
        index <- index+1
      }
    }
    Q[i,i] <- -sum(Q[i,])
   }
   # Compute P matrix:
   P <- P%*%MatrixExp(mat=Q,t=h)
   # Integrand:
   e[t-1,] <- as.vector(t(P[1:(nstates-1),1:(nstates-1)]))
  }
  # Expectancy in yrs using <method> for numerical approx to the integral:
  # Note the role of scale:
  LE <- rep(NA,(nstates-1)^2)
  # Simple grid approximation of integral:
  if(method=="step"){
       for(i in 1:length(LE)){LE[i] <- scale*sum(h*e[,i])}
  }
  # Middle Riemann sum (grid has been shifted already):
  if(method=="MiddleRiemann"){
      for(i in 1:length(LE)){LE[i] <- scale*sum(h*e[,i])}
  }
  # Simpson's rule:
  if(method=="Simpson"){
      # Work with even number of intervals (so with uneven number of nnodes):
      L <- length(grid)-1
      nnodes <- ifelse(round(L/2)!=L/2,L,L-1)
      # Do Simpson:
      for(i in 1:length(LE)){
        n <- nnodes-1
        adapt <- 1
        SS <- e[0+adapt,i]
        for(j in seq(1,n/2-1,by=1)){  SS <- SS + 2 * e[2*j+adapt,i]}
        for(j in seq(1,n/2,by=1)){    SS <- SS + 4 * e[2*j-1+adapt,i]}
        S <- S+e[n+adapt,i]
        LE[i] <- scale*h*SS/3
      }
  }
  # Return:
  LE
}

# Estimate LEs given beta:
LE <- msm.le.internal(beta=beta)


####################################
# No marginal LEs if nstates > 100 (unreasonable number?):
if(nstates>100){
   LE.pnt <- c(LE=LE)
   LEs <- NA
}


####################################
# Marginal LEs
if(nstates<=100){ # cap number of states to a reasonable number

# Distribution of states (if needed for marginal LEs):
nstatesBaseline <- length(table(statedistdata$state))
# Use model for baseline states when more than one state:
if(nstatesBaseline>1){
  # Data:
  y  <- statedistdata$state-1
  # Build formula for prevalence model:
  formula.txt <- "y~"
  formula.txt <- paste(formula.txt,statedist.covariates[1],sep="")
  n.sd.covars <- length(statedist.covariates)
  if(n.sd.covars>1){
    for(i in 2:n.sd.covars){
     add.txt     <- paste("+", statedist.covariates[i],sep="")
     formula.txt <- paste(formula.txt,add.txt,sep="")
    }
  }
   # Multinomial regression model:
  sd.model <- multinom(formula=as.formula(formula.txt),
                       data=statedistdata,trace=F)

  # Define new data for prediction with prevalence model <bmodel>:
  newdata <- statedist.covariates[1:n.sd.covars]

  # Specify values of covariates for prevalence model:
  sd.covars.values <- c(1, rep(NA,n.sd.covars) )
  for(i in 1:n.sd.covars){
    sd.covars.values[i+1] <- b.covariates[[i]]
  }

  # State distribution model: definitions for  later use:
  gamma  <- c(summary(sd.model)$coeff)
  L2     <- length(gamma)
  SigmaG <- vcov(sd.model)

}else{
 sd.model <- NA
 sd.model.inits <- NA
 sd.covars.values <- NA
}

####################################
# For simulation: point estimate and covariance for msm model:
# Remember: msm maximises 2LL. This is taken into account by the factor 1/2:
# If MC is fitted, remove the parameters for the MC model:
MCfitted <- any(names(model$opt$par)=="p")
if(!MCfitted){
  Sigma  <- solve(1/2*model$opt$hessian)
  mu     <- model$opt$par
  fixedpars <- model$fixedpars
}else{
  nparLE <-  max(which(names(model$opt$par)=="qcov"))
  Sigma  <- solve(1/2*model$opt$hessian)[1:nparLE,1:nparLE]
  mu     <- model$opt$par[1:nparLE]
  fixedpars <- 0
  if(!is.null(model$call$fixedpars) & "qcov"%in%names(model$fixedpars)){
    nfixed <- max(which(names(model$fixedpars)=="qcov"))
    fixedpars <- model$fixedpars[1:nfixed]
  }
}
L <- length(mu)

###########################
###########################
# K States:

# Marginal LEs:
mLE <- rep(NA,nstates-1)
if(nstatesBaseline>1){
   # Derive inits probs:
   gamma.index <- seq(1,(nstates-2)*length(sd.covars.values),by=nstates-2)
   # loop across states:
   inits <- rep(NA,nstates-1)
   inits.index <- c(1)
   exp.lp <- 0
   for (K in 1:(nstates-2)){
    assign(paste("lp",K,sep=""),gamma[gamma.index+(K-1)]%*%sd.covars.values)
    exp.lp <- exp.lp + exp(get(paste("lp",K,sep="")))
    inits.index <- c(inits.index,K*(nstates-1)+1)
   }
   inits[1] <- 1/(1+exp.lp)
   for (K in 1:(nstates-2)){
     inits[K+1] <- exp(get(paste("lp",K,sep="")))*inits[1]
   }
}else{
   inits <- c(1,rep(0,nstates-2))
   inits.index <- c(1)
   for (K in 1:(nstates-2)){
     inits.index <- c(inits.index,K*(nstates-1)+1)
   }
}

for (K in 1:(nstates-1)){
  mLE[K] <- inits%*%LE[inits.index+(K-1)]
}
tLE    <- sum(mLE)

# Save the initial distribution:
sd.model.inits <- inits

# Point estimate LEs:
LE.pnt <- c(LE=LE,mLE=mLE,tLE=tLE)

# No simulation:
if(S==0){LEs<-NA}

#######################################
# Estimated variance via ML simulation:
if(S>0){

# Set seed:
set.seed(setseed)

# Prelim:
LEs <- matrix(NA,S,length(LE.pnt))

# Simulation:
for(s in 1:S){
  # For msm: Sampling beta from a Multivariate Normal
  # Drawing univariate:
  z <- matrix(0,L,1)
  for(j in 1:L){z[j] <- rnorm(1,0,1)};
  # Multivariate draw:
  R <- chol(Sigma)
  p <- as.vector(mu+t(R)%*%z)
  # No restrictions:
  p0 <- p
  # Deal with fixed pars:
  if(length(fixedpars)!=0 & is.null(RestrAndConst)){
   p0 <- rep(NA,nbeta)
   p0[fixedpars] <- 0
   p0[is.na(p0)] <- p
  }
  # Deal with contraints AND fixed pars:
  if(!is.null(RestrAndConst)){
      p0 <- length(p)
      for(i in 1:length(RestrAndConst)){
         if(RestrAndConst[i]==0){
           p0[i] <-0
         }else{
           p0[i] <- p[RestrAndConst[i]]
         }
      }
    }

  # Parameters:
  betab <- matrix(p0,ntrans,ncovs,byrow=FALSE)

  if(nstatesBaseline>1){
   # For sd.model: Sampling gamma from a Multivariate Normal
   # Drawing univariate:
   z <- matrix(0,L2,1)
   for(j in 1:L2){ z[j] <- rnorm(1,0,1) }
   # Multivariate draw:
   R <- chol(SigmaG)
   gammag <- as.vector(gamma+t(R)%*%z)
   # Derive inits probs:
   # loop across states
   inits <- rep(NA,nstates-1)
   exp.lp<-0
   for (K in 1:(nstates-2)){
     assign(paste("lp",K,sep=""),gammag[gamma.index+(K-1)]%*%sd.covars.values)
     exp.lp<-exp.lp+exp(get(paste("lp",K,sep="")))
   }
   inits[1] <- 1/(1+exp.lp)
   for (K in 1:(nstates-2)){
     inits[K+1] <- exp(get(paste("lp",K,sep="")))*inits[1]
   }
}

  # LEs:
  e <- msm.le.internal(beta=betab)
  # loop across states
  c.eK<-c(e)
  sum.eK<-0
  for (K in 1:(nstates-1)){
    assign(paste("e",K,sep=""),inits%*%e[inits.index+(K-1)])
    c.eK<-c(c.eK,get(paste("e",K,sep="")))
    sum.eK<-sum.eK+get(paste("e",K,sep=""))
  }
  LEs[s,] <- c(c.eK,sum.eK)
}
}
}


# Names:
names.LE1<-c(0)
names.LE2<-c(0)
for (K1 in 1:(nstates-1)){
  for (K2 in 1:(nstates-1)){
  names.LE1<-c(names.LE1,paste("e",K1,K2,sep=""))
  }
  names.LE2<-c(names.LE2,paste("e.",K1,sep=""))
  }
names(LE.pnt)<-c(names.LE1[-1],names.LE2[-1],"e")


# Function returns list with R objects:
LEs <- list(pnt=LE.pnt,sim=LEs,h=h,covars=b.covariates,S=S,
       model=model,sd.model=sd.model,
       sd.model.inits=sd.model.inits,
       sd.covars.values = sd.covars.values)

# Class:
LEs$call <- match.call()
class(LEs) <- "elect"

# Return:
LEs
}


