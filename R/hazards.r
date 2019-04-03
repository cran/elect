# Function to plot age-dependent hazards for
# a model fitted with 'msm'

# Ardo van den Hout, Cambridge 2010 - UCL 2019
# Mei Sum Chan, UCL 2018

# With thanks to Ying Lou for discussing functionality

#################################################################
hazards <- function(x, b.covariates ,
                    no.years, trans = NULL,
                    max.haz = .5, min.haz =0, CI = FALSE,
                    col = NULL, lty = NULL, lwd = NULL,
                    LEGEND = TRUE, location="topleft",
                    age.shift = 0){

# Explanation:
# <x>  = model fitted with 'msm'. Should contain <age>
#            as time-dependent covariate
# <trans>  = transition hazard(s) to be plotted. Default to all
# <col>    = colour of curve(s)
# <type>   = line type of curve(s)


####################
# Redefine:
model <- x

# Redefine if needed:
if(is.vector(trans)){
    trans <- matrix(trans,1,2)
}

##############################
# Determine trans if needed:
if(!is.null(trans)){
  ntrans <- nrow(trans)
}else{
  Q <- model$qmodel$imatrix
  ntrans <- sum(Q)
  trans <- matrix(NA,ntrans,2)
  index <- 1
  for(i in 1:nrow(Q)){
    for(j in 1:nrow(Q)){
      if(Q[i,j]){
        trans[index,]<- c(i,j)
        index <- index + 1
      }
    }
  }
}

##########################
# Input checks:
if(!is.null(col)){
  if(length(col)!=nrow(trans)){
    stop("\nSpecify <col> for rows in <trans> only.\n\n")
  }
}
if(!is.null(lty)){
  if(length(lty)!=nrow(trans)){
    stop("\nSpecify <lty> for rows in <trans> only.\n\n")
  }
}
if(!is.null(lwd)){
  if(length(lwd)!=nrow(trans)){
    stop("\nSpecify <lwd> for rows in <trans> only.\n\n")
  }
}
# Check no. covariates (same check as in elect.r):
nstates <- nrow(model$Qmatrices$baseline)
Q.null  <- matrix(as.numeric(model$Qmatrices$baseline!=0),nstates,nstates)
diag(Q.null) <- 0
ntrans <- sum(Q.null)
ncovs <- 1+length(b.covariates)
nbeta <- max(which(names(model$estimates)=="qcov"))
# Check:
if(nbeta!=(ntrans*ncovs)){
  stop("\nAll covariates in msm model should be specified.\n\n")
}

# Define age grid:
b.age <- b.covariates$age
age.grid <- seq(b.age, b.age+no.years, by=1/12)

# Colours and types:
if(is.null(col)){col <- rep(1,ntrans)}
if(is.null(lty)){lty <- 1:ntrans}
if(is.null(lwd)){lwd <- rep(2,ntrans)}

# Scale in plot:
plot.age.grid <- age.grid - age.shift


# Plot transition hazards:
L <- length(age.grid)
plot(x=c(plot.age.grid[1],plot.age.grid[L]),y=c(min.haz,max.haz),
     type="n", xlab = "Age", ylab = "Hazards",las=1)
# Loop over the hazards:
for(i in 1:ntrans){
  haz   <- rep(NA,L)
  # Lower and upperbound for CI (if needed):
  hazLB <- rep(NA,L)
  hazUB <- rep(NA,L)
  for(j in 1:L){
    if(ncovs==2){
      covariates <- list(age = age.grid[j])
    }else{
      covariates <- c(age = age.grid[j], b.covariates[2:length(b.covariates)])
    }
    haz[j] <- qmatrix.msm(model,
                    covariates = covariates,
                    ci="none")[trans[i,1],trans[i,2]]
    # Compute CIs if asked for:
    if(CI){
      Q.CI <- qmatrix.msm(model, covariates = covariates,
                          ci="delta")
      hazLB[j] <- Q.CI$L[trans[i,1],trans[i,2]]
      hazUB[j] <- Q.CI$U[trans[i,1],trans[i,2]]
     }
    }
  lines(plot.age.grid,haz,lwd=lwd[i],col=col[i],lty=lty[i])
  # Plot CIs if asked for:
  if(CI){
    lines(plot.age.grid,hazLB,lwd=lwd[i]/2,col=col[i],lty=lty[i])
    lines(plot.age.grid,hazUB,lwd=lwd[i]/2,col=col[i],lty=lty[i])
    }
  }


# Add legend or not:
if(LEGEND){
  legend <- rep(NA,ntrans)
  for(i in 1:ntrans){
    legend[i] <- paste(as.character(trans[i,1]), "to",
                       as.character(trans[i,2]))
  }
  legend(x= location, legend = as.character(legend) , col = col,
         text.col = 1, lty = lty, lwd=lwd,
         merge = TRUE, bg = "white")
}


}




