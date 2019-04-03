# Function to explore data for use with 'msm' and 'elect'

# Ardo van den Hout, Cambridge 2010 - UCL 2019
# Mei Sum Chan, UCL 2018

# With thanks to Ying Lou for discussing functionality


###########################################################
# This function explores longitudinal data in
# long format with <age> as time scale.
explore <- function(data = NULL, id = NULL, state = NULL,
             age = NULL, digits = 3, HIST = TRUE,
             hist.col = c("green","red","blue"),
             INFO = FALSE){

# Explanation:
# <data>   = data frame with <id>, <state>, and <age>
# <id>     = identifier. Specify if no <data>
# <state>  = state variable. Specify if no <data>
# <age>    = age or transformed age. Specify if
#            no <data>
# <digits> = number of decimal digits in ouput
# <HIST>   = TRUE for histogram. FALSE otherwise
# <hist.col> = three colours for three histograms
# <INFO>   = TRUE for return object with <id> info on
#            length time intervals


##########################
# Input checks:
if(!is.null(data)){
  if(!is.data.frame(data)){
   stop("\nData <data> should be a data frame.\n\n")
  }
  if(!("id"%in%names(data))){
    stop("\nData should contain variable <id>.\n\n")
  }
  if(!("state"%in%names(data))){
    stop("\nData should contain variable <state>.\n\n")
  }
  if(!("age"%in%names(data))){
    stop("\nData should contain variable <age>.\n\n")
  }
}
if(is.null(data)){
  if(is.null(id)){
   stop("\nSpecify <id> or used data with <id>.\n\n")
  }
  if(is.null(state)){
    stop("\nSpecify <state> or used data with <state>.\n\n")
  }
  if(is.null(age)){
    stop("\nSpecify <age> or used data with <age>.\n\n")
  }
}

###########################################
# Define data with <id> etc. if needed:
if(is.null(data)){
  data <- as.data.frame(cbind(id=id, state=state, age=age))
}

#########################
# Additional definitions:

# Function for state table (code taken from 'msm'):
state.tab <- function(data){
  n <- length(data$state)
  subject <- match(data$id, unique(data$id))
  prevsubj <- c(NA, subject[1:(n - 1)])
  previous <- c(NA, data$state[1:(n - 1)])
  previous[prevsubj != subject] <- NA
  ntrans <- table(previous, data$state)
  names(dimnames(ntrans)) <- c("From state", "To state")
  ntrans
}

# Individuals:
subjects <- unique(data$id)
N <- length(subjects)

# Add baseline indicator:
data$bsline <- 0
for(i in 1:N){
  select <- which(data$id==subjects[i])
  data$bsline[select[1]] <- 1
}

##########################
# Basic stats:

# Sample size:
cat("\n---------------------------------\n")
cat("Size of the data:\n")
cat("\nNumber of individuals =",N,"\n")
cat("Number of records =",nrow(data),"\n")
cat("Number of variables =",ncol(data),"\n")
cat("-----------------------------------\n")

# Frequencies number of observation per individual:
cat("\n-----------------------------------------\n")
cat("Observed states during follow-up:\n")
freq <- t(c(table(table(data$id))))
labs <- as.numeric(names(table(table(data$id))))
tab <- as.data.frame(matrix(c(labs,freq),2,length(labs),byrow = TRUE),
                     row.names=c("Number of records per individual:",
                                 "Corresponding frequencies:"))
names(tab) <- NULL
print(tab)

# Table states:
cat("\nFrequencies of observed states:")
print(table(data$state))

# State table:
cat("\nState table:\n")
print(state.tab(data))
cat("------------------------------------------\n")

# Age:
cat("\n-------------------------------\n")
cat("Distribution of age and age intervals:\n")
minage <- min(data$age)
maxage <- max(data$age)
cat("\nMinimum age = ",round(minage,digits),"\n")
cat("Maximum age = ",round(maxage,digits),"\n")

# Mean length of age intervals:
len <- rep(NA,nrow(data))
ids <- rep(NA,nrow(data))
ind <- 1
for(i in 2:nrow(data)){
  if(data$id[i]==data$id[i-1]){
    int <- data$age[i]-data$age[i-1]
    #if(int>0){
      len[ind] <- int
      ids[ind] <- data$id[i]
      ind <- ind+1
    #  }
  }
}
len <- len[!is.na(len)]
ids <- ids[!is.na(ids)]
cat("\nMin length of intervals = ",round(min(len),digits),"\n")
cat("Median length of intervals = ",round(median(len),digits),"\n")
cat("Max length of intervals = ",round(max(len),digits),"\n")
if(HIST){
  cat("\n(See also histograms for age,\n")
  cat("and length of intervals)\n")
}
cat("-------------------------------\n")

##########################
# Basic histograms (if asked for):
if(HIST){
  opar <- par(mfrow=c(1,3), mex=0.8,mar=c(5,5,2,1)+.1)
  hist(data$age[data$bsline==1], col=hist.col[1],
     xlab="Age at baseline", main = "")
  hist(data$age, col=hist.col[2],
     xlab="Age during follow-up", main ="")
  hist(len, col=hist.col[3],
     xlab="Length of intervals between observations", main ="")
  # Tidy up:
  opar <- par(mfrow=c(1,1))
}


##########################
# Return info (if asked for)
if(INFO){
  intervals <- as.data.frame(cbind(interval.length=len,id=ids) )
  return(intervals)
}
}
