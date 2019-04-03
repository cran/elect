\name{check.RestrAndConst}
\alias{check.RestrAndConst}
\alias{check.RestrAndConst.elect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checking the definition of \code{RestrAndConst} in \code{elect}
}
\description{
Function to check the definition of \code{RestrAndConst} in the \code{elect} call.
}
\usage{
check.RestrAndConst(x, RestrAndConst, PRINT = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 Fitted \code{msm} model
}
  \item{RestrAndConst}{
Vector which indexes the independent parameters in \code{model\$opt\$par} w.r.t. to the model parameters
}
  \item{PRINT}{
TRUE for printing the comparison
}
}
\details{
Function to help defining the parameter constraints in the \code{elect} call\\
when the \code{constraint} option is used in \code{msm}.}
\value{
TRUE  when \code{RestrAndConst} is well-defined, \code{FALSE} otherwise.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Ardo van den Hout
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{elect}}
}
\examples{
# Fit model with constraints:
Q     <- rbind(c(0,0.01,0.01), c(0,0,0.01), c(0,0,0))
model <- msm(state~age, subject = id, data = electData,
             center = FALSE, qmatrix = Q, deathexact = TRUE,
             fixedpars = c(7,8), covariates = ~age+x)
RestrAndConst <- c(1:6,0,0,7)
check.RestrAndConst(model, RestrAndConst, PRINT= TRUE)

# Usage: elect(x = model,..., RestrAndConst = RestrAndConst)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
