\name{plusmin.elect}
\alias{plusmin.elect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Computation of a function of life expectancies
}
\description{
Compute an additive function of a series of estimated life expectancies
}
\usage{
plusmin.elect(x, index = NA, func = "plus", probs = c(0.025, 0.5, 0.975), digits = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Life expectancies estimated by \code{elect}
}
  \item{index}{
Selection of the life expectancies in the function (following the order in \code{summary})
}
  \item{func}{
Required series of "plus" and "minus". Example: \code{func = c("plus","plus")}
}
  \item{probs}{
Probabilities for the quantiles
}
  \item{digits}{
 Number of digits in the output
}
}
\details{
This function requires that uncertainty of LEs is estimated; that is, \code{S} is not 0.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
# Fit msm model:
Q     <- rbind(c(0,0.01,0.01), c(0,0,0.01), c(0,0,0))
model <- msm(state~age, subject = id, data = ELECTData,
             center = FALSE, qmatrix = Q, death = TRUE,
             covariates = ~age+x)

# Estimate function of life expectancies:
sddata <- ELECTData[ELECTData$state < 3,]
LEs    <- elect(x = model, b.covariates = list(age = 0, x = 1),
                statedistdata=sddata, h = 0.5, age.max = 50, S = 25)
plusmin.elect(LEs, index = c(2,4), func = "plus")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
