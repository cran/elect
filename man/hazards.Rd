\name{hazards}
\alias{hazards}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Age-dependent hazards based on a 'msm' model
}
\description{
Graph with transition-specific hazard functions derived from an age-dependent model fitted using 'msm'
}
\usage{
hazards(x, b.covariates, no.years, trans = NULL,
        max.haz = .5, min.haz = 0, CI = FALSE, col = NULL,
        lty = NULL, lwd = NULL, LEGEND = TRUE,
        location = "topleft", age.shift = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
Fitted \code{msm} model with age as the Gompertz time scale
}
\item{b.covariates}{
List with specified covariates values for the prediction (ignore intercept)
}
\item{no.years}{
 Number of years for the prediction
 }
\item{trans}{
 Matrix with rows (r,s) for hazard of going from state r to state s. Default to all (r,s)-hazards that are modelled in \code{x}
 }
\item{max.haz}{
 Upperbound hazard-axis
 }
\item{min.haz}{
 Lowerbound hazard-axis
}
\item{CI}{
 TRUE for plotting 95\% confidence bands. FALSE otherwise
}
\item{col}{
 Colour for each hazard curve. Example for two curves: \code{col = c("red","green")}
}
\item{lty}{
 Line type for each hazard curve. Example for two curves:  \code{lty = c(1,2)}
}
\item{lwd}{
 Width of line for each hazard curve. Example for two curves:  \code{lwd = c(2,2)}
}
\item{LEGEND}{
 TRUE for adding a legend. FALSE otherwise
}
\item{location}{
 Location for legend. Default to \code{"topleft"}. See help file for \code{legend} for furhter details
}
\item{age.shift}{
 Value to shift the age scale in the graph. Useful when age in the model is on a shifted scale. Default to 0.
}
}
\details{
This function is an add-on to the functionality in the 'msm' package. A Gompertz model with age as the time scale can be fitted in 'msm' piecewise-constantly by adding age as a covariate. The function \code{qmatrix.msm} is used by \code{hazards} repeatedly to take into account the age dependence when calculating the hazards.
}
%\value{
%\item{intervals }{
%  ...
%}
%}
\references{
  Jackson, C.H. (2011). Multi-State Models for Panel Data: The msm Package
  for R., Journal of Statistical Software, 38(8), 1-29.
}
\author{
Ardo van den Hout. With thanks to Ying Lou.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{
%%....} ~

%\seealso{
%\code{\link{elect}}
%}
\examples{
# Fit msm model:
Q     <- rbind(c(0,0.01,0.01), c(0,0,0.01), c(0,0,0))
model <- msm(state~age, subject = id, data = electData,
             center = FALSE, qmatrix = Q, deathexact = TRUE,
             covariates = ~age+x)
# Plot hazards on original age scale:
hazards(model, b.covariates = list(age = 0, x = 0),
        no.years = 20, max.haz = 0.25,  age.shift = -70)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
