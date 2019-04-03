\name{summary.elect}
\alias{summary.elect}
\alias{summary}
\title{Summarises the estimation of the life expectancies}
\description{
Summary of estimated state-specific and marginal life expectancies
}
\usage{
\method{summary}{elect}(object, probs = c(0.025, 0.5, 0.975),
               digits = 3, StartStateTotals = FALSE,
               print = TRUE, sd.model = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Life expectancies estimated by \code{elect}
}
  \item{probs}{
Numeric vector of probabilities with values in [0,1] for quantiles
}
  \item{digits}{
 Number of decimal places in output
}
  \item{StartStateTotals}{
 TRUE for output on start-state totals e_{r.} (for \code{S}>0)
}
  \item{print}{
TRUE for printing output to screen, FALSE otherwise
}
  \item{sd.model}{
TRUE for printing a summary of the fitted model for the prevalence
}
 \item{...}{
 Other arguments (not yet implemented)
}
}
\details{
A summary for the state-specific and marginal life expectancies as derived in
\code{elect}. Quantiles are derived from simulation based on the maximum likelihood estimation.
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
\code{\link{elect}}, \code{\link{plusmin}}, \code{\link{plot.elect}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
