\name{ELECTData}
\alias{ELECTData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Three-state data for exploring \code{elect}
}
\description{
Longitudinal toy data for a three-state illness-death process.
}
\usage{
ELECTData
}
%- maybe also 'usage' for other objects documented here.
\format{
Data frame with 487 rows, grouped by 100 individuals. Interval-censored data for living states 1 and 2, and exact times for death state 3.
Variables: \code{id}  = {identification number},
\code{state} = {state}, \code{age} = age in years on a shifted scale,
\code{x} = binary time-independent covariate.
}
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
# Sample size:
print(length(unique(ELECTData$id)))
# Frequencies number of observation per individual:
print(table(table(ELECTData$id)))
# State table
print(statetable.msm(state, id, data=ELECTData))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
