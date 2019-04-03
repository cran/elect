\name{explore}
\alias{explore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data statistics for an age-dependent model
}
\description{
Data statistics that are important for fitting an age-dependent multi-state model. Basic summaries and additional information on the age distribution in the data
}
\usage{
explore(data = NULL, id = NULL, state = NULL, age = NULL,
        digits = 3, HIST = TRUE, hist.col = c("green","red",
        "blue"), INFO = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
Data frame with variables \code{id}, \code{state}, and \code{age}
}
\item{id}{
Identifier. Specify if \code{data} is not provided
}
 \item{state}{
State variable. Specify if \code{data} is not provided
}
\item{age}{
Age or transformed age. Specify if \code{data} is not provided
}
\item{digits}{
 Number of digits in the output
}
\item{HIST}{
 TRUE for histograms of the age distribution. FALSE otherwise
}
\item{hist.col}{
 Colours for the three histograms of the age distribution
}
\item{INFO}{
 TRUE for returning a list which links \code{id} with time intervals between the subsequent records. FALSE otherwise
}
}
%%\details{
%%This function requires that uncertainty of LEs is estimated; %%that is, \code{S} is not 0.
%%}
\value{
\item{intervals }{
  Data frame with variables \code{interval.length} and corresponding \code{id} (when \code{INFO = TRUE})
}
}
%\references{
%% ~put references to the literature/web site here ~
%}
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
explore(electData)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
