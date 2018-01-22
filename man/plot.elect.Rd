\name{plot.elect}
\alias{plot.elect}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting of life expectancies
}
\description{
Graphical representation by smoothed densities of the life expectancies as estimated by \code{elect}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{plot}{elect}(x, which = NULL, kernel = "gaussian", col = "red", lwd = 2, cex.lab = 1,...)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Life expectancies estimated by \code{elect}
}
  \item{which}{
 Subselection for plotting (following the order in \code{summary}). Example: \code{which = c(1,3,5)}
}
  \item{kernel}{ Character string for smoothing kernel (\code{"gaussian"}, \code{"rectangular"}, \code{"triangular"}, \code{"epanechnikov"}, \code{"biweight"}, \code{"cosine"}, or \code{"optcosine"})
}
  \item{col}{
  Colour of curves
}
  \item{lwd}{
  Line width of curves
}
  \item{cex.lab}{
Magnification to be used for axis-labels
}
  \item{...}{
  Other arguments (not yet implemented)
 }
}
\details{
Presents distributions of the estimated life expectancies derived from
 the maximum likelihood estimate of the model parameters.
The smoothing is undertaken using the \code{R} function \code{density}.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
