\name{elect}
\alias{elect}
\title{
State-specific and marginal life expectancies
}
\description{
Estimation of state-specific and marginal life expectancies given
a multi-state survival model fitted using the \code{msm}-package
}
\usage{
elect(x, b.covariates, statedistdata, time.scale.msm = "years",
        h, age.max, S = 0, setseed = NULL, RestrAndConst = NULL,
        statedist.covariates = "age", method = "step")
}
\arguments{
 \item{x}{ Fitted \code{msm} model }
 \item{b.covariates}{ List with specified covariates values (ignore intercept) }
 \item{statedistdata}{ Data used to estimate prevalence distribution of living states }
 \item{time.scale.msm}{ Time scale in multi-state model: \code{"years"},
                       \code{"months"},\code{"weeks"}, or a  value  in (0,1]}
 \item{h}{ Grid parameter for integration where scale is \code{time.scale.msm}}
 \item{age.max}{Assumed maximum age in same time 	scale as in fitted model}
 \item{S}{ Number of replications for estimation of uncertainty (\code{S}=0 for no estimation)}
 \item{setseed}{Seed for the random number generation in the simulation}
 \item{RestrAndConst}{Vector which indexes the independent model parameters. Only
              needed when \code{constraint} is used in \code{msm} call}
 \item{statedist.covariates}{Names of covariates for model for prevalence distribution of living states}
 \item{method}{Approximation of integral: \code{"step"} for simple step function,   \code{"MiddleRiemann"} or \code{"Simpson"}}
}

\details{
The \code{elect}-package is an add-on to the \code{msm}-package. In the \code{msm} call for fitting the model use \code{center=FALSE} and \code{death=TRUE}, and
names \code{state} and \code{age}. Do not use variables encoded as factor by \code{R}. Covariate \code{age} should be the first entry in \code{ b.covariates}. The other covariates in \code{b.covariates} should follow the order
 in the \code{msm} call. The life expectancies are computed by approximating the
 integral numerically with a grid defined by \code{h}. The specification of \code{statedist.covariates} should
 be a subset of  \code{b.covariates}.}

\value{A list containing the following components:
 \item{pnt }{Life expectancies derived from MLE of model parameters}
 \item{sim }{Simulated life expectancies using the MLE of model parameters}
 \item{h}{As specified in \code{elect} call}
 \item{covars}{ Covariates as specified in \code{elect} call}
 \item{S}{ \code{S} as specified in \code{elect} call}
 \item{sd.model}{Fitted model for the prevalence distribution of living states}
}
\references{
  Jackson, C.H. (2011). Multi-State Models for Panel Data: The msm Package
  for R., Journal of Statistical Software, 38(8), 1-29.

  Van den Hout, A. (2017). Multi-State Survival Models for Interval-Censored Data. Boca Raton: CRC/Chapman & Hall.
}
\author{
Ardo van den Hout and Mei Sum Chan
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{summary.elect}}, \code{\link{plot.elect}}
}
\examples{
# Fit msm model:
Q     <- rbind(c(0,0.01,0.01), c(0,0,0.01), c(0,0,0))
model <- msm(state~age, subject = id, data = ELECTData,
             center = FALSE, qmatrix = Q, death = TRUE,
             covariates = ~age+x)

# Estimate life expectancies:
sddata <- ELECTData[ELECTData$bsline == 1,]
LEs    <- elect(x = model, b.covariates = list(age = 0, x = 0),
                statedistdata = sddata, h = 0.5, age.max = 50, S = 25)
summary(LEs)
plot(LEs)     # For smoother graphs, increase S
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
