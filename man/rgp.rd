\name{rgp}
\alias{rgp}
\title{ Generalized Pareto random numbers}
\description{
This function generates pseudo random numbers from a generalized Pareto distribution (gPd).
}
\usage{
rgp(n,shape,scale)
}
\arguments{
  \item{n}{ sample size.}
  \item{shape}{ shape parameter.}
  \item{scale}{ scale parameter. Default \code{scale=1}.}
}
\details{
 The distribution function  of the gPd with \code{shape} and  \code{scale}  parameters  \eqn{\gamma}{gamma} and \eqn{\sigma}{sigma} is

  \deqn{F(x) = 1 - \left[ 1 + \frac{\gamma x}{ \sigma } \right] ^ { - 1 /\gamma}}{ F(x) = 1 - [ 1 + gamma x  / sigma ]^(-1/gamma)}

  where   \eqn{\gamma}{gamma}  is a real number, \eqn{\sigma > 0}{sigma > 0} and \eqn{1 + \gamma  x  / \sigma > 0}{1 + gamma x / sigma > 0}. When \eqn{\gamma = 0}{gamma =
    0}, we have the exponential distribution with \code{scale} parameter \eqn{\sigma}{sigma}.
}
\value{
  A vector of length \code{n}.
}
\author{ Elizabeth Gonzalez Estrada, Jose A. Villasenor Alva }
\seealso{  \code{\link{gpd.test} for testing the gPd hypothesis}}
\examples{
rgp(30,shape=1.5)    ## Generates 30 random numbers from a gPd with shape parameter 1.5.
}
\keyword{ distribution }
