\name{gPd.fit}
\alias{gPd.fit}
\title{ Fitting the generalized Pareto distribution to data}
\description{
This function fits a generalized Pareto distribution (gPd) to a data set using either the asymptotic maximum likelihood method (amle)  or the combined  method proposed by Villasenor-Alva and Gonzalez-Estrada (2009).
}
\usage{
gPd.fit(x,method)
}
\arguments{
  \item{x}{ numeric data vector containing a random sample from a distribution function with support on the positive real numbers.}
  \item{method}{ a character string giving the name of the parameter estimation method to be used. There are two available methods:  \code{"combined"} and \code{"amle"}. Use \code{"combined"} for fitting a gPd with shape parameter <0. Use \code{"amle"}
                for fitting a gPd with shape parameter >= 0.  }
}
\details{
The distribution function  of the gPd is given in the details section  of the function \code{\link{gPd.test}}.
}
\value{
  The parameter estimates.
}
\references{ Villasenor-Alva, J.A. and Gonzalez-Estrada, E. (2009). A bootstrap goodness of fit test for the generalized Pareto distribution. \emph{Computational Statistics and Data Analysis},\bold{53},11,3835-3841.  }
\author{ Elizabeth Gonzalez Estrada, Jose A. Villasenor Alva }
\seealso{ \code{\link{gPd.test}} for testing the gPd hypothesis, \code{\link{rgPd}} for generating gPd random numbers.}
\examples{
x <-  rgPd(20,shape = 1)   	## Random sample of size 20
gPd.fit(x,"amle")           ## Fitting a gPd to  x  using the "amle" method
}
\keyword{ htest }
\concept{ fitting the generalized Pareto distribution to data }
