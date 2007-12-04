\name{nls2}
\alias{nls2}
\title{Nonlinear Least Squares with Brute Force}
\description{
  Determine the nonlinear least-squares estimates of the
  parameters of a nonlinear model.
}
\usage{
nls2(formula, data = parent.frame(), start, control = nls.control(), 
        algorithm = c("default", "plinear", "port", "brute-force"), \dots, 
	all = FALSE)
}
\arguments{
  \item{formula}{same as \code{formula} parameter in \code{nls}.}
  \item{data}{same as \code{data} parameter in \code{nls}.}
  \item{start}{same as \code{start} parameter in \code{nls} except that
    it may alternately be a data frame in which case an optimization will
    be run with the starting value for each row successively.}
  \item{control}{same as \code{control} parameter in \code{nls}.}
  \item{algorithm}{same as \code{algorithm} parameter in \code{nls}
    with the addition of the \code{"brute-force"} method.}
  \item{\dots}{other arguments passed to \code{nls}.}
  \item{all}{if \code{all} is true then a list of \code{nls} objects is
    returned, one for each row in \code{start}; otherwise, only the one
    with least residual sum of squares is returned.}
}
\details{
  Similar to \code{nls} except that \code{start} may be a data frame
  whose rows are starting values.  If \code{all} is false (the default)
  it will run \code{nls} from the starting values defined by each
  row and then return the \code{nls} object corresponding to the one
  with least residual sum of squares.  The \code{"brute-force"} algorithm
  just returns the nls object corresponding to the starting values.

  By using \code{algorithm="brute-force"} and a data frame of starting
  values one can get a brute force solution.  This can then be fed into
  a subsequent invocation of \code{nls2} (or \code{nls}) as starting
  values.

  By using a data frame of starting values one can perform mulitple
  optimizations in one call returning the best one (or all of them
  if \code{all = TRUE}.
}
\seealso{
  \code{\link{nls}}.
}
\examples{

richness <- c(44,36,31,39,38,26,37,33,34,48,25,22,44,5,9,13,17,15,21,10,16,22,
13,20,9,15,14,21,23,23,32,29,20,26,31,4,20,25,24,32,23,33,34,23,28,30,10,29,
40,10,8,12,13,14,56,47,44,37,27,17,32,31,26,23,31,34,37,32,26,37,28,38,35,27,
34,35,32,27,22,23,13,28,13,22,45,33,46,37,21,28,38,21,18,21,18,24,18,23,22,
38,40,52,31,38,15,21)

area <- c(26.22,20.45,128.68,117.24,19.61,295.21,31.83,30.36,13.57,60.47,
205.30,40.21,7.99,1.18,5.40,13.37,4.51,36.61,7.56,10.30,7.29,9.54,6.93,12.60,
2.43,18.89,15.03,14.49,28.46,36.03,38.52,45.16,58.27,67.13,92.33,1.17,
29.52,84.38,87.57,109.08,72.28,66.15,142.27,76.41,105.76,73.47,1.71,305.75,
325.78,3.71,6.48,19.26,3.69,6.27,1689.67,95.23,13.47,8.60,96.00,436.97,
472.78,441.01,467.24,1169.11,1309.10,1905.16,135.92,438.25,526.68,88.88,31.43,
21.22,640.88,14.09,28.91,103.38,178.99,120.76,161.15,137.38,158.31,179.36,
214.36,187.05,140.92,258.42,85.86,47.70,44.09,18.04,127.84,1694.32,34.27,
75.19,54.39,79.88,63.84,82.24,88.23,202.66,148.93,641.76,20.45,145.31,
27.52,30.70)

st <- expand.grid(Const = seq(-100, 100, 25), B = seq(-100, 100, 25), 
	A = seq(-1, 1, .25))
fo <- richness ~ Const + B * (area ^ A)
# use brute force to select best row in st
mod0 <- nls2(fo, st = st, algorithm = "brute-force")
mod0
# now use mod0 as starting value for optimization
nls2(fo, st = coef(mod0)) 

}
\keyword{nonlinear}
\keyword{regression}
\keyword{models}
