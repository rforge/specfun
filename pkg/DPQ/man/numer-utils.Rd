\name{numer-utils}
\title{Numerical Utilities - Functions, Constants}
%% Numeric Constants .. ---------------------------
\alias{M.LN2}
\alias{M_cutoff}
\alias{G_half}
%% Functions ... ----------------------------------
\alias{all_mpfr}
\alias{any_mpfr}
\alias{logr}
%--------------------------------------------------
\description{
  The \pkg{DPQ} package provides some numeric constants used in some of its
  distribution computations.

  \code{all_mpfr()} and \code{any_mpfr()} return \code{\link{TRUE}}
  iff all (or \sQuote{any}, respectively) of their arguments inherit from
  \code{\link{class}} \code{"mpfr"} (from package \CRANpkg{Rmpfr}).

  \code{logr(x,a)} computes \code{\link{log}(x / (x + a))} in a numerically stable way.
}
\usage{
## Numeric Constants :
M.LN2        # = log(2) = 0.693....
M_cutoff     # = .......=  3196577161300663808 ~= 3.2e+18
G_half       # = sqrt(pi) = Gamma( 1/2 )

## Functions :
all_mpfr(...)
any_mpfr(...)
logr(x, a)    # == log(x / (x + a)) -- but numerically smart; x > 0, a >= 0 > -x
}
\arguments{
  \item{\dots}{numeric or \code{"mpfr"} numeric vectors.}
  \item{x, a}{number-like .. }
}
\value{
  The numeric constant in the first case; a numeric (or "mpfr") vector of
  appropriate size in the 2nd case.
}
\details{
  \describe{
    \item{\code{all_mpfr()},}{}
    \item{\code{all_mpfr()} :}{test if \code{\link{all}} or
      \code{\link{any}} of their arguments or of class \code{"mpfr"} (from
      package \CRANpkg{Rmpfr}).  The arguments are evaluated only until
      the result is determined, see the example.}

    \item{\code{logr()}}{computes \eqn{\log( x / (x+a) )} in a numerically
      stable way.}%% which means .... (TODO)

    %% \item{\code{..()}}{ ..}
    %% \item{\code{..()}}{ ..}
  }
}
%% \references{
%%   ...
%% }
\author{Martin Maechler}
\seealso{
  \code{\link{.Machine}}
}
\examples{
(Ms <- ls("package:DPQ", pattern = "^M"))
lapply(Ms, function(nm) { cat(nm,": "); print(get(nm)) }) -> .tmp
}
\keyword{math}