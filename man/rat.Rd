\name{rat}
\Rdversion{1.1}
\alias{rat}

\docType{data}

\title{
Rat Data
}

\description{
A comparison of three chemotherapy drugs applied to 30 rats that had an induced leukemic condition.
}

\usage{data(rat)}
\format{
  A data frame with 120 observations on the following 5 variables.
  \describe{
    \item{\code{Drug}}{Between-rat covariate: drug applied.}
    \item{\code{RedBloodCells}}{Within-rat covariate: red blood cell count.}
    \item{\code{Subject}}{Number of a rat.}
    \item{\code{WhiteBloodCells}}{Within-rat covariate: white red blood cell count.}
    \item{\code{Y}}{Number of cancer cell colonies.}
  }
}



\source{
Myers, P.H., Montgomery, D.C., and Vining, G.G. (2002). \emph{Generalized Linear Models with Applications in Engineering and the Sciences}. New York: John Wiley & Sons
}

\references{
Lee, Y., Nelder, J.A., and Pawitan, Y. (2006). \emph{Generalized Linear Models with Random Effects}. Boca Raton: Chapman & Hall/CRC.
}

\examples{
data(rat)

}
\keyword{datasets}
