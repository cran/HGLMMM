\name{salamander}
\Rdversion{1.1}
\alias{salamander}

\docType{data}

\title{
Salamander Data
}
\description{
Experiment on salamander mating. Three experiments were conducted: two were done with the same salamanders in the summer and autumn and another one in the autumn of the same year
using different salamanders.
}
\usage{data(salamander)}
\format{
  A data frame with 360 observations on the following 8 variables.
  \describe{
    \item{\code{Cross}}{Interaction term of \code{TypeF} and \code{typeM}: \code{RR} \code{RW} \code{WR} \code{WW}.}
    \item{\code{Experiment}}{Number of the experiment.}
    \item{\code{Female}}{Number of the female salamander.}
    \item{\code{Male}}{Number of the male salamander.}
    \item{\code{Mate}}{Response variable indicating success of mating.}
    \item{\code{Season}}{a factor with levels \code{Fall} \code{Summer}.}
    \item{\code{TypeF}}{Type of the female salamander \code{R} (Roughbutt) \code{W} (Whiteside).}
    \item{\code{TypeM}}{Type of the male salamander \code{R} (Roughbutt) \code{W} (Whiteside).}
  }
}


\source{
McCullagh, P. and Nelder, J.A. (1989). \emph{Generalized Linear Models}, 2nd edition. London: Chapman & Hall.
}

\references{
Lee, Y., Nelder, J.A., and Pawitan, Y. (2006). \emph{Generalized Linear Models with Random Effects}. Boca Raton: Chapman & Hall/CRC.
}

\examples{
data(salamander)

}
\keyword{datasets}
