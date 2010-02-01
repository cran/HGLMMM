\name{cake}
\Rdversion{1.1}
\alias{cake}

\docType{data}

\title{
Cake Data
}

\description{
In an experiment on the preparation of chocolate cakes, conducted at Iowa State College, 3 reciepes for preparing the batter
were compared.
}

\usage{data(cake)}

\format{
  A data frame with 270 observations on the following 5 variables.
  \describe{
    \item{\code{Angle}}{Breaking angle.}
    \item{\code{Batch}}{For each mix, enough batter was prepared for 6 cakes, each of which was baked at different temperature.}
    \item{\code{Recipe}}{Reciepe 1 and 2 differed in that the chocolate was added at 40 degrees and 60 degrees respectively,
			  while reciepe 3 contained extra sugar.}
    \item{\code{Replicate}}{There were 15 replications for each reciepe and temperature.}
    \item{\code{Temperature}}{Temperature at which cakes were baked.}
  }
}



\source{
Cochran, W.G. and Cox, G.M. (1957). \emph{Experimental Designs}, 2nd edition. New York: John Wiley and Sons. 
}

\references{
Lee, Y., Nelder, J.A., and Pawitan, Y. (2006). \emph{Generalized Linear Models with Random Effects}. Boca Raton: Chapman & Hall/CRC.
}

\examples{
data(cake)
}

\keyword{datasets}
