\name{summary.HGLM}
\Rdversion{1.1}
\alias{summary.HGLM}

\title{
Summary Method for HGLM Objects
}

\description{
Summarizes the fit of the generalized linear model with random effects via H-likelihood.
}

\usage{
\method{summary}{HGLM}(object,..., V = FALSE)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{an object inheriting from class \code{HGLMObject}.}
  \item{V}{If V=TRUE prints summaries for random effects.}
  \item{...}{other arguments}
}

\author{Marek Molas \email{m.molas@erasmusmc.nl}}


\seealso{\code{\link{HGLMObject}},\code{\link{HGLMfit}}}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

