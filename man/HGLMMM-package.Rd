\name{HGLMMM-package}
\Rdversion{1.1}
\alias{HGLMMM-package}
\alias{HGLMMM}
\docType{package}

\title{
Hierarchical Generalized Linear Models (H-likelihood)
}

\description{
This package fits generalized linear models with random effects via H-likelihood. Complex multilevel designs or multi-membership designs
can be handled. Modelling of dispersion components as a function of covariates is allowed.
}

\details{
\tabular{ll}{
Package: \tab HGLMMM\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2009-10-08\cr
License: \tab \acronym{GPL} Version >=2 \cr
LazyLoad: \tab yes\cr
}
The package main fitting function is \code{\link{HGLMfit}}. It allowes the distribution of the data to follow one of the four distributions from 
exponential family: Normal, Binomial, Poisson and Gamma. Multiple independent random effects can be specified in the model, which are allowed to 
follow Normal, Beta, Gamma or Inverse-Gamma random effects. Overdispersion (residual dispersion) parameter can be fixed or estimated and is allowed 
to depend on covariates. Similarly a structure can be specified for each dispersion component of random effects. The estimation of the dispersion
parameters is done using adjusted profile likelihood corresponding to the restricted likelihood. The mean structure parameters are estimated either
using H-likelihood or first order Laplace approximation to the marginal likelihood. 
}

\author{
Marek Molas

Maintainer: Marek Molas \email{m.molas@erasmusmc.nl}
}
\note{
	Part of the code of the user interface function \code{HGLMfit} (and functions it invokes) is based on package \code{lme4} by Douglas Bates and Martin Maechler
}

\references{
Lee, Y. and Nelder, J.A. (1996). Hierarchical generalized linear models (with discussion). \emph{Journal of the Royal Statistical Society. Series B (Methological)}
\bold{58}, 619-678.

Lee, Y. and Nelder, J.A. (2001). Hierarchical generalised linear models: A synthesis of generalised linear models, random-effect models and structured dispersions.
\emph{Biometrika} \bold{88}, 987-1006. 

Lee, Y., Nelder, J.A., and Pawitan, Y. (2006). \emph{Generalized Linear Models with Random Effects}. Boca Raton: Chapman & Hall/CRC.

Noh, M. and Lee, Y. (2007). REML estimation for binary data in GLMMs. \emph{Journal of Multivariate Analysis} \bold{98}, 896-915.

}
\keyword{ package }


\seealso{\code{\link{HGLMfit}}}

