\name{HGLMObject}

\alias{HGLMObject}

\title{Fitted HGLM Object}

\description{
An object returned by the \code{HGLMfit} function, inheriting from class \code{HGLM} and representing a fitted
HGLM model. Objects of this class have methods for the generic functions 
\code{print} and \code{summary}.
}

\value{
  The \code{HGLMObject} consists of the following components: 
  \item{Results}{a list with the results of the fit:
    \describe{
        \item{\code{Beta}}{the vector of fixed effects in the mean structure.}
        \item{\code{Vs}}{the vector of random effects in the mean structure.}
	\item{\code{Dispersion}}{the vecotr of fixed effects in the dispersion structure of random effects.}
	\item{\code{ResidualDispersion}}{the vector or fixed effects in the overdispersion (residual dispersion) structure.}
	\item{\code{GradientDisp}}{the gradient of dispersion parameters}
	\item{\code{HessianDisp}}{the hessian of dispersion parameters}
	\item{\code{SEDisp}}{standard errors of dispersion parameters}
	\item{\code{GradientFix}}{the gradient of fixed effects in the mean structure}
	\item{\code{HessianFix}}{the hessian of fixed effects in the mean structure}
	\item{\code{SEFix}}{standard errors of fixed effects in the mean structure}
	\item{\code{GradientResDisp}}{the gradient of overdispersion (residual dispersion) parameters}
	\item{\code{HessianResDisp}}{the hessian of overdispersion (residual dispersion) parameters}
	\item{\code{SEResDisp}}{standard errors of overdispersion (residual dispersion) parameters}
	\item{\code{HLikelihood}}{the value of the h-likelihood at maximum}
	\item{\code{APFix}}{the value of the adjusted profile likelihood (Laplace approximation to the marginal likelihood) at maximum}
	\item{\code{APDis}}{the value of the adjusted profile likelihood (REML likelihood) at maximum}
	\item{\code{CLikelihood}}{the value of the conditional part of the h-likelihood \code{f(y|v)}}
	\item{\code{SEVs}}{standard errors of random effects in the mean structure}
    }
  }
  \item{Details}{a list with the detailed result of the fit:
    \describe{
      	\item{\code{Beta}}{the vector of fixed effects in the mean structure.}
	   \item{\code{Vs}}{the vector of random effects in the mean structure.}
	   \item{\code{YDispParms}}{the vecotr of fixed effects in the dispersion structure of random effects.}
	   \item{\code{RDispParms}}{the vector or fixed effects in the overdispersion (residual dispersion) structure.}
	   \item{\code{DevianceResidualY}}{a matrix of deviance residuals for the response}
	   \item{\code{DevianceResidualR}}{list of vectors of the deviance residuals from the random components, each component gives a vector of residuals.}
	   \item{\code{DevianceResidualODisp}}{a matrix of the deviance residuals from the gamma model of the overdispersion structre}
	   \item{\code{DevianceResidualDisp}}{a list of deviance residuals from the gamma models for the dispersion paramters, each random component yields a vector of residuals}
	   \item{\code{StdDevianceResidualY}}{a matrix of standardized \code{DevianceResidualY}.}
	   \item{\code{StdDevianceResidualR}}{a list of vectors of standardized \code{DevianceResidualR}.}
	   \item{\code{StdDevianceResidualODisp}}{a matrix of standardized \code{DevianceResidualODisp}.}
	   \item{\code{StdDevianceResidualDisp}}{a list of vectors of standardized \code{DevianceResidualDisp}.}
    }
  }
  \item{NAMES}{a list of names}
  \item{CALL}{original CALL used to obtain this object}
}  

\author{Marek Molas \email{m.molas@erasmusmc.nl}}


\seealso{
    \code{\link{HGLMfit}}
}

\keyword{regression}
