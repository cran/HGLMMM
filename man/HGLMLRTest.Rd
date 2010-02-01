\name{HGLMLRTest}
\Rdversion{1.1}
\alias{HGLMLRTest}


\title{
 Likelihood Ratio Tests of Two Nested HGLM Models
}

\description{
 Compares two HGLM models using profile likelihoods and h-likelihood. First profile likelihood is a Laplace approximation to the marginal likelihood, while second 
 profile likelihood is REML type likelihood. Also values of h-likelihood are compared. 
}

\usage{
HGLMLRTest(x1, x2)
}

\arguments{
  \item{x1}{an object inheriting from class \code{HGLMObject}.}
  \item{x2}{an object inheriting from class \code{HGLMObject}.}
}

\author{Marek Molas <m.molas@erasmusmc.nl>}

