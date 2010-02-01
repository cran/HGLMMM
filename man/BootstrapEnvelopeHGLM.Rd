\name{BootstrapEnvelopeHGLM}
\Rdversion{1.1}
\alias{BootstrapEnvelopeHGLM}

\title{
QQ-plot for deviance residuals of the response with bootstrap envelopes
}

\description{
This function creates qq-plot of the deviance residuals for the response together with bootstraped envelopes.
}

\usage{
BootstrapEnvelopeHGLM(model, numBoot, seed)
}

\arguments{
  \item{model}{
    Specify a HGLM model for which qq-plot of deviance residuals together with bootstrap envelopes needs to be created.
}
  \item{numBoot}{
    Number of bootstraps. For 95\% envelopes specify 19. 
}
  \item{seed}{
    Random seed to use for sampling. Argument to the function ``set.seed''.
}
}
\details{
    Under the assumptions of the correctness of the final model, data is generated and the model is refitted. Deviance residuals for the
    response are obtained and qq-plot is computed. Function plots the qq-plot for the original model, together with samples envelopes 
    (maximum and minimum of qq-plots of all the generated qq-plots).
}

\author{
Marek Molas \email{m.molas@erasmusmc.nl}.
}


