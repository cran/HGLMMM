\name{HGLMfit}
\Rdversion{1.1}
\alias{HGLMfit}

\title{
Hierarchical Generalized Linear Models via H-likelihood
}

\description{
This function fits generalized linear models with random effects using h-likelihood approach. 
}

\usage{
HGLMfit(DistResp = "Normal", DistRand = NULL, Link = NULL, 
            LapFix = FALSE, ODEst = NULL, ODEstVal = 0, 
	        formulaMain, formulaOD, formulaRand, DataMain, 
            DataRand, Offset = NULL, BinomialDen = NULL, 
            StartBeta = NULL, StartVs = NULL, StartRGamma = NULL, 
            INFO = TRUE, DEBUG = FALSE, na.action, contrasts = NULL, 
            CONV = 1e-04)
}

\arguments{
  \item{DistResp}{Specify which distribution the response follows: "Normal", "Binomial", "Poisson" and "Gamma"are allowed.}
  \item{DistRand}{Specify which distribution your random effects follow: "Normal", "Beta", "Gamma" and "IGamma" are allowed, 
	the argument is a vector of characters, each entry corresponds to subsequent random component.}
  \item{Link}{Specify which link is used, the following links are currently available: "Identity" for the "Normal" response; 
	"Logit" for the "Binomial" response; "Log" for "Poisson" response; "Inverse" or "Log" for the "Gamma" response.}
  \item{LapFix}{Specify whether fixed effect in the mean structure are estimated by H-likelihood (FALSE) or by Laplace approximation to the
 	marginal likelihood (TRUE).}
  \item{ODEst}{Specify whether overdispersion (residual dispersion) parameters are fixed (FALSE) or estimated (TRUE).}
  \item{ODEstVal}{Provide (starting) values for the overdispersion structure - vector of length equal to the number of parameters in the overdispersion structure.}
  \item{formulaMain}{Specify the formula for the mean structure of the model: \code{ Response} \code{~ Fixed.Effect.1+ ... +Fixed.Effect.N}
                    \code{+(Random.Effect.1|Subject.1)+ ... +(Random.Effect.N|Subject.N).}}
  \item{formulaOD}{Specify the formula for the overdispersion (residual dispersion) structure: \code{~ Effect.1+ ... +Effect.N.}}
  \item{formulaRand}{Specify the list of formulas for dispersion strucutre of each random effects:
            \code{list(one=~Effect.1.1+...,two=~Effect.2.1+ ...,three=...,...).}}
  \item{DataMain}{Specify dataframe for \code{formulaMain} and \code{formulaOD}.}
  \item{DataRand}{List of dataframes for \code{formulaRand}: \code{list(one=data.R.1,two=data.R.2,...).}}
  \item{Offset}{In case "Poisson" response is used specify the OFFSET variable, variable specifies factor \code{t} in the following: \code{log(mu/t)=eta.}}
  \item{BinomialDen}{In case "Binomial" response is used specify the number of trials for each observation.}
  \item{StartBeta}{Specify optional starting values for fixed effects in the mean structure (one vector of numeric values).}
  \item{StartVs}{Specify optional starting values for random effects in the mean structure (one vector of numeric values).}
  \item{StartRGamma}{Specify optional starting values for parameters of dispersion components of random effects (one vector of numeric values).}
  \item{INFO}{Request display of iteration information if TRUE.}
  \item{DEBUG}{Request details display of iteration mechanism progress if TRUE.}
  \item{na.action}{See \code{\link{lm}}}
  \item{contrasts}{See \code{\link{lm}}, caution as it is currently not fully developed.}
  \item{CONV}{Specify the convergence criterion, the default is 1e-4, for models with many random effects could be set less strict.}
}

%%\details{
%%}

\value{See \code{\link{HGLMObject}} for the components of the fit.}

\references{
Lee, Y. and Nelder, J.A. (1996). Hierarchical generalized linear models (with discussion). \emph{Journal of the Royal Statistical Society. Series B (Methological)}
\bold{58}, 619-678.

Lee, Y. and Nelder, J.A. (2001). Hierarchical generalised linear models: A synthesis of generalised linear models, random-effect models and structured dispersions.
\emph{Biometrika} \bold{88}, 987-1006. 

Lee, Y., Nelder, J.A., and Pawitan, Y. (2006). \emph{Generalized Linear Models with Random Effects}. Boca Raton: Chapman & Hall/CRC.

Noh, M. and Lee, Y. (2007). REML estimation for binary data in GLMMs. \emph{Journal of Multivariate Analysis} \bold{98}, 896-915.
}

\author{
Marek Molas \email{m.molas@erasmusmc.nl}.
}

\note{
	Part of the code of the user interface function \code{HGLMfit} (and functions it invokes) is based on package \code{lme4} by Douglas Bates and Martin Maechler
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

\examples{
#################################
######## Salamander Data ########
#################################

# Bernoulli model with crossed-random effects and estimated overdispersion #
RSal<-data.frame(int=rep(1,60))

modBin<-HGLMfit(DistResp="Binomial",DistRand=c("Normal","Normal"),
            Link="Logit",LapFix=TRUE,ODEst=FALSE, ODEstVal=c(0),
            formulaMain=Mate~TypeF+TypeM+TypeF*TypeM+(1|Female)+(1|Male), 
            formulaOD=~1,formulaRand=list(one=~1,two=~1),
            DataMain=salamander,DataRand=list(RSal,RSal),
            Offset=NULL,BinomialDen=rep(1,360),INFO=TRUE,DEBUG=FALSE) 
modBin
summary(modBin)

###########################
######## Cake Data ########
###########################

# Linear mixed model with nested random effects #
cake$repbatch<-100*cake$Replicate+cake$Batch
R1Cake<-data.frame(int=rep(1,15))
R2Cake<-data.frame(int=rep(1,45))
modCake1<-HGLMfit(DistResp="Normal",DistRand=c("Normal","Normal"),
            Link="Identity",LapFix=FALSE,ODEst=TRUE,ODEstVal=c(0),
            formulaMain=Angle~as.factor(Recipe)+as.factor(Temperature)
            +as.factor(Recipe)*as.factor(Temperature)+
            (1|Replicate)+(1|repbatch),formulaOD=~1,
            formulaRand=list(one=~1,two=~1),DataMain=cake,
            DataRand=list(R1Cake,R2Cake),Offset=NULL,BinomialDen=rep(1,360),
            ,INFO=TRUE,DEBUG=FALSE)     
modCake1
summary(modCake1)
HGLMLikeDeriv(modCake1)
modCake2<-HGLMfit(DistResp="Normal",DistRand=c("Normal","Normal"),
            Link="Identity",LapFix=FALSE,ODEst=TRUE,ODEstVal=c(0),
            formulaMain=Angle~as.factor(Recipe)+as.factor(Temperature)+
            (1|Replicate)+(1|repbatch),formulaOD=~1,
            formulaRand=list(one=~1,two=~1),DataMain=cake,
            DataRand=list(R1Cake,R2Cake),Offset=NULL,BinomialDen=rep(1,360),
            ,INFO=TRUE,DEBUG=FALSE)  
HGLMLikeDeriv(modCake2)
HGLMLRTest(modCake2,modCake1)

##########################
######## Rat Data ########
##########################

Rrat<-data.frame(WBC=tapply(rat$WhiteBloodCells,rat$Subject,mean),
            RBC=tapply(rat$RedBloodCells,rat$Subject,mean))
            
# Poisson-Normal Model #
modRat1<-HGLMfit(DistResp="Poisson",DistRand=c("Normal"),Link="Log",
            LapFix=FALSE,ODEst=FALSE,ODEstVal=c(0),formulaMain=Y~WhiteBloodCells
            +RedBloodCells+as.factor(Drug)+(1|Subject),formulaOD=~1,
            formulaRand=list(one=~1),DataMain=rat,DataRand=list(Rrat),
            ,INFO=TRUE,DEBUG=FALSE)  
BootstrapEnvelopeHGLM(modRat1,9,9999)

# Quasi-Poisson-Normal Model #
modRat2<-HGLMfit(DistResp="Poisson",DistRand=c("Normal"),Link="Log",
            LapFix=FALSE,ODEst=TRUE,ODEstVal=c(0),formulaMain=Y~WhiteBloodCells+
            RedBloodCells+as.factor(Drug)+(1|Subject),formulaOD=~1,
            formulaRand=list(one=~WBC+I(WBC^2)),DataMain=rat,DataRand=list(Rrat),
            ,INFO=TRUE,DEBUG=FALSE)    
# Gamma-Normal Model # 
modRat3<-HGLMfit(DistResp="Gamma",DistRand=c("Normal"),Link="Log",
            LapFix=FALSE,ODEst=TRUE,ODEstVal=c(0),formulaMain=Y~WhiteBloodCells+
            RedBloodCells+as.factor(Drug)+(1|Subject),formulaOD=~1,
            formulaRand=list(one=~WBC+I(WBC^2)),DataMain=rat,DataRand=list(Rrat),
            ,INFO=TRUE,DEBUG=FALSE)                   
# Normal-Normal Model #
modRat4<-HGLMfit(DistResp="Normal",DistRand=c("Normal"),Link="Identity",
            LapFix=FALSE,ODEst=TRUE,ODEstVal=c(0),formulaMain=Y~WhiteBloodCells+
            RedBloodCells+as.factor(Drug)+(1|Subject),formulaOD=~1,
            formulaRand=list(one=~WBC+I(WBC^2)),DataMain=rat,DataRand=list(Rrat),
            ,INFO=TRUE,DEBUG=FALSE)   
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

 
