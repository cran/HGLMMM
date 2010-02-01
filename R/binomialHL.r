

################################
##### Copied from BinHlike #####
################################


HLikeBinom<-function(Y=NULL,B=NULL,X=NULL,ZZ=NULL,SS=NULL,Beta=NULL,Vstart=NULL,EstimateOverDisp=FALSE,LaplaceFixed=FALSE,OFFSET=NULL,Link="Logit",
                        RandDist=c("Normal"),DDR=NULL,DDY=NULL,DYgamma=NULL,DRgamma=NULL,Info=FALSE,DEBUG=FALSE){
  
    # Creator of U the inverse of the link for V #
    LinkR<-function(x,RandDist){
        if (RandDist=="Normal") out<-x
        if (RandDist=="Gamma") out<-exp(x)
        if (RandDist=="IGamma") out<-(-1/x)
        if (RandDist=="Beta") out<-exp(x)/(1+exp(x))
        out
    }
    
    # Random effects W vector creator - this takes as an argument u vector#
    WRVC<-function(x,RandDist){
        if (RandDist=="Normal") out<-rep(1,length(x))
        if (RandDist=="Gamma") out<-x
        if (RandDist=="IGamma") out<-x^2
        if (RandDist=="Beta") out<-x*(1-x)
        out
    }
    # x- vscale y- uscale - computes deviances for the estimation of the lambda paramters #
    DevRand<-function(x,y,RandDist){
        if (RandDist=="Normal") out<-y^2
        if (RandDist=="Gamma") out<-2*(y-x-1)
        if (RandDist=="IGamma") out<-2*(log(y)-x-1)
        if (RandDist=="Beta") out<--log(4*y*(1-y))
        out
    }
    DWRDU<-function(x,RandDist){
        if (RandDist=="Normal") out<-rep(0,length(x))
        if (RandDist=="Gamma") out<-rep(1,length(x))
        if (RandDist=="IGamma") out<-2*x
        if (RandDist=="Beta") out<-1-2*x
        out
    }
    # link of the main distribution part - choice between canonical inverse and logarithm #
    LinkY<-function(mu,Link){
        if (Link=="Inverse")    eta<--(1/mu)
        if (Link=="Log")        eta<-log(mu)
        if (Link=="Identity")   eta<-mu
        if (Link=="Logit")      eta<-log(mu/(1-mu))
        if (Link=="Probit")     eta<-qnorm(mu)
        if (Link=="CLogLog")    eta<-log(-log(1-mu))
        eta
    }
    # Inverse of the link #
    InvLinkY<-function(eta,Link){
        if (Link=="Inverse")    mu<--(1/eta)
        if (Link=="Log")        mu<-exp(eta)
        if (Link=="Identity")   mu<-eta
        if (Link=="Logit")      mu<-exp(eta)/(1+exp(eta))
        if (Link=="Probit")     mu<-pnorm(eta)
        if (Link=="CLogLog")    mu<-1-exp(-exp(eta))
        mu
    }
    # Generation of the weight matrix W#
    Wmatgen<-function(mu,Link){
        if (Link=="Inverse")    Wvec<-((B-mu)*(B/mu^3))^(-1)
        if (Link=="Log")        Wvec<-((B-mu)*(1/(mu*B)))^(-1)
        if (Link=="Identity")   Wvec<-((B-mu)*(mu/B^3))^(-1)
        if (Link=="Logit")      Wvec<-(B-mu)*(mu/B)
        if (Link=="Probit")     Wvec<-((B-mu)*(mu/B)*(1/((B^2)*dnorm(qnorm(mu/B))^2)))^(-1)
        if (Link=="CLogLog")    Wvec<-(B-mu)*(B/mu)*(log(1-(mu/B))^2)
        Wvec
    }
    # Generation of the adjustment to the weight matrix - factor (y-mu)/phi is left outside #
    #Amatgen<-function(mu,Link){
    #    if (Link=="Inverse")    Avec<-rep(0,length(mu))
    #    if (Link=="Log")        Avec<--(1/mu)
    #    Avec
    #}
    # Generation of the derivative dWdmu #
    dWdmugen<-function(mu,Link){
        if (Link=="Inverse")    dWdmu<-(3*mu^2*(B-mu)*B+mu^3*B)/(B*(B-mu))^2
        if (Link=="Log")        dWdmu<-(B^2)/(B-mu)^2
        if (Link=="Identity")   dWdmu<-(-(B^3)*(B-2*mu))/((B-mu)*mu)^2
        if (Link=="Logit")      dWdmu<-(1-2*(mu/B))
        if (Link=="Probit")     dWdmu<-(-(2*B^2)*dnorm(qnorm(mu/B))*qnorm(mu/B))/((B-mu)*mu)-((B^3)*(dnorm(qnorm(mu/B))^2)*(B-2*mu))/((B-mu)*mu)^2
        if (Link=="CLogLog")    dWdmu<-(-(B/mu)*log(1-(mu/B))^2)-((B-mu)*B*(1/mu^2)*log(1-(mu/B))^2)-((B/mu)*2*log(1-(mu/B)))
        dWdmu
    }
    # Generation of the derivative dmudeta #
    dmudetagen<-function(mu,Link){
        if (Link=="Inverse")    dmudeta<-(mu^2)/B
        if (Link=="Log")        dmudeta<-mu
        if (Link=="Identity")   dmudeta<-B
        if (Link=="Logit")      dmudeta<-(mu/B)*(B-mu)
        if (Link=="Probit")     dmudeta<-B*dnorm(qnorm(mu/B))
        if (Link=="CLogLog")    dmudeta<--B*log(1-(mu/B))*(1-(mu/B))
        dmudeta
    }
    # Generation of the derivative dAdmu (y-mu)/Phi is outside#
    #dAdmugen<-function(mu,Link){
    #    if (Link=="Inverse")    dAdmu<-rep(0,length(mu))
    #    if (Link=="Log")        dAdmu<-(1/mu^2)
    #    dAdmu
    #}
    
    n<-nrow(Y)
    p<-ncol(X)
    nrand<-length(ZZ)
    q<-rep(0,nrand)
    if (length(RandDist)!=nrand) stop("Random effects distribution not properly defined")
    
    # Creating the design matrix for random effects used in the program #
    for (i in 1:nrand){
        if (i==1) { Z<-ZZ[[1]]
                    q[i]<-dim(ZZ[[1]])[2]}
        else {      Z<-cbind(Z,ZZ[[i]])
                    q[i]<-dim(ZZ[[i]])[2]}
    }
    
    # Defining cumulative dimensions #
    qcum<-cumsum(c(0,q))
    if (is.null(Vstart)) Vstart<-rep(0,sum(q))
    V<-list(0)
    U<-list(0)
    TT<-cbind(X,Z)
    PsiM<-rep(0,sum(q))
    for (i in 1:nrand){
        TT<-rbind(TT,cbind(matrix(0,q[i],p+qcum[i]),diag(q[i]),matrix(0,q[i],qcum[nrand+1]-qcum[i+1])))
        V[[i]]<-as.matrix(Vstart[(qcum[i]+1):(qcum[i+1])])
        if (i==1) VT<-V[[1]]
        else VT<-c(VT,list(V[[i]]))
        if (RandDist[i]=="Normal") PsiM[(qcum[i]+1):qcum[i+1]]<-0
        if (RandDist[i]=="Beta")  PsiM[(qcum[i]+1):qcum[i+1]]<-0.5
    }
    
            # Defining GammaM #        
        Lambda <- exp(DDR%*%DRgamma)
        Phi <- exp(DDY%*%DYgamma)

        GammaMvec<-c(Phi,Lambda)
        
        # Mean values #
        eta<-TT[1:n,]%*%as.matrix(c(Beta,unlist(VT)))
        mu<-B*InvLinkY(eta,Link)
        
        #second derivatives matricies #
        Wvec<-Wmatgen(mu,Link)
        WR<-list(0)
        UT<-0
        WRT<-list(0)
        for (i in 1:nrand) {
            U[[i]]<-LinkR(V[[i]],RandDist[i])
            if (i==1) UT<-U[[1]]
            else UT<-c(UT,U[[i]])
            WR[[i]]<-WRVC(U[[i]],RandDist[i])
            if (i==1) WRT<-WR[[1]]
            else WRT<-c(WRT,WR[[i]])
        }

        # Truncated Computations #
        #MTheta<-1-exp(-mu)
        #M1Theta<-exp(-mu)*mu
        #M2Theta<-exp(-mu)*mu*(1-mu)
        #M3Theta<-M2Theta*(1-mu)-mu*M1Theta
        #WTildevec<-as.vector(Wvec+((M2Theta/MTheta)-(M1Theta/MTheta)^2))
        WTotvec<-c(Wvec,WRT)

        
        ISIGMAMvec<-as.vector((1/GammaMvec)*WTotvec)
        
        HLikelihood<-0
        HLikelihood<-sum(Y*log(mu/(B-mu))+B*log(1-(mu/B)))
        CondLikelihood<-HLikelihood
    
        for (i in 1:nrand){
        if (RandDist[i]=="Normal") HLikelihood<-HLikelihood+sum((1/Lambda[(1+qcum[i]):qcum[i+1]])*(-(V[[i]]^2)/2)-0.5*log(Lambda[(1+qcum[i]):qcum[i+1]]*2*pi))
        if (RandDist[i]=="Gamma") HLikelihood<-HLikelihood+sum((1/Lambda[(1+qcum[i]):qcum[i+1]])*(V[[i]]-exp(V[[i]]))-
                                        (1/Lambda[(1+qcum[i]):qcum[i+1]])*log(Lambda[(1+qcum[i]):qcum[i+1]])-lgamma(1/Lambda[(1+qcum[i]):qcum[i+1]]))
        if (RandDist[i]=="IGamma"){HLikelihood<-HLikelihood+sum((1/Lambda[(1+qcum[i]):qcum[i+1]])*(V[[i]]+log(-V[[i]])))+sum(1+(1/Lambda[(1+qcum[i]):qcum[i+1]])*log((1/Lambda[(1+qcum[i]):qcum[i+1]])))-sum(lgamma(1+(1/Lambda[(1+qcum[i]):qcum[i+1]])))}
        if (RandDist[i]=="Beta"){HLikelihood<-HLikelihood+sum(lgamma(1/Lambda[(1+qcum[i]):qcum[i+1]]))-2*sum(lgamma(1/(2*Lambda[(1+qcum[i]):qcum[i+1]])))+sum((1/(2*Lambda[(1+qcum[i]):qcum[i+1]]))*log(U[[i]]*(1-U[[i]])))}        
    }
    
    # Computing Hessian when joint estimation also for the value of pd(h) #
    
    HessianD<-(t(TT)*rep(ISIGMAMvec,each=ncol(TT)))%*%TT
    HessianF<-HessianD[1:p,1:p]
    SES<-sqrt(diag(solve(HessianD)))
    SEBeta<-SES[1:p]
    SEVs<-SES[(p+1):(p+sum(q))]
    
    
    # Computing value of pd(h) #
    logdet<-determinant(HessianD/(2*pi),log=T)$modulus
    attributes(logdet)<-NULL
    APDispLikelihood<-HLikelihood-0.5*logdet
    
    TT2<-TT[,(p+1):(p+sum(q))]
    # Computing Hessian with respect to v #
    HessianV<-(t(TT2)*rep(ISIGMAMvec,each=ncol(TT2)))%*%TT2
    # Computing Value of pv(h)
    logdet<-determinant(HessianV/(2*pi),log=T)$modulus
    attributes(logdet)<-NULL
    APFixdLikelihood<-HLikelihood-0.5*logdet
    
    if (LaplaceFixed==TRUE){
          
        # Computing Standard Errors #
        SEVs<-sqrt(diag(solve(HessianV)))
        SEBeta<-rep(NA,p)
        

        
    }
    dmudeta<-dmudetagen(mu,Link)
    # Computing gradient of fixed effects #
    if (LaplaceFixed==FALSE){
        gradientFixed<-t(X)%*%(as.vector((Y-mu)/Phi)*as.vector(Wvec)*as.vector(1/dmudeta))
    }
    else gradientFixed<-rep(NA,length(Beta))
    
    PROC_OUTPUT<-list(HLikelihood=HLikelihood, APFixdLikelihood=APFixdLikelihood, APDispLikelihood=APDispLikelihood, CLikelihood=CondLikelihood, 
                GradientFixed=gradientFixed,ASEBeta=SEBeta, SEVs=SEVs,HessianF=HessianF)
    PROC_OUTPUT                             
}
