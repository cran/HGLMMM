################################################
##### Normal Distribution - FixedDispGrad ######
################################################


ComputeGradHess<-function(ObjectiveFunction=NULL,...,d=0.1){
    check<-require(numDeriv)
    if (!check) stop("Load numDeriv first please")
    g<-grad(ObjectiveFunction,...)
    print("######################################## Hessian ######################################")
    h<-hessian(ObjectiveFunction,...,method.args=list(d=d))
    se<-sqrt(diag(solve(-h)))
    output<-list(funcvalue=ObjectiveFunction(...),gradient=g,hessian=h,se=se)
    output
}


NormalAdjProfV<-function(Beta=NULL,Y=NULL,B=NULL,X=NULL,ZZ=NULL,SS=NULL,OFFSET=NULL,Vstart=NULL,EstimateOverDisp=FALSE,LaplaceFixed=FALSE,RandDist=c("Normal"),Link=c("Inverse"),
                            DDR=NULL,DDY=NULL,DYgamma=NULL,DRgamma=NULL,Info=FALSE,DEBUG=FALSE,RespDist=c("Normal"),RespLink=NULL){

# Adjustment of the standard estimation procedure to compute vhat given delta and beta #
                      
    # Creator of U the inverse of the link for V #
    # Creator of U the inverse of the link for V #
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
        eta
    }
    # Inverse of the link #
    InvLinkY<-function(eta,Link){
        if (Link=="Inverse")    mu<--(1/eta)
        if (Link=="Log")        mu<-exp(eta)
        if (Link=="Identity")   mu<-eta
        mu
    }
    # Generation of the weight matrix W#
    Wmatgen<-function(mu,Link){
        if (Link=="Inverse")    Wvec<-mu^4
        if (Link=="Log")        Wvec<-mu^2
        if (Link=="Identity")   Wvec<-rep(1,length(mu))
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
        if (Link=="Inverse")    dWdmu<-4*(mu^3)
        if (Link=="Log")        dWdmu<-2*mu
        if (Link=="Identity")   dWdmu<-rep(0,length(mu))
        dWdmu
    }
    # Generation of the derivative dmudeta #
    dmudetagen<-function(mu,Link){
        if (Link=="Inverse")    dmudeta<-mu^2
        if (Link=="Log")        dmudeta<-mu
        if (Link=="Identity")   dmudeta<-rep(1,length(mu))
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
    if (is.null(OFFSET)) OFFSET<-1
    V<-list(0)
    U<-list(0)
    TT<-cbind(X,Z)
    PsiM<-rep(0,sum(q))
    for (i in 1:nrand){
        TT<-rbind(TT,cbind(matrix(0,q[i],p+qcum[i]),diag(q[i]),matrix(0,q[i],qcum[nrand+1]-qcum[i+1])))
        V[[i]]<-as.matrix(Vstart[(qcum[i]+1):(qcum[i+1])])
        if (i==1) VT<-V[[1]]
        else VT<-c(VT,list(V[[i]]))
        if (RandDist[i]=="Normal")  PsiM[(qcum[i]+1):qcum[i+1]]<-0
        if (RandDist[i]=="Gamma")   PsiM[(qcum[i]+1):qcum[i+1]]<-1
        if (RandDist[i]=="IGamma")  PsiM[(qcum[i]+1):qcum[i+1]]<-1
        if (RandDist[i]=="Beta")    PsiM[(qcum[i]+1):qcum[i+1]]<-0.5
    }
          
    TT2<-TT[,(p+1):(p+sum(q))]
    XX2<-as.matrix(TT[,1:p])
    Iteration<-0
    Convergence<-10
    
    while (Convergence>1e-4){
    #for (iii in 1:3){
        Iteration<-Iteration+1
        if (Info) cat("\n Iteration: ",Iteration,"     Convergence: ",Convergence,"\n")
        MeanParmsLast<-c(Beta,unlist(VT))
  
        # Defining GammaM #        
        Lambda <- exp(DDR%*%DRgamma)
        Phi <- exp(DDY%*%DYgamma)
        
        GammaMvec<-c(Phi,Lambda)

        # Mean values #
        eta<-TT[1:n,]%*%as.matrix(c(Beta,unlist(VT)))
        mu<-InvLinkY(eta,Link)
        
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

        WTotvec<-c(Wvec,WRT)      
        ISIGMAMvec<-as.vector((1/GammaMvec)*WTotvec)
        
        # Adjustment computation for the Laplace Approximation to the mean #
        Scorr<-rep(0,n)
#        if (LaplaceFixed==TRUE){
#            DWTildeDthetavec<-WTildevec-(M2Theta/MTheta)+(M3Theta/MTheta)-3*((M1Theta*M2Theta)/(MTheta^2))+((M1Theta^2)/(MTheta^2))+2*((M1Theta^3)/(MTheta^3))
#            TT2<-TT[,(p+1):(p+sum(q))]
#            print("DimTT2");print(dim(TT2));print("Dim2");print(length(rep(ISIGMAMvec,each=nrow(t(TT2)))))
#            INV1<-solve((t(TT2)*rep(ISIGMAMvec,each=nrow(t(TT2))))%*%TT2)
#            PP2<-TT2%*%INV1%*%(t(TT2)*rep(ISIGMAMvec,each=nrow(t(TT2))))
#            MOD<-tcrossprod(INV1,Z)
#            ADJDER1<-list(0)
#            ADJDER2<-list(0)
#            WTildemuvec<-as.vector(WTildevec/mu)
#            for (i in 1:nrand){
#                ADJDER1[[i]]<--ZZ[[i]]%*%(MOD[(qcum[i]+1):(qcum[i+1]),]*rep(WTildemuvec,each=nrow(MOD[(qcum[i]+1):(qcum[i+1]),])))
#                ADJDER2[[i]]<--t(t(MOD[(qcum[i]+1):(qcum[i+1]),])*WTildemuvec)
#            }
          
            # Computing correction quantities for the Laplace Approximation of fixed effects #
#            CorrTerms<-list(0)
#            CorrTerms[[1]]<-diag(PP2[1:n,1:n])*as.vector(1/WTildevec)*as.vector(1/Wvec)*as.vector(DWTildeDthetavec)
#            for (i in 1:nrand){
#                ADJ1<-rep(0,n)
#                ADJ2<-rep(0,n)
#                ADJ1<-t(ADJDER1[[i]])%*%(as.vector(diag(PP2[1:n,1:n]))*as.vector(1/WTildevec)*as.vector(DWTildeDthetavec))
#                if (RandDist[i]=="Gamma") ADJ2<-t(ADJDER2[[i]])%*%as.vector(diag(PP2[(n+qcum[i]+1):(n+qcum[1+i]),(n+qcum[i]+1):(n+qcum[1+i])]))
#                CorrTerms<-c(CorrTerms,list(ADJ1,ADJ2))
#            }
   
#            CorrTermsLength<-length(CorrTerms)
#            CorrTerms<-as.matrix(unlist(CorrTerms))
          
#            dim(CorrTerms)<-c(n,CorrTermsLength)
#            print("Corrterms");print(CorrTerms)
                                 
#            Scorr<-0.5*mu*apply(CorrTerms,1,sum)
#            print("Scorr");print(Scorr)   
#         }
         Ystar<-Y
         dmudeta<-dmudetagen(mu,Link)
         zmain<-eta+(Ystar-mu)/dmudeta
         
         PsiMstar<-PsiM
         
         zrand<-list(0)
         for (i in 1:nrand){
            zrand[[i]]<-V[[i]]+(PsiMstar[(qcum[i]+1):qcum[i+1]]-U[[i]])/WR[[i]]
         }
         zrand<-as.matrix(unlist(zrand))
         
         zTot<-as.matrix(c(zmain,zrand))
         zTotV<-zTot-XX2%*%Beta
         # Updating Equations #
         VParmsLast<-c(unlist(VT))
         CPT2ISIGMAM<-t(TT2)*rep(ISIGMAMvec,each=nrow(t(TT2)))

         VParms<-solve(CPT2ISIGMAM%*%TT2,CPT2ISIGMAM%*%zTotV)
         
#         Beta<-MeanParms[1:p]
         for (i in 1:nrand){
            V[[i]]<-VParms[(qcum[i]+1):(qcum[i+1])]
            if (i==1) VT<-V[[i]]
            else VT<-c(VT,V[[i]])
         }
         Convergence<-sum(abs(VParms-VParmsLast))

    }

    eta<-TT[1:n,]%*%as.matrix(c(Beta,unlist(VT)))
    mu<-InvLinkY(eta,Link)
         
    Wvec<-Wmatgen(mu,Link)
    WR<-list(0)
    UT<-0
    WRT<-0
    for (i in 1:nrand) {
        U[[i]]<-LinkR(V[[i]],RandDist[i])
        if (i==1) UT<-U[[1]]
        else UT<-c(UT,U[[i]])
        WR[[i]]<-WRVC(U[[i]],RandDist[i])
        if (i==1) WRT<-WR[[1]]
        else WRT<-c(WRT,WR[[i]])
        }
 
    WTotvec<-c(Wvec,unlist(WRT))
        
    ISIGMAMvec<-(1/GammaMvec)*WTotvec                    
                            
    # Computing likelihood #
    
    Likelihood<-0
    
    # Response contribution to h #
    
    Likelihood<-Likelihood+sum(-0.5*log(2*pi*Phi))+sum(-0.5*((Y-mu)^2)/Phi)
    #print("Likelihood 1");print(Likelihood)
    # Random effects contribution to the likelihood #
    #print("U");print(U[[1]])
    for (i in 1:nrand){
        if (RandDist[i]=="Normal"){Likelihood<-Likelihood-0.5*sum(log(2*pi*Lambda[(1+qcum[i]):qcum[i+1]]))-0.5*sum(((VT[(1+qcum[i]):qcum[i+1]])^2)/Lambda[(1+qcum[i]):qcum[i+1]])}
        if (RandDist[i]=="Gamma"){Likelihood<-Likelihood+sum((1/Lambda[(1+qcum[i]):qcum[i+1]])*(V[[i]]-exp(V[[i]]))-
                                        (1/Lambda[(1+qcum[i]):qcum[i+1]])*log(Lambda[(1+qcum[i]):qcum[i+1]])-lgamma(1/Lambda[(1+qcum[i]):qcum[i+1]]))}
        if (RandDist[i]=="IGamma"){Likelihood<-Likelihood+sum((1/Lambda[(1+qcum[i]):qcum[i+1]])*(V[[i]]+log(-V[[i]])))+sum(1+(1/Lambda[(1+qcum[i]):qcum[i+1]])*log((1/Lambda[(1+qcum[i]):qcum[i+1]])))-sum(lgamma(1+(1/Lambda[(1+qcum[i]):qcum[i+1]])))}
        if (RandDist[i]=="Beta"){Likelihood<-Likelihood+sum(lgamma(1/Lambda[(1+qcum[i]):qcum[i+1]]))-2*sum(lgamma(1/(2*Lambda[(1+qcum[i]):qcum[i+1]])))+sum((1/(2*Lambda[(1+qcum[i]):qcum[i+1]]))*log(U[[i]]*(1-U[[i]])))}
 
    }
    #print("Likelihood 2");print(Likelihood)
    # Modification of the likelihood by the trace #
                
        
    Hmat<-(t(TT2)*rep(ISIGMAMvec,each=ncol(TT2)))%*%TT2
    Hmat<-Hmat/(2*pi)
    
    Likelihood<-Likelihood-0.5*determinant(Hmat,log=T)$modulus
    Likelihood
   
}   
