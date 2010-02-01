####################################
##### Copied from BinGradHess2 #####
####################################


ProfileDispBinom<-function(DRgamma=NULL,Y=NULL,B=NULL,X=NULL,ZZ=NULL,SS=NULL,OFFSET=NULL,Beta=NULL,VStart=NULL,EstimateOverDisp=FALSE,LaplaceFixed=FALSE,RandDist=c("Normal"),Link="Logit",
                            DDR=NULL,DDY=NULL,DYgamma=NULL,Info=FALSE,DEBUG=FALSE){
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
    if (is.null(VStart)) VStart<-rep(0,sum(q))
    V<-list(0)
    U<-list(0)
    TT<-cbind(X,Z)
    PsiM<-rep(0,sum(q))
    for (i in 1:nrand){
        TT<-rbind(TT,cbind(matrix(0,q[i],p+qcum[i]),diag(q[i]),matrix(0,q[i],qcum[nrand+1]-qcum[i+1])))
        V[[i]]<-as.matrix(VStart[(qcum[i]+1):(qcum[i+1])])
        if (i==1) VT<-V[[1]]
        else VT<-c(VT,list(V[[i]]))
        if (RandDist[i]=="Normal") PsiM[(qcum[i]+1):qcum[i+1]]<-0
        if (RandDist[i]=="Gamma")  PsiM[(qcum[i]+1):qcum[i+1]]<-1
        if (RandDist[i]=="IGamma") PsiM[(qcum[i]+1):qcum[i+1]]<-1
        if (RandDist[i]=="Beta")   PsiM[(qcum[i]+1):qcum[i+1]]<-0.5
    }
          
    Iteration<-0
    Convergence<-10
    
    while (Convergence>1e-4){
    #for (iii in 1:2){
        Iteration<-Iteration+1
        if (Info) cat("\n Iteration: ",Iteration,"     Convergence: ",Convergence,"\n")
        MeanParmsLast<-c(Beta,unlist(VT))
  
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

        
        dWdmu<-dWdmugen(mu,Link)
        dmudeta<-dmudetagen(mu,Link)
        
        ISIGMAMvec<-as.vector((1/GammaMvec)*WTotvec)

        # Adjustment computation for the Laplace Approximation to the mean - SKIP SKIP SKIP for the moment SKIP SKIP SKIP#
        Scorr<-rep(0,n)
        if (LaplaceFixed==TRUE){
            #DWTildeDthetavec<-WTildevec-(M2Theta/MTheta)+(M3Theta/MTheta)-3*((M1Theta*M2Theta)/(MTheta^2))+((M1Theta^2)/(MTheta^2))+2*((M1Theta^3)/(MTheta^3))
            TT2<-TT[,(p+1):(p+sum(q))]

            INV1<-solve((t(TT2)*rep(ISIGMAMvec,each=nrow(t(TT2))))%*%TT2)
            PP2<-TT2%*%INV1%*%(t(TT2)*rep(ISIGMAMvec,each=nrow(t(TT2))))
            MOD<-INV1%*%(t(Z)*rep((1/Phi),each=ncol(Z)))
            ADJDER1<-list(0)
            ADJDER2<-list(0)
            #WTildemuvec<-as.vector(WTildevec/mu)
            for (i in 1:nrand){
                ADJDER1[[i]]<--ZZ[[i]]%*%(MOD[(qcum[i]+1):(qcum[i+1]),])
                ADJDER2[[i]]<--(MOD[(qcum[i]+1):(qcum[i+1]),])
            }
          
            # Computing correction quantities for the Laplace Approximation of fixed effects #
            CorrTerms<-list(0)
            CorrTerms[[1]]<-diag(PP2[1:n,1:n])*as.vector(1/Wvec)*as.vector(1/Wvec)*(dWdmu)*dmudeta
            for (i in 1:nrand){
                ADJ1<-rep(0,n)
                ADJ2<-rep(0,n)
                ADJ1<-t(ADJDER1[[i]])%*%(as.vector(diag(PP2[1:n,1:n]))*as.vector(1/Wvec)*as.vector(dWdmu)*as.vector(dmudeta))
                if (RandDist[i]=="Gamma") ADJ2<-t(ADJDER2[[i]])%*%as.vector(diag(PP2[(n+qcum[i]+1):(n+qcum[1+i]),(n+qcum[i]+1):(n+qcum[1+i])]))
                if (RandDist[i]=="IGamma") ADJ2<-t(ADJDER2[[i]])%*%as.vector(diag(PP2[(n+qcum[i]+1):(n+qcum[1+i]),(n+qcum[i]+1):(n+qcum[1+i])])*as.vector(2*U[[i]]))
                if (RandDist[i]=="Beta") ADJ2<-t(ADJDER2[[i]])%*%as.vector(diag(PP2[(n+qcum[i]+1):(n+qcum[1+i]),(n+qcum[i]+1):(n+qcum[1+i])])*as.vector(1-2*U[[i]]))
                CorrTerms<-c(CorrTerms,list(ADJ1,ADJ2))
            }
   
            CorrTermsLength<-length(CorrTerms)
            CorrTerms<-as.matrix(unlist(CorrTerms))
          
            dim(CorrTerms)<-c(n,CorrTermsLength)

                                 
            Scorr<-0.5*Phi*dmudeta*apply(CorrTerms,1,sum) #dmudeta - change with distribution

         }
         Ystar<-Y-Scorr  
         zmain<-eta+(Ystar-mu)/dmudeta
         
         PsiMstar<-PsiM+(Lambda*crossprod(Z,as.matrix(Wvec*Scorr*(1/Phi)/dmudeta))) ##### Change with distribution /mu*(1-mu) #####
         
         zrand<-list(0)
         for (i in 1:nrand){
            zrand[[i]]<-V[[i]]+(PsiMstar[(qcum[i]+1):qcum[i+1]]-U[[i]])/WR[[i]] ##### KEEP IT NOT CHANGED FOR A WHILE #####
         }
         zrand<-as.matrix(unlist(zrand))
         
         zTot<-as.matrix(c(zmain,zrand))
         # Updating Equations #
         MeanParmsLast<-c(Beta,unlist(VT))
         CPTISIGMAM<-t(TT)*rep(ISIGMAMvec,each=nrow(t(TT)))

         MeanParms<-solve(CPTISIGMAM%*%TT,CPTISIGMAM%*%zTot)
         Beta<-MeanParms[1:p]
         for (i in 1:nrand){
            V[[i]]<-MeanParms[(p+qcum[i]+1):(p+qcum[i+1])]
            if (i==1) VT<-V[[i]]
            else VT<-c(VT,V[[i]])
         }
         Convergence<-sum(abs(MeanParms-MeanParmsLast))
         
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
        #M4Theta<-M3Theta*(log(1-MTheta)+1)+2*M2Theta*log(1-MTheta)+M1Theta*log(1-MTheta)
        #WTildevec<-as.vector(Wvec+((M2Theta/MTheta)-(M1Theta/MTheta)^2))
        WTotvec<-c(Wvec,WRT)
        #print("WTotvec");print(length(WTotvec))
        #print("GammaMvec");print(length(GammaMvec))
        
        ISIGMAMvec<-as.vector((1/GammaMvec)*WTotvec)
        #print("Dim1");print(length(ISIGMAMvec))
        # Adjustment computation for the Laplace Approximation to the mean #
        #DWTildeDthetavec<-WTildevec-(M2Theta/MTheta)+(M3Theta/MTheta)-3*((M1Theta*M2Theta)/(MTheta^2))+((M1Theta^2)/(MTheta^2))+2*((M1Theta^3)/(MTheta^3))
        
                   
        
         
        # Computing the derivative of vhat with respect to lambda #
        
        DVhatDlambda<-list(0)
                    
        for (i in 1:nrand){
            DVhatDlambda[[i]]<--solve((t(ZZ[[i]])*rep(Wvec/Phi,each=ncol(ZZ[[i]])))%*%ZZ[[i]]+diag(WR[[i]]/Lambda[(1+qcum[i]):qcum[i+1]]))%*%(PsiM[(qcum[i]+1):(qcum[i+1])]-U[[i]])/(Lambda[(1+qcum[i]):qcum[i+1]]^2)
        }
        INV1<-solve((t(TT)*rep(ISIGMAMvec,each=ncol(TT)))%*%TT)
        Derivative<-0
        
        # Reading in from DDR how many covariates are assigned to each lambda #
        #ncovdisp<-rep(0,nrand)
        #for (i in 1:nrand){
        #    mattemp<-DDR[(1+qcum[i]):qcum[i+1],]
        #    index<-0
        #    for (j in 1:ncol(mattemp)){
        #        index[j]<-sum(mattemp[,j]!=0)
        #    }
        #    mattemp<-mattemp[,(index!=0)]
        #    ncovdisp[i]<-ncol(mattemp)
        #}
        #print("ncovdisp");print(ncovdisp)
        
        # Computing derivative #
        dhhatdlambda<-matrix(0,sum(q),1)
        for (i in 1:nrand){
            if (RandDist[i]=="Normal") dhhatdlambda[(1+qcum[i]):qcum[i+1]]<-(-1/Lambda[(1+qcum[i]):qcum[i+1]]^2)*(-(V[[i]]^2)/2)-0.5*(1/Lambda[(1+qcum[i]):qcum[i+1]])
            if (RandDist[i]=="Gamma")  dhhatdlambda[(1+qcum[i]):qcum[i+1]]<-(-1/Lambda[(1+qcum[i]):qcum[i+1]]^2)*(V[[i]]-exp(V[[i]])-digamma(1/Lambda[(1+qcum[i]):qcum[i+1]])+1-log(Lambda[(1+qcum[i]):qcum[i+1]]))
            if (RandDist[i]=="Beta")   dhhatdlambda[(1+qcum[i]):qcum[i+1]]<-(-1/Lambda[(1+qcum[i]):qcum[i+1]]^2)*(0.5*V[[i]]-log(1+exp(V[[i]]))+digamma(1/Lambda[(1+qcum[i]):qcum[i+1]])-digamma(1/(2*Lambda[(1+qcum[i]):qcum[i+1]])))
            if (RandDist[i]=="IGamma") dhhatdlambda[(1+qcum[i]):qcum[i+1]]<-(-1/Lambda[(1+qcum[i]):qcum[i+1]]^2)*(V[[i]]-(-log(-V[[i]]))+(Lambda[(1+qcum[i]):qcum[i+1]]+1)+log(1/Lambda[(1+qcum[i]):qcum[i+1]])-digamma(1+1/Lambda[(1+qcum[i]):qcum[i+1]]))
        }
        
        
        Derivtemp<-as.vector(dhhatdlambda)*as.vector(Lambda)*DDR
        Derivative<-apply(Derivtemp,2,sum)
        #print("Result:")
        #print(Derivative)
        
        # Computing derivative of the log det #
        # Compute the derivative of the whole diagnoal multiply by lambda and DDR - this should do the job #
        DiagDerivative<-matrix(0,nrand,ncol(DDR))
        AdjustmentDerivative<-matrix(0,nrand,ncol(DDR))
        for (i in 1:nrand){
        DiagDeriv<-matrix(0,n+sum(q),ncol(DDR))

        DiagDeriv[1:n,]<-as.vector(1/Phi)*as.vector(dWdmu)*as.vector(dmudeta)*(ZZ[[i]]%*%(as.vector(DVhatDlambda[[i]])*(Lambda[(1+qcum[i]):qcum[i+1]]*DDR[(1+qcum[i]):qcum[i+1],])))       # This is already multiplied by lambda
        DiagDeriv[(n+1+qcum[i]):(n+qcum[i+1]),]<-(DWRDU(U[[i]],RandDist[i])*WR[[i]]*as.vector(DVhatDlambda[[i]])- # This is already multiplied by lambda
                                                                                WR[[i]]*(1/Lambda[(1+qcum[i]):qcum[i+1]]))*DDR[(1+qcum[i]):qcum[i+1],]
        #DiagDeriv has already derivatives with respect to all gammas - some are zeros where gamma=0 for all observations #
                
        # Computing adjustment for each gamma #
        
        for (j in 1:ncol(DDR)){
                AdjustmentDerivative[i,j]<--0.5*sum(diag(INV1%*%(t(TT)*rep(DiagDeriv[,j],each=ncol(TT)))%*%TT))
            }
        }
        Derivative<-Derivative+apply(AdjustmentDerivative,2,sum)
        PROC_OUTPUT<-list(Derivative=Derivative)
  
        Derivative      
  }
