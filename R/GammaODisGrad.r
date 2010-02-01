# In this file we develope the OD inference for the Gamma model #
ProfileODispGamma<-function(DYgamma=NULL,Y=NULL,X=NULL,ZZ=NULL,SS=NULL,Beta=NULL,Vstart=NULL,OFFSET=NULL,Link=c("Inverse"),
            EstimateOverDisp=FALSE,LaplaceFixed=FALSE,RandDist=c("Normal"),DDR=NULL,DDY=NULL,DRgamma=NULL,Info=FALSE,DEBUG=FALSE){

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
    LinkY<-function(x,Link){
        if (Link=="Inverse")    out<--(1/x)
        if (Link=="Log")        out<-log(x)
        if (Link=="Identity")   out<-x
        out
    }
    # Inverse of the link #
    InvLinkY<-function(x,Link){
        if (Link=="Inverse")    out<--(1/x)
        if (Link=="Log")        out<-exp(x)
        if (Link=="Identity")   out<-x
        out
    }
    # Generation of the weight matrix W#
    Wmatgen<-function(x,Link){
        if (Link=="Inverse")    Wvec<-x^2
        if (Link=="Log")        Wvec<-rep(1,length(x))
        if (Link=="Identity")   Wvec<-1/(x^2)
        Wvec
    }
    # Generation of the adjustment to the weight matrix - factor (y-mu)/phi is left outside #
    #Amatgen<-function(mu,Link){
    #    if (Link=="Inverse")    Avec<-rep(0,length(mu))
    #    if (Link=="Log")        Avec<--(1/mu)
    #    Avec
    #}
    # Generation of the derivative dWdmu #
    dWdmugen<-function(x,Link){
        if (Link=="Inverse")    dWdmu<-2*x
        if (Link=="Log")        dWdmu<-rep(0,length(x))
        if (Link=="Identity")   dWdmu<--2*(1/x^3)
        dWdmu
    }
    # Generation of the derivative dmudeta #
    dmudetagen<-function(x,Link){
        if (Link=="Inverse")    dmudeta<-x^2
        if (Link=="Log")        dmudeta<-x
        if (Link=="Identity")   dmudeta<-rep(1,length(x))
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
        if (RandDist[i]=="Normal") PsiM[(qcum[i]+1):qcum[i+1]]<-0
        if (RandDist[i]=="Gamma")  PsiM[(qcum[i]+1):qcum[i+1]]<-1
        if (RandDist[i]=="IGamma") PsiM[(qcum[i]+1):qcum[i+1]]<-1
        if (RandDist[i]=="Beta")   PsiM[(qcum[i]+1):qcum[i+1]]<-0.5
    }
          
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

        # Create an expected version of the ISIGMAMvec #
        #EISIGMAMvec<-ISIGMAMvec
        # Here add a correction for a non-canonical link #
        #ISIGMAMvec[1:n]<-ISIGMAMvec[1:n]-((Y-mu)/Phi)*Amatgen(mu,Link)

        #Bvec<-(Wvec)-(Y-mu)*Amatgen(mu,Link) # This is the entry of ISIGMAMvec[1:n] with difference in Phi #
        dWdmu<-dWdmugen(mu,Link)
        dmudeta<-dmudetagen(mu,Link)
        #dAdmu<-dAdmugen(mu,Link)
        # Adjustment computation for the Laplace Approximation to the mean #
        
        Scorr<-rep(0,n)
        if (LaplaceFixed==TRUE){
            
            TT2<-TT[,(p+1):(p+sum(q))]
            
            INV1<-solve((t(TT2)*rep(ISIGMAMvec,each=ncol(TT2)))%*%TT2)
            PP2<-TT2%*%INV1%*%(t(TT2)*rep(ISIGMAMvec,each=ncol(TT2)))
            MOD<-INV1%*%(t(Z)*rep(1/Phi,each=ncol(Z)))
            ADJDER1<-list(0)
            ADJDER2<-list(0)
     
            for (i in 1:nrand){
                ADJDER1[[i]]<--ZZ[[i]]%*%(MOD[(qcum[i]+1):(qcum[i+1]),])#*rep(Wvec,each=nrow(MOD[(qcum[i]+1):(qcum[i+1]),]))) 
                ADJDER2[[i]]<--t(t(MOD[(qcum[i]+1):(qcum[i+1]),]))
            }        

            # Computing correction quantities for the Laplace Approximation of fixed effects #
            CorrTerms<-list(0)
            CorrTerms[[1]]<-diag(PP2[1:n,1:n])*as.vector(1/Wvec)*as.vector(1/Wvec)*as.vector(dWdmu)*as.vector(dmudeta)
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

                                 
            Scorr<-0.5*Phi*dmudeta*apply(CorrTerms,1,sum)

         }

         Ystar<-Y-Scorr  
         zmain<-eta+(Ystar-mu)/dmudeta       
         PsiMstar<-PsiM+(Lambda*crossprod(Z,as.matrix(Wvec*Scorr*(1/Phi)/dmudeta)))
         
         zrand<-list(0)
         for (i in 1:nrand){
            zrand[[i]]<-V[[i]]+(PsiMstar[(qcum[i]+1):qcum[i+1]]-U[[i]])/WR[[i]]
         }
         zrand<-as.matrix(unlist(zrand))
         
         zTot<-as.matrix(c(zmain,zrand))
         
         # Updating Equations #
         MeanParmsLast<-c(Beta,unlist(VT))
         CPTISIGMAM<-t(TT)*rep(ISIGMAMvec,each=nrow(t(TT)))


         aaa<-solve(CPTISIGMAM%*%TT)

         MeanParms<-solve(CPTISIGMAM%*%TT,CPTISIGMAM%*%zTot)

         Beta<-MeanParms[1:p]
         for (i in 1:nrand){
            V[[i]]<-MeanParms[(p+qcum[i]+1):(p+qcum[i+1])]
            if (i==1) VT<-V[[i]]
            else VT<-c(VT,V[[i]])
         }
         Convergence<-sum(abs(MeanParms-MeanParmsLast))
   }
   
         # Reevaluation of mean and u #
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
        
        # Create an expected version of the ISIGMAMvec #
        # EISIGMAMvec<-ISIGMAMvec
        # Here add a correction for a non-canonical link #
        #ISIGMAMvec[1:n]<-ISIGMAMvec[1:n]-((Y-mu)/Phi)*Amatgen(mu,Link)
        #Bvec<-(Wvec)-(Y-mu)*Amatgen(mu,Link) # This is the entry of ISIGMAMvec[1:n] with difference in Phi #
        dWdmu<-dWdmugen(mu,Link)
        dmudeta<-dmudetagen(mu,Link)
        #dAdmu<-dAdmugen(mu,Link)
        #Amat<-Amatgen(mu,Link)
        
        CPTISIGMAM<-t(TT)*rep(ISIGMAMvec,each=ncol(TT))
        PMAT<-TT%*%solve(CPTISIGMAM%*%TT)%*%CPTISIGMAM
        
        # Computing the derivative of vhat with respect to lambda #
        DVhatDlambda<-list(0)   
        
        INV1<-solve((t(TT)*rep(ISIGMAMvec,each=ncol(TT)))%*%TT)
        Derivative<-0
        
        dhhatdphi<-matrix(0,n,1)
        #Change below for different distribution #
        dhhatdphi<-(((Y/mu)+log(mu/Y))/Phi^2)-(1/Phi^2)+(log(Phi)/Phi^2)+(digamma(1/Phi)/Phi^2)
        
        Derivtemp<-as.vector(dhhatdphi)*as.vector(Phi)*DDY
        Derivative<-apply(Derivtemp,2,sum)
        
        # Computing derivative of the log det #
        # Compute the derivative of the whole diagnoal multiply by lambda and DDR - this should do the job #
        DiagDerivative<-matrix(0,1,ncol(DDY))
        AdjustmentDerivative<-matrix(0,1,ncol(DDY))
        for (i in 1:1){
        DiagDeriv<-matrix(0,n+sum(q),ncol(DDY))


        DiagDeriv[1:n,]<-as.vector(Wvec)*as.vector(-1/(Phi^2))*as.vector(Phi)*DDY
        DiagDeriv[(n+1+qcum[i]):(n+qcum[i+1]),]<-0
        #DiagDeriv has already derivatives with respect to all gammas - some are zeros where gamma=0 for all observations #
                
        # Computing adjustment for each gamma #


        for (j in 1:ncol(DDY)){
                AdjustmentDerivative[i,j]<--0.5*sum(diag(INV1%*%(t(TT)*rep(DiagDeriv[,j],each=ncol(TT)))%*%TT))
            }
        }
        Derivative<-Derivative+apply(AdjustmentDerivative,2,sum)
        PROC_OUTPUT<-list(Derivative=Derivative)

        Derivative

}        
