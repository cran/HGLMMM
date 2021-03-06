##########################################
##### Normal Distribution Estimation #####
##########################################
#####################################################
##### Estimation algorithm - gamma distribution #####
#####################################################

# We introduce link here: Inverse or Logarithm #

IWLS_Normal<-function(Y=NULL,X=NULL,ZZ=NULL,SS=NULL,Beta=NULL,Vstart=NULL,OFFSET=NULL,Link=c("Identity"),
            EstimateOverDisp=FALSE,LaplaceFixed=FALSE,RandDist=c("Normal"),DDR=NULL,DDY=NULL,DYgamma=NULL,DRgamma=NULL,Info=FALSE,DEBUG=FALSE,CONV=CONV){

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
        if (RandDist[i]=="Normal") PsiM[(qcum[i]+1):qcum[i+1]]<-0
        if (RandDist[i]=="Gamma")  PsiM[(qcum[i]+1):qcum[i+1]]<-1
        if (RandDist[i]=="IGamma") PsiM[(qcum[i]+1):qcum[i+1]]<-1
        if (RandDist[i]=="Beta")   PsiM[(qcum[i]+1):qcum[i+1]]<-0.5
    }
          
    Iteration<-0
    Convergence<-10
    
    while (Convergence>CONV){
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
        # Here add a correction for a non-canonical link #
        
        #print("ISIGMAMvec2");print(ISIGMAMvec)
        #Bvec<-(Wvec)-(Y-mu)*Amatgen(mu,Link) # This is the entry of ISIGMAMvec[1:n] with difference in Phi #
        dWdmu<-dWdmugen(mu,Link)
        dmudeta<-dmudetagen(mu,Link)
        #dAdmu<-dAdmugen(mu,Link)
        #Amat<-Amatgen(mu,Link)
        # Adjustment computation for the Laplace Approximation to the mean #
        
        Scorr<-rep(0,n)
        if (LaplaceFixed==TRUE){
            
            TT2<-TT[,(p+1):(p+sum(q))]
            
            INV1<-solve((t(TT2)*rep(ISIGMAMvec,each=ncol(TT2)))%*%TT2)
            PP2<-TT2%*%INV1%*%(t(TT2)*rep(ISIGMAMvec,each=ncol(TT2)))
            MOD<-INV1%*%(t(Z)*rep((1/Phi),each=ncol(Z)))
            ADJDER1<-list(0)
            ADJDER2<-list(0)
     
            for (i in 1:nrand){
                ADJDER1[[i]]<--ZZ[[i]]%*%(MOD[(qcum[i]+1):(qcum[i+1]),])#*rep(Wvec,each=nrow(MOD[(qcum[i]+1):(qcum[i+1]),]))) 
                ADJDER2[[i]]<--t(t(MOD[(qcum[i]+1):(qcum[i+1]),]))
            }        

            # Computing correction quantities for the Laplace Approximation of fixed effects #
            CorrTerms<-list(0)
            CorrTerms[[1]]<-diag(PP2[1:n,1:n])*as.vector(1/Wvec)*as.vector(1/Wvec)*(dWdmu)*dmudeta
            for (i in 1:nrand){
                ADJ1<-rep(0,n)
                ADJ2<-rep(0,n)
                ADJ1<-t(ADJDER1[[i]])%*%(as.vector(diag(PP2[1:n,1:n]))*as.vector(1/Wvec)*as.vector(dWdmu)*as.vector(dmudeta)) # Check this one
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
         #print("Block 2");print(CPTISIGMAM);#print("EISIGMAMvec");#print(EISIGMAMvec);print("zTot");print(zTot)
         

         MeanParms<-solve(CPTISIGMAM%*%TT,CPTISIGMAM%*%zTot)

         Beta<-MeanParms[1:p]
         for (i in 1:nrand){
            V[[i]]<-MeanParms[(p+qcum[i]+1):(p+qcum[i+1])]
            if (i==1) VT<-V[[i]]
            else VT<-c(VT,V[[i]])
         }

         Convergence<-sum(abs(MeanParms-MeanParmsLast))
        
         
         ###############################
         ##### Variance Components #####
         ###############################
        
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

        #Computing Deviances #
        DevianceRand<-list(0)
        DevianceRandCpy<-list(0)
        DevianceRandQCpy<-list(0)
        for (i in 1:nrand){
            DevianceRand[[i]]<-DevRand(V[[i]],U[[i]],RandDist[i])
        }            
    
        WTotvec<-c(Wvec,unlist(WRT))
        ISIGMAMvec<-(1/GammaMvec)*WTotvec
        # Here an adjustment for non-canonical links #
        
        # Create an expected version of the ISIGMAMvec #
        
        # Here add a correction for a non-canonical link #
        #ISIGMAMvec[1:n]<-ISIGMAMvec[1:n]-((Y-mu)/Phi)*Amatgen(mu,Link)
        #Bvec<-(Wvec)-(Y-mu)*Amatgen(mu,Link) # This is the entry of ISIGMAMvec[1:n] with difference in Phi #
        dWdmu<-dWdmugen(mu,Link)
        dmudeta<-dmudetagen(mu,Link)
        #dAdmu<-dAdmugen(mu,Link)
        #Amat<-Amatgen(mu,Link)
        
        CPTISIGMAM<-t(TT)*rep(ISIGMAMvec,each=ncol(TT))
        PMAT<-TT%*%solve(CPTISIGMAM%*%TT)%*%CPTISIGMAM
        DVhatDlambda<-list(0)    
    
        for (i in 1:nrand){
            DVhatDlambda[[i]]<--solve((t(ZZ[[i]])*rep(Wvec/Phi,each=ncol(ZZ[[i]])))%*%ZZ[[i]]+diag(WR[[i]]/Lambda[(1+qcum[i]):qcum[i+1]]))%*%(PsiM[(qcum[i]+1):(qcum[i+1])]-U[[i]])/(Lambda[(1+qcum[i]):qcum[i+1]]^2)
        }    

        qmod<-list(0)
        qCur<-list(0)
        qrr<-list(0)
        
        for (i in 1:nrand){
            SSCur<-rep(0,n+sum(q))
            SSCur[1:n]<-SS[[i]]
            SSCur[(n+1+qcum[i]):(n+qcum[i+1])]<-(1:(qcum[i+1]-qcum[i]))
            #print("SSCur");print(SSCur)
            DlogWDloglambda<-matrix(0,n+sum(q),n+sum(q))
            DlogWDloglambda[1:n,1:n]<-diag(as.vector(1/Wvec)*as.vector(dWdmu)*as.vector(dmudeta)*(as.vector(ZZ[[i]]%*%(DVhatDlambda[[i]]*Lambda[(1+qcum[i]):qcum[i+1]])))) # RRRRR
            DlogWDloglambda[(n+1+qcum[i]):(n+qcum[i+1]),(n+1+qcum[i]):(n+qcum[i+1])]<-diag(DWRDU(U[[i]],RandDist[i])*as.vector(DVhatDlambda[[i]]))*as.vector(Lambda[(1+qcum[i]):qcum[i+1]])  # RRRRR
            qmod[[i]]<-diag(PMAT%*%DlogWDloglambda)
            qCur[[i]]<-cbind(qmod[[i]],SSCur)

            qCur[[i]]<-tapply(qCur[[i]][,1],qCur[[i]][,2],sum)
            #print("qCur1");print(qCur[[i]]);print(sum(qCur[[i]]))
            qCur[[i]]<-qCur[[i]][row.names(qCur[[i]])!="0"]
            qrr[[i]]<-diag(PMAT[(n+1+qcum[i]):(n+qcum[i+1]),(n+1+qcum[i]):(n+qcum[i+1])])
            qrr[[i]]<-qrr[[i]]-qCur[[i]]

            # Correction to estimate the true likelihood instead of EQL #
            if (RandDist[i]=="Gamma") qrr[[i]]<-qrr[[i]]+1+2*(log(Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])+2*(digamma(1/Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])
            if (RandDist[i]=="IGamma") qrr[[i]]<-qrr[[i]]+1+(2/Lambda[(1+qcum[i]):(qcum[i+1])])-2*(log(1/Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])-2*((1+Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])+2*(digamma(1+(1/Lambda[(1+qcum[i]):(qcum[i+1])]))/Lambda[(1+qcum[i]):(qcum[i+1])])
            if (RandDist[i]=="Beta")  qrr[[i]]<-qrr[[i]]+1-2*(digamma(1/Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])+2*(digamma(1/(2*Lambda[(1+qcum[i]):(qcum[i+1])]))/Lambda[(1+qcum[i]):(qcum[i+1])])+log(4)/Lambda[(1+qcum[i]):(qcum[i+1])]
                                                
            # Applying the correction for the deviances #
            #DevianceRand[[i]][qrr[[i]]>0.999999]<-DevianceRand[[i]][qrr[[i]]>0.999999]+0.01*Lambda[(1+qcum[i]):(qcum[i+1])][qrr[[i]]>0.999999]
            #qrr[[i]][qrr[[i]]>0.999999]<-0.99
            # If qrr equal to one substract 0.01 and add to deviance #

            DevianceRandCpy[[i]]<-DevianceRand[[i]]
            DevianceRand[[i]]<-DevianceRand[[i]]/(1-qrr[[i]]) 
            DevianceRandQCpy<-DevianceRand

        }            

        ######################################
        ##### Gamma model for dispersion #####
        ######################################

        invSigmaGammaR<-list(0)
        for (i in 1:nrand){ 
            invSigmaGammaR[[i]]<-((1-qrr[[i]])/4)
        }
        muGammaR<-Lambda
        
        oldDRgamma<-DRgamma
        
        ksiR<-DDR%*%DRgamma
        
        DevianceRand<-unlist(DevianceRand)
        
        ZRresp<-ksiR+(DevianceRand-muGammaR)/muGammaR
        
        invSigmaGammaR<-diag(unlist(invSigmaGammaR))

        DRgamma<-solve(crossprod(DDR,invSigmaGammaR)%*%DDR,crossprod(DDR,invSigmaGammaR)%*%ZRresp)

        ##########################
        ##### Overdispersion #####
        ##########################
        
        if (EstimateOverDisp==TRUE) {
        Lambda <- exp(DDR%*%DRgamma)
        Phi <- exp(DDY%*%DYgamma)

        WTotvec<-c(Wvec,unlist(WRT))
        ISIGMAMvec<-(1/GammaMvec)*WTotvec
        #ISIGMAMvec[1:n]<-ISIGMAMvec[1:n]-((Y-mu)/Phi)*Amatgen(mu,Link)
        
        CPTISIGMAM<-t(TT)*rep(ISIGMAMvec,each=ncol(TT))
        PMAT<-TT%*%solve(CPTISIGMAM%*%TT)%*%CPTISIGMAM     
              
        qmodO<-0
        qrrO<-rep(0,n)
        
        qrrO<-diag(PMAT[1:n,1:n])
        qrrO<-qrrO

        # Applying the correction for the deviances #
        DevianceResp<-rep(0,n)

        DevianceResp<-(Y-mu)^2
        
        DevianceRespCpy<-DevianceResp
        DevianceResp<-DevianceResp/(1-qrrO) 

        # Algorithm for Gamma model #
        invSigmaGammaO<-((1-qrrO)/4)
        
        muGammaO<-Phi
        
        oldDYgamma<-DYgamma
        
        ksiO<-DDY%*%DYgamma
             
        ZOresp<-ksiO+(DevianceResp-muGammaO)/muGammaO
        if (DEBUG) {print("ZOresp");print(ZOresp)}
        invSigmaGammaO<-diag(invSigmaGammaO)

        DYgamma<-solve(crossprod(DDY,invSigmaGammaO)%*%DDY,crossprod(DDY,invSigmaGammaO)%*%ZOresp)
        if (DEBUG) {print("DYgamma");print(DYgamma)            }
        Convergence<-Convergence+sum(abs(DYgamma-oldDYgamma))
        if (DEBUG) {print("Convergence overdisp");print(Convergence)}
        }
                
        
        
        #print("Lambda");print(Lambda);
        #print("Beta");print(Beta);
   
        
        Convergence<-Convergence+sum(abs(DRgamma-oldDRgamma))        
        if (DEBUG) {
        print("DRConvergence");print(Convergence)
        print("NewDR");print(DRgamma)
        print("OldDR");print(oldDRgamma)    
        }
    }
    
    DevianceRand<-list(0)
    for (i in 1:nrand){
        DevianceRand[[i]]<-DevRand(V[[i]],U[[i]],RandDist[i])
    }        
    
    
    # Computation of Deviance residuals for Gamma models #
    DevianceRespCpy<-(Y-mu)^2
    DevianceRespOut<-sign(Y-mu)*sqrt(DevianceRespCpy)
    DevianceResp<-DevianceRespCpy/(1-qrrO)
    
    DevianceODisp<-2*(-log(DevianceResp/Phi)+(DevianceResp-Phi)/Phi)
    DevianceODisp<-sign(DevianceResp-Phi)*sqrt(DevianceODisp)
    DevianceDisp<-list(0)
    DevianceRandCpyOut<-list(0)
    DevianceRandQCpy<-list(0)
    
    DevianceRand<-list(0)
    for (i in 1:nrand){
        DevianceRand[[i]]<-DevRand(V[[i]],U[[i]],RandDist[i])
        DevianceRandCpy[[i]]<-DevianceRand[[i]]
        DevianceRandCpyOut[[i]]<-sign(PsiMstar[(qcum[i]+1):qcum[i+1]]-U[[i]])*sqrt(DevianceRandCpy[[i]])
        DevianceRandQCpy[[i]]<-DevianceRand[[i]]/(1-qrr[[i]]) 
        DevianceDisp[[i]]<-2*(-log(DevianceRandQCpy[[i]]/Lambda[(1+qcum[i]):(qcum[i+1])])+(DevianceRandQCpy[[i]]-Lambda[(1+qcum[i]):(qcum[i+1])])/Lambda[(1+qcum[i]):(qcum[i+1])])   
        DevianceDisp[[i]]<-sign(DevianceRandQCpy[[i]]-Lambda[(1+qcum[i]):(qcum[i+1])])*sqrt(DevianceDisp[[i]])
    }
    # Computation of standardized residuals #
    # Standardized residuals from the main model #
    HAT1<-diag(diag(n+sum(q))-TT%*%solve(CPTISIGMAM%*%TT)%*%CPTISIGMAM)

    StdDevianceRespOut<-DevianceRespOut/sqrt(Phi*HAT1[1:n])
    # Standardized residuals from the random effects model #
    StdDevianceRandCpyOut<-list(0)
    for (i in 1:nrand){
        StdDevianceRandCpyOut[[i]]<-DevianceRandCpyOut[[i]]/sqrt(Lambda[(1+qcum[i]):(qcum[i+1])]*HAT1[(n+1+qcum[i]):(n+qcum[i+1])])
    }
    # Standardized residuals from the overdispersion model #
    HAT2<-diag(diag(n)-DDY%*%solve(t(DDY)%*%invSigmaGammaO%*%DDY)%*%t(DDY)%*%invSigmaGammaO)

    StdDevianceODisp<-DevianceODisp/sqrt((2/(1-qrrO))*HAT2)
    # Standardized residuals from the dispersion model#
    HAT3<-diag(diag(sum(q))-DDR%*%solve(t(DDR)%*%invSigmaGammaR%*%DDR)%*%t(DDR)%*%invSigmaGammaR)

    StdDevianceDisp<-list(0)
    for (i in 1:nrand){
        StdDevianceDisp[[i]]<-DevianceDisp[[i]]/sqrt((2/(1-qrr[[i]]))*HAT3[(1+qcum[i]):(qcum[i+1])])
    }
    
    ##### Transformed scale fitted values #####
    # Still to do for each dispersion component and for the mean + overdispersion #
    
    
    PROC_OUTPUT<-list(Beta=Beta,Vs=VT,YDispParms=DYgamma,RDispParms=DRgamma,DevianceResidualY=DevianceRespOut,DevianceResidualR=DevianceRandCpyOut,
                        DevianceResidualODisp=DevianceODisp,DevianceResidualDisp=DevianceDisp,
                        StdDevianceResidualY=StdDevianceRespOut,StdDevianceResidualR=StdDevianceRandCpyOut,
                        StdDevianceResidualODisp=StdDevianceODisp,StdDevianceResidualDisp=StdDevianceDisp)
    PROC_OUTPUT
}
