BootstrapEnvelopeHGLM<-function(model,numBoot,seed){
    callit<-model$CALL
    # Check the distribution for the response #
    DistResp<-callit[[match("DistResp",names(callit))]]
    ODEst<-callit[[match("ODEst",names(callit))]]
    ODEstVal<-callit[[match("ODEstVal",names(callit))]]
    if (DistResp=="Poisson" | DistResp=="Binomial") {
        if (ODEst==TRUE) stop("Distribution Poisson/Binomial with estimated overdispersion \n Sampling not possible")
        ODEstVal<-eval(ODEstVal)
        # Below we leave the room for improvement - it can happen there are two parameters -2 2 and their sum is equal to zero #
        if (ODEst==FALSE & sum(ODEstVal)!=0) stop("Distribution Poisson/Binomial with fixed overdispersion different from one \n Sampling not possible")
    }
       
    # Extract design matrix #
    formulaMain<-as.formula(callit[[match("formulaMain",names(callit))]])
    fr<-HGLMFrames(callit,formulaMain,contrasts)
    Xmat<-fr$X
    # Extract the overdispersion design matrix #
    formulaOD<-as.formula(callit[[match("formulaOD",names(callit))]])
    frOD<-HGLMFramesOD(callit,formulaOD)$X
    # Define the Frame for each Dispersion modelling #
    frRD<-list(0)
    frRDX<-list(0)
    namesRD<-list(0)
    DistRand<-eval(callit[[match("DistRand",names(callit))]])
    formulaRand<-list(0)
    formulaRand<-eval(callit[[match("formulaRand",names(callit))]])   
    for (i in 1:length(DistRand)){
        frRD[[i]]<-HGLMFramesRD(callit,formulaRand[[i]],whichDispComp=i)
        frRDX[[i]]<-frRD[[i]]$X
        namesRD[[i]]<-names(frRD[[i]]$fixef)
    }  
    FL<-HGLMFactorList(formulaMain,fr,0L,0L) 
    
    Y<-matrix(fr$Y,length(fr$Y),1)
    X<-fr$X
    ZZ<-FL$Design
    SS<-FL$Subject    
    q<-0
    
    for (i in 1:length(ZZ)){
        if (i==1) { Z<-ZZ[[1]]
                    q[i]<-dim(ZZ[[1]])[2]}
        else {      Z<-cbind(Z,ZZ[[i]])
                    q[i]<-dim(ZZ[[i]])[2]}
    }
    
    # Bootstraping the random effects #
    # Count the dimension of the parameters for the design matrix #
    numpar<-unlist(lapply(frRDX,ncol))
    numpar<-c(0,cumsum(numpar))
    
    CurDes<-list(0)
    for (i in 1:length(DistRand)){
        # Create the design matrix for the dispersion components #
        Disp<-model$Results$Dispersion[(numpar[i]+1):numpar[i+1]]
        CurDes[[i]]<-exp(frRDX[[i]]%*%Disp)
    }
    
    # Determine link #
    Link<-callit[[match("Link",names(callit))]]
    ###########################
    # Perform actual sampling #
    ###########################
    
    SampleResponse<-matrix(0,nrow(Y),numBoot)
    set.seed(seed)
    for (iii in 1:numBoot){
        # Sample Random Effects #
        uuu<-list(0)
        vvv<-list(0)
        for (i in 1:length(DistRand)){
            if (DistRand[i]=="Normal") uuu[[i]]<-rnorm(nrow(CurDes[[i]]),mean=0,sd=sqrt(CurDes[[i]]))
            if (DistRand[i]=="Gamma")  uuu[[i]]<-rgamma(nrow(CurDes[[i]]),shape=(1/CurDes[[i]]),scale=CurDes[[i]])
            if (DistRand[i]=="Beta")   uuu[[i]]<-rbeta(nrow(CurDes[[i]]),shape1=(2*CurDes[[i]]),shape2=(2*CurDes[[i]]))
            if (DistRand[i]=="IGamma") uuu[[i]]<-(1/rgamma(nrow(CurDes[[i]]),shape=(1+(1/CurDes[[i]])),scale=CurDes[[i]]))   
            if (DistRand[i]=="Normal") vvv[[i]]<-uuu[[i]]
            if (DistRand[i]=="Gamma")  vvv[[i]]<-log(uuu[[i]])
            if (DistRand[i]=="IGamma") vvv[[i]]<-(-1/uuu[[i]])
            if (DistRand[i]=="Beta")   vvv[[i]]<-log(uuu[[i]]/(1-uuu[[i]]))
        }
        vvv<-unlist(vvv)
        # Create the linear predictor #
        eta<-X%*%model$Results$Beta+Z%*%vvv
        
        if (DistResp=="Binomial") {
            BinomialDen<-eval(callit[[match("BinomialDen",names(callit))]])
            if (is.na(BinomialDen)) stop("(BDWBD) Binomial Data Without Binomial Denominators")
            if (Link!="Logit") stop("(LINK) For binomial data only logit link is currently used")
            if (Link=="Logit") pSuccess<-(exp(eta)/(1+exp(eta)))
            SampleResponse[,iii]<-rbinom(nrow(Y),BinomialDen,pSuccess)
        }
        if (DistResp=="Poisson"){
            Offset<-eval(callit[[match("Offset",names(callit))]])
            if (is.null(Offset)) Offset<-1
            if (Link!="Log") stop("(LINK) For Poisson data only log link is currently used")
            if (Link=="Log") mu<-Offset*exp(eta)
            SampleResponse[,iii]<-rpois(nrow(Y),mu)
        }
        if (DistResp=="Normal"){
            ResidualDispersion<-exp(frOD%*%model$Results$ResidualDispersion)
            if (Link!="Identity") stop("(LINK) For normal data only identity link is currently used")
            if (Link=="Identity") mu<-eta
            SampleResponse[,iii]<-rnorm(nrow(Y),mean=mu,sd=sqrt(ResidualDispersion))
        }
        if (DistResp=="Gamma"){
            ResidualDispersion<-exp(frOD%*%model$Results$ResidualDispersion)
            if (Link!="Log" | Link!="Inverse") stop("(LINK) For gamma data only log and inverse links are currently used")
            if (Link=="Log") mu<-exp(eta)
            if (Link=="Inverse") mu<--1/eta
            SampleResponse[,iii]<-rgamma(nrow(Y),shape=1/ResidualDispersion,scale=mu*ResidualDispersion)
        }
        
    }
    # Now fit model to each generated data and retain residuals #
    # In the final step create the qq-plot - GOOD LUCK #
    FinalResiduals<-matrix(0,nrow(Y),numBoot)
    Finalqq<-matrix(0,nrow(Y),numBoot)  
    # Select MainData #
    MD<-eval(callit[[match("DataMain",names(callit))]])
    RESP<-formulaMain[[2]]
    selmat<-matrix(0,length(names(MD)),2)
    selmat[,1]<-names(MD)
    selmat[,2]<-1:nrow(selmat)
    selcol<-as.numeric(selmat[selmat[,1]==RESP,2])
    for (iii in 1:numBoot){
        MD[,selcol]<-SampleResponse[,iii]
        callitfit<-callit
        callitfit$DataMain<-as.name("MD")
        modeltemporary<-eval(callitfit)
        FinalResiduals[,iii]<-modeltemporary$Details$StdDevianceResidualY
        tempqq<-qqnorm(FinalResiduals[,iii],plot.it=FALSE)
        tempqq<-data.frame(x=tempqq$x,y=tempqq$y)
        tempqq<-tempqq[order(tempqq$x),]
        Finalqq[,iii]<-tempqq$y
    }
    # Now create a plot #
    a<-model$Details$StdDevianceResidualY
    a<-qqnorm(a,plot.it=FALSE)
    a<-data.frame(x=a$x,y=a$y)
    a<-a[order(a$x),]
    mmin<-apply(Finalqq,1,min)
    mmax<-apply(Finalqq,1,max)
    plot(a$x,a$y,type="n",xlim=c(-4,4),ylim=c(-4,4),xlab="Theoretical Quantiles",ylab="Sample Quantiles")
    lines(a$x,a$y,lty=1,col="black",lwd=1);abline(a=0,b=1,col="black",lty=3)
    lines(a$x,mmin,lty=4,col="black",lwd=1)
    lines(a$x,mmax,lty=4,col="black",lwd=1)   
}


