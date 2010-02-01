###################################################
##### Main Program for the Gamma distribution #####
###################################################

GammaHGLM.fit<-function(YP,XP,ZZP,SSP,BetaP=NULL,VstartP=NULL,OFFSET=NULL,EstimOD=FALSE,LapFixP=FALSE,RDP=c("Normal"),
                        Link=c("Inverse"),DDRP=NULL,DDYP=NULL,DYgammaP=NULL,DRgammaP=NULL,DEBUG=FALSE,Info=FALSE,CONV=1e-4){
                 
    check<-require(numDeriv)
    if (!check) stop("I require numDeriv package to be installed ...")                     
    
    PROC_OUTPUT_TP<-0
  
    #################
    # Poisson Model #
    #################
    
    # Fit the model #
    fitmodelTP<-IWLS_Gamma(Y=YP,X=XP,ZZ=ZZP,SS=SSP,Beta=BetaP,Vstart=VstartP,OFFSET=OFFSET,EstimateOverDisp=EstimOD,LaplaceFixed=LapFixP,RandDist=RDP,
                            Link=Link,DDR=DDRP,DDY=DDYP,DYgamma=DYgammaP,DRgamma=DRgammaP,DEBUG=DEBUG,Info=Info,CONV=CONV)

    # Compute standard errors of dispersion components #
    
    gradDispTP<-ProfileDispGamma(Y=YP,X=XP,ZZ=ZZP,SS=SSP,OFFSET=OFFSET,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,LaplaceFixed=LapFixP,
                    RandDist=RDP,Link=Link,DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms,DRgamma=fitmodelTP$RDispParms)
                    
    hessDispTP<-jacobian(ProfileDispGamma,fitmodelTP$RDispParms,Y=YP,X=XP,ZZ=ZZP,SS=SSP,OFFSET=OFFSET,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,LaplaceFixed=LapFixP,
                    RandDist=RDP,Link=Link,DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms)
    hessDispTP<-(t(hessDispTP)+hessDispTP)/2
    SEDispTP<-sqrt(diag(solve(-hessDispTP)))

    # Overdispersion part #
    if (EstimOD==TRUE) {
    gradODispTP<-ProfileODispGamma(Y=YP,X=XP,ZZ=ZZP,SS=SSP,OFFSET=OFFSET,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,LaplaceFixed=LapFixP,
                    RandDist=RDP,Link=Link,DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms,DRgamma=fitmodelTP$RDispParms)


    hessODispTP<-jacobian(ProfileODispGamma,fitmodelTP$YDispParms,Y=YP,X=XP,ZZ=ZZP,SS=SSP,OFFSET=OFFSET,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,LaplaceFixed=LapFixP,
                    RandDist=RDP,Link=Link,DDR=DDRP,DDY=DDYP,DRgamma=fitmodelTP$RDispParms)

    hessODispTP<-(t(hessODispTP)+hessODispTP)/2
    SEODispTP<-sqrt(diag(solve(-hessODispTP)))

    }
    else {
    gradODispTP<-"Not Estimated"
    hessODispTP<-"Not Estimated"
    SEODispTP<-"Not Estimated"
    }

    print("StdErr Fixed")
    # Compute standard errors of fixed components if laplace approximation is used #
    
    if (LapFixP) {
        gradhessFixTP<-ComputeGradHess(GammaAdjProfV,fitmodelTP$Beta,OFFSET=OFFSET,Y=YP,X=XP,ZZ=ZZP,SS=SSP,Vstart=fitmodelTP$Vs,LaplaceFixed=LapFixP,RandDist=RDP,Link=Link,
                            DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms,DRgamma=fitmodelTP$RDispParms)
    
        HLL<-HlikeGamma(Y=YP,X=XP,ZZ=ZZP,SS=SSP,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,OFFSET=OFFSET,LaplaceFixed=LapFixP,
                    Link=Link,RandDist=RDP,DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms,DRgamma=fitmodelTP$RDispParms)
    
        PROC_OUTPUT_TP<-list(Beta=fitmodelTP$Beta,Vs=fitmodelTP$Vs,Dispersion=fitmodelTP$RDispParms,ResidualDispersion=fitmodelTP$YDispParms,GradientDisp=gradDispTP,HessianDisp=hessDispTP,
                                    SEDisp=SEDispTP,GradientFix=gradhessFixTP$gradient,HessianFix=gradhessFixTP$hessian,SEFix=gradhessFixTP$se,
                                    GradientResDisp=gradODispTP,HessianResDisp=hessODispTP,SEResDisp=SEODispTP,
                                    HLikelihood=HLL$HLikelihood,APFix=HLL$APFixdLikelihood,APDis=HLL$APDispLikelihood,CLikelihood=HLL$CLikelihood,SEVs=HLL$SEVs)
    
    }
    else {
        HLL<-HlikeGamma(Y=YP,X=XP,ZZ=ZZP,SS=SSP,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,OFFSET=OFFSET,LaplaceFixed=LapFixP,
                    Link=Link,RandDist=RDP,DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms,DRgamma=fitmodelTP$RDispParms)
    
        PROC_OUTPUT_TP<-list(Beta=fitmodelTP$Beta,Vs=fitmodelTP$Vs,Dispersion=fitmodelTP$RDispParms,ResidualDispersion=fitmodelTP$YDispParms,GradientDisp=gradDispTP,HessianDisp=hessDispTP,
                                    SEDisp=SEDispTP,GradientFix=HLL$GradientFixed,HessianFix=HLL$HessianF,SEFix=HLL$ASEBeta,
                                    GradientResDisp=gradODispTP,HessianResDisp=hessODispTP,SEResDisp=SEODispTP,
                                    HLikelihood=HLL$HLikelihood,APFix=HLL$APFixdLikelihood,APDis=HLL$APDispLikelihood,CLikelihood=HLL$CLikelihood,SEVs=HLL$SEVs)    
    }
   
       PROC_OUTPUT_TP<-list(Results=PROC_OUTPUT_TP,Details=fitmodelTP)
       PROC_OUTPUT_TP
}
