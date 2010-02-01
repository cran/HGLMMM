###################################################
##### Main Program for the Gamma distribution #####
###################################################

NormalHGLM.fit<-function(YP,XP,ZZP,SSP,BetaP=NULL,VstartP=NULL,OFFSET=NULL,EstimOD=FALSE,LapFixP=FALSE,RDP=c("Normal"),
                        Link=c("Identity"),DDRP=NULL,DDYP=NULL,DYgammaP=NULL,DRgammaP=NULL,Info=FALSE,DEBUG=FALSE,CONV=1e-4){
                         
    check<-require(numDeriv)
    if (!check) stop("I require numDeriv package to be installed ...")                     
    
    PROC_OUTPUT_TP<-0
  
    #################
    # Poisson Model #
    #################
    
    # Fit the model #
    fitmodelTP<-IWLS_Normal(Y=YP,X=XP,ZZ=ZZP,SS=SSP,Beta=BetaP,Vstart=VstartP,OFFSET=OFFSET,EstimateOverDisp=EstimOD,LaplaceFixed=LapFixP,RandDist=RDP,
                            Link=Link,DDR=DDRP,DDY=DDYP,DYgamma=DYgammaP,DRgamma=DRgammaP,DEBUG=DEBUG,Info=Info,CONV=CONV)

    # Compute standard errors of dispersion components #
    
    gradDispTP<-ProfileDispNormal(Y=YP,X=XP,ZZ=ZZP,SS=SSP,OFFSET=OFFSET,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,LaplaceFixed=LapFixP,
                    RandDist=RDP,Link=Link,DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms,DRgamma=fitmodelTP$RDispParms)       
    hessDispTP<-jacobian(ProfileDispNormal,fitmodelTP$RDispParms,Y=YP,X=XP,ZZ=ZZP,SS=SSP,OFFSET=OFFSET,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,LaplaceFixed=LapFixP,
                    RandDist=RDP,Link=Link,DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms)
    hessDispTP<-(t(hessDispTP)+hessDispTP)/2                 
    SEDispTP<-sqrt(diag(solve(-hessDispTP)))
    gradODispTP<-ProfileODispNormal(Y=YP,X=XP,ZZ=ZZP,SS=SSP,OFFSET=OFFSET,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,LaplaceFixed=LapFixP,
                    RandDist=RDP,Link=Link,DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms,DRgamma=fitmodelTP$RDispParms)
    hessODispTP<-jacobian(ProfileODispNormal,fitmodelTP$YDispParms,Y=YP,X=XP,ZZ=ZZP,SS=SSP,OFFSET=OFFSET,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,LaplaceFixed=LapFixP,
                    RandDist=RDP,Link=Link,DDR=DDRP,DDY=DDYP,DRgamma=fitmodelTP$RDispParms)
    hessODispTP<-(t(hessODispTP)+hessODispTP)/2
    SEODispTP<-sqrt(diag(solve(-hessODispTP)))
    # Compute standard errors of fixed components if laplace approximation is used #
    
    if (LapFixP) {
        gradhessFixTP<-ComputeGradHess(NormalAdjProfV,fitmodelTP$Beta,OFFSET=OFFSET,Y=YP,X=XP,ZZ=ZZP,SS=SSP,Vstart=fitmodelTP$Vs,LaplaceFixed=LapFixP,RandDist=RDP,Link=Link,
                            DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms,DRgamma=fitmodelTP$RDispParms)
    
        HLL<-HlikeNormal(Y=YP,X=XP,ZZ=ZZP,SS=SSP,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,OFFSET=OFFSET,LaplaceFixed=LapFixP,
                    Link=Link,RandDist=RDP,DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms,DRgamma=fitmodelTP$RDispParms)
    
        PROC_OUTPUT_TP<-list(Beta=fitmodelTP$Beta,Vs=fitmodelTP$Vs,Dispersion=fitmodelTP$RDispParms,ResidualDispersion=fitmodelTP$YDispParms,GradientDisp=gradDispTP,HessianDisp=hessDispTP,
                                    SEDisp=SEDispTP,GradientFix=gradhessFixTP$gradient,HessianFix=gradhessFixTP$hessian,SEFix=gradhessFixTP$se,
                                    GradientResDisp=gradODispTP,HessianResDisp=hessODispTP,SEResDisp=SEODispTP,
                                    HLikelihood=HLL$HLikelihood,APFix=HLL$APFixdLikelihood,APDis=HLL$APDispLikelihood,CLikelihood=HLL$CLikelihood,SEVs=HLL$SEVs,modelAll=fitmodelTP)
    
    }
    else {
        HLL<-HlikeNormal(Y=YP,X=XP,ZZ=ZZP,SS=SSP,Beta=fitmodelTP$Beta,Vstart=fitmodelTP$Vs,OFFSET=OFFSET,LaplaceFixed=LapFixP,
                    Link=Link,RandDist=RDP,DDR=DDRP,DDY=DDYP,DYgamma=fitmodelTP$YDispParms,DRgamma=fitmodelTP$RDispParms)
        PROC_OUTPUT_TP<-list(Beta=fitmodelTP$Beta,Vs=fitmodelTP$Vs,Dispersion=fitmodelTP$RDispParms,ResidualDispersion=fitmodelTP$YDispParms,GradientDisp=gradDispTP,HessianDisp=hessDispTP,
                                    SEDisp=SEDispTP,GradientFix=HLL$GradientFixed,HessianFix=HLL$HessianF,SEFix=HLL$ASEBeta,
                                    GradientResDisp=gradODispTP,HessianResDisp=hessODispTP,SEResDisp=SEODispTP,
                                    HLikelihood=HLL$HLikelihood,APFix=HLL$APFixdLikelihood,APDis=HLL$APDispLikelihood,CLikelihood=HLL$CLikelihood,SEVs=HLL$SEVs)
            
    }
        PROC_OUTPUT_TP<-list(Results=PROC_OUTPUT_TP,Details=fitmodelTP)
        PROC_OUTPUT_TP
}
