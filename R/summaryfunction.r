# Testing space: we create the class object #

# This function cleans the object gives it the desired shape/dimensions and assigns names #



summary.HGLM<-function(object,...,V=FALSE){
    assignNames<-function(object){
        # Assign names to Beta #
        mc<-object$CALL
        names(object$Results$Beta)<-object$NAMES$namesX
        # Extract which dispersion componentes belong to which names and assign names #
        ddims<-lapply(object$NAMES$namesRD,length)
        namesRD<-unlist(object$NAMES$namesRD)
        row.names(object$Results$Dispersion)<-namesRD
        if (!is.matrix(object$Results$ResidualDispersion)){
            object$Results$ResidualDispersion<-as.matrix(object$Results$ResidualDispersion,length(object$Results$ResidualDispersion),1)
        }
        row.names(object$Results$ResidualDispersion)<-object$NAMES$namesOD
        # Determine the lengths of the random effects vectors #
        vdims<-unlist(lapply(object$NAMES$indexRE,length))
        vdims<-c(0,cumsum(vdims))
        indexRE<-unlist(object$NAMES$indexRE)
        for (i in 1:(length(vdims)-1)){
            names(object$Results$Vs)[(vdims[i]+1):vdims[i+1]]<-paste(object$NAMES$namesRE[i],indexRE[(vdims[i]+1):vdims[i+1]])
        }    
        if (!is.matrix(object$Results$SEDisp)){
            object$Results$SEDisp<-as.matrix(object$Results$SEDisp,length(object$Results$SEDisp),1)
        }     
        object
    }

    object<-assignNames(object)
    numbercomp<-length(object$NAMES$namesRE)
    displength<-unlist(lapply(object$NAMES$namesRD,length))
    displength<-c(0,cumsum(displength))
    
    # print("block 1")
    # Should we display random effects summary #
    ans<-list(0)
    ans[[1]]<-V
    names(ans)[1]<-"VsDisplay"
    ans$DispNames<-object$NAMES$namesRE
    ans$call<-object$CALL
    object<-object$Results
    # Fixed effects in the mean structure calculation #
    ans$fixedcoeff<- matrix(NA, 0L, 4L)
    est<-object$Beta
    ses<-object$SEFix
    zstat<-est/ses
    pval<-2*pnorm(abs(zstat),lower.tail=FALSE)
    grad<-object$GradientFix
    ans$fixedcoeff<-cbind(est,ses,zstat,pval)
    colnames(ans$fixedcoeff) <- c("Estimate", 
            "Std. Error", "Z value", "Pr(>|Z|)")
            
    # Random effects in the mean structure calculation #
    ans$randomcoeff<-matrix(NA,0L,4L)
    est<-object$Vs
    ses<-object$SEVs
    zstat<-est/ses
    pval<-2*pnorm(abs(zstat),lower.tail=FALSE)
    ans$randomcoeff<-cbind(est,ses,zstat,pval)
    colnames(ans$randomcoeff) <- c("Estimate", 
            "Std. Error", "Z value", "Pr(>|Z|)")
    
    # Overdispersion #
    if (object$GradientResDisp[1]=="Not Estimated"){
        ans$ODcoeff<-object$ResidualDispersion
    }
    else {
        ans$ODcoeff<-matrix(NA,0L,4L)
        est<-object$ResidualDispersion
        ses<-object$SEResDisp
        zstat<-est/ses
        pval<-2*pnorm(abs(zstat),lower.tail=FALSE)
        ans$ODcoeff<-cbind(est,ses,zstat,pval)
        colnames(ans$ODcoeff) <- c("Estimate", 
            "Std. Error", "Z value", "Pr(>|Z|)")
    }
    # Dispersion #
    ans$Discoeff<-list(0)
    for (i in 1:numbercomp){
        est<-object$Dispersion[(displength[i]+1):displength[i+1],1]
        ses<-object$SEDisp[(displength[i]+1):displength[i+1],1]
       # print("block 2")
        zstat<-est/ses
        pval<-2*pnorm(abs(zstat),lower.tail=FALSE)
        ans$Discoeff[[i]]<-cbind(est,ses,zstat,pval)
        colnames(ans$Discoeff[[i]]) <- c("Estimate", 
            "Std. Error", "Z value", "Pr(>|Z|)")
    }
    # Likelihood values #
    ans$HLikelihood<-object$HLikelihood
    ans$APFix<-object$APFix
    ans$APDis<-object$APDis
    ans$CLikelihood<-object$CLikelihood
    class(ans)<-"summary.HGLM"
    ans        
}



print.summary.HGLM<-function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) {
    
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")

    
    
    
    cat("===== Fixed Coefficients - Mean Structure =====\n\n")
    coefs <- x$fixedcoeff
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
    na.print = "NA", ...)  
    cat("\n\n")
    if (x$VsDisplay) {
        cat("===== Random Coefficients - Mean Structure =====\n")
        coefs<-x$randomcoeff
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)  
        cat("\n\n")
    }
    # Remember Overdispersion can be fixed #
    
    if (ncol(x$ODcoeff)==1){
        cat("===== Overdispersion Parameters Fixed =====\n")
        coefs <- x$ODcoeff
        colnames(coefs)<-"Fixed Value"
        print(coefs)
        cat("\n\n")
    }
    else {
        cat("===== Overdispersion Parameters Estimated =====\n\n")
        coefs <- x$ODcoeff
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
        na.print = "NA", ...)  
        cat("\n\n")  
    }
    
    # Printing Dispersion #
    cat("===== Dispersion Parameters Estimated =====\n\n")
    for (i in 1:length(x$Discoeff)){
        cat(paste("Dispersion Component:",x$DispNames[i],"\n"))
        coefs <- x$Discoeff[[i]]
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
        na.print = "NA", ...)  
        cat("\n")
    }
    cat("\n\n")
    # Likelihood Values #
    cat("===== Likelihood Functions Value =====\n")
    cat("H-likelihood       :",x$HLikelihood,"\n")
    cat("Marginal likelihood:",x$APFix,"\n")
    cat("REML likelihood    :",x$APDis,"\n")
    cat("C-likelihood       :",x$CLikelihood,"\n")
}    
            


# What is needed is a function for diagnostic plots and for the numerical convergence check #
print.HGLM<-function(x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
   signif.stars = getOption("show.signif.stars"),plain=FALSE, ...){
   if (plain==TRUE) {
        class(x)<-NULL
        print(x)
        class(x)<-"HGLM"
   }
   else {
        mc<-x$CALL
        temp<-match("DistResp",names(mc))
        DistResp<-mc[[temp]]
        cat("===== HGLM Model Information =====\n")
        cat(paste("Response Distribution:",DistResp,"\n"))  
        namesRE<-x$NAMES$namesRE
        temp<-match("DistRand",names(mc))
        DistRand<-eval(mc[[temp]])
        for (i in 1:length(DistRand)){
            DistRandtemp<-DistRand[i]
            DistRandname<-namesRE[i]
            cat("Random Effect",i,":",DistRandtemp,"/",DistRandname,"\n")
        }
        temp<-match("Link",names(mc))
        LinkTemp<-mc[[temp]]
        cat("Link:",LinkTemp,"\n")
        if (x$Results$GradientResDisp[1]=="Not Estimated") ODTemp<-"Fixed"
        else ODTemp<-"Estimated"
        cat("Overdispersion Structure is",ODTemp,"\n")
        temp<-match("LapFix",names(mc))
        temp<-mc[[temp]]
        if (temp==TRUE) LapFix="Laplace Approximation"
        else LapFix="H-likelihood"
        cat("Estimation of fixed effects:",LapFix,"\n")
        temp<-match("DataMain",names(mc))
        temp<-mc[[temp]]
        cat("Dataset used:",temp,"\n")
        
        temp<-match("formulaMain",names(mc))
        temp<-mc[[temp]]
        cat("Model Equation:\n")
        cat("\t");print(temp);cat("\n")
        temp<-match("formulaOD",names(mc))
        temp<-mc[[temp]]
        cat("Overdispersion Equation:\n")
        cat("\t");print(temp);cat("\n")
        cat("Dispersion Equation(s):\n")
        temp<-match("formulaRand",names(mc))
        temp<-eval(mc[[temp]])
        for (i in 1:length(DistRand)){
            cat("Component",i,":");print(temp[[i]]);cat("\n")
        }       
   }
}

# Convergence Diagnostics #

HGLMLikeDeriv<-function(x){
    if (class(x)!="HGLM") stop("This function is for HGLM objects")
    cat("Gradient Fixed Structure:\n")
    print(x$Results$GradientFix)
    cat("Gradient Overdispersion Structure:\n")
    print(x$Results$GradientResDisp)
    cat("Gradient Dispersion:\n")
    print(x$Results$GradientDisp)
}


HGLMLRTest<-function(x1,x2){
    if (class(x1)!="HGLM" | class(x2)!="HGLM") stop("This function is for HGLM objects")
    x1<-x1$Results
    x2<-x2$Results
    
    LH1<-x1$HLikelihood
    LH2<-x2$HLikelihood
    
    LM1<-x1$APFix
    LM2<-x2$APFix
    
    LD1<-x1$APDis
    LD2<-x2$APDis
    
    if (LH1<LH2) cat("H-likelihood of model 2 is higher\n")
    if (LH1>LH2) cat("H-likelihood of model 1 is higher\n")
    if (LH1==LH2) cat("H-likelihoods are equal\n")
    
    DFM1<-length(x1$Beta)
    DFM2<-length(x2$Beta)
    
    DFMDIFF<-abs(DFM1-DFM2)
    LRM<-2*abs(LM1-LM2)
    pLRM<-1-pchisq(LRM,DFMDIFF)
    if (DFMDIFF==0) pLRM<-NA    

    DFD1<-nrow(x1$Dispersion)
    DFD2<-nrow(x2$Dispersion)
    
    DFDDIFF<-abs(DFD1-DFD2)
    LRD<-2*abs(LD1-LD2)
    pLRD<-1-pchisq(LRD,DFDDIFF)
    if (DFDDIFF==0) pLRD<-NA
	    
    cat("Marginal likelihood comparison:\n")
    cat("LR test p-value:",pLRM,"\n")
    cat("LR test statistics:",LRM,"\n")
    cat("LR difference df:",DFMDIFF,"\n")
    
    cat("REML likelihood comparison:\n")
    cat("LR test p-value:",pLRD,"\n")
    cat("LR test statistics:",LRD,"\n")
    cat("LR difference df:",DFDDIFF,"\n")
    
}
    

    
