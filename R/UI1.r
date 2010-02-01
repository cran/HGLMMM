################################################
##### First approach to a general function #####
################################################

# Description of the parameters : DistResp - Distribution of the Response
#                               : DistRand - List containing the distributions of the random components
#                               : Link     - Link for the Response
#                               : LapFix   - Specify whether to use the Laplace Approximation for the estimation of the fixed paramters
#                               : ODEst    - Should the Overdispersion parameter of the response be estimated
#                               : ODEstVal - (Starting) Value for the overdispersion parameter
#                               : formulaMain - formula specifying response and fixed effects 
#                               : formulaOD   - formula specifying the structure of the overdispersion
#                               : formulaRand - list of formulas specifying the structure of dispersion of each random component (length of the list should be the same as 
#                               :               length of the list DistRand
#                               : DataMain    - specify the data for the response dependence on covariates and overdispersion modelling
#                               : DataRand    - list of datas for modelling the dependence of dispersion parameters on covariates (random effects dependence) 
#                               : Offset      - specify Offset in Poisson model only
#                               : BinomialDen - Specify number of trials for bernoulli distributed response (assume 1 for binary)
#                               : StartBeta   - starting values for fixed effects
#                               : StartVs     - starting values for Vs
#                               : StartYgamma - starting values for overdispersion dependence parameters
#                               : StartRgamma - starting values for dispersion dependence parameters
#                               : INFO/DEGUBG - display additional information
#                               : na.action   - what to do with missing values as in lme4
#                               : contrasts   - works as in lme4 (figure it out what it does)

# At the moment make it work for a Normal HGLM as an example #
HGLMfit<-function(DistResp="Normal",DistRand=NULL,Link=NULL,LapFix=FALSE,ODEst=NULL,ODEstVal=0,formulaMain,formulaOD,formulaRand,DataMain,DataRand,Offset=NULL,BinomialDen=NULL,
            StartBeta=NULL,StartVs=NULL,StartRGamma=NULL,INFO=TRUE,DEBUG=FALSE,na.action,contrasts=NULL,CONV=1e-4){
    # take care of a binomial denominator #
    require(Matrix)
    require(numDeriv)
    if (Link=="Probit"|Link=="CLogLog") stop ("This link not implemented yet")
    if (DistResp=="Binomial" & is.null(BinomialDen)) stop("Binomial Denominator not specified")
    if (!is.null(BinomialDen)) BinomialDen<-as.vector(BinomialDen)
    # get number of random effects #
    if (is.null(DistRand)) stop("At least one random effect must be specified for the analysis - otherwise use function glm") 
    nrand<-length(DistRand)
    if (length(formulaRand)!=nrand) stop("Number of formulas for random effects is wrong")
    if (length(DataRand)!=nrand) stop("Number of datasets for random effects dispersion is wrong")
       
    mc <- match.call()
    stopifnot(length(formula <- as.formula(formulaMain)) == 3)
    for(i in 1:nrand) stopifnot(length(formula <- as.formula(formulaRand[[i]])) == 2)
    
    if ((DistResp=="Normal" | DistResp=="Gamma") & is.null(ODEst)) ODEst<-TRUE
    if ((DistResp=="Poisson" | DistResp=="Binomial") & is.null(ODEst)) ODEst<-FALSE
    
    if (ODEst==TRUE){
        stopifnot(length(formula <- as.formula(formulaOD)) == 2)
    }
    # Frame for the model of the response #
    fr<-HGLMFrames(mc,formulaMain,contrasts)
    namesX<-names(fr$fixef)
    namesY<-names(fr$mf)[1]

   # Define the Frame for the Overdispersion modelling #
    frOD<-HGLMFramesOD(mc,formulaOD)
    namesOD<-names(frOD$fixef)
    
    # Define the Frame for each Dispersion modelling #
    frRD<-list(0)
    frRDX<-list(0)
    namesRD<-list(0)
    for (i in 1:length(DistRand)){
        frRD[[i]]<-HGLMFramesRD(mc,formulaRand[[i]],whichDispComp=i)
        frRDX[[i]]<-frRD[[i]]$X
        namesRD[[i]]<-names(frRD[[i]]$fixef)
    }  
     
    # Creating the design for random effects = this returns design for random effects and subject indexes#
    FL<-HGLMFactorList(formulaMain,fr,0L,0L) 

    namesRE<-FL$namesRE
    # Now apply the above to the normal HGLM #
    # Extract Y #
    Y<-matrix(fr$Y,length(fr$Y),1)
    X<-fr$X
    ZZ<-FL$Design
    SS<-FL$Subject
    # Generate starting values for Beta #
    # Maybe needs an adjustment for gamma models #
    if (!is.null(StartBeta)) Beta<-StartBeta
    else {
        Beta<-rep(0,ncol(X))
        for (i in 1:ncol(X)){
            if (sum(!(X[,i]==1))==0) {
                if (is.null(Link)) {
                    if (DistResp=="Normal") Beta[i]<-mean(Y)
                    if (DistResp=="Gamma") Beta[i]<--1/mean(Y)
                    if (DistResp=="Poisson") Beta[i]<-log(mean(Y))
                    if (DistResp=="Binomial") Beta[i]<-log(mean(Y/BinomialDen)/(1-mean(Y/BinomialDen)))
                }
                if (DistResp!="Binomial"){
                    if (Link=="Identity") Beta[i]<-mean(Y)
                    if (Link=="Log") Beta[i]<-log(mean(Y))
                    if (Link=="Inverse") Beta[i]<--1/mean(Y)
                    if (Link=="Logit") Beta[i]<-log(mean(Y)/(1-mean(Y)))
                    if (Link=="Probit") Beta[i]<-qnorm(mean(Y))
                    if (Link=="CLogLog") Beta[i]<-log(-log(1-mean(Y)))
                }
                if (DistResp=="Binomial"){
                    if (Link=="Identity") Beta[i]<-mean(Y/BinomialDen)
                    if (Link=="Log") Beta[i]<-log(mean(Y/BinomialDen))
                    if (Link=="Inverse") Beta[i]<--1/mean(Y/BinomialDen)
                    if (Link=="Logit") Beta[i]<-log(mean(Y/BinomialDen)/(1-mean(Y/BinomialDen)))
                    if (Link=="Probit") Beta[i]<-qnorm(mean(Y/BinomialDen))
                    if (Link=="CLogLog") Beta[i]<-log(-log(1-mean(Y/BinomialDen)))
                }                
                
                
            }    
        }
    }
    # Generate starting values for Vs #
    # Compute dimension of random effects #
    nrand<-length(ZZ)
    q<-rep(0,nrand)
     for (i in 1:nrand) q[i]<-dim(ZZ[[i]])[2]
    
    # Defining cumulative dimensions #
    qcum<-cumsum(c(0,q))
    if (!is.null(StartVs)) Vs<-StartVs
    else {
        Vs<-rep(0,qcum[nrand+1])
        for (i in 1:nrand){
            if (DistRand[i]=="IGamma") Vs[(qcum[i]+1):qcum[i+1]]<--0.01
        }
    }
    
    # keep Offset and Link and LapFix and DistRand as it is in the argument #
    # estimate overdispersion is defined above ODEst #
    # Overdispersion specification #
    DDY<-frOD$X
    
    # Starting value for the OD or fixed value #
    DYgamma<-ODEstVal

    
    # Create DDR matrix #
    DDR<-as.matrix(bdiag(frRDX))
    if (!is.null(StartRGamma)) DRgamma<-StartRGamma
    else DRgamma<-rep(0,ncol(DDR))
    
    # Define default link #
    if (is.null(Link)){
        if (DistResp=="Normal")     Link<-"Identity"
        if (DistResp=="Poisson")    Link<-"Log"
        if (DistResp=="Gamma")      Link<-"Inverse"
        if (DistResp=="Binomial")   Link<-"Logit"
    }
    
    
    # Select a function to handle response #
    #print("frOD");print(frOD)
    #print("Y");print(Y)
    #print("X");print(X)
    #print("ZZ");print(ZZ)
    #print("SS");print(SS)
    #print("DDR");print(DDR)
    #print("DDY");print(DDY)
    #print("Beta");print(Beta)
    #print("Vs");print(Vs)
    #print("Offset");print(Offset)
    #print("Link");print(Link)
    #print("EstimateOverDisp");print(ODEst)
    #print("LaplaceFixed");print(LapFix)
    #print("RandDist");print(DistRand)
    #print("DYgamma");print(DYgamma)
    #print("DRgamma");print(DRgamma)
    #print("frRD");print(frRD)
    
    if (DistResp=="Normal"){          
        OUTPUT<-NormalHGLM.fit(YP=Y,XP=X,ZZP=ZZ,SSP=SS,BetaP=Beta,VstartP=Vs,OFFSET=Offset,EstimOD=ODEst,LapFixP=LapFix,RDP=DistRand,Link=Link,
                        DDRP=DDR,DDYP=DDY,DYgammaP=DYgamma,DRgammaP=DRgamma,Info=INFO,DEBUG=DEBUG,CONV=CONV)
    }  
    if (DistResp=="Poisson"){
        OUTPUT<-PoissonHGLM.fit(YP=Y,XP=X,ZZP=ZZ,SSP=SS,BetaP=Beta,VstartP=Vs,OFFSET=Offset,EstimOD=ODEst,LapFixP=LapFix,RDP=DistRand,Link=Link,
                        DDRP=DDR,DDYP=DDY,DYgammaP=DYgamma,DRgammaP=DRgamma,Info=INFO,DEBUG=DEBUG,CONV=CONV)
    }
    if (DistResp=="Binomial"){
        OUTPUT<-BinomialHGLM.fit(YP=Y,XP=X,ZZP=ZZ,SSP=SS,B=BinomialDen,BetaP=Beta,VstartP=Vs,OFFSET=Offset,EstimOD=ODEst,LapFixP=LapFix,RDP=DistRand,Link=Link,
                        DDRP=DDR,DDYP=DDY,DYgammaP=DYgamma,DRgammaP=DRgamma,Info=INFO,DEBUG=DEBUG,CONV=CONV)
    }
    if (DistResp=="Gamma"){
        OUTPUT<-GammaHGLM.fit(YP=Y,XP=X,ZZP=ZZ,SSP=SS,BetaP=Beta,VstartP=Vs,OFFSET=Offset,EstimOD=ODEst,LapFixP=LapFix,RDP=DistRand,Link=Link,
                        DDRP=DDR,DDYP=DDY,DYgammaP=DYgamma,DRgammaP=DRgamma,Info=INFO,DEBUG=DEBUG,CONV=CONV)
    }
    class(OUTPUT)<-"HGLM"
    OUTPUT$NAMES<-list(namesY=namesY,namesX=namesX,namesOD=namesOD,namesRD=namesRD,namesRE=namesRE,indexRE=lapply(FL$Subject,levels))
    OUTPUT$CALL<-mc
    OUTPUT
    
}



HGLMFactorList<-function (formula, fr, rmInt, drop) 
{
    mf <- fr$mf
    bars <- expandSlash(findbars(formula[[3]]))
    for (i in 1:length(bars)){
        checkcorr<-findplus(bars[[i]])
        if (checkcorr==1) stop("Correlated random effects are not currently allowed in the HGLM routines")
        if (checkcorr==-1) stop("You do not need to specify '-1' for no intercept it is done be default")
    }
    

    if (!length(bars)) 
        stop("No random effects terms specified in formula")

    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))

       
    fl <- lapply(bars, function(x) {
        ff <- eval(substitute(as.factor(fac)[, drop = TRUE], 
            list(fac = x[[3]])), mf)
        im <- as(ff, "sparseMatrix")
        if (!isTRUE(validObject(im, test = TRUE))) 
            stop("invalid conditioning factor in random effect: ", 
                format(x[[3]]))
        
        if (is.name(x[[2]])){
        tempexp<-paste("~",as.character(x[[2]]),"-1")
        tempexp<-as.formula(tempexp)[[2]]
        }
        else tempexp<-x[[2]]
        mm <- model.matrix(eval(substitute(~expr, list(expr = tempexp))), 
            mf)    
            
        if (rmInt) {
            if (is.na(icol <- match("(Intercept)", colnames(mm)))) 
                break
            if (ncol(mm) < 2) 
                stop("lhs of a random-effects term cannot be an intercept only")
            mm <- mm[, -icol, drop = FALSE]
        }
        ans <- list(f = ff, A = do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) im)), Zt = do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) {
                im@x <- mm[, j]
                im
            })), ST = matrix(0, ncol(mm), ncol(mm), dimnames = list(colnames(mm), 
            colnames(mm))))
        if (drop) {
            ans$A@x <- rep(0, length(ans$A@x))
            ans$Zt <- drop0(ans$Zt)
        }
        ans
    })
   
 
##### Skip the below code for now #####   
#    dd <- VecFromNames(dimsNames, "integer", c(list(n = nrow(mf), 
#        p = ncol(fr$X), nt = length(fl), q = sum(sapply(fl, function(el) nrow(el$Zt)))), 
#        dimsDefault))
         
       
#    nlev <- sapply(fl, function(el) length(levels(el$f)))
#    if (any(diff(nlev)) > 0) 
#        fl <- fl[rev(order(nlev))]
#    trms <- lapply(fl, "[", -1)
#    names(trms) <- NULL
#    fl <- lapply(fl, "[[", "f")
#    attr(fl, "assign") <- seq_along(fl)
#    fnms <- names(fl)
#    if (length(fnms) > length(ufn <- unique(fnms))) {
#        fl <- fl[match(ufn, fnms)]
#        attr(fl, "assign") <- match(fnms, ufn)
#    }
#    names(fl) <- ufn
#    dd["nest"] <- all(sapply(seq_along(fl)[-1], function(i) isNested(fl[[i - 
#        1]], fl[[i]])))
#    list(trms = trms, fl = fl, dims = dd)
Design<-list(0)
Subject<-list(0)
for (i in 1:length(fl)){
    Subject[[i]]<-as.factor(fl[[i]]$f)
    tempmat<-fl[[i]]$Zt
    tempmat<-as.matrix(t(tempmat))
    Design[[i]]<-tempmat
    }
list(Design=Design,Subject=Subject,namesRE=names(bars))
}

dimsNames<-c("nt","n","p","q","s","np","LMM","REML","fTyp","lTyp","vTyp","nest","useSc","nAGQ","verb","mxit","mxfn","cvg")
dimsDefault<-structure(list(s = 1L, mxit = 300L, mxfn = 900L, verb = 0L, np = 0L, 
    LMM = 0L, REML = 0L, fTyp = 2L, lTyp = 5L, vTyp = 1L, useSc = 1L, 
    nAGQ = 1L, cvg = 0L), .Names = c("s", "mxit", "mxfn", "verb", 
"np", "LMM", "REML", "fTyp", "lTyp", "vTyp", "useSc", "nAGQ", 
"cvg"))


expandSlash<-function (bb) 
{
    if (!is.list(bb)) 
        return(expandSlash(list(bb)))
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]]))) 
            return(lapply(unlist(makeInteraction(trms)), function(trm) substitute(foo | 
                bar, list(foo = x[[2]], bar = trm))))
        x
    }))
}

findbars<-function(term) 
{
    if (is.name(term) || !is.language(term)) 
        return(NULL)
    if (term[[1]] == as.name("(")) 
        return(findbars(term[[2]]))
    if (!is.call(term)) 
        stop("term must be of class call")
    if (term[[1]] == as.name("|")) 
        return(term)
    if (length(term) == 2) 
        return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

##### Check if there were two terms entered in a correlated fashion #####
##### Current software does not allow for correlated random effects #####
findplus<-function(term)
{
    if (is.numeric(term))
        return(0)        
    if (!is.language(term)) 
        return(NULL)
    if (length(term)==1) return(0)
    if (term[[1]] == as.name("|"))
        return(findplus(term[[2]]))
    if (!is.call(term))
        stop("term must be of class call")
    if (term[[1]] == as.name("+"))
        return(1)
    if (term[[1]] == as.name("-"))
        return(-1)
}
    
###########################################################################
###########################################################################
    
slashTerms<-function (x) 
{
    if (!("/" %in% all.names(x))) 
        return(x)
    if (x[[1]] != as.name("/")) 
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}


HGLMFrames<-function (mc, formula, contrasts, vnms = character(0)) 
{
    mf <- mc
    m <- match(c("DataMain", "weights", "na.action", "offset"), 
        names(mf), 0)

    mf <- mf[c(1, m)] 

    frame.form <- subbars(formula)
    if (length(vnms) > 0) 
        frame.form[[3]] <- substitute(foo + bar, list(foo = parse(text = paste(vnms, 
            collapse = " + "))[[1]], bar = frame.form[[3]]))
    

    fixed.form <- nobars(formula)
    
    if (inherits(fixed.form, "name")) 
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))
    environment(fixed.form) <- environment(frame.form) <- environment(formula)
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE

    mf[[1]] <- as.name("model.frame")
    names(mf)[2]<-"data"
    fe <- mf
    mf <- eval(mf, parent.frame(2))
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2))
    fe
    Y <- model.response(mf, "any")  
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    mt <- attr(fe, "terms")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf,contrasts) # work out how the contrasts work 
    else matrix(, NROW(Y), 0)
    storage.mode(X) <- "double"

    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL
    wts <- model.weights(mf)
    if (is.null(wts)) 
        wts <- numeric(0)
    off <- model.offset(mf)
    if (is.null(off)) 
        off <- numeric(0)
    if (any(wts <= 0)) 
        stop(gettextf("negative weights or weights of zero are not allowed"))
    if (length(off) && length(off) != NROW(Y)) 
        stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
            length(off), NROW(Y)))
    attr(mf, "terms") <- mt
    list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), 
        mf = mf, fixef = fixef)
}

HGLMFramesOD<-function (mc, formula) 
{
    mf <- mc
    m <- match(c("DataMain","na.action"), 
        names(mf), 0)

    mf <- mf[c(1, m)] 

    #frame.form <- subbars(formula)
    #if (length(vnms) > 0) 
    #    frame.form[[3]] <- substitute(foo + bar, list(foo = parse(text = paste(vnms, 
    #        collapse = " + "))[[1]], bar = frame.form[[3]]))
    

    fixed.form <- nobars(formula)
    fixed.form
    if (inherits(fixed.form, "name")) 
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))
    environment(fixed.form) <- environment(formula)
    mf$formula <- fixed.form
    mf$drop.unused.levels <- TRUE

    mf[[1]] <- as.name("model.frame")
    names(mf)[2]<-"data"
    fe <- mf
    mf <- eval(mf, parent.frame(2))
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2))
    # Y <- model.response(mf, "any")  
    #if (length(dim(Y)) == 1) {
    #    nm <- rownames(Y)
    #    dim(Y) <- NULL
    #    if (!is.null(nm)) 
    #        names(Y) <- nm
    #}
    mt <- attr(fe, "terms")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf) # work out how the contrasts work 
    else matrix(1,NROW(Y), 1)
    storage.mode(X) <- "double"

    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL
 
    attr(mf, "terms") <- mt
    list(X = X, mf = mf, fixef = fixef)
}

HGLMFramesRD<-function (mc, formula, whichDispComp) 
{
    mf <- mc
    m <- match(c("DataRand","na.action"), 
        names(mf), 0)
    mf <- mf[c(1, m)] 

    mf[[2]]<-mf[[2]][[whichDispComp+1]]
    
    #frame.form <- subbars(formula)
    #if (length(vnms) > 0) 
    #    frame.form[[3]] <- substitute(foo + bar, list(foo = parse(text = paste(vnms, 
    #        collapse = " + "))[[1]], bar = frame.form[[3]]))
    

    fixed.form <- nobars(formula)
    if (inherits(fixed.form, "name")) 
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))
    environment(fixed.form) <- environment(formula)
    mf$formula <- fixed.form
    mf$drop.unused.levels <- TRUE

    mf[[1]] <- as.name("model.frame")
    names(mf)[2]<-"data"
    fe <- mf
    mf <- eval(mf, parent.frame(2))
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2))
    # Y <- model.response(mf, "any")  
    #if (length(dim(Y)) == 1) {
    #    nm <- rownames(Y)
    #    dim(Y) <- NULL
    #    if (!is.null(nm)) 
    #        names(Y) <- nm
    #}
    mt <- attr(fe, "terms")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf) # work out how the contrasts work 
    else matrix(1,NROW(Y), 1)
    storage.mode(X) <- "double"

    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL
 
    attr(mf, "terms") <- mt
    list(X = X, mf = mf, fixef = fixef)
}

subbars<-function (term) 
{
    if (is.name(term) || !is.language(term)) 
        return(term)
    if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name("|")) 
        term[[1]] <- as.name("+")
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

nobars<-function (term) 
{
    if (!("|" %in% all.names(term))) 
        return(term)
    if (is.call(term) && term[[1]] == as.name("|")) 
        return(NULL)
    if (length(term) == 2) {
        nb <- nobars(term[[2]])
        if (is.null(nb)) 
            return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) 
        return(nb3)
    if (is.null(nb3)) 
        return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

VecFromNames<-function (nms, mode = "numeric", defaults = list()) 
{
    ans <- vector(mode = mode, length = length(nms))
    names(ans) <- nms
    ans[] <- NA
    if ((nd <- length(defaults <- as.list(defaults))) > 0) {
        if (length(dnms <- names(defaults)) < nd) 
            stop("defaults must be a named list")
        stopifnot(all(dnms %in% nms))
        ans[dnms] <- as(unlist(defaults), mode)
    }
    ans
}

isNested<-function (f1, f2) 
{
    f1 <- as.factor(f1)
    f2 <- as.factor(f2)
    stopifnot(length(f1) == length(f2))
    sm <- as(new("ngTMatrix", i = as.integer(f2) - 1L, j = as.integer(f1) - 
        1L, Dim = c(length(levels(f2)), length(levels(f1)))), 
        "CsparseMatrix")
    all(diff(sm@p) < 2)
}
