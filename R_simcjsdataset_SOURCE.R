#
require(mvtnorm)

sim.cov.dat <- function(T,n,N.cont=5, N.factor=5, time.as.factor=TRUE,method_control=list(method="sparse",betanorm=list(p=3,s=3), taper_decay=list(p=c(a=0.4,0.77),s=c(a=0.4,0.77)), sparse_Ntruecov=list(p=3,s=3)),nl_control=list(nonlinear=FALSE,n.time.interactions=0,cor_strength=8.3),print_=FALSE,add.t.at.zero=TRUE,intercept_range=list(p=c(min=0.4,max=0.7),s=c(min=0.55,max=0.85))){
#    T=7;n=200;N.cont=5; N.factor=5; time.as.factor=TRUE;method_control=list(method="sparse",betanorm=list(p=3,s=3), taper_decay=list(p=c(a=0.4,0.77),s=c(a=0.4,0.77)), sparse_Ntruecov=list(p=3,s=3));nl_control=list(nonlinear=FALSE,interactions=FALSE,n.time.interactions=0,cor_strength=8.3)# setting cor_strength to 2 means no correlation
    N.tot =N.cont+N.factor # total number of covariates
    nlatent=1.4;rtdf=2;
    make_sigma=exp(-as.matrix(dist(matrix(rt(round(nlatent*N.tot)*N.tot,2),N.tot,round(N.tot*nlatent))))/nl_control$cor_strength) # median ~ 0.3,max~0.55,min~0.01-0.1
    xx =apply(rmvnorm(n=n,mean=rep(0,N.tot),sigma=make_sigma),2,function(x) x-mean(x)) # simulated correlated covariates
        print(summary(make_sigma[lower.tri(N.tot,diag=FALSE)]))
    xx.fact <- apply(xx[,1:N.factor,drop=FALSE],2, function(z){nf=sample(c(2,3,4),1); do.call(what="cut",args=list(x=z,breaks=quantile(z,seq(0,1,1/nf)),labels=1:nf,include.lowest=TRUE))})
    # make a LONG (x(T-1)) data.frame
    xx.wide<-cbind(xx[,-(1:N.factor)],as.data.frame(xx.fact))
    xxx=as.data.frame(cbind(data.frame(time=factor(rep(2:T,n),levels=2:T)),apply(xx[,-(1:N.factor)], 2,rep,each=(T-1)),apply(xx.fact, 2,rep,each=(T-1))))
    names(xxx) <- c("time",paste0(letters[1:N.tot]))
    names(xx.wide) <- paste0(letters[1:N.tot])
    intercepts=logit(c(p=do.call("runif",c(list(n=1),intercept_range$p)),s=do.call("runif",c(list(n=1),intercept_range$s))))
    beta<-list(p=rep(0,N.tot+1),s=rep(0,N.tot+1))
       if(method_control$method=="sparse"){ # sparsity
           beta$p[sample(1:(N.tot+1),method_control$sparse_Ntruecov$p)]<-method_control$betanorm$p/method_control$sparse_Ntruecov$p*sample(c(-1,1),method_control$sparse_Ntruecov$p,replace=TRUE)
           beta$s[sample(1:(N.tot+1),method_control$sparse_Ntruecov$s)]<-method_control$betanorm$s/method_control$sparse_Ntruecov$s*sample(c(-1,1),method_control$sparse_Ntruecov$s,replace=TRUE) # end sparse
       } else if(method_control$method=="tapering"){
           beta$p[sample(1:(N.tot+1))]<- method_control$betanorm$p*do.call(what=function(x) x/sum(x),args=list(x=method_control$taper_decay$p[1]*method_control$taper_decay$p[2]^(0:N.tot)))
           beta$s[sample(1:(N.tot+1))]<- method_control$betanorm$s*do.call(what=function(x) x/sum(x),args=list(x=method_control$taper_decay$s[1]*method_control$taper_decay$s[2]^(0:N.tot)))
       } # end effect tapering
       # begin making model matrice
       mm <-vector(length=N.tot+1,mode="list")
       mm[[1]]<- model.matrix(~time-1,data=xxx) # mm for time variable
       beta.mm <- list(p=vector(length=N.tot+1,mode="list"),s=vector(length=N.tot+1,mode="list"))
       beta.mm[["p"]][[1]]= intercepts["p"]+2*qnorm(p=c(0.70))*(0.5*seq(-1*beta$p[[1]],beta$p[[1]],length.out=ncol(mm[[1]])+1))[-ceiling(ncol(mm[[1]])/2)]
       beta.mm[["s"]][[1]]= intercepts["s"]+2*qnorm(p=c(0.70))*(0.5*seq(-1*beta$s[[1]],beta$s[[1]],length.out=ncol(mm[[1]])+1))[-ceiling(ncol(mm[[1]])/2)]
    names(beta.mm[["s"]][[1]])<-names(beta.mm[["p"]][[1]]) <- colnames(mm[[1]])    
       if(!nl_control$nonlinear){ # linear effects
           for(i_ in 2:(N.tot+1)){
               x=xxx[,i_,drop=FALSE]
               if(is.factor(x[,1])){ # factors 
                   mm[[i_]]=model.matrix(~.,data=xxx[,i_,drop=FALSE])[,-1,drop=FALSE]
                   beta.mm[["p"]][[i_]]= 2*qnorm(p=c(0.05))*(0.5*seq(-1*beta$p[[i_]],beta$p[[i_]],length.out=ncol(mm[[i_]])+1))[-ceiling(ncol(mm[[i_]])/2)]*2^(ncol(mm[[i_]])==1)
                   beta.mm[["s"]][[i_]]= 2*qnorm(p=c(0.05))*(0.5*seq(-1*beta$s[[i_]],beta$s[[i_]],length.out=ncol(mm[[i_]])+1))[-ceiling(ncol(mm[[i_]])/2)]*2^(ncol(mm[[i_]])==1)
               }else if(is.numeric(x[,1])){
                   mm[[i_]]=model.matrix(~.,data=xxx[,i_,drop=FALSE])[,-1,drop=FALSE]
                   beta.mm[["p"]][[i_]]<-beta$p[[i_]]
                   beta.mm[["s"]][[i_]]<-beta$s[[i_]]
               }
               names(beta.mm[["s"]][[i_]])<-names(beta.mm[["p"]][[i_]]) <- colnames(mm[[i_]])
           } # i_
       } else if(nl_control$nonlinear){ # polynomials
           for(i_ in 2:(N.tot+1)){
               x=xxx[,i_,drop=FALSE]
               if(is.factor(x[,1])){
                   mm[[i_]]=model.matrix(~.,data=xxx[,i_,drop=FALSE])[,-1,drop=FALSE]
                   beta.mm[["p"]][[i_]]= (0.5*seq(-1*beta$p[[i_]],beta$p[[i_]],length.out=ncol(mm[[i_]])+1))[-ceiling(ncol(mm[[i_]])/2)]*2^(ncol(mm[[i_]])==1)
                   beta.mm[["s"]][[i_]]= (0.5*seq(-1*beta$s[[i_]],beta$s[[i_]],length.out=ncol(mm[[i_]])+1))[-ceiling(ncol(mm[[i_]])/2)]*2^(ncol(mm[[i_]])==1)
               }else if(is.numeric(x[,1])){
                   x[,paste0(names(x),"^2")]=as.numeric(scale((x[,1]-runif(1,min(x[,1]),max(x[,1])))^2,center=TRUE))
                   mm[[i_]]=model.matrix(~.,data=x)[,-1,drop=FALSE]
                   beta.mm[["p"]][[i_]]<-(c(0,1)+c(1,-1)*runif(1,0.25,0.75))*beta$p[[i_]]
                   beta.mm[["s"]][[i_]]<-(c(0,1)+c(1,-1)*runif(1,0.25,0.75))*beta$s[[i_]]
               }
               names(beta.mm[["s"]][[i_]])<-names(beta.mm[["p"]][[i_]]) <- colnames(mm[[i_]])               
           } # i_
           #for(j in 2:6){plot(mm[[j]][,1], mm[[j]]%*%c(0.7,0.3))}; plot(1,1,typ="n")
       } # end if non-linear
    # resemble design.matrix and betas
    modmat <- do.call("cbind",mm)
    beta.true <- list(p=do.call("c",beta.mm[["p"]]),s=do.call("c",beta.mm[["s"]]))
    true.long <- list(p=modmat%*%beta.true[["p"]],s=modmat%*%beta.true[["s"]])
    if(print_){
        par(mfrow=c(3,4))
        for(i in 1:(N.tot+1)){
            plot(xxx[rank(xxx[,i],ties.method="first"),i],true.long[["p"]][rank(xxx[,i],ties.method="first")])
        }
        for(i in 1:(N.tot+1)){
            plot(xxx[rank(xxx[,i],ties.method="first"),i],true.long[["s"]][rank(xxx[,i],ties.method="first")])
        }
    }
    if(add.t.at.zero){
        t1_ix <-rep(2:T,n)+T*rep(0:(n-1),each=(T-1))
        xxx.2 <- rbind(xxx,xxx)[1:(nrow(xxx)+n),]
        nt1_ix <-which(1:nrow(xxx.2)%in%t1_ix==FALSE)
        t1_ix <-rep(2:T,n)+T*rep(0:(n-1),each=(T-1))
        nt1_ix <-which(1:nrow(xxx.2)%in%t1_ix==FALSE)
        xxx.2[t1_ix,]<-xxx; xxx.2[nt1_ix,]<-NA; xxx=xxx.2
    #    xxx.2[nt1_ix,"time"]<-1
        modmat.2 <- rbind(modmat,modmat)[1:(nrow(modmat)+n),]
        modmat.2[t1_ix,]<-modmat; modmat.2[nt1_ix,]<-NA;modmat=modmat.2
        true.long2=true.long;true.long2[[2]]<-true.long2[[1]]<-matrix(rep(NA,nrow(xxx.2))); true.long2[[1]][t1_ix,1]<-true.long[[1]]; true.long2[[2]][t1_ix,1]<-true.long[[2]];
        true.long=true.long2
    }
    return(list(truelong=true.long,mm=modmat,beta=beta.true,datalong=xxx,datawide=xx.wide))}

# done simulating FAKE COVARIATE DATA
# SIMULATION CHAPTURE-RECAPTURE RESPONSE
sim.ch <- function(T,p,s,n.target,entry.probs=c(0.4,0.6*rep(1/(T-2),T-2),0)){
  n.start <- length(p)/T;
  if(any(is.na(p))) p[is.na(p)]<-1 # if they enter at time 1, they are seen
  if(any(is.na(s))) s[is.na(s)]<-1
  attempts <- 0
  n.realized <- 0
  first_ <- rep(1:T,n.start)==1
  id_ <- rep(1:n.start,each=T)
  while((n.realized < n.target) & (attempts < 11)){
      deaths <- 1*((runif(length(s))>s))
      foo<-do.call("rbind",tapply(deaths,INDEX=id_,FUN=function(x,probs) {
          histor <- numeric(T)
          enter=sample(1:T,1,prob=entry.probs)
          # enter+1: this means we condition on entry
          histor[enter:T] <- c(1,1*(cumsum(x[(enter+1):T])==0))
          return(histor)
      },probs=entry.probs))
      # make sure they are seen on first capture/entry
      first(foo)
      ch <- foo*matrix(runif(length(p))<=p,n.start,T,byrow=TRUE)
      obsvec <-apply(ch,1,function(x) any(x[-length(x)]==1))
      n.realized <- sum(obsvec)
      attempts <- attempts +1
  }
if(attempts>10){ print(paste("oops, tried 10 times to meet at least",n.target,"and failed"))}  
  if(n.realized>n.target){
      discard<-which((1:n.realized)%in% sort(sample(which(obsvec),n.target))==FALSE)
      obsvec[discard]<-FALSE
  }
return(list(ch=ch,obs=obsvec))}

## wrapper function to simulate BOTH covariate and capture histoires
sim.cjs <- function(N.cont, # number of continuous covariates
                    N.factor,# number of categorical covariates
                    n.start, # initial N individuals, for random generation (latter kulled to n.start)
                    n.target, # final number of individuals (n.start > n.target)
                    T, # number of capture periods
                    entry.probs=c(0.4,0.6*rep(1/(T-2),T-2),0), # probs, to randomly generate first-captures periods
                    time.as.factor=TRUE,# include time-as-factor
                    method_control=list(method="sparse",betanorm=list(p=3,s=3),taper_decay=list(p=c(a=0.4,0.77),s=c(a=0.4,0.77)), sparse_Ntruecov=list(p=3,s=3)),# 'sparse' is just a few important covariates; taper is they all have a small effect. Controls how many are impactedful
                    nl_control=list(nonlinear=FALSE,n.time.interactions=0,cor_strength=8.3),
                    add.t.at.zero=TRUE,
                    intercept_range=list(p=c(min=0.4,max=0.7),s=c(min=0.55,max=0.85)),                    
                    print_=FALSE # print/plot the randomly generated covariates
                    ){
    if(n.target>n.start){ stop("n.target must be <= n.start")}
    n.target <- sample(c(200:300),1)# 300 ORIGINALLY I USED 300 individuals
  # simulate a high dimensional dataset
    dd<- sim.cov.dat(T=T,n=n.start,N.cont=N.cont, N.factor=N.factor, time.as.factor=time.as.factor,method_control=method_control,nl_control=nl_control,print_=print_,add.t.at.zero=add.t.at.zero,intercept_range=intercept_range)
    # find the significant effects
    sig.effs<-lapply(dd$beta, function(x) {
            effect.names= unique(gsub("[[:digit:]]{1,}","",names(x)))
            sig.eff= matrix(0,length(effect.names));row.names(sig.eff)=effect.names
            for(i_ in effect.names){wheff<- which(names(x)%in%paste0(i_,c("",1:100)))
                if(all(x[wheff]==0) | (length(unique(x[wheff]))==1 & length(x[wheff])!=1)){
                    sig.eff[i_,]=0
                } else { sig.eff[i_,]=1}
                                }
            return(sig.eff)})
    # simulate the capture histories
    simdat <- sim.ch(T=T,p=1/(1+exp(-dd$truelong$p)),s=1/(1+exp(-dd$truelong$s)),n.target=n.target,entry.probs=entry.probs)
    ch.wide.obs <- simdat$ch[which(simdat$obs),] # a simple matrix of width T, length n, only animals observed at least once
    colnames(ch.wide.obs)<-paste(1:T)
    cov.dat.wide <- dd$datawide[which(simdat$obs),];cov.dat.wide$id <- row.names(cov.dat.wide)
    cov.dat.wide.means <- lapply(cov.dat.wide,function(x) if(is.numeric(x)){mean(x)} else {0})
                                        # mean center the covariates!!!
    cov.dat.wide[,which(unlist(lapply(cov.dat.wide,is.numeric)))] <- apply(cov.dat.wide[,which(unlist(lapply(cov.dat.wide,is.numeric)))],2,function(x) x-mean(x))
    simdatobs.long <- which(rep(simdat$obs,each=T))
    cov.dat.long <-dd$datalong[simdatobs.long,]
    cov.dat.long$id <- rep(row.names(ch.wide.obs),each=T)
    p.true.obs <- dd$truelong$p[simdatobs.long]
    s.true.obs <- dd$truelong$s[simdatobs.long]
    first_ <- first(ch.wide.obs)
    true.wide = list(s = matrix(s.true.obs,nrow(ch.wide.obs),T,byrow=TRUE),p = matrix(p.true.obs,nrow(ch.wide.obs),T,byrow=TRUE))
    # delete all covariate information before first capture
    for(i in 1:length(first_)){
        for(par_ in 1:length(true.wide)){
            true.wide[[par_]][i,1:first_[i]]<-NA
        }
    }
    # remove intercept from betas
    beta=dd$beta
    intercepts = lapply(beta, function(x) mean(x[grep("time",names(x))]))
    beta = mapply(beta,intercepts, FUN=function(b_,i_){ c(b_[grep("time",names(b_))]-i_, b_[grep("time",names(b_),invert=TRUE)],interc=i_)},SIMPLIFY=FALSE)
    return(list(ch.data = ch.wide.obs,beta=beta,mm=dd$mm,cov.data.wide=cov.dat.wide,true.wide=true.wide,sig.effs=sig.effs))
} # end function sim.cjs
