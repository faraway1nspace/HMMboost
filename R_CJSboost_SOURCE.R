# I developed and tested these algorithms in the document R_CJS_likelihood_demo.R
# there is a CJS likelihood and a two.slice.marginal whatever
require(mboost)
require(Rcpp)
require(RcppArmadillo)
require(inline)

# find first capture
first <- function(wide) apply(wide,1,function(x) min(which(as.logical(x))))

# convert a Mark ch into a data.frame 
ch2wide <- function(ch,col=1,col.stem="t"){
    ncol_ <- nchar(ch[1,col])
    structure(.Data=apply(do.call("rbind",strsplit(ch[,col],split="")),2,as.numeric),dim=c(nrow(ch),ncol_),dimnames=list(1:nrow(ch),paste0(col.stem,1:ncol_)))
}

# prediction to new data
boost_predict <- function(mod, newdata, agg=TRUE, mstop=NULL, offsets = NULL){
      if(is.null(mstop)) mstop <- mod$mstop
      if(is.null(offsets)) offsets <- mod$offsets
      preds <- list()
	# RETURN AGGREGATE (SINGLE VECTOR) RESULTS
	if(agg) {
            for(comp in 1:length(mod$formula)){
                bl <- mod$bl[[comp]]; xselect <- mod$xselect[[comp]][1:mstop]; ens <- mod$ens[[comp]][1:mstop]; nu <- mod$nu[[comp]]; n_iter <- length(xselect)
                pfun <- function(w) {
                    ix <- which(xselect == w & (1:n_iter) <= mstop)
                    ret <- nu * bl[[w]]$predict(ens[ix], newdata = newdata, agg="sum")
                    ret
                }
                pr <- lapply(unique(xselect), pfun)
                for (i in 1:length(pr)) pr[[i]] <- as(pr[[i]], "matrix")
                pr[[length(pr)+1]] <- matrix( offsets[[comp]] , nrow=nrow(newdata),ncol=1)
                preds[[comp]] <- Reduce("+",pr)
            }
            names(preds) <- names(mod$formula)
            return(preds)
	# RETURN UN-AGGREGATED (per m) RESULTS
	} else {
            for(comp in 1:length(mod$formula)){
                bl <- mod$bl[[comp]]; xselect <- mod$xselect[[comp]][1:mstop]; ens <- mod$ens[[comp]][1:mstop]; nu <- mod$nu[[comp]]; n_iter <- length(xselect)
                pfun <- function(w) {
                    ix <- which(xselect == w & (1:n_iter) <= mstop)
                    ret <- nu * bl[[w]]$predict(ens[ix], newdata = newdata, agg="none")
                    m <- Matrix(0, nrow = nrow(newdata), ncol = sum((1:n_iter)<=mstop))
                    m[, ix] <- ret
                    m}
                pr <- lapply(unique(xselect), pfun)
                for (i in 1:length(pr)) pr[[i]] <- as(pr[[i]], "matrix")
                pr[[length(pr)+1]] <- cbind(matrix(offsets[[comp]],nrow=nrow(newdata),ncol=1), matrix(0,nrow=nrow(newdata), ncol=(mstop-1)))
                preds[[comp]] <- Reduce("+",pr)
            }
            names(preds) <- names(mod$formula)
            return(preds)
 } # if agg
 } # END boost_predict

# coeficients
boost_coef <- function(mod, mstop=NULL){
      if(is.null(mstop)) mstop <- mod$mstop
      coefs <- vector(mode="list",length=length(mod$formula))
      names(coefs) = names(mod$formula)
      for(comp in 1:length(mod$formula)){
          bl <- mod$bl[[comp]]; xselect <- mod$xselect[[comp]][1:mstop]; ens <- mod$ens[[comp]][1:mstop]; nu <- mod$nu[[comp]]; n_iter <- length(xselect)
          # only calculate for bols, bsplines,
          any.bl_lin <- unlist(lapply(bl,function(x) any(class(x) == "bl_lin")))
          if(any(any.bl_lin)){
              for(b_ in 1:length(bl)){ # loop through baselearners
                  if(b_ %in% xselect & any.bl_lin[b_]){ # check if a) selected, b) is a linear baselearner
                      coefs[[comp]] <- c(coefs[[comp]], list(nu*Reduce("+",lapply(ens[which((xselect ==b_) & (1:n_iter <=mstop))],function(x) x$model)))) # append coefficients
                  } else if(any.bl_lin[b_]){
                      fakedat = get("weights",envir=environment(bl[[1]]$fit))*0
                      fakefit = bl[[b_]]$fit(fakedat)
                      coefs[[comp]] <- c(coefs[[comp]], list(fakefit$model)) # append 0 coefficents
                  } else{
                      coefs[[comp]] <- c(coefs[[comp]], list(NA)) # append NA for btrees
                  }
              }
          }
          names(coefs[[comp]]) <- mod$learnernames[[comp]]
      }
      names(coefs) <- names(mod$formula)      
      return(coefs)
}

# weighted versions
cjs.loglike.include2='
double loglike_CH_i(arma::ivec yrow, arma::vec dSurv, arma::vec dP){
   int ifirst = (int) arma::min(arma::find(yrow==1));
   int ilastsight = arma::max(arma::find(yrow==1));
   int T = yrow.size();
   double loglike;
   if(ilastsight!=(T-1)){
        arma::vec chi_p_neverseen(T-ilastsight);
        chi_p_neverseen[chi_p_neverseen.size()-1] = 1;
        for(int t_ =(chi_p_neverseen.size()-2), c_ = T-1; t_>-1;t_--,c_--){
           chi_p_neverseen[t_] = (1-dSurv[c_])+(1-dP[c_])*dSurv[c_]*chi_p_neverseen[t_+1];
        }
        loglike = log(chi_p_neverseen[0]);
   } else {
       loglike = 0; 
   }
   for(int t_=(ifirst+1);t_<(ilastsight+1);t_++){
       loglike += log(dSurv[t_]*pow(dP[t_],(yrow[t_]))*pow((1-dP[t_]),(1-yrow[t_])));
   }
   return loglike;
}'
#  WORKS
cjs.loglike.vec = cxxfunction(sig = c(Y="integer",logitp="numeric",logits="numeric",first="integer",weights="integer"),'
arma::imat yvec = Rcpp::as<arma::ivec>(Y);
arma::vec s = 1/(1+exp(-1*Rcpp::as<arma::vec>(logits)));
arma::vec p = 1/(1+exp(-1*Rcpp::as<arma::vec>(logitp)));
arma::vec w = Rcpp::as<arma::vec>(weights);
int nyx =yvec.size();
arma::uvec switch_ = arma::find((Rcpp::as<arma::vec>(first))==1);
arma::uvec::iterator i = switch_.begin(), j=switch_.begin(), jend = switch_.end();
j++;
double ll=0;
for(;j!=jend;j++){
   ll += loglike_CH_i(yvec.rows(*i,(*j)-1), s.rows(*i,(*j)-1), p.rows(*i,(*j)-1)) * w(*j-1);
   i++;
}
ll += loglike_CH_i(yvec.rows(*i,nyx-1), s.rows(*i,nyx-1), p.rows(*i,nyx-1))* w(nyx-1);
return wrap(ll);',plugin="RcppArmadillo",includes=cjs.loglike.include2)




# This samples from the two-slice marginals. INPUTS:  Y vector of counts; p vector of logit p , s vector of logit s (same length as Y), and first : vector with 1 scoring if it is a new sequence of ch (i.e., the 'first' capture)
#OUTPUTS: three columns representing: z1z2 = {1,1}, {1,0}, {0,0}
twoslicemarginal.include='
arma::mat alphabeta_gamma(arma::ivec yrow, arma::vec dSurv, arma::vec dP, arma::ivec switcher){
   int T = yrow.size();
   arma::mat alph= arma::ones(2,T), beta = arma::ones(2,T), gamma = arma::zeros(3,T);
   alph(1,0)=0;
   arma::vec alph_unnorm(2), gam_unnorm(3);
   arma::vec::iterator pi = dP.begin(), si = dSurv.begin();
   arma::ivec::iterator yi = yrow.begin(), sw = switcher.begin();
   for(int i=1; i<T; i++){
     yi++;
     si++;
     pi++;
     sw++;
     alph_unnorm(0)= (1-*sw)*(alph(0,i-1)*(*si)*pow(*pi,*yi)*pow(1-*pi,1-*yi)) + (*sw);
     alph_unnorm(1)= (1-*sw)*(alph(0,i-1)*(1-*si)*pow(1,1-*yi)*pow(0,*yi) + alph(1,i-1)*pow(1,1-*yi)*pow(0,*yi));
     alph.col(i) = alph_unnorm/(arma::accu(alph_unnorm));
   }
   for(int j=T-2; j>-1;j--){
     beta(0,j) = (1-*sw)*(beta(0,j+1)*(*si)*pow(*pi,*yi)*pow(1-*pi,1-*yi) + beta(1,j+1)*(1-*si)*pow(1,1-*yi)*pow(0,*yi)) + (*sw);
     beta(1,j) = (1-*sw)*(beta(1,j+1)*pow(1,1-*yi)*pow(0,*yi)) + (*sw);
     yi--;
     si--;
     pi--;
     sw--;
   }
   gamma(0,0) = 1;
   for(int i=1; i<T; i++){
     yi++;
     si++;
     pi++;
     gam_unnorm(0) = alph(0,i-1)*pow(*pi,*yi)*pow(1-*pi,(1-*yi))*beta(0,i)*(*si);
     gam_unnorm(1) = alph(0,i-1)*pow(0,*yi)*beta(1,i)*(1-*si);
     gam_unnorm(2) = alph(1,i-1)*pow(0,*yi)*beta(1,i);
     gamma.col(i) = gam_unnorm/arma::accu(gam_unnorm);
}
return gamma;
}'
twoslicemarginal = cxxfunction(sig = c(y="integer",s="numeric",p="numeric",first="integer"),'
arma::ivec yvec = Rcpp::as<arma::ivec>(y);
arma::vec ls = 1/(1+exp(-1*Rcpp::as<arma::vec>(s)));
arma::vec lp = 1/(1+exp(-1*Rcpp::as<arma::vec>(p)));
int nyx =yvec.size();
arma::ivec switch_ = Rcpp::as<arma::ivec>(first);
arma::mat z11_z10_z00 = alphabeta_gamma(yvec,ls,lp,switch_);
return wrap(arma::trans(z11_z10_z00));',plugin="RcppArmadillo",includes=twoslicemarginal.include)

# R version of the above
two.slice.marginal<-function(y,s,p,id){
    # row 1 is alive; row 2 is dead
  pvec <-1/(1+exp(-p)); n_ <- length(pvec)
  svec <- 1/(1+exp(-s))
  bet <- alph <- matrix(1,2,n_)
  switch_ <- 1*(diff(c(0,id))!=0)#switch id's
  alph[,1] <- c(1,0)#pvec[1]*c(1,0) # by definition, the first case = 1 (alive)
  # forward pass
  for(i in 2:n_){
      alph[1,i] <- (1-switch_[i])*(alph[1,i-1]*svec[i]*(pvec[i]^(y[i])*(1-pvec[i])^(1-y[i])) + alph[2,i-1]*0)+pvec[i]*switch_[i]
      alph[2,i] <- (1-switch_[i])*(alph[1,i-1]*(1-svec[i])*1^(1-y[i])*0^(y[i])+ alph[2,i-1]*1^(1-y[i])*0^(y[i]))
      alph[,i] <- alph[,i]/sum(alph[,i])
  }
  bet[,n_]<-1
  for(i in (n_-1):1){
      bet[1,i] <- (1-switch_[i+1])*(bet[1,i+1]*(svec[i+1])*(pvec[i+1]^y[i+1])*(1-pvec[i+1])^(1-y[i+1]) +  bet[2,i+1]*(1-svec[i+1])*1^(1-y[i+1])*0^(y[i+1])) + 1*switch_[i+1]
      bet[2,i] <- (1-switch_[i+1])*(bet[1,i+1]*0 + bet[2,i+1]*1^(1-y[i+1])*0^(y[i+1]))+ 1*switch_[i+1]
  }
   # norma <- alph*bet; norma <- norma/(c(1,1)%x%t(colSums(norma)))
   #alphbeta.vec(y,pvec,svec,id)
  # backward pass
  z3 <- matrix(0,n_,3)# {t-1,t} (1,1) (1,0) (0,0)
  z3[1,]<-c(1,0,0)
  for(i in 2:n_){
      z3[i,1] <- alph[1,i-1]*(pvec[i]^y[i]*(1-pvec[i])^(1-y[i]))*bet[1,i]*svec[i] # 1->1
      z3[i,2] <- alph[1,i-1]*(0^y[i])*bet[2,i]*(1-svec[i]) # 1->0
      z3[i,3] <- alph[2,i-1]*(0^y[i])*bet[2,i]*1 # 1->0
  }
  normaz <- z3/(rowSums(z3)%x%matrix(1,1,3))
  return(normaz)}

##########################################################################
# This samples from the posterior of z: 0 is alive 1 is dead

forwfilter.backsamp.include='
 int sample_category(arma::vec prob){
   GetRNGstate();
   arma::vec randp=arma::randu(1);
   PutRNGstate();
   arma::uvec whichk = arma::find(arma::cumsum(prob) > randp(0));
   return ((int) arma::min(whichk));
 }
 arma::vec concat_normalize(double elem1, double elem2, arma::vec vec12){
    vec12(0) *= elem1;
    vec12(1) *= elem2;
    return (vec12/arma::accu(vec12));
 }
 arma::vec alpha_backsamp(arma::ivec yrow, arma::vec dSurv, arma::vec dP, arma::ivec switcher){
    int T = yrow.size();
    arma::mat alph= arma::ones(2,T);
    alph(1,0)=0;
    arma::vec alph_unnorm(2), z(T);
    arma::vec::iterator pi = dP.begin(), si = dSurv.begin();
    arma::ivec::iterator yi = yrow.begin(), sw = switcher.begin();
    for(int i=1; i<T; i++){
      yi++;
      si++;
      pi++;
      sw++;
      alph_unnorm(0)= (1-*sw)*(alph(0,i-1)*(*si)*pow(*pi,*yi)*pow(1-*pi,1-*yi)) + (*sw);
      alph_unnorm(1)= (1-*sw)*(alph(0,i-1)*(1-*si)*pow(1,1-*yi)*pow(0,*yi) + alph(1,i-1)*pow(1,1-*yi)*pow(0,*yi));
      alph.col(i) = alph_unnorm/(arma::accu(alph_unnorm));
    }
    z(T-1) = sample_category(alph.col(T-1));
    for(int j=T-2; j>-1;j--){
       if(*sw == 1) {
         z(j) = sample_category(alph.col(j));
         sw--;
       } else {
         sw--;
         z(j) = sample_category(concat_normalize( ((1-z(j+1))*(*si)*pow(*pi,(*yi))*pow(1-*pi,1-*yi) + (z(j+1))*(1-*si)),  pow(0,1-z(j+1)), alph.col(j)));
       }
       yi--;
       si--;
       pi--;
    }
    return z;
 }'
forwfilter.backsamp = cxxfunction(sig = c(y="integer",s="numeric",p="numeric",first="integer"),'
 arma::ivec yvec = Rcpp::as<arma::ivec>(y);
 arma::vec ls = 1/(1+exp(-1*Rcpp::as<arma::vec>(s)));
 arma::vec lp = 1/(1+exp(-1*Rcpp::as<arma::vec>(p)));
 int nyx =yvec.size();
 arma::ivec switch_ = Rcpp::as<arma::ivec>(first);
 arma::vec latentseq = alpha_backsamp(yvec,ls,lp,switch_);
 return wrap(latentseq);',plugin="RcppArmadillo",includes=forwfilter.backsamp.include)


#############################################3
# Boosting helper function: function to go from 1 entry per individual to long vectorized input
ind.elem_to_long.vector <- function(x,T, pre.first.mask,first.mask.long, replace.first="next"){
    x.long = rep(x,each=T)[!pre.first.mask]
    if(replace.first == "next"){
       x.long[which(first.mask.long==1)] <- x.long[which(first.mask.long==1)+1]
    } else {
        x.long[which(first.mask.long==1)]<- replace.first
    }
    return(x.long)
}

# BLACK-BOX OPENER: function to tally which predictors were selected by a btree model
gettreeIDs <- function(treem){
    ret <- numeric(0)
    curtrees <- list(treem)
    while(length(curtrees)>0){
        nexttrees <- list()
        for(tt in 1:length(curtrees)){
            curtt_varid <- c(-999999,curtrees[[tt]][[5]][[1]])
            if(any(curtt_varid>0)){ ret <- c(ret,curtt_varid[-1])}
            if(!curtrees[[tt]][[4]]){
                nexttrees[[length(nexttrees)+1]] <- curtrees[[tt]][[8]];
                nexttrees[[length(nexttrees)+1]] <- curtrees[[tt]][[9]]
            }
        } # end greedy trees
        curtrees <- nexttrees
    } # end while
    ret
}

# BLACK-BOX OPENER: function to tally which predictors were selected by a btree model
selected_in_ensemble <- function(mod, mstop= mod$mstop){
    varids_all <- list()
    for(comp in 1:length(mod$formula)){		
        varids <- vector(mstop, mode="list") # container for the variables used
        xselect <- mod$xselect[[comp]]; ens <- mod$ens[[comp]]; baseln <- mod$learnernames[[comp]]
        uxselect <- unique(xselect); uxselect <- uxselect[which(uxselect!=0)]
        for(w in uxselect){
            bm_vars <- all.vars(formula(paste("~",baseln[w])))	# litte macro to grab the variable names in the baselearner formula
            ix <- which(xselect == w & (1:mod$mstop) <= mstop)
            if(length(grep("btree", baseln[w])) > 0){	# test if this wth baselearner is a btree
                varids[ix] <- lapply(ix, function(bm_i) {   bm_vars[gettreeIDs(ens[[bm_i]]$model)] })
            } else {		# simple to just count the variables in static baselearners
                varids[ix] <- list(c(bm_vars))	
            }
        }
        varids_all[[comp]] <- varids
    }
names(varids_all) <- names(mod$formula)
return(varids_all)}


# a baselearner to avoid overfitting: returns coefficients of zero
# bnofit should return an object of class "blg" which has entries "dpp", etc, get_data, all this unnecessary stuff
# blg$dpp() is a function that takes WEIGHTS and returns a list fit, hatvalues, predict, etc., and is a class of "bl_lin","bl"
bnofit <- function(...,by = NULL, index = NULL, intercept = FALSE, df = NULL, lambda = 0){
    ret<-list(
        dpp=function(weights){
             if(is.matrix(weights) | is.data.frame(weights)){ n_ <- nrow(weights) } else { n_ <- length(weights)}
             force(n_)
             dppret <- list(
                 fit=function(y){
                     n_n <- length(y); force(n_n)
                     fitret <- list(model=matrix(0),fitted=function(){rep(0,n_n)})
                     class(fitret) <- c("bm_lin","bm")
                     return(fitret)
                 },
                 hatvalues=0,
                 predict=function(bm,newdata=NULL,aggregate="sum"){
                 if(!is.null(newdata)){ pred <- rep(0,nrow(newdata))
                               } else {
                                   pred<-rep(0,n_)
                               }
                 return(pred)},
                 df=0,
                 Xnames="nofit")
             class(dppret)<-c("bl","bl_lin")
             return(dppret)}
        )
    class(ret)<-"blg"    
    return(ret)}


###############################################
# CJSboost function
cjsboost <-function(formula, # named list of R formula's (response not necessary)
                    ch.data,  # matrix for WIDE format capture-recapture
                    cov.data.wide=NULL,  # data.frame WIDE format for individual level covariates
                    cov.data.array=NULL, # data.frame LONG format for time+individual varying cova
                    mstop = 3000, # stopping criteria, either single integer (for all components) or named list for different criteria per component (named the same as in formula)
                    m_estw=3,# how often to perform E-step? (every m_estw'th of boosting iteratn)
                    nu=lapply(formula,function(x){r<-0.01; r}), # named list of the shrinkage rate, for different shrinkage per component (named the same as in formula)
                    offsets=NA, # named list of start values per component (named the same as in formula
                    weights=1, # weights should be nrow=nrow(ch.data)
                    oobag_risk = FALSE, # estimate risk descent as empirical risk (FALSE) or out-of-bag risk (TRUE); default is FALSE; generally only used in the cvrisk routine
                    id = NULL, # optional vector of IDs to identify rows in data with ch data
                    timefactor.constraint = FALSE, # whether to enforce a p_{T-1} = p_{T} constraint
                    add.intercept=TRUE, # option to automatically add an intercept variable called 'interc'    
                    allblg=NULL, # optional patch-in for prior defined base-learner (e.g.,for CV)
                    oldother=NULL,# optional patch-in for prior defined base-learner(e.g.,for CV)
                    return_blg=TRUE, # option to return baselearners (e.g., for CV and prediction)
                    time.spacing=NULL # not implimented
    ){
    first = apply(ch.data,1,function(y.row) min(which(y.row!=0))) # find first capture
    T=ncol(ch.data); nrowch = nrow(ch.data);
    if(any(first==T)){ stop("error: please remove capture histories that only have data in the final capture period") }
    y.vec = as.vector(t(ch.data)) # vectorise the data
    rm(ch.data);
    if(is.null(id)) {id =1:nrowch};
    id.vec = rep(id,each=T) # identify individuals 
    if(length(weights)==1 & is.numeric(weights)){weights=rep(weights,nrowch)}
    pre.first.mask = as.vector(sapply(first,function(ist,T) {c(1-numeric(ist-1),rep(0,T-ist+1))},T)) # find pre-first capture
    first.mask.long = as.vector(sapply(first,function(ist,T) {c(rep(1,ist),rep(0,T-ist))},T))[!pre.first.mask] # only first+ captures
    inv.first.mask.long = 1-first.mask.long # post-first capture
    y.long = y.vec[!pre.first.mask] # only first+ captures: used internally
    Ny = length(y.long)
    id.long = id.vec[!pre.first.mask] # only first+ id's
    weights.long  = ind.elem_to_long.vector(weights,T,pre.first.mask,first.mask.long,replace.first=0)# zero-weights on the firstcapt
    # for cross-validation: insample and oob weights
    if(oobag_risk & any(weights == 0)){ # if estimating risk on an out-of-bag sample (implied by zero entries in weights argument)
        assess.weights = 1*((1:nrowch) %in% which(weights==0)) # weights on each INDIVIDUAL
        assess.weights.long = ind.elem_to_long.vector(assess.weights,T,pre.first.mask,first.mask.long,replace.first=0)
    } else { # if no oob_risk, then just assess on the insample training weights
        assess.weights.long = weights.long
    }
    # base learner indices
    allvars <- unique(unlist(lapply(formula, function(x) all.vars(formula(x)))))        
    learnernames <- lapply(formula, function(x) attr(terms(x),"term.labels")) # names of base-learners in formula
    lrners.uniq <- unique(unlist(learnernames))  # unique baselearners
    lrners.uniq.formu <- lapply(learnernames,match,lrners.uniq) # map unique baselearners to their posit,ion in the formula
    # process the covariate data
    if(is.null(allblg)){ # notice we can shunt-in allblg so we don't need to re-make the blg (expensive)
        timefactor.long = rep(1:T,nrowch)[!pre.first.mask] #
        if(timefactor.constraint){ timefactor.long[which(timefactor.long == T)] <- T-1 } # enforcing a p(T-1)=p(T) constraint
        timefactor.long <- factor(timefactor.long,levels=1:( (T-1)*timefactor.constraint + T*!timefactor.constraint ))
        dat = data.frame(id=factor(id.long), timefactor=timefactor.long)
        timefactor.levels = levels(dat$timefactor)
        if(add.intercept) dat$interc<-1    
        dat$time <- as.numeric(scale(rep(c(2,2:T),nrowch)[!pre.first.mask],scale=FALSE))
        time.scaling<-unique(data.frame(t= rep(c(2,2:T),nrowch)[!pre.first.mask],covariate= dat$time,stringsAsFactors=FALSE));time.scaling<-time.scaling[order(time.scaling$t),]
        if(!is.null(cov.data.wide)){
            for(nam_ in names(cov.data.wide)[names(cov.data.wide)%in%allvars]){
                dat[[nam_]]<-ind.elem_to_long.vector(cov.data.wide[,nam_],T,pre.first.mask,first.mask.long,replace.first="next")
            }}
        if(!is.null(cov.data.array)){
            for(nam_ in dimnames(cov.data.array)[[3]][dimnames(cov.data.array)[[3]]%in%allvars]){
                dat[[nam_]]<-as.vector(t(cov.data.array[,,nam_]))[!pre.first.mask]
            }}
        if(!all(allvars%in%names(dat))) stop(paste("missing data called:",paste(allvars[which(allvars%in%names(dat)==FALSE)],collapse=",")))
         # INITIALIZE BASE-LEARNERS:
        allblg <- lapply(lrners.uniq, function(bltext) {ret <- eval(parse(text = paste("with(data=dat,expr=",bltext,")"))); return(ret)})
        attr(allblg,"learnernames") <- learnernames
        attr(allblg,"lrners_uniq") <- lrners.uniq
        attr(allblg,"lrners_uniq_formu") <- lrners.uniq.formu
        rm(cov.data.wide,cov.data.array)
    } else { # done if(is.null(allblg))
        time.scaling = NULL
        timefactor.levels=NULL
    }
    allbl <- lapply(allblg, function(curblg,w.inbag) {curblg$dpp(w.inbag)},w.inbag=weights.long)       
    blg <- lapply(lrners.uniq.formu, function(ix) allblg[ix])# reorganize base-learner progenitors (ordered by formula
    bl <- lapply(lrners.uniq.formu, function(ix) allbl[ix])
    names(blg) <- names(bl)  <- names(formula)
    # optimize offsets: maximum likelihood on the intercept model
    if(any(is.na(offsets))){
        opt.offsets = lapply(1:10, function(x) { optim(par=runif(2,-0.5,0.5), fn=function(x,Y,first,w){ -1*cjs.loglike.vec(Y,rep(x[2],length(Y)),rep(x[1],length(Y)),first,w)}, gr=NULL, Y=y.long,first=first.mask.long,w=weights.long)[c("par","value","convergence")]});
        offsets = as.list(opt.offsets[[ which.min(unlist(lapply(opt.offsets, function(x) (x[["convergence"]] ==0)*x[["value"]] + (x[["convergence"]] !=0)*10^12)))]]$par); names(offsets) = names(formula)
    }
    # calcalculate expections of z1->z2: smooth two-slice-marginasl:{1->1}, {1->0},{0->0}
    #w_z <- two.slice.marginal(y.vec,rep(offsets$s,length(y.vec)),rep(offsets$p,length(y.vec)),id.vec) # smooth two-slice marginals
    #w_z.long1 <- w_z[!pre.first.mask,]
    w_z.long <- twoslicemarginal(y.vec[!pre.first.mask],rep(offsets$s,sum(!pre.first.mask)),rep(offsets$p,sum(!pre.first.mask)),first.mask.long) # smooth two-slice marginals
    fit = lapply(offsets, function(o){ rep(o,Ny)}) # starting fit vector
    # gradient functions
    n.gradient = list(s = function(f,p,y,z11,z10,pre.first.mask) (z11-z10*exp(f))/(1+exp(f))*pre.first.mask,
                  p = function(f,s,y,z11,z10,pre.first.mask) ((z11*exp(f)+z11)*y-z11*exp(f))/(exp(f)+1)*pre.first.mask
                  )
    Qfunc = function(y,p,s,z11,z10,z00,weights.long) (-z11*(log(1/(1+exp(-s)))+ y*log(1/(1+exp(-p)))+(1-y)*log(1-1/(1+exp(-p)))) - z10*log(1-1/(1+exp(-s))) - z00)*weights.long
    # initial negative gradients
    u = list(s = n.gradient[["s"]](fit[[1]],fit[[-1]],y=y.long, z11=w_z.long[,1], z10=w_z.long[,2],inv.first.mask.long),
         p = n.gradient[["p"]](fit[[2]],fit[[-2]],y=y.long, z11=w_z.long[,1], z10=w_z.long[,2],inv.first.mask.long))
    # processes for the boosting: ens, xselect, mrisk, ss, ...
    ens <- lapply(formula, function(x) vector(mstop, mode="list")) # container for selected ensemble
    xselect <- lapply(formula, function(x) rep(NA,mstop)) # container for selected variables
    #startrisk <- 
    mrisk <- mqfunc <- numeric(mstop) # track the risk and the qfunction
    ss <- lapply(blg, function(x) vector(length(x), mode="list")) # container for iterative sums of squares per base-learner
    tsums <- lapply(blg, function(x) rep(NA, length(x)))
    # pick baselearners by minimizing the target loglikelkhood
    altselectfunct = list(s = function(tmpf,f,nu,altf,Y,first,weights) -1*cjs.loglike.vec(Y=Y,logitp=altf[[1]],logits=f+nu*tmpf$fitted(), first=first,weights=weights), p = function(tmpf,f,nu,altf,Y,first,weights) -1*cjs.loglike.vec(Y=Y,logitp=f+nu*tmpf$fitted(),logits=altf[[1]], first=first,weights=weights))
    # START GRADIENT DESCENT
    totalcomps <- length(formula)
    for(m in 1:mstop){ # booting iterations
      for(cp in 1:totalcomps){# loop through components
        # recalculate thie gradient
        u[[cp]] <- n.gradient[[cp]](f=fit[[cp]],fit[[-cp]],y=y.long,z11=w_z.long[,1],z10=w_z.long[,2],inv.first.mask.long)
       #tmpfitted <- lapply(bl[[cp]], function(blfit,u) try(blfit$fit(y=u),silent=TRUE),u=u[[cp]]) # fit to neg gradient
        tmpfitted <- lapply(bl[[cp]], function(blfit,u)  blfit$fit(y=u),u=u[[cp]]) # fit to neg gradient
        #xselect[[cp]][m]<-bestbl<-which.min(lapply(tmpfitted, FUN=altselectfunct[[cp]], f=fit[[cp]],nu=nu[[cp]],altf=fit[-cp],Y=y.long,first=first.mask.long,weights=weights.long)) # selection by minimizing likelihood
        xselect[[cp]][m]<-bestbl<-which.min(lapply(tmpfitted, function(fitd,w,u) { sum(w*(fitd$fitted()-u)^2)/sum(w)},w=weights.long,u=u[[cp]])) # selection but best fit to negative gradient
        # update fit
        fit[[cp]]<-fit[[cp]]+nu[[cp]]*tmpfitted[[bestbl]]$fitted()
        # ensure proper class labels for mboost/modeltools
        ens[[cp]][[m]]<-list(model=tmpfitted[[bestbl]]$model); class(ens[[cp]][[m]]) <- class(tmpfitted[[bestbl]])        
      } # cp
      # calculate the Q.func
      mqfunc[m] = sum(Qfunc(y=y.long,p=fit[["p"]],s=fit[["s"]],z11=w_z.long[,1],z10=w_z.long[,2],z00=w_z.long[,3],weights.long))
      mrisk[m] = -1*cjs.loglike.vec(Y=y.long,logitp=fit[["p"]],logits=fit[["s"]], first=first.mask.long,weights=assess.weights.long)
      # redo E step
      if(m %% m_estw ==0){
        w_z.long <- twoslicemarginal(y.long,fit[["s"]],fit[["p"]],first.mask.long) # 
      } # E step
    } # m
    summary_ <- lapply(xselect,table) # lapply(fit,function(x) inv.logit(unique(x)))u
    for(cp in 1:totalcomps){names(summary_[[cp]])<-learnernames[[cp]][as.numeric(names(summary_[[cp]]))]}
    # return the fit vector as a matrix (padded with NA's for pre-first capture)
    fit.matrix <- lapply(fit, function(x, frst,frst.long,nrowch,T){ ret=matrix(NA,nrowch,T)
        wfirst = which(first.mask.long==1);
        for(i in 1:(length(wfirst)-1)){  ret[i,(frst[i]+1):T]<- x[(wfirst[i]+1):(wfirst[i+1]-1)] }
        ret[nrow(ret),(frst[length(frst)]+1):T] <- x[(wfirst[length(wfirst)]+1):length(x)]
        return(ret)
    },frst=first, frst.long=first.mask.long, nrowch=nrowch,T=T)
    res <- list(formula=formula, learnernames=learnernames, nu=nu, mstop=mstop, ens = ens, m_estw=m_estw, xselect=xselect,fit=fit.matrix, qfunc=mqfunc, risk=mrisk, offsets=offsets, bestm = which.min(mrisk),id=id,summary=summary_, add.intercept=add.intercept, allblg=NULL, bl=bl, time.scaling=time.scaling, timefactor.levels=timefactor.levels)
    if(return_blg){ # for CV; option to return prior defined baselearners
        res$allblg=allblg
    }
    return(res)
} # end HMMboost


#################################################
# cross-validation / subsampling / bootstrap function
subsampF<-function(labels, # vector of grouping variable
    ntimes=15, # number of times to repeat
    oobfraction=1/3, # oobag fraction (doesn't apply for bootstrap because it always returns about ~ 63% of the data (maybe max ~77)
    method="none",# none, stratified, bootstrap (NOTE: conditional on labels...0
    preW=rep(1,length(labels))){ # preweights (e.g., if there is a greater bootstrap going on
  N <- length(labels); Ncl <- length(unique(labels)) 
  if(method=="none"){ # just random subsampling (by Weights)
    ret <- lapply(1:ntimes, function(run,N,preW,oobf){
      sel<-sample(1:N, round(sum(preW>0)*(1-oobf)), replace=FALSE, prob=1*(preW>0))
      retw <- ((1:N %in% sel)*1)*preW; return(list(inbag=retw))}, N=N,preW=preW,oobf=oobfraction) # end default
  } else if(method=="bootstrap"){ # bootstrapping (by Weights)
    ret<-lapply(1:ntimes, function(run,N,preW){
          sel <- sample(1:N,sum(preW),replace=TRUE,prob=preW/sum(preW)) # sample based on preweights
          retw <- as.numeric(table(sel)[as.character(1:N)])
          retw[is.na(retw)]<-0; return(list(inbag=retw))},N=N, preW=preW) # end bootstrap
  } else if(method=="stratify"){ # subsample, but ensure each class in labels is in inbag and oobag
    ret <- lapply(1:ntimes, function(run,labels,N,preW,oobf) {
      classes <- unique(labels[preW>0]); ncuts <- round(1/oobf)
      sing.perm <- lapply(classes, function(cl) {
        tmpcl <- which(labels==cl & preW>0); lencl <- length(tmpcl); rep(tmpcl[sample(x=1:lencl, size=lencl)], length.out=max(ncuts,lencl)) })
      permut <- unlist(sing.perm) 
      indices <- vector(ncuts, mode="list")
      for (i in 1:length(permut)) {
        k = ((i-1) %% ncuts) +1
        indices[[k]] <- c(indices[[k]], permut[i])} # permutation is done by reordering and stepping through by 1/oobfaction
      retw <- ((1:N %in% unique(unlist(indices[1:(ncuts-1)])))*1)*preW
      return(list(inbag=retw))},labels=labels,N=N,preW=preW,oobf=oobfraction) # end stratified
  } else if(method=="strataboot"){ # bootstrap by entire groups, multiplying weights
    ret <- lapply(1:ntimes, function(run,labels,N,preW){
      classes <- unique(labels[which(preW>0)]); Ncl<-length(classes)
      selcl<-sample(classes,Ncl,replace=TRUE)
      retw <- as.numeric(table(selcl)[as.character(labels)])
      retw[is.na(retw)] <- 0
      retw <- retw * preW; return(list(inbag=retw)) },labels=labels,N=N,preW=preW)
  } else if(method=="stratifyboot"){ # bootstrap within groups...
    ret <- lapply(1:ntimes, function(run,labels,N,preW){
      classes <- unique(labels[which(preW>0)]); Ncl<-length(classes)
      selind <- lapply(classes,function(cl){
        ixcl <- which(labels==cl)
        if(length(ixcl)>1){ sample(x=ixcl, size=sum(preW[ixcl]), replace=TRUE,prob=preW[ixcl]/sum(preW[ixcl]))
                          } else { return(ixcl)}
      })
      retw <-as.numeric(table(unlist(selind))[as.character(1:N)])
      retw[is.na(retw)]<-0; return(list(inbag=retw))},labels=labels, N=N, preW=preW) # end bootstrap
  } else if(method=="constrainedboot"){ # bootstrap ensure all (sub) members are represented: the probabilities are not exact, but for just a few constrains, they are approximately correct. CONSTRAINTS: for each column of matrix labels, ensure each level is represented at least once.
      if(!(class(labels)=="matrix" | class(labels)=="data.frame")){
          labels=matrix(labels)
      }
      selN<-nrow(labels)
      ret <- lapply(1:ntimes, function(run,labels){
          selind<-numeric(0)
          for(col_ in 1:ncol(labels)){
              if(!all(levels(labels[,col_])%in%unique(labels[selind,col_]))){
                  for(cl in levels(labels[,col_])[which(levels(labels[,col_])%in%unique(labels[selind,col_])==FALSE)]){
                      selind[length(selind)+1]<-sample(which(labels[,col_]==cl),1)
                  }
              }
          }
          if(is.na(oobfraction)){
              retw=tabulate(c(selind,sample(1:selN,size=selN-length(selind),replace=TRUE,prob=(rep(2,selN)*(!((1:selN)%in%selind))+rep(1,selN)*((1:selN)%in%selind))/(2*(selN-length(selind))+length(selind)))),nbins=selN)
          } else {
              oob_remove=sample((1:selN)[which((1:selN)%in%selind==FALSE)],size=round(selN*oobfraction),replace=FALSE) # discard
              retw=tabulate(c((1:selN)[which((1:selN)%in%oob_remove==FALSE)],sample((1:selN)[which((1:selN)%in%oob_remove==FALSE)],size=round(oobfraction*selN),replace=TRUE)),nbins=selN)
          }
      return(list(inbag=retw))},labels=labels)
  } else { stop(paste("you specified an unknown method",method,"in function subsampF"))}
  return(ret)} # end subsampF


###################################
# cross-validation sub-function: function to parallelize; used internally in cjsboost_cvrisk
cjsboost_cv_fit <- function(X,formula,ch.data,mstop,m_estw,nu,offsets,add.intercept,id,allblg,oldother){
    cvmod=NULL
    try(cvmod<-cjsboost(formula = formula,ch.data=ch.data,cov.data.wide =NULL,cov.data.array=NULL,mstop=mstop,m_estw = m_estw,nu = nu,offsets=offsets,weights=X,oobag_risk=TRUE,add.intercept=add.intercept,id=id,allblg=allblg,return_blg=FALSE,oldother=oldother),silent=FALSE) # run the Cross-validation model
    # weights=X; ch.data=NULL;cov.data.wide =NULL;cov.data.array=NULL;oobag_risk=TRUE
    ret<-list()
    if(is.null(cvmod)){ ret$success=0
                        ret$message=geterrmessage()
                        ret$risk=Inf
                        ret$newoffset=lapply(offsets,function(x) return(NA))
                        ret$bestm=0
                        ret$stabilsel=lapply(formula,function(x) return(NA))
    } else { # cvmod success
        allvars=lapply(formula,all.vars) # all variable names
        sel_ens= selected_in_ensemble(cvmod,mstop=mstop)
        ret$success=1
        ret$risk=cvmod$risk
        ret$bestm=cvmod$bestm
        ret$message=""
        ret$stabilsel=mapply(iform=formula,selv=sel_ens,FUN=function(iform,selv){ tvarz= all.vars(iform); tvarz=tvarz[which(tvarz!=paste(iform)[2])]; res=matrix(Inf,nrow=length(tvarz),ncol=1,dimnames=list(tvarz,"first"));
            for(tt in tvarz){ res[tt,"first"]<-min(c(Inf,which(unlist(lapply(selv, function(selvi,ttt){ttt%in%selvi},ttt=tt)))))}
            return(res)},SIMPLIFY=FALSE)
    }
    return(ret)}
# done CV-sub function


cjsboost_cvrisk <- function(            #
    formula, # named list of R formula's (response not necessary)
    ch.data,  # matrix for WIDE format capture-recapture
    cov.data.wide=NULL,  # data.frame WIDE format for individual level covariates
    cov.data.array=NULL, # data.frame LONG format for time+individual varying cova
    mstop = 3000, # stopping criteria, either single integer (for all components) or named list for different criteria per component (named the same as in formula)
    m_estw=2,# how often to perform E-step? (every m_estw'th of boosting iteratn)
    nu=lapply(formula,function(x){return(0.1)}), # named list of the shrinkage rate, for different shrinkage per component (named the same as in formula)
    offsets=NA, # named list of start values per component (named the same as in formula
    add.intercept=TRUE, # option to automatically add an intercept variable called 'interc'
    id = NULL, # optional vector of IDs to identify rows in data with ch data 
    allblg=NULL, # optional patch-in for prior defined base-learner
    time.spacing=NULL, # not implimented
    N_bootstrap=30,# number of bootstrap iteration
    mc.cores=1, # package(parallel): parallelize the bootstrap runs, number of cores 
    rerun_failures=TRUE,
    bootstrap_weights=NULL, # option to patch in your own weights,
    bootstrap_method="constrainedboot", # if no bootstrap weights supplied, then this is passed to subsampF to generate weights
    oobfraction=0.36, # if no bootstrap weights supplied, then this is passed to subsampF to generate weights
    best_mstop_mean_trim = 0.1, #  when estimating the optimal mstop, trim the mean by this fraction (on both ends)
    plot=FALSE # plot the gradient descent
    ) {
    first_ = first(ch.data); T=ncol(ch.data);
    # check first-captures != T
    if(any(first_==T)){ stop("please remove capture histories whose first observation is the final capture period T") }
    if(!all(names(formula)==c("s","p"))){ stop("please supply argument 'formula' as a named list, with the names 's' for survival, and 'p' for capture probability, in that order") } # check named list
    # base-model, just to set up the base-learners; check for errors
    basemod=cjsboost(formula = formula,ch.data = ch.data,cov.data.wide = cov.data.wide,cov.data.array = cov.data.array,mstop = 2,m_estw = 1,nu = nu,offsets = offsets,weights = rep(1,nrow(ch.data)),oobag_risk = FALSE,add.intercept = add.intercept,id=id,allblg=allblg,return_blg =TRUE)
    allvars=unique(unlist(lapply(formula,all.vars)))
    findfactors <- which(unlist(lapply(cov.data.wide[,which(names(cov.data.wide)%in%allvars),drop=FALSE],class))%in%c("factor","character"))
    bootstrap_labels= apply(cbind(data.frame(first=first(ch.data)),cov.data.wide[,findfactors,drop=FALSE]),1,function(x_) paste0(x_,collapse=""))
    # need to bootstarp WITHIN factors to ensure numerical stability
    if(is.null(bootstrap_weights)){ # option to enter your own
        bootstrap_weights <- lapply(subsampF(bootstrap_labels, ntimes=N_bootstrap, method=bootstrap_method,oobfraction=oobfraction), function(x) x$inbag)
    }
    # RUN THE MODELS in parallel
    allmods=parallel::mclapply(X=bootstrap_weights, FUN=cjsboost_cv_fit, formula,ch.data,mstop,m_estw,nu,offsets=NA,add.intercept,id,allblg=basemod$allblg, mc.cores=mc.cores)
    cv_successes=unlist(lapply(allmods,function(x) x$success))
    messages=''
    if(any(!cv_successes)){
        print(paste("the following", sum(!cv_successes)," errors occurred:"))
        messages=lapply(allmods[which(!cv_successes)],function(err_) err_$message)
        for(err_ in which(!cv_successes)){ # print collected errors
            print(allmods[[err_]]$message)
        } 
        if(rerun_failures){ # try to rerun with new weights
           bootstrap_weights[which(!cv_successes)]<-lapply(subsampF(bootstrap_labels, ntimes=sum(!cv_successes), method=bootstrap_method,oobfraction=oobfraction), function(x) x$inbag)
       }
    }
    cv_successes=unlist(lapply(allmods,function(x) x$success))
    which_suc <-which(cv_successes==1)    
    if(any(cv_successes==1)){
        print(paste0(length(which_suc)," runs succeeded"))
    } else {
        stop("every run failed. Perhaps you have too many categorical variables?")
    }
    #N_bootstrap <- sum(cv_successes==1)
    cvrisk = matrix(NA,mstop,N_bootstrap)
    for(j in which_suc){ # gather all CV Risk descent profiles
            cvrisk[,j]=allmods[[j]]$risk
    } # done gather cv risk profiles
    # rescale the cvrisk: Divide by Nj, the number of post-first-capture observations
    Nj.holdout <- unlist(lapply(bootstrap_weights,function(w_,frst,T){ sum((T-frst)*(w_==0))},frst=first(ch.data),T=ncol(ch.data))) # number of post-capture observations
    cvrisk.rescale <- cvrisk * 1/(t(Nj.holdout) %x% (rep(1,nrow(cvrisk)))) # divid neg.loglike by number of observations
    #cvrisk.rescale <- cvrisk -  t(cvrisk[1,]%x%t(rep(1,nrow(cvrisk))))
    # stability selection
    stabilsel=vector(mode="list",length=length(formula)); names(stabilsel) <- names(formula)# container for stability selection
    for(cp_ in 1:length(formula)){ # gather results of stability selection
        cv_first_sel <- do.call("cbind", lapply(allmods[which_suc], function(rs,cp_) {rs[["stabilsel"]][[cp_]]},cp_=cp_))
        stabilsel[[cp_]]<-apply(cv_first_sel,1, function(firstsel_){
            sapply(1:mstop, function(thres,regv2) {
                sum(regv2<=thres)/length(regv2) }, regv2=firstsel_)
        }) # apply
    } # end cp_ through stability selection
    # find best m
    bestms=apply(cvrisk,2,function(x) which.min(x)) # best m
    if(plot==TRUE){
        plot(c(1,mstop),range(cvrisk[,which_suc]),typ="n",ylab="holdout risk",xlab="boosting iteration");
        for(j in which_suc){lines(cvrisk[,j],col=rpois(1,10))}
    }
    return(list(cvrisk=cvrisk, cvrisk.rescale=cvrisk.rescale, stabilsel=stabilsel,bestm=list(mean=mean(bestms,trim=best_mstop_mean_trim),median=median(bestms)), messages=messages,bootstrap_weights=bootstrap_weights,which_succeeded=which_suc))}
# DONE CV stability selection
# cvrisk: is the matrix (mstop x N_boot) with the holdout (raw) risk estimated per m, per subsample/cv run
# cvrisk.rescale: same as cvrisk, but the risk is divided by Nj, the number of post-firstcapture observations, a truer defn of risk
# stabilsel: named list, with an entry for each parameters; each entry is a matrix showing each covariates stability select path over mstop iterations. The mean of the stability selection probabilities, for each covariate, would be an approximate posterio inclusion probability
#

cjsboost_hyperparam <- function(            #]
    formula, # named list of R formula's (response not necessary)
    ch.data,  # matrix for WIDE format capture-recapture
    cov.data.wide=NULL,  # data.frame WIDE format for individual level covariates
    cov.data.array=NULL, # data.frame LONG format for time+individual varying cova
    mstop = 3000, # stopping criteria, either single integer (for all components) or named list for different criteria per component (named the same as in formula)
    m_estw=2,# how often to perform E-step? (every m_estw'th of boosting iteratn)
    nu.start=lapply(formula,function(x){return(0.1)}), # named list of the shrinkage rate, for different shrinkage per component (named the same as in formula)
    nu_search_steps = 7, # how many steps to take to find an optimal nu
    offsets=NA, # named list of start values per component (named the same as in formula
    add.intercept=TRUE, # option to automatically add an intercept variable called 'interc'
    id = NULL, # optional vector of IDs to identify rows in data with ch data 
    allblg=NULL, # optional patch-in for prior defined base-learner
    time.spacing=NULL, # not implimented
    N_bootstrap=30,# number of bootstrap iteration
    mc.cores=1, # package(parallel): parallelize the bootstrap runs, number of cores 
    rerun_failures=TRUE, 
    bootstrap_weights=NULL, # option to patch in your own weights,
    bootstrap_method="bootstrap", # if no bootstrap weights supplied, then this is passed to subsampF to generate weights
    oobfraction=0.36, # if no bootstrap weights supplied, then this is passed to subsampF to geinerate weights
    best_mstop_mean_trim = 0.1, #  when estimating the optimal mstop, trim the mean by this fraction (on both ends)
    return_all_cv_models=FALSE #  return ALL the cvmods (not a good idea); otherwise, just return the CV runs for the optimal nu (default)
    ){
    #     formula, # named list of R formula's (response not necessary)
#    formula=formu.PLS;ch.data=y; cov.data.wide=cov.dat.indiv; cov.data.array=cov.dat.array;mstop=3000; m_estw=2; nu_search_steps = 7; offsets=NA; add.intercept=TRUE; id = NULL; allblg=NULL; time.spacing=NULL; bootstrap_method = "bootstrap",rerun_failures=FALSE; best_mstop_mean_trim = 0.1; plot=FALSE
    # find all the variables in formula    
    allvars=unique(unlist(lapply(formula,all.vars)))
    # check bootstrap weights
    if(is.null(bootstrap_weights)){
        # check which method of subsampling, and whether it requires labels
        if((bootstrap_method%in%c("stratify","strataboot","stratifyboot","constrainedboot")) & !is.null(cov.data.wide)){ #
            findfactors <- names(which(unlist(lapply(cov.data.wide[,which(names(cov.data.wide)%in%allvars),drop=FALSE],class))%in%c("factor","character"))) # find the factors                 
            bootstrap_labels= apply(cbind(data.frame(first=first(ch.data)),cov.data.wide[,findfactors,drop=FALSE]),1,function(x_) paste0(x_,collapse=""))
            bootstrap_weights <- lapply(subsampF(bootstrap_labels, ntimes=N_bootstrap, method=bootstrap_method,oobfraction=oobfraction), function(x) x$inbag) # subsampling with some sort of grouping variable(s)
        } else {
            bootstrap_weights <- lapply(subsampF(bootstrap_labels=rep(1,nrow(ch.data)), ntimes=N_bootstrap, method=bootstrap_method,oobfraction=oobfraction), function(x) x$inbag) # basic subsampling
        } # type of bootstrap
    } # if no user-supplied bootstrap
    nu.mean <- mean(unlist(nu.start))
    cur.nu.ratio <- nu.ratio <- nu.start$p/nu.start$s # ratio of nu(p)/nu(s)
    cur.p = uniroot(function(x,nu.mu,nu.ratio){ (mean(c(x,x/nu.ratio))-nu.mu)},nu.mu=nu.mean,nu.ratio=cur.nu.ratio,interval=c(0.00000001,1))$root    
    nu.cur <- list(s=cur.p/cur.nu.ratio, p = cur.p)
    # run the base-model, just to generate the requisite base-learners 
    basemod=cjsboost(formula = formula,ch.data = ch.data,cov.data.wide = cov.data.wide,cov.data.array = cov.data.array,mstop = 2,m_estw = 1,nu = nu.cur,offsets = offsets,weights = rep(1,nrow(ch.data)),oobag_risk = FALSE,add.intercept = add.intercept,id=id,allblg =allblg,return_blg =TRUE)
    # run first BOOTSTRAP-VALIDATION 
    cvmod <- cjsboost_cvrisk(formula=formula, ch.data=ch.data, cov.data.wide=cov.data.wide, cov.data.array=cov.data.array, mstop = mstop, m_estw=m_estw, nu=nu.cur, offsets=offsets, add.intercept=add.intercept, id = id, allblg=basemod$allblg, time.spacing=NULL,N_bootstrap=length(bootstrap_weights), mc.cores=mc.cores, rerun_failures=rerun_failures, bootstrap_weights=bootstrap_weights, bootstrap_method=bootstrap_method, oobfraction=oobfraction, best_mstop_mean_trim = best_mstop_mean_trim, plot=FALSE)
    bestm = round(cvmod$bestm$median) # best m (mean is likely unreliable)    
    cvmods <- cvres <-list() # storage container for cv-risk results summary and models
    cvres[[1]]<-c(nu.s = nu.cur$s, nu.p = nu.cur$p, nu.ratio = cur.nu.ratio, bestm = bestm, cvrisk = mean(cvmod$cvrisk[round(cvmod$bestm$median),]))
    cvmods[[length(cvmods)+1]] <- cvmod # store the cross-validation model for later
    # next nu to search
    next.nu.ratio <-nu.ratio*c(0.5,2)
    # loop through the search proceedure
    while(nu_search_steps >0){ 
        for(subi in 1:length(next.nu.ratio)){
            cur.nu.ratio <- next.nu.ratio[subi] # update the nu-ratio
            cur.p = uniroot(function(x,nu.mu,nu.ratio){ (mean(c(x,x/nu.ratio))-nu.mu)},nu.mu=nu.mean,nu.ratio=cur.nu.ratio,interval=c(0.00000001,1))$root
            nu.cur <- list(s=cur.p/cur.nu.ratio, p = cur.p)            
            # run next cv
            cvmod <- cjsboost_cvrisk(formula=formula, ch.data=ch.data, cov.data.wide=cov.data.wide, cov.data.array=cov.data.array, mstop = mstop, m_estw=m_estw, nu=nu.cur, offsets=offsets, add.intercept=add.intercept, id = id, allblg=basemod$allblg, time.spacing=NULL,N_bootstrap=length(bootstrap_weights), mc.cores=mc.cores, rerun_failures=rerun_failures, bootstrap_weights=bootstrap_weights, bootstrap_method=bootstrap_method, oobfraction=oobfraction, best_mstop_mean_trim = best_mstop_mean_trim, plot=FALSE)             
            bestm = round(cvmod$bestm$median) # best m (mean is likely unreliable)
            if(bestm >= (mstop-1)){ print(paste0("warning: best-m at ",bestm," has hit the boundary 'mstop'=",mstop,". Results may be unreliable. Please restart and supply a larger mstop, or larger initial values for nu")) }
            cvres[[length(cvres)+1]] = c(nu.s = nu.cur$s, nu.p = nu.cur$p, nu.ratio = cur.nu.ratio, bestm = bestm, cvrisk = mean(cvmod$cvrisk[bestm,]))
            cvmods[[length(cvmods)+1]] <- cvmod
        }
        cvres = cvres[order(unlist(lapply(cvres,function(x) x[["nu.ratio"]])))]
        cvmods <- cvmods[order(unlist(lapply(cvres,function(x) x[["nu.ratio"]])))]
                                        # find next best nu
        bestcv.cur = min(do.call("rbind",cvres)[,"cvrisk"]); bestcv.cur.ix = which.min(do.call("rbind",cvres)[,"cvrisk"])
        range.ratio.cur = range(do.call("rbind",cvres)[,"nu.ratio"])
        if(any(cvres[[bestcv.cur.ix]][["nu.ratio"]] == range.ratio.cur)){ # test if the ratio is on a boundary
            next.nu.ratio = cvres[[bestcv.cur.ix]][["nu.ratio"]]*( 2^(cvres[[bestcv.cur.ix]][["nu.ratio"]] == max(range.ratio.cur))*0.5^(cvres[[bestcv.cur.ix]][["nu.ratio"]] == min(range.ratio.cur)) ) # double or half nu.p, if it is on a boundary
        } else { # else, best nu is caught inbetween two run: take average
            next.nu.ratio = 0.5* cvres[[bestcv.cur.ix]][["nu.ratio"]] + 0.5*(do.call("rbind",cvres)[bestcv.cur.ix+c(-1,1),"nu.ratio"][which.min(do.call("rbind",cvres)[bestcv.cur.ix+c(-1,1),"cvrisk"])]) # weight inbetween neighbours with lowest risk
        }
        nu_search_steps = nu_search_steps-1
    }
    bestcv.mod = which.min(do.call("rbind",cvres)[,"cvrisk"])
    cvmod = cvmods[[bestcv.mod]]
    bestm = round(cvmod$bestm$median) # the mean is too unreliable here
    nu = list(s = cvres[[bestcv.mod]][["nu.s"]],p = cvres[[bestcv.mod]][["nu.p"]])
    print(paste0("Best cvrisk=",round(cvres[[bestcv.mod]][["cvrisk"]],5)," at nu.s = ",round(nu$s,5),", and nu.p=",round(nu$p,5),", with a nu.ratio of ",round(cvres[[bestcv.mod]][["nu.ratio"]],5)))
    if(return_all_cv_models){
        return(list(cvmods = cvmods, bestcv.mod = bestcv.mod, bestm=bestm, nu.optimization.summary = do.call("rbind",cvres),nu.start=nu.start,bootstrap_weights=bootstrap_weights))
    } else {
        cvmod <- c(cvmod, list(optimal.nu = nu, nu.optimization.summary=do.call("rbind",cvres)))
        return(cvmod)
    }
} # DONE cjsboost_hyperpar function

##############################################################################
# stochastic gradient descent:
cjsboost.stochastic <-function(formula, # named list of R formula's (response not necessary)
                    ch.data,  # matrix for WIDE format capture-recapture
                    cov.data.wide=NULL,  # data.frame WIDE format for individual level covariates
                    cov.data.array=NULL, # data.frame LONG format for time+individual varying cova
                    mstop = 3000, # stopping criteria, either single integer (for all components) or named list for different criteria per component (named the same as in formula)
                    m_estw=20,# how many draws of the posterior to do before updating FIT and recalculating
                    nu=lapply(formula,function(x){r<-0.001; r}), # named list of the shrinkage rate, for different shrinkage per component (named the same as in formula)
                    offsets=NULL, # named list of start values per component (named the same as in formula
                    weights=1, # weights should be nrow=nrow(ch.data)
                    oobag_risk = FALSE, # estimate risk descent as empirical risk (FALSE) or out-of-bag risk (TRUE); default is FALSE; generally only used in the cvrisk routine
                    id = NULL, # optional vector of IDs to identify rows in data with ch data
                    timefactor.constraint = FALSE, # whether to enforce a p_{T-1} = p_{T} constraint
                    add.intercept=TRUE, # option to automatically add an intercept variable called 'interc'    
                    allblg=NULL, # optional patch-in for prior defined base-learner (e.g.,for CV)
                    oldother=NULL,# optional patch-in for prior defined base-learner(e.g.,for CV)
                    return_blg=TRUE, # option to return baselearners (e.g., for CV and prediction)
                    time.spacing=NULL # not implimented
    ){
    if((mstop %% m_estw != 0)){ print("mstop should be a multiple of m_estw")}
    first = apply(ch.data,1,function(y.row) min(which(y.row!=0))) # find first capture
    T=ncol(ch.data); nrowch = nrow(ch.data);
    if(any(first==T)){ stop("error: please remove capture histories that only have data in the final capture period") }
    y.vec = as.vector(t(ch.data)) # vectorise the data
    rm(ch.data);
    if(is.null(id)) {id =1:nrowch};
    id.vec = rep(id,each=T) # identify individuals 
    if(length(weights)==1 & is.numeric(weights)){weights=rep(weights,nrowch)}
    pre.first.mask = as.vector(sapply(first,function(ist,T) {c(1-numeric(ist-1),rep(0,T-ist+1))},T)) # find pre-first capture
    first.mask.long = as.vector(sapply(first,function(ist,T) {c(rep(1,ist),rep(0,T-ist))},T))[!pre.first.mask] # only first+ captures
    inv.first.mask.long = 1-first.mask.long # post-first capture
    y.long = y.vec[!pre.first.mask] # only first+ captures: used internally
    Ny = length(y.long)
    id.long = id.vec[!pre.first.mask] # only first+ id's
    weights.long  = ind.elem_to_long.vector(weights,T,pre.first.mask,first.mask.long,replace.first=0)# zero-weights on the firstcapt
    # for cross-validation: insample and oob weights
    if(oobag_risk & any(weights == 0)){ # if estimating risk on an out-of-bag sample (implied by zero entries in weights argument)
        assess.weights = 1*((1:nrowch) %in% which(weights==0)) # weights on each INDIVIDUAL
        assess.weights.long = ind.elem_to_long.vector(assess.weights,T,pre.first.mask,first.mask.long,replace.first=0)
    } else { # if no oob_risk, then just assess on the insample training weights
        assess.weights.long = weights.long
    }
    # base learner indices
    allvars <- unique(unlist(lapply(formula, function(x) all.vars(formula(x)))))        
    learnernames <- lapply(formula, function(x) attr(terms(x),"term.labels")) # names of base-learners in formula
    lrners.uniq <- unique(unlist(learnernames))  # unique baselearners
    lrners.uniq.formu <- lapply(learnernames,match,lrners.uniq) # map unique baselearners to their posit,ion in the formula
    # process the covariate data
    if(is.null(allblg)){ # notice we can shunt-in allblg so we don't need to re-make the blg (expensive)
        timefactor.long = rep(1:T,nrowch)[!pre.first.mask] #
        if(timefactor.constraint){ timefactor.long[which(timefactor.long == T)] <- T-1 } # enforcing a p(T-1)=p(T) constraint
        timefactor.long <- factor(timefactor.long,levels=1:( (T-1)*timefactor.constraint + T*!timefactor.constraint ))
        dat = data.frame(id=id.long, timefactor=timefactor.long)
        timefactor.levels = levels(dat$timefactor)        
        if(add.intercept) dat$interc<-1    
        dat$time <- as.numeric(scale(rep(c(2,2:T),nrowch)[!pre.first.mask],scale=FALSE))
        time.scaling<-unique(data.frame(t= rep(c(2,2:T),nrowch)[!pre.first.mask],covariate= dat$time,stringsAsFactors=FALSE));time.scaling<-time.scaling[order(time.scaling$t),]        
        if(!is.null(cov.data.wide)){
            for(nam_ in names(cov.data.wide)[names(cov.data.wide)%in%allvars]){
                dat[[nam_]]<-ind.elem_to_long.vector(cov.data.wide[,nam_],T,pre.first.mask,first.mask.long,replace.first="next")
            }}
        if(!is.null(cov.data.array)){
            for(nam_ in dimnames(cov.data.array)[[3]][dimnames(cov.data.array)[[3]]%in%allvars]){            
                dat[[nam_]]<-as.vector(t(cov.data.array[,,nam_]))[!pre.first.mask]
            }}
        if(!all(allvars%in%names(dat))) stop(paste("missing data called:",paste(allvars[which(allvars%in%names(dat)==FALSE)],collapse=",")))
         # INITIALIZE BASE-LEARNERS:
        allblg <- lapply(lrners.uniq, function(bltext) {ret <- eval(parse(text = paste("with(data=dat,expr=",bltext,")"))); return(ret)})
        attr(allblg,"learnernames") <- learnernames
        attr(allblg,"lrners_uniq") <- lrners.uniq
        attr(allblg,"lrners_uniq_formu") <- lrners.uniq.formu
        rm(cov.data.wide,cov.data.array)
    } else { # done if(is.null(allg))
        time.scaling = NULL
        timefactor.levels=NULL
    } 
    allbl <- lapply(allblg, function(curblg,w.inbag) {curblg$dpp(w.inbag)},w.inbag=weights.long)       
    blg <- lapply(lrners.uniq.formu, function(ix) allblg[ix])# reorganize base-learner progenitors (ordered by formula
    bl <- lapply(lrners.uniq.formu, function(ix) allbl[ix])
    names(blg) <- names(bl)  <- names(formula)
    # optimize offsets: maximum likelihood on the intercept model
    if(any(is.na(offsets))){
        opt.offsets = lapply(1:10, function(x) { optim(par=runif(2,-0.5,0.5), fn=function(x,Y,first,w){ -1*cjs.loglike.vec(Y,rep(x[2],length(Y)),rep(x[1],length(Y)),first,w)}, gr=NULL, Y=y.long,first=first.mask.long,w=weights.long)[c("par","value","convergence")]});
        offsets = as.list(opt.offsets[[ which.min(unlist(lapply(opt.offsets, function(x) (x[["convergence"]] ==0)*x[["value"]] + (x[["convergence"]] !=0)*10^12)))]]$par); names(offsets) = names(formula)
    }
    # make an initial posterior draw from 
    w_z.long <- twoslicemarginal(y.vec[!pre.first.mask],rep(offsets$s,sum(!pre.first.mask)),rep(offsets$p,sum(!pre.first.mask)),first.mask.long) # smooth two-slice marginals
    fit = lapply(offsets, function(o){ rep(o,Ny)}) # starting fit vector
    # gradient functions: z0 is one time step back; z1 is current
    n.gradient = list(
        s = function(f,p,y,z0,z1,inv.first){ inv.first*(z0*z1/(1+exp(f)) - z0*(1-z1)*exp(f)/(exp(f)+1))},
        p = function(f,s,y,z0,z1,inv.first){inv.first*(z0*z1*((1+exp(f))*y-exp(f))/(exp(f)+1))})    
    # complete data likelihood: z0 is one time step back; z1 is current
    Qfunc <- function(y,s,p,z0,z1,weights.long){ (-z0*z1*( y*log(1/(1+exp(-p))) + (1-y)*log(1/(1+exp(p)))+log(1/(1+exp(-s)))) - (1-z1)*z0*(log(1/(1+exp(s)))))*weights.long}
    # initial negative gradients
    u = vector(mode="list",length=length(formula));names(u)=names(formula)
    # processes for the boosting: ens, xselect, mrisk, ss, ...
    ens <- lapply(formula, function(x) vector(mstop, mode="list")) # container for selected ensemble
    xselect <- lapply(formula, function(x) rep(NA,mstop)) # container for selected variables
    #startrisk <- 
    mrisk <- mqfunc <- numeric(mstop) # track the risk and the qfunction
    ss <- lapply(blg, function(x) vector(length(x), mode="list")) # container for iterative sums of squares per base-learner
    tsums <- lapply(blg, function(x) rep(NA, length(x)))
    # pick baselearners by minimizing the target loglikelkhood
    altselectfunct = list(s = function(tmpf,f,nu,altf,Y,first,weights) -1*cjs.loglike.vec(Y=Y,logitp=altf[[1]],logits=f+nu*tmpf$fitted(), first=first,weights=weights), p = function(tmpf,f,nu,altf,Y,first,weights) -1*cjs.loglike.vec(Y=Y,logitp=f+nu*tmpf$fitted(),logits=altf[[1]], first=first,weights=weights))
    # START GRADIENT DESCENT
    totalcomps <- length(formula)
    for(bigm in 1:(mstop/m_estw)){ # booting
        parallel.fit<-fit # temporary storage for each (parallel) boosting run
        for(subm in 1:m_estw){ # subloop
            # draw from the posterior of z: 0 = alive, 1=dead (but I want to switch it, so 1: alive)
            m <- subm + m_estw*(bigm-1)
            post.z = 1-forwfilter.backsamp(y=y.long,s=fit[["s"]],p=fit[["p"]],first=first.mask.long)
            for(cp in 1:totalcomps){# loop through components
                # recalculate thie gradient
                u[[cp]] = n.gradient[[cp]](fit[[cp]],fit[[-cp]],y=y.long, z0=c(0,post.z[-Ny]),z1=post.z,inv.first=inv.first.mask.long)
#               tmpfitted <- lapply(bl[[cp]], function(blfit,u) try(blfit$fit(y=u),silent=TRUE),u=u[[cp]]) # fit to neg gradient
                tmpfitted <- lapply(bl[[cp]], function(blfit,u)  blfit$fit(y=u),u=u[[cp]]) # fit to neg gradient
                xselect[[cp]][m]<-bestbl<-which.min(lapply(tmpfitted, FUN=altselectfunct[[cp]], f=fit[[cp]],nu=nu[[cp]],altf=fit[-cp],Y=y.long,first=first.mask.long,weights=weights.long)) # selection by minimizing likelihood
                #xselect[[cp]][m]<-bestbl<-which.min(lapply(tmpfitted, function(fitd,w,u) { sum(w*(fitd$fitted()-u)^2)/sum(w)},w=weights.long,u=u[[cp]])) # selection but best fit to negative gradient
                parallel.fit[[cp]] = parallel.fit[[cp]]+nu[[cp]]*tmpfitted[[bestbl]]$fitted()
                # ensure proper class labels for mboost/modeltools
                ens[[cp]][[m]]<-list(model=tmpfitted[[bestbl]]$model); class(ens[[cp]][[m]]) <- class(tmpfitted[[bestbl]])        
            } # different components
        } # done parallel runs
        # update fit
        fit = parallel.fit
        # calculate the Q.func
        mqfunc[(m-m_estw+1):m] = sum(Qfunc(y=y.long, s=fit[["s"]], p=fit[["p"]], z0=c(0,post.z[-Ny]), z1=post.z,weights.long))
        mrisk[(m-m_estw+1):m] = -1*cjs.loglike.vec(Y=y.long,logitp=fit[["p"]],logits=fit[["s"]], first=first.mask.long,weights=assess.weights.long)
      # redo E step
    } # boosting m
    summary_ <- lapply(xselect,table) # lapply(fit,function(x) inv.logit(unique(x)))u
    for(cp in 1:totalcomps){names(summary_[[cp]])<-learnernames[[cp]][as.numeric(names(summary_[[cp]]))]}
    # make the fit vector into a matrix (for easy comparisons, later)
        fit.matrix <- lapply(fit, function(x, frst,frst.long,nrowch,T){ ret=matrix(NA,nrowch,T)
        wfirst = which(first.mask.long==1);
        for(i in 1:(length(wfirst)-1)){  ret[i,(frst[i]+1):T]<- x[(wfirst[i]+1):(wfirst[i+1]-1)] }
        ret[nrow(ret),(frst[length(frst)]+1):T] <- x[(wfirst[length(wfirst)]+1):length(x)]
        return(ret)
    },frst=first, frst.long=first.mask.long, nrowch=nrowch,T=T)        
    res <- list(formula=formula, learnernames=learnernames, nu=nu, mstop=mstop, ens = ens, m_estw=m_estw, xselect=xselect, fit=fit.matrix, qfunc=mqfunc, risk=mrisk, offsets=offsets, bestm = which.min(mrisk),id=id,summary=summary_, add.intercept=add.intercept, allblg=NULL, bl=bl,time.scaling=time.scaling, timefactor.levels=timefactor.levels)
    if(return_blg){ # for CV; option to return prior defined baselearners
        res$allblg=allblg
    }
    return(res)
} # end stochastic cjsboost

cjsboost_cv_fit.stochastic <- function(X,formula,ch.data,mstop,m_estw,nu,offsets,add.intercept,id,allblg,oldother){
    cvmod=NULL
    try(cvmod<-cjsboost.stochastic(formula = formula,ch.data=ch.data,cov.data.wide =NULL,cov.data.array=NULL,mstop=mstop,m_estw = m_estw,nu = nu,offsets=offsets,weights=X,oobag_risk=TRUE,add.intercept=add.intercept,id=id,allblg=allblg,return_blg=FALSE,oldother=oldother),silent=FALSE) # run the Cross-validation model
    # weights=X; ch.data=NULL;cov.data.wide =NULL;cov.data.array=NULL;oobag_risk=TRUE
    ret<-list()
    if(is.null(cvmod)){ ret$success=0
                        ret$message=geterrmessage()
                        ret$risk=Inf
                        ret$newoffset=lapply(offsets,function(x) return(NA))
                        ret$bestm=0
                        ret$stabilsel=lapply(formula,function(x) return(NA))
    } else { # cvmod success
        allvars=lapply(formula,all.vars) # all variable names
        sel_ens= selected_in_ensemble(cvmod,mstop=mstop)
        ret$success=1
        ret$risk=cvmod$risk
        ret$bestm=cvmod$bestm
        ret$message=""
        ret$stabilsel=mapply(iform=formula,selv=sel_ens,FUN=function(iform,selv){ tvarz= all.vars(iform); tvarz=tvarz[which(tvarz!=paste(iform)[2])]; res=matrix(Inf,nrow=length(tvarz),ncol=1,dimnames=list(tvarz,"first"));
            for(tt in tvarz){ res[tt,"first"]<-min(c(Inf,which(unlist(lapply(selv, function(selvi,ttt){ttt%in%selvi},ttt=tt)))))}
            return(res)},SIMPLIFY=FALSE)
    }
    return(ret)}
# done CV-sub function


cjsboost_cvrisk.stochastic <- function(            #
    formula, # named list of R formula's (response not necessary)
    ch.data,  # matrix for WIDE format capture-recapture
    cov.data.wide=NULL,  # data.frame WIDE format for individual level covariates
    cov.data.array=NULL, # data.frame LONG format for time+individual varying cova
    mstop = 3000, # stopping criteria, either single integer (for all components) or named list for different criteria per component (named the same as in formula)
    m_estw=3,# how often to perform E-step? (every m_estw'th of boosting iteratn)
    nu=lapply(formula,function(x){r=0.01; r}), # named list of the shrinkage rate, for different shrinkage per component (named the same as in formula)
    offsets=NA, # named list of start values per component (named the same as in formula
    add.intercept=TRUE, # option to automatically add an intercept variable called 'interc'
    id = NULL, # optional vector of IDs to identify rows in data with ch data 
    allblg=NULL, # optional patch-in for prior defined base-learner
    time.spacing=NULL, # not implimented
    N_bootstrap=30,# number of bootstrap iteration
    mc.cores=1, # package(parallel): parallelize the bootstrap runs, number of cores 
    rerun_failures=TRUE,
    bootstrap_weights=NULL, # option to patch in your own weights,
    bootstrap_method="constrainedboot", # if no bootstrap weights supplied, then this is passed to subsampF to generate weights
    oobfraction=0.25, # if no bootstrap weights supplied, then this is passed to subsampF to generate weights
    best_mstop_mean_trim = 0.1, #  when estimating the optimal mstop, trim the mean by this fraction (on both ends)
    plot=TRUE # plot the gradient descent
    ) {
    first_ = first(ch.data); T=ncol(ch.data);
    # check first-captures != T
    if(any(first_==T)){ stop("please remove capture histories whose first observation is the final capture period T") }
    if(!all(names(formula)==c("s","p"))){ stop("please supply argument 'formula' as a named list, with the names 's' for survival, and 'p' for capture probability, in that order") } # check named list
    # basic model
    basemod=cjsboost.stochastic(formula = formula,ch.data = ch.data,cov.data.wide = cov.data.wide,cov.data.array = cov.data.array,mstop = 2,m_estw = 1,nu = nu,offsets = offsets,weights = rep(1,nrow(ch.data)),oobag_risk = FALSE,add.intercept = add.intercept,id=id,allblg =NULL,return_blg =TRUE)
    # if no boo
    allvars=unique(unlist(lapply(formula,all.vars)))
    findfactors <- which(unlist(lapply(cov.data.wide[,which(names(cov.data.wide)%in%allvars),drop=FALSE],class))%in%c("factor","character"))
    bootstrap_labels= apply(cbind(data.frame(first=first(ch.data)),cov.data.wide[,findfactors,drop=FALSE]),1,function(x_) paste0(x_,collapse=""))
    # need to bootstarp WITHIN factors to ensure numerical stability
    if(is.null(bootstrap_weights)){ # option to enter your own
        bootstrap_weights <- lapply(subsampF(bootstrap_labels, ntimes=N_bootstrap, method=bootstrap_method,oobfraction=oobfraction), function(x) x$inbag)
    }
    # RUN THE MODELS in parallel
    allmods=parallel::mclapply(X=bootstrap_weights, FUN=cjsboost_cv_fit.stochastic, formula,ch.data,mstop,m_estw,nu,offsets=NA,add.intercept,id,allblg=basemod$allblg, mc.cores=mc.cores)
    cv_successes=unlist(lapply(allmods,function(x) x$success))
    messages=''
    if(any(!cv_successes)){
        print(paste("the following", sum(!cv_successes)," errors occurred:"))
        messages=lapply(allmods[which(!cv_successes)],function(err_) err_$message)
        for(err_ in which(!cv_successes)){ # print collected errors
            print(allmods[[err_]]$message)
        } 
        if(rerun_failures){ # try to rerun with new weights
           bootstrap_weights[which(!cv_successes)]<-lapply(subsampF(bootstrap_labels, ntimes=sum(!cv_successes), method=bootstrap_method,oobfraction=oobfraction), function(x) x$inbag)
       }
    }
    cv_successes=unlist(lapply(allmods,function(x) x$success))
    which_suc <-which(cv_successes==1)    
    if(any(cv_successes==1)){
        print(paste0(length(which_suc)," runs succeeded"))
    } else {
        stop("every run failed. Perhaps you have too many categorical variables?")
    }
    #N_bootstrap <- sum(cv_successes==1)
    cvrisk = matrix(NA,mstop,N_bootstrap)
    for(j in which_suc){ # gather all CV Risk descent profiles
            cvrisk[,j]=allmods[[j]]$risk
    } # done gather cv risk profiles
    # rescale the cvrisk estimates, for a better visualization
    cvrisk.rescale = cvrisk -  t(cvrisk[1,]%x%t(rep(1,nrow(cvrisk))))
    stabilsel=vector(mode="list",length=length(formula)); names(stabilsel) <- names(formula)# container for stability selection
    for(cp_ in 1:length(formula)){ # gather results of stability selection
        cv_first_sel <- do.call("cbind", lapply(allmods[which_suc], function(rs,cp_) {rs[["stabilsel"]][[cp_]]},cp_=cp_))
        stabilsel[[cp_]]<-apply(cv_first_sel,1, function(firstsel_){
            sapply(1:mstop, function(thres,regv2) {
                sum(regv2<=thres)/length(regv2) }, regv2=firstsel_)
        }) # apply
    } # end cp_ through stability selection
    bestms=apply(cvrisk,2,function(x) which.min(x)) # best m
    if(plot==TRUE){
        plot(c(1,mstop),range(cvrisk[,which_suc]),typ="n",ylab="CV empirical risk",xlab="boosting iteration");
        for(j in which_suc){lines(cvrisk[,j],col=rpois(1,10))}
    }
    return(list(cvrisk=cvrisk, cvrisk.rescale=cvrisk.rescale, stabilsel=stabilsel,bestm=list(mean=mean(bestms,trim=best_mstop_mean_trim),median=median(bestms)), messages=messages,bootstrap_weights=bootstrap_weights,which_succeeded=which_suc))}
