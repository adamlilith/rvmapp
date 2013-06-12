vmapp <-
function(d,pred,x_vars=NULL,sim_n=1000,predict_delta=TRUE,give_p=TRUE,lifn_F1=li_F1,lifn_F2=li_F2,pars_f1=c(0.1,1,0.2,0.1,1,0.2),pars_f2=c(1,1,1,1,0),F1fn=F1,F2fn=F2)
{
    if(is.null(dim(pred))) ##MLE, so expand into a big matrix
    {
      pred <- array(rep(pred,sim_n),dim=c(length(pred),sim_n))
      pred <- t(pred)
    }
    return_val=list()
    use_x <- !is.null(x_vars)

    n <- dim(pred)[2]
    sim_n <- dim(pred)[1]

    x_<-pred
    if(use_x)
        x_<-array(rep(x_vars,sim_n),dim=c(length(x_vars),sim_n))


    total_inv_d<-sum(d)
    if(total_inv_d<=1) ##If only zero or ones in d
    {
        print('Validation data must be heterogeneous (ie not all 0s or all 1s)')
         return(NA)
    }

    ## Simulate outcomes from the prediction vectors
    S<-apply(pred,2,bs)

    discr <- t(apply(S,1,function(x){x-d}))
    discr_abs <- abs(discr)
    
    discr <- (discr + 1)/2 # code S>d as 1, S<d as 0
    discr[discr==0.5] <- NA
    
    ## Keep input parameters in the returned object (for plotting, prediction, etc.)
    return_val$d <- d
    return_val$pred <- pred
    return_val$x_vars <- x_vars
    return_val$F1fn <- F1fn
    return_val$F2fn <- F2fn
    return_val$S <- S
    return_val$discrepencies <- discr_abs
    return_val$direction_discrepencies <- discr

    ## P-values ##
    return_val$p_val_slope <- NULL
    return_val$p_val_intercept <- NULL
    if(give_p)
    {
        ## Slope of F1
        m1 <- cbind(x_,discr)
        f1_coeffs<-apply(m1,1,function(x){
            yy <- x[(n+1):(2*n)]
            xx <- x[1:n]
            fit <- glm(yy ~ xx,family=binomial(logit))
            return(fit$coefficients)
        })
        p_slope<-mean(f1_coeffs[2,]>=0)
        mean_discr<-apply(discr,1,mean,na.rm=TRUE)
        p_intercept<-mean(mean_discr>=0.5)        

        return_val$p_val_slope=p_slope
        return_val$p_val_intercept=p_intercept
    }

    ## Delta prediction ##
    return_val$delta <- NULL
    return_val$f1_pars <- NULL
    return_val$f2_pars <- NULL

    if(predict_delta)
    {
        ## FIT F1 ##
        m1 <- cbind(x_,discr)        
        n_par_f1 <- length(pars_f1)    
        f1_pars <- t(apply(m1,1,function(x){ 
            fit_f1<- optim(par=pars_f1, fn=lifn_F1,x=x[1:n],y=x[(n+1):(2*n)],control = list(maxit = 500,reltol=1e-3))
            return(fit_f1$par)
        }))

        ## FIT F2 ##
        m2 <- cbind(x_,discr_abs) ## Paste prediction and discrepencies together for fast apply of the fitting.
        n_par_f2 <- length(pars_f2)
        f2_pars<-t(apply(m2,1,function(x){
            fit_f2 <- optim(par=pars_f2, fn=lifn_F2,x=x[1:n],y=x[(n+1):(2*n)], control = list(maxit = 500,reltol=1e-3))
            return(fit_f2$par)
        }))

        
        ## Delta ##
        f1f2_pars<-cbind(x_,f1_pars,f2_pars)
        delta<-t(apply(f1f2_pars,1,function(x){
            f1 <- F1fn(x=x[1:n],pars=x[(n+1):(n+n_par_f1)])
            f2 <- F2fn(x=x[1:n],pars=x[(n+1+n_par_f1):(n+n_par_f1+n_par_f2)])
            return(2 * (f1 - 0.5) * f2)
        }))

        return_val$delta=delta
        return_val$f1_pars=f1_pars
        return_val$f2_pars=f2_pars
    }

    class(return_val) <- "vmapp"
    return(return_val)
}


plot.vmapp<-function(x,...)
{
    plot(x$delta[1,])
    for(i in 2:nrow(x$delta))
        points(x$delta[i,])
}

