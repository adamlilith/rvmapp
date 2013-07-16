vmapp <-
function(d,
    pred,
    x_vars=NULL,
    sim_n=1000,
    predict_delta=TRUE,
    give_p=TRUE,
    lifn_F1=li_F1,
    lifn_F2=li_F2,
    pars_f1=c(0.1,1,0.2437,0.1,1,0.2437),
    pars_f2=c(1,1,1,1,0),
    F1fn=F1,
    F2fn=F2)
{
    if(is.null(dim(pred))) ##MLE, so expand into a big matrix
    {
        pred <- array(rep(pred,sim_n),dim=c(length(pred),sim_n))
        pred <- t(pred)
    }else{
        ## If a bootstrapped/bayesian prediction matrix
        ## Sample with replacement sim_n times.
        reps_index <- sample(1:dim(pred)[1],sim_n,replace=TRUE)
        pred <- pred[reps_index,]
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
    return_val$p_val_overall <- NULL
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
        return_val$p_val_overall=p_intercept
    }

    ## Delta prediction ##
    return_val$delta <- NULL
    return_val$f1_pars <- NULL
    return_val$f2_pars <- NULL

    if(predict_delta)
    {
        ## FIT F2 ##
        m2 <- cbind(x_,discr_abs) ## Paste prediction and discrepencies together for fast apply of the fitting.
        n_par_f2 <- length(pars_f2)
        f2_pars<-t(apply(m2,1,function(x){
            fit_f2 <- optim(par=pars_f2, 
                fn=lifn_F2,x=x[1:n],
                y=x[(n+1):(2*n)], 
                control = list(maxit = 500,reltol=1e-6))
            return(fit_f2$par)
        }))

        ## Get F2 preds for use in restricting F1 ##
        mp <- cbind(pred,f2_pars)
        f2_preds<-t(apply(mp,1,function(x){
            F2fn(x[1:n],x[(n+1):n_par_f2])
        }))

        ## FIT F1 ##
        m1 <- cbind(x_,discr,f2_preds,pred)        
        n_par_f1 <- length(pars_f1)    
        f1_pars <- t(apply(m1,1,function(x){ 
            fit_f1<- optim(par=pars_f1, 
                fn=lifn_F1,x=x[1:n],
                y=x[(n+1):(2*n)],
                f2_preds=x[(2*n+1):(3*n)],
                pred=x[(3*n+1):(4*n)],
                control = list(maxit = 500,reltol=1e-3))
            return(fit_f1$par)
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
        return_val$f2_preds=f2_preds
    }

    class(return_val) <- "vmapp"
    return(return_val)
}


.deltafn <- function(xx,f1,f2,f1pars,f2pars) 
{
    2 * (f1(xx,f1pars)-0.5) * f2(xx,f2pars)
}

plot.vmapp<-function(x,...)
{
    pars <- x$f1_pars
    n_levels <- 30
    range <- range(x$pred)
    xx <- seq(range[1],range[2],length.out=n_levels)
    mu_y <- numeric(n_levels)
    upper_y <- numeric(n_levels)
    lower_y <- numeric(n_levels)
    y_dist <- array(dim=c(nrow(pars),n_levels))

    for(i in 1:nrow(pars))
        y_dist[i,] <- .deltafn(xx,f1=x$F1fn,
                f2=x$F2fn,
                f1pars=x$f1_pars[i,],
                f2pars=x$f2_pars[i,])

    mu_y <- apply(y_dist,2,mean)
    upper_y <- apply(y_dist,2,quantile,probs=0.975) 
    lower_y <- apply(y_dist,2,quantile,probs=0.025)

    plot(xx,mu_y,
        xlab=expression(hat(p)),
        ylab=expression(delta),
        type='l',
        lwd=2,
        col='red',
        ylim=c(-0.5,0.5),
        ...)
    lines(xx,upper_y,lty=2,lwd=2)
    lines(xx,lower_y,lty=2,lwd=2)
    abline(h=0,lty=3,col='grey')

    legend('topleft',
        legend=c('Mean delta','95% CI'),
        lty=1:2,
        lwd=2,
        col=c('red','black'))
}

print.vmapp <- function(x,...)
{
    cat('###########################################################\n')
    cat('Validation Metric Applied to Probabilistic Predictions\n\n')
    cat(paste(length(x$d),'validation data points used.\n'))
    cat('  Test for overall bias:\n')
    cat('    Average difference between predicted\n')
    cat('    and actual probabilities: ')
    cat(signif(mean(x$pred)-mean(x$d),3) ) 
    cat('.\n')
    cat('    (p < 0.05: Under-estimation,\n')
    cat('     p > 0.95: Over-estimation)\n')   
    cat('    P-value: ')
    cat(signif(x$p_val_overall,3))
    cat('\n')
    cat('  Test for direction change in bias:\n')
    cat('    (p < 0.05: Over/Under,\n')
    cat('     p > 0.95: Under/Over)\n')
    cat('    P-value: ')
    cat(signif(x$p_val_slope,3))
    cat('\n')
    cat('\n###########################################################\n')
}

predict.vmapp <- function(object,x,CI=0.95,rawdist=FALSE,rawprob=FALSE,...)
{
    ## If rawprob == FALSE, x will be interpreted as  ##
    ## the index of a case in the validation set for  ##
    ## which to predict delta.                        ##
    pars <- object$f1_pars
    y_dist <- numeric(nrow(pars))
    if(!rawprob)
    {
        for(i in 1:nrow(pars))
            y_dist[i] <- .deltafn(object$pred[i,x],f1=object$F1fn,
                    f2=object$F2fn,
                    f1pars=object$f1_pars[i,],
                    f2pars=object$f2_pars[i,])
    } else {
    ## If rawprob == TRUE, x will be interpreted as  ##
    ## a predicted probability for which to predict delta##
    ## WARNING: Should only be used for x in range p_hat ##    
        for(i in 1:nrow(pars))
            y_dist[i] <- .deltafn(x,f1=object$F1fn,
                    f2=object$F2fn,
                    f1pars=object$f1_pars[i,],
                    f2pars=object$f2_pars[i,])
    }
    ## If rawdist == TRUE, the samples from the predicted delta
    ## will be returned instead of mean and CIs.
    if(!rawdist) {
        lower_prob <- (1 - CI) / 2
        upper_prob <- 1 - lower_prob
        cat(paste('Using a',CI,'% CI\n'))
        lowerCI <- quantile(y_dist,probs=lower_prob)
        upperCI <- quantile(y_dist,probs=upper_prob)
        result <- data.frame(mean=mean(y_dist),lowerCI=lowerCI,upperCI=upperCI)
        rownames(result) <- NULL
        return(result)
    } else {
        return(y_dist)
    }
}

