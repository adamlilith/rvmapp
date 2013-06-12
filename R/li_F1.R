li_F1 <-
function(pars,x,y)
{    
    x<-x[!is.na(y)]
    y<-y[!is.na(y)]

    n<-length(x)
    check <- (F1(x=0:1,pars=pars))
    if(sum( check > 0 ) != 2 | sum( check < 1 ) != 2 | sum(is.finite(check)) != 2 ) ## Constrain 0 < F1 < 1
        return(1000000)

    ## Deviation at p_hat = 0 cannot be < 0, nor > 0 at p_hat=1
    if(F1(x=0,pars=pars) < 0.5 | F1(x=1,pars=pars) > 0.5 )
        return(1000000)

    ll_pos <- log(F1(x=x[y==1],pars=pars))
    ll_neg <- log(1-F1(x=x[y!=1],pars=pars))

    NLL <- -sum(c(ll_neg,ll_pos))
    return(NLL)
}
