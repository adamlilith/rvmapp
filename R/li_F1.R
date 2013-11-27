li_F1 <-
function(pars,x,y)
{    
    n<-length(x)
    if( sum(pars>0) !=5 ) ## All pars must be positive
        return(1000000)

    ll_pos <- log(F1(x=x[y==1],pars=pars))
    ll_neg <- log(1-F1(x=x[y!=1],pars=pars))

    NLL <- -sum(c(ll_neg,ll_pos))
    return(NLL)
}
