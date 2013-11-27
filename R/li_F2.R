li_F2 <-
function(pars,x,y,pred,f1_preds)
{    
    x<-x[!is.na(y)]
    y<-y[!is.na(y)]

    n<-length(x)

    ## Check if any would-be predicted delta's are impossible
    f2_preds <- F2(x=pred,pars=pars)
    delta <- 2 * (f2_preds  - 0.5) * f1_preds
    if( any( (delta < -(1-pred)) | (delta > pred) | (f2_preds < 0) | (f2_preds > 1) ) | any(F2(x=0:1,pars=pars) > 1 | F2(x=0:1,pars=pars) < 0 ) ) 
        return(1000000)

    ll_pos <- log(F2(x=x[y==1],pars=pars))
    ll_neg <- log(1-F2(x=x[y!=1],pars=pars))

    NLL <- -sum(c(ll_neg,ll_pos))
    return(NLL)
}
