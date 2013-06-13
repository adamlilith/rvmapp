li_F1 <-
function(pars,x,y,pred,f2_preds)
{    
    x<-x[!is.na(y)]
    y<-y[!is.na(y)]

    n<-length(x)

    ## Check if any would-be predicted delta's are impossible
    f1_preds <- F1(x=pred,pars=pars)
    delta <- 2 * (f1_preds  - 0.5) * f2_preds
    if( any( (delta < -(1-pred)) | (delta > pred) | (f1_preds < 0) | (f1_preds > 1) ) | any(F1(x=0:1,pars=pars) > 1 | F1(x=0:1,pars=pars) < 0 ) ) 
        return(1000000)

    ll_pos <- log(F1(x=x[y==1],pars=pars))
    ll_neg <- log(1-F1(x=x[y!=1],pars=pars))

    NLL <- -sum(c(ll_neg,ll_pos))
    return(NLL)
}
