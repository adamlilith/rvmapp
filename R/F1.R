F1 <-
function(x,pars)
{
    a1=pars[1];b1=pars[2];c1=pars[3];a2=pars[4];b2=pars[5];c2=pars[6]
    y <- (exp(-a1*x^b1)-c1)*(exp(-a2*(1-x)^b2)-c2)
    return(y)
}
