F1 <-
function(x,pars)
{
    a1=pars[1];b1=pars[2];a2=pars[3];b2=pars[4];gamma=pars[5]
    y= (1-exp(-a1*(x)^b1)) * (1-exp(-a2*(1-x)^b2)) + gamma
    return(y)
}
