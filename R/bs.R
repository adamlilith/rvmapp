bs <-
function(p)
{
   U<-runif(length(p),0,1)
   outcomes<-U<p
   return(outcomes)
}
