\name{predict.vmapp}
\alias{predict.vmapp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict Method for Class 'vmapp'
}
\description{
Provides predicted deviations (\eqn{\delta}) between actuall probabilities and 
those predicted by the model. 
}
\usage{
\method{predict}{vmapp}(object, x, CI = 0.95, rawdist = FALSE, rawprob = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of type \code{"vmapp"}.
}
  \item{x}{
If \code{rawprob} is set to \code{FALSE}, then \code{x} will be interpreted as
the index of a case in the validation set for which to predict \eqn{\delta}. 
Alternatively, if \code{rawprob} is set to \code{TRUE}, x will be interpreted 
as a predicted probability for which to predict \eqn{\delta}.
WARNING: The latter should only be used for values of \code{x} in range of
predicted probabilities from the model.
}
  \item{CI}{
Width of the confidence interval to return. Defaults to 0.95, meaning a 95\% 
confidence interval.
}
  \item{rawdist}{
Boolean: If set to \code{FALSE} (default), a data frame with the mean, and lower,
and upper confidence boundaries will be returned.
If set to \code{TRUE}, then instead of a confidence interval, raw data
points representing the distribution of predicted \eqn{\delta} will be returned. 
Summary statistics can then be calculated on this vector.
}
  \item{rawprob}{
Boolean: If set to \code{TRUE}, then predictions will be returned for the data
point at index \code{x}. If set to \code{FALSE}, the predictive distribution of 
\eqn{\delta} will be returned which corresponds to predicted probability \code{x}.
}
  \item{\dots}{
Additional arguments will be ignored. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
As outlined in the argument descriptions, the default behaviour is to return 
the mean, and confidence boundaries of the predicted \eqn{\delta}. If \code{rawdist} 
is set to \code{TRUE}, then the returned value will be a vector of samples from
the predictive distribution of \eqn{\delta}.
}
\references{
    Chivers, C., Leung, B. Yan, N. (2013) Probabilistic predictions in ecology and their validation. In prep.
}
\author{
Corey Chivers <corey.chivers@mail.mcgill.ca>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.vmapp}},\code{\link{vmapp}},\code{\link{print.vmapp}}
}
\examples{
    ## Generate outcomes
    p<-runif(100)
    outcomes<-bs(p)

    ## We have a perfect model 
    p_hat<-p

    ## Run vmapp
    #VMAPP <- vmapp(d=outcomes, pred=p_hat)
    #predict(VMAPP,1)   #Predicted delta for the first datapoint.
    #predict(VMAPP,0.2,rawprob=TRUE) #Predicted delta p_hat=0.2
    #predict(VMAPP,1,rawdist=TRUE)   #Predicted delta distribution. This returns the
                                     #each datapoint in the monte carlo simulation of delta.
    
}

