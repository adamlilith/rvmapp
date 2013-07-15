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
predict.vmapp(object, x, CI = 0.95, rawdist = FALSE, rawprob = FALSE, ...)
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, x, CI = 0.95, rawdist = FALSE, rawprob = FALSE, 
    ...) 
{
    pars <- object$f1_pars
    y_dist <- numeric(nrow(pars))
    if (!rawprob) {
        for (i in 1:nrow(pars)) y_dist[i] <- .deltafn(object$pred[i, 
            x], f1 = object$F1fn, f2 = object$F2fn, f1pars = object$f1_pars[i, 
            ], f2pars = object$f2_pars[i, ])
    }
    else {
        for (i in 1:nrow(pars)) y_dist[i] <- .deltafn(x, f1 = object$F1fn, 
            f2 = object$F2fn, f1pars = object$f1_pars[i, ], f2pars = object$f2_pars[i, 
                ])
    }
    if (!rawdist) {
        lower_prob <- (1 - CI)/2
        upper_prob <- 1 - lower_prob
        cat(paste("Using a", CI, "\% CI\n"))
        lowerCI <- quantile(y_dist, probs = lower_prob)
        upperCI <- quantile(y_dist, probs = upper_prob)
        result <- data.frame(mean = mean(y_dist), lowerCI = lowerCI, 
            upperCI = upperCI)
        rownames(result) <- NULL
        return(result)
    }
    else {
        return(y_dist)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
