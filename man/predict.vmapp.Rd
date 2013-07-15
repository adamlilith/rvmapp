\name{predict.vmapp}
\alias{predict.vmapp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.vmapp(object, x, CI = 0.95, rawdist = FALSE, rawprob = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{CI}{
%%     ~~Describe \code{CI} here~~
}
  \item{rawdist}{
%%     ~~Describe \code{rawdist} here~~
}
  \item{rawprob}{
%%     ~~Describe \code{rawprob} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
