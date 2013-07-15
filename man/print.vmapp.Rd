\name{print.vmapp}
\alias{print.vmapp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.vmapp(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, ...) 
{
    cat("###########################################################\n")
    cat("Validation Metric Applied to Probabilistic Predictions\n\n")
    cat(paste(length(x$d), "validation data points used.\n"))
    cat("  Test for overall bias:\n")
    cat("    Average difference between predicted\n")
    cat("    and actual probabilities: ")
    cat(signif(mean(x$pred) - mean(x$d), 3))
    cat(".\n")
    cat("    (p < 0.05: Under-estimation,\n")
    cat("     p > 0.95: Over-estimation)\n")
    cat("    P-value: ")
    cat(signif(x$p_val_intercept, 3))
    cat("\n")
    cat("  Test for direction change in bias:\n")
    cat("    (p < 0.05: Over/Under,\n")
    cat("     p > 0.95: Under/Over)\n")
    cat("    P-value: ")
    cat(signif(x$p_val_slope, 3))
    cat("\n")
    cat("\n###########################################################\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
