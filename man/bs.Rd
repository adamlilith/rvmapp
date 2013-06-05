\name{bs}
\alias{bs}

\title{
Fast Bernoulli Sampler
}

\description{
Draw Bernoulli \eqn{n} outcomes \eqn{ k \in \{0,1\} } for an \eqn{n} length vector of probability values \eqn{ 0 \leq p \leq 1 }.
}

\usage{
bs(p)
}

\arguments{
  \item{p}{
a vector containing probability values from which to generate Bernoulli outcomes. All elements must be \eqn{ 0 \leq p \leq 1 }.
}
}

\value{
A logical vector of Bernoulli outcomes of the same length as \code{p}.
}

\references{
See \href{http://bayesianbiologist.com/2012/11/01/variable-probability-bernoulli-outcomes-fast-and-slow/}{bayesianbiologist} for more details.
}

\author{
Corey Chivers <corey.chivers@mail.mcgill.ca>
}

\examples{
    ## Generate a vector of probabilities ##
    probs <- runif(100)
    ## Generate outcomes ##
    outcomes <- bs(probs)
}
