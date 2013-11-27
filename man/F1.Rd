\name{F1}
\alias{F1}

\title{
Functional form of the probability that a predicted outcome is not equal to the observed data.
}

\description{
This is the default functional form used to fit the probability of deviation between predicted outcomes and observed data. It models \eqn{ P(S \ne D ) }, where \eqn{S} is the simulated outcome from the model predicted probabilities \eqn{\hat{p}} and \eqn{D} is the observed outcome (data).
}

\usage{
F1(x, pars)
}

\arguments{
  \item{x}{
    a vector of locations at which to evaluate \code{F1}.
}
  \item{pars}{
    a vector of parameter values of length 5.
}
}
\details{
This functional form is the default passed to \code{vmapp}. It is defined as:
%% 
    \eqn{ F_1 = ( 1 - e^{-a_1 (x)^{b_1}}) (1 - e^{-a_2 (1 - x)^{b_2}}) + 
        \gamma }
}
\value{
    A vector of the same length as \code{x}.
}

\author{
    Corey Chivers <corey.chivers@mail.mcgill.ca>
}

\examples{
    ## Default starting parameters for \code{vmapp} fitting.
    pars <- c(1,1,1,1,0)
    
    ## Plot F1 over the range 0-1
    curve(F1(x,pars=pars))
}
