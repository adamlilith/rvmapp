\name{vmapp}
\alias{vmapp}

\title{
Validation Metric Applied to Probabilistic Predictions (VMAPP).
}
\description{
Provides an inferential metric of goodness of fit of a binary prediction model using independent validation data. Also provides an estimate of the differences between probabilities predicted by a model and the underlying probabilities generating the observations.
}
\usage{
vmapp(d, pred, x_vars = NULL, sim_n = 1000, predict_delta = TRUE, 
    give_p = TRUE, lifn_F1 = li_F1, lifn_F2 = li_F2, 
    pars_f1 = c(0.1, 1, 0.2437, 0.1, 1, 0.2437),
    pars_f2 = c(1, 1, 1, 1, 0), F1fn = F1, F2fn = F2)
}

\arguments{
  \item{d}{
    A vector of binary outcomes. This is validation data NOT used to fit the model.
}
  \item{pred}{
    A vector of predicted probabilities \eqn{\hat{p}}. Alternatively, a \eqn{m} by \eqn{n} matrix of prediction vectors obtained via bootstrapping or Bayesian methods. \eqn{n} is the number of validation data points, \eqn{m} is the number of prediction vectors provided by bootstrapping or Bayesian methods. 
}
  \item{x_vars}{
    A vector of quantities on which to assess deviations. If not given, the default is to use \eqn{\hat{p}}.
}
  \item{sim_n}{
    Integer number of simulated outcomes from the prediction vector to use.
}
  \item{predict_delta}{
    Boolean value indicating whether or not to return estimated deviation between model predictions (\eqn{\hat{p}}) and generating probabilities (\eqn{p}). This process is computationally expensive and should be set to \code{FALSE} if not needed.
}
  \item{give_p}{
    Boolean value indicating whether or not to return p-values for overall and directional biases.
}
  \item{lifn_F1}{
    A function which returns the negative log-likelihood of \eqn{ P(S > D \mid S \ne  D) }, where \eqn{S} is the simulated outcome from \eqn{\hat{p}} and \eqn{D} is the observed outcome (data). See \code{\link{li_F1}}, which is the default.
}
  \item{lifn_F2}{
    A function which returns the negative log-likelihood of \eqn{P(S \ne  D)}. See \code{\link{li_F2}}, which is the default.
}
  \item{pars_f1}{
    A vector of starting parameters for \code{F1} to be passed to \code{optim}. 
}
  \item{pars_f2}{
    A vector of starting parameters for \code{F2} to be passed to \code{optim}. 
}
  \item{F1fn}{
    A function which returns \eqn{ P(S > D \mid S \ne  D) }, where \eqn{S} is the simulated outcome from \eqn{\hat{p}} and \eqn{D} is the observed outcome (data). See \code{\link{F1}}, which is the default.
}
  \item{F2fn}{
    A function which returns \eqn{P(S \ne  D)}. See \code{\link{F2}}, which is the default.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{p_val_slope}{
    The proportion of simulated outcomes for which deviation direction was possitively related to \eqn{\hat p} (or \code{x_vars},if given). This is a two-tailed test, where p-value \eqn{ \le \alpha/2 } indicates a significant negative relationship. Alternatively, p-value \eqn{ \ge 1 - \alpha/2 } indicates a significant positive relationship.
    }
  \item{p_val_intercept }{
    The proportion of simulated outcomes for which the mean deviation direction greater than zero (zeror indicating no overall bias). This is also a two-tailed test, where p-value \eqn{ \le \alpha/2 } indicates a significant negative bias (\eqn{\hat p < p}). Alternatively, p-value \eqn{ \ge 1 - \alpha/2 } indicates a significant positive (\eqn{\hat p > p}).
    }
  \item{delta }{
    A matrix of estimated differences between model predictions and the true underlying probabilities (\eqn{\hat p - p}). This matrix constitutes a distribution of differences at each data point. 
    }
  \item{f1_pars }{
    A matrix of the fitted parameters of \eqn{F_1}. 
    }
  \item{f2_pars }{
    A matrix of the fitted parameters of \eqn{F_2}. 
    }

}
\references{
    Chivers, C., Leung, B. Yan, N. (2013) Probabilistic predictions in ecology and their validation. In prep.
}
\author{
Corey Chivers <corey.chivers@mail.mcgill.ca>
}

\examples{

## Simulate data on which to build a model ##
x <- runif(100,-2,2)        ## Independent vars
p <- 1/(1+exp(-x))          ## Simple logistic probabilities
outcomes<-bs(p)             ## Binary outcomes

par(mfrow=c(2,2))
plot(x,p,ylim=c(0,1))



################################################################################
## First VMAPP for a model with only stochastic uncertainty (MLE of parameters)
model.fit <- glm(outcomes ~ x, family=binomial(logit)) ## Simple logistic Regression

## Simulate some validation data
x_val <- runif(100,-2,2)    ## New independent vars
p_val <- 1/(1+exp(-x_val))  ## Probabilities from the same model as before
outcomes_val <- bs(p_val)   ## Binary outcomes (validation data)

## Predict the probabilities associated with our validation data 
p_hat <- predict(model.fit, newdata = data.frame(x=x_val),type="response")

## Plot predictions on top of actual probabilities
points(x_val,p_hat,col='red')
legend('topleft',
    legend=c('Actual Probability','Predicted Probability'),
    col=c('black','red'),
    pch=1)

## Run VMAPP ##
VMAPP <- vmapp(d=outcomes_val, pred=p_hat)
VMAPP

## Plot VMAPP
plot(VMAPP,main='VMAPP -- No Parametric Uncertainty')


################################################################################
## Now VMAPP on bootstrapped predictions which encorporate parameter uncertainty

## Set up the same plot as before 
plot(x,p,ylim=c(0,1))
legend('topleft',
    legend=c('Actual Probability','Predicted Probability'),
    col=c('black','red'),
    pch=1)

n_boot <- 1000 ## Number of bootstraps
p_hat <- array(dim=c(n_boot,100)) ## An array to store predictions
for(i in 1:n_boot)
{
    boot_index <- sample(1:100,replace=TRUE) ## Sample original data with replacement
    x_boot <- x[boot_index]
    model.fit <- glm(outcomes[boot_index] ~ x_boot, family=binomial(logit))
    p_hat[i,] <- predict(model.fit, newdata = data.frame(x_boot=x_val),type="response")
    points(x_val,p_hat[i,],col=rgb(1,0,0,0.1),pch=20) ## Add bootstrapped preds to plot
}
points(x,p) ## Redraw actual probabilities on top

## Run VMAPP again, this time with parameter uncertainty included
VMAPP <- vmapp(d=outcomes_val, pred=p_hat)
VMAPP

plot(VMAPP)


}
