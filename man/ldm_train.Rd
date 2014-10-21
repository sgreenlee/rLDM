\name{ldm_train}
\alias{ldm_train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ldm_train(y, X, gamma, lambda_1, lambda_2, cost)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{lambda_1}{
%%     ~~Describe \code{lambda_1} here~~
}
  \item{lambda_2}{
%%     ~~Describe \code{lambda_2} here~~
}
  \item{cost}{
%%     ~~Describe \code{cost} here~~
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
function (y, X, gamma, lambda_1, lambda_2, cost) 
{
    m <- nrow(X)
    K <- matrix(nrow = m, ncol = m)
    for (i in 1:m) {
        for (j in 1:m) {
            if (i == j) 
                K[i, j] <- 1
            else if (i > j) 
                K[i, j] <- K[j, i]
            else K[i, j] <- exp(-gamma * sum((X[i, ] - X[j, ])^2))
        }
    }
    beta <- rep(0, m)
    Q <- 4 * lambda_1 * ((m * t(K) \%*\% K) - (K \%*\% y \%*\% t(y) \%*\% 
        t(K)))/(m^2) + K
    alpha <- lambda_2/m * solve(Q) \%*\% K \%*\% y
    Y <- diag(y)
    A <- solve(Q) \%*\% K \%*\% Y
    h <- diag(Y \%*\% K \%*\% solve(Q) \%*\% K \%*\% Y)
    BETA_EPSILON <- 10^-6
    n <- 0
    while (TRUE) {
        beta_old <- vector(length = m)
        del_beta <- vector(length = m)
        for (i in 1:m) {
            del_beta[i] <- (Y \%*\% K \%*\% alpha)[i] - 1
            beta_old[i] <- beta[i]
            beta[i] <- min(max(beta[i] - del_beta[i]/h[i], 0), 
                cost)
            alpha <- alpha + (beta[i] - beta_old[i]) * A[, i]
        }
        if (sum((beta - beta_old)^2) < BETA_EPSILON) 
            break
        if (1000 < (n <- n + 1)) 
            break
    }
    list(alpha = alpha, K = K)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
