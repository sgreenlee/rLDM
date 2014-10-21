\name{ldm.default}
\alias{ldm.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
<filler>%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ldm.default(y, X, kernel = "radial", cost = 10, gamma = 1, degree = 2, lambda_1 = 1/4, lambda_2 = 1/4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{kernel}{
%%     ~~Describe \code{kernel} here~~
}
  \item{cost}{
%%     ~~Describe \code{cost} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{lambda_1}{
%%     ~~Describe \code{lambda_1} here~~
}
  \item{lambda_2}{
%%     ~~Describe \code{lambda_2} here~~
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
function (y, X, kernel = "radial", cost = 10, gamma = 1, degree = 2, 
    lambda_1 = 1/4, lambda_2 = 1/4) 
{
    kernel_type <- pmatch(kernel, c("linear", "polynomial", "radial", 
        "sigmoid")) - 1
    if (!is.factor(y)) 
        stop("y must be a factor")
    if (length(levels(y)) != 2) 
        stop("y must be a two-level factor: multi-class not implemented yet")
    else y <- c(-1, 1)[unclass(y)]
    X <- scale(X)
    tuning_params <- c(cost, lambda_1, lambda_2, degree, gamma)
    n <- dim(X)[1]
    m <- dim(X)[2]
    dyn.load("rldm_train.so")
    gmat <- .C("makeGramMatrix", as.double(X), as.integer(n), 
        as.integer(m), G = double(n^2), as.integer(kernel_type), 
        as.double(tuning_params))
    G <- matrix(gmat$G, nrow = n)
    Gy <- G \%*\% y
    Y <- diag(y)
    GY <- G \%*\% Y
    Q <- 4 * lambda_2 * ((t(G) \%*\% G)/n - Gy \%*\% t(Gy)/n^2) + 
        G
    invQ <- .Internal(La_solve(Q, diag(n), .Machine$double.eps))
    A <- invQ \%*\% GY
    h <- diag(t(GY) \%*\% A)
    alpha <- lambda_2/n * invQ \%*\% Gy
    beta <- double(n)
    beta_old <- double(n)
    del_beta <- double(n)
    cret <- .C("coordDescent", as.integer(n), as.double(A), as.double(GY), 
        alpha = as.double(alpha), beta = as.double(beta), as.double(beta_old), 
        as.double(del_beta), as.double(h), as.double(cost))
    dyn.unload("rldm_train.so")
    ret <- list(alpha = cret$alpha)
    ret$G <- G
    ret$fitted <- sign(G \%*\% ret$alpha)
    ret$cost <- cost
    ret$lambda_1 <- lambda_1
    ret$lambda_2 <- lambda_2
    ret$gamma <- gamma
    ret$degree <- degree
    class(ret) <- "ldm"
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
