#' Calculates mean(log(1+exp(xs))) of a vector xs.
#'
#' Calculates \code{mean(log(1+exp(xs)))} of a vector \code{xs}.
#'
#' @param xs A vector of random numbers.
#' @details Numerically stable calculation of \code{mean(log(1+exp(xs)))}.
#' @return A number.
#' @examples
#' logsumexp_vec(10^(-5:5))
#' logsumexp_vec(rep(0, 5))
#' @export
logsumexp_vec <- function(xs){
  return (mean(pmax(xs, 0) + log(1 + exp(-abs(xs)))))
}

#' Calculates the mean log probability of a vector of i.i.d. 1-d Hurdle random variables under the abk parametrization.
#'
#' Calculates the mean log probability of a vector of i.i.d. 1-d Hurdle random variables under the abk parametrization.
#'
#' @param V A logical vector, indicating if each entry in \code{Y} is non-zero, i.e. \code{V = (Y != 0)}.
#' @param Y A numerical vector of i.i.d. 1-d Hurdle random variables.
#' @param a A number, the \code{a} parameter.
#' @param b A number, the \code{b} parameter.
#' @param k A number, the \code{k} parameter.
#' @details
#' The log probability of a sample \code{y} from the 1-d Hurdle model with abk parametrization with respect to the sum of the Lebesgue measure and a point mass at 0 is
#' \eqn{av+by-y^2k/2-\log(1+\sqrt(2\pi/k)\exp(a+b^2/(2k)))}{a*v+b*y-y^2*k/2-log(1+sqrt(2pi/k)*exp(a+b^2/(2k)))}
#' where \code{v = (y != 0)}. The function calculates the mean of log probabilities over a vector of i.i.d samples.
#' @return A number, the mean log probability.
#' @examples
#' y <- rhurdle1d_abk(n=100, a=1, b=2, k=3)
#' log_dhurdle1d_abk(y != 0, y, a=1, b=2, k=3)
#' @export
log_dhurdle1d_abk <- function(V, Y, a, b, k){
  if (length(V) != length(Y))
    stop("V, Y must have the same length.")
  if (k <= 0) stop("k must be positive.")
  mean(a*V+b*Y-Y^2*k/2) - logsumexp_vec(a+b^2/2/k - log(k/2/pi)/2)
}

#' Calculates the mean log probability of a vector of i.i.d. 1-d Hurdle random variables under the pms parametrization.
#'
#' Calculates the mean log probability of a vector of i.i.d. 1-d Hurdle random variables under the pms parametrization.
#'
#' @param V A logical vector, indicating if each entry in \code{Y} is non-zero, i.e. \code{V = (Y != 0)}.
#' @param Y A numerical vector of i.i.d. 1-d Hurdle random variables.
#' @param p A number, the \code{p} parameter.
#' @param mu A number, the \code{mu} parameter.
#' @param sigmasq A number, the \code{sigmasq} parameter.
#' @details
#' The log probability of a sample \code{y} from the 1-d Hurdle model with pms parametrization with respect to the sum of the Lebesgue measure and a point mass at 0 is
#' \eqn{\log(1-p)}{log(1-p)} if it is equal to 0, or \eqn{\log(p)-(y-mu)^2/2/sigmasq}{log(p)-(y-mu)^2/2/sigmasq} otherwise.
#' The function calculates the mean of log probabilities over a vector of i.i.d samples.
#' @return A number, the mean log probability.
#' @examples
#' y <- rhurdle1d_pms(n=100, log_odds=1.2, mu=2.3, sigmasq=3.4)
#' log_dhurdle1d_pms(y != 0, y, p=exp(1.2)/(exp(1.2)+1), mu=2.3, sigmasq=3.4)
#' @export
log_dhurdle1d_pms <- function(V, Y, p, mu, sigmasq){ # 1-d hurdle
  # V, Y: vectors
  V_on <- mean(V)
  return (log(1-p)*(1-V_on)+log(p)*V_on-
            sum((Y[V]-mu)^2)/2/sigmasq/length(V)
          -log(2*pi*sigmasq)/2*V_on)
}


#' Calculates the mean log probability of a vector of 1-d Hurdle random variables under the abk parametrization, each variable with different a and b and shared k.
#'
#' Calculates the mean log probability of a vector of 1-d Hurdle random variables under the abk parametrization, each variable with different \code{a} and \code{b} and shared \code{k}.
#'
#' @param V A logical vector, indicating if each entry in \code{Y} is non-zero, i.e. \code{V = (Y != 0)}.
#' @param Y A numerical vector of i.i.d. 1-d Hurdle random variables.
#' @param A A numerical vector of the same size as \code{V} and \code{Y}, the \code{a} parameters for each sample from \code{Y}.
#' @param B A numerical vector of the same size as \code{V} and \code{Y}, the \code{b} parameters for each sample from \code{Y}.
#' @param k A positive number, the common \code{k} parameter for all samples from \code{Y}.
#' @details
#' The log probability of a sample \code{y} from the 1-d Hurdle model with abk parametrization with respect to the sum of the Lebesgue measure and a point mass at 0 is
#' \eqn{av+by-y^2k/2-\log(1+\sqrt(2\pi/k)\exp(a+b^2/(2k)))}{a*v+b*y-y^2*k/2-log(1+sqrt(2pi/k)*exp(a+b^2/(2k)))}
#' where \code{v = (y != 0)}. The function calculates the mean of the log probabilities over a vector of samples with different \code{a} and \code{b} parameters and the same \code{k}.
#' This is used for the multivariate case where in the conditional 1-d Hurdle model the \code{a} and \code{b} are functions in other (parent) nodes/variables.
#' @return A number, the mean log probability.
#' @examples
#' set.seed(1)
#' A <- rnorm(100)
#' B <- rnorm(100)
#' y <- generate_one_each_a_b_k(A=A, B=B, k=2)
#' log_dhurdle_vec_abk(y != 0, y, A=A, B=B, k=2)
#' @export
log_dhurdle_vec_abk <- function(V, Y, A, B, k){
  if (length(unique(length(V), length(Y), length(A), length(B))) > 1)
    stop("V, Y, A, B must have the same length.")
  if (k <= 0) stop("k must be positive.")
  c(crossprod(V, A) + crossprod(Y, B)) / length(V) - mean(Y^2)*k/2 - logsumexp_vec(A+B^2/2/k-log(k/2/pi)/2)
}


#' Calculate the a parameter in the conditional Hurdle density of the child node given the parents for one sample, assuming a is linear in the parents and their indicators.
#'
#' Calculate the \code{a} parameter in the conditional Hurdle density of the child node given the parents for one sample, assuming \code{a} is linear in the parents and their indicators.
#'
#' @param aa A numerical vector of length \code{2*length(yo)+1}. See details.
#' @param vo A numerical vector of the same length as \code{yo} indicating if each entry in \code{yo} is non-zero, i.e. \code{vo = (yo != 0)}.
#' @param yo A numerical vector, a sample for the parent nodes (regressors).
#' @details
#' Assuming \code{aa = c(a1, a2, a3)} where \code{a1} is a scalar, \code{a2} and \code{a3}, the function returns \code{a1 + sum(a2*vo) + sum(a3*yo)}.
#' This is used for calculating the \code{a} parameter of the conditional Hurdle density of the child node given the parent nodes, assuming the \code{a} parameter is linear in the parents and their indicators.
#' @return A number, the mean log probability.
#' @examples
#' set.seed(1)
#' vo <- stats::rbinom(100, 1, 0.8)
#' yo <- stats::rnorm(100) * vo
#' sum_a(stats::rnorm(201), vo, yo)
#' @export
sum_a <- function(aa, vo, yo){
  #aa : c(g11, g1, h1o) ## 1, m-1, m-1
  c(crossprod(aa, c(1, vo, yo)))
}

#' Calculate the b parameter in the conditional Hurdle density of the child node given the parents for one sample, assuming b is linear in the parents and their indicators.
#'
#' Calculate the \code{b} parameter in the conditional Hurdle density of the child node given the parents for one sample, assuming \code{b} is linear in the parents and their indicators.
#'
#' @param bb A numerical vector of length \code{2*length(yo)+1}. See details.
#' @param vo A numerical vector of the same length as \code{yo} indicating if each entry in \code{yo} is non-zero, i.e. \code{vo = (yo != 0)}.
#' @param yo A numerical vector, a sample for the parent nodes (regressors).
#' @param minus_Y A logical, see details.
#' @details
#' Assuming \code{bb = c(b1, b2, b3)} where \code{b1} is a scalar, \code{b2} and \code{b3}, the function returns \code{b1 + sum(b2*vo) - sum(b3*yo)} if \code{minus_Y == TRUE}, otherwise \code{b1 + sum(b2*vo) + sum(b3*yo)}.
#' This is used for calculating the \code{a} parameter of the conditional Hurdle density of the child node given the parent nodes, assuming the \code{b} parameter is linear in the parents and their indicators.
#' @return A number, the mean log probability.
#' @examples
#' set.seed(1)
#' vo <- stats::rbinom(100, 1, 0.8)
#' yo <- stats::rnorm(100) * vo
#' sum_b(stats::rnorm(201), vo, yo, minus_Y=TRUE)
#' @export
sum_b <- function(bb, vo, yo, minus_Y=TRUE){
  if (minus_Y)
    c(crossprod(bb, c(1, vo, -yo)))
  else
    c(crossprod(bb, c(1, vo, yo)))
}

#' Calculate the a parameter in the conditional Hurdle density of the child node given the parents for multiple samples, assuming a is linear in the parents and their indicators.
#'
#' Calculate the \code{a} parameter in the conditional Hurdle density of the child node given the parents for multiple samples, assuming \code{a} is linear in the parents and their indicators.
#'
#' @param aa A numerical vector of length \code{2*ncol(Yo)+1}. See details.
#' @param Vo A numerical vector of the same dimension as \code{Yo} indicating if each entry in \code{Yo} is non-zero, i.e. \code{Vo = (Yo != 0)}.
#' @param Yo A numerical vector, a sample for the parent nodes (regressors).
#' @details
#' Matrix version of \code{sum_a()}. See examples.
#' @return A number, the mean log probability.
#' @examples
#' set.seed(1)
#' n <- 100; p_others <- 10
#' Vo <- matrix(stats::rbinom(n*p_others, 1, 0.8), nrow=n, ncol=p_others)
#' Yo <- matrix(stats::rnorm(n*p_others) * Vo, nrow=n, ncol=p_others)
#' aa <- rnorm(2*p_others+1)
#' sum_A_mat(aa, Vo, Yo) - sapply(1:n, function(i){sum_a(aa, Vo[i,], Yo[i,])})
#' @export
sum_A_mat <- function(aa, Vo, Yo){
  c(tcrossprod(aa, cbind(1, Vo, Yo)))
}

#' Calculate the b parameter in the conditional Hurdle density of the child node given the parents for multiple samples, assuming b is linear in the parents and their indicators.
#'
#' Calculate the \code{b} parameter in the conditional Hurdle density of the child node given the parents for multiple samples, assuming \code{b} is linear in the parents and their indicators.
#'
#' @param bb A numerical vector of length \code{2*ncol(Yo)+1}. See details.
#' @param Vo A numerical vector of the same dimension as \code{Yo} indicating if each entry in \code{Yo} is non-zero, i.e. \code{Vo = (Yo != 0)}.
#' @param Yo A numerical vector, a sample for the parent nodes (regressors).
#' @param minus_Y A logical, see details in \code{sum_b()}.
#' @details
#' Matrix version of \code{sum_b()}. See examples.
#' @return A number, the mean log probability.
#' @examples
#' set.seed(1)
#' n <- 100; p_others <- 10
#' Vo <- matrix(stats::rbinom(n*p_others, 1, 0.8), nrow=n, ncol=p_others)
#' Yo <- matrix(stats::rnorm(n*p_others) * Vo, nrow=n, ncol=p_others)
#' bb <- rnorm(2*p_others+1)
#' sum_B_mat(bb, Vo, Yo) - sapply(1:n, function(i){sum_b(bb, Vo[i,], Yo[i,])})
#' @export
sum_B_mat <- function(bb, Vo, Yo, minus_Y=TRUE){
  if (minus_Y)
    c(tcrossprod(bb, cbind(1, Vo, -Yo)))
  else
    c(tcrossprod(bb, cbind(1, Vo, Yo)))
}

#' Convert parameters in the abk parametrization to pms parametrization or vice versa.
#'
#' Gradient of \code{log_dhurdle1d_abk()} (with scalar \code{v} and \code{y}) with respect to \code{a}, \code{b}, \code{k}.
#'
#' @param paras A vector of 3 numbers. If \code{pms_to_abk == TRUE}, \code{paras} should contain the log odds and the conditional mean and variance for the Gaussian part. If \code{pms_to_abk == FALSE}, it should contain the \code{a}, \code{b} and \code{k} parameters.
#' @param pms_to_abk A logical.
#' @details
#' The log density of a sample \code{y} from the 1-d Hurdle model with abk parametrization with respect to the sum of the Lebesgue measure and a point mass at 0 is
#' \eqn{av+by-y^2k/2-\log(1+\sqrt(2\pi/k)\exp(a+b^2/(2k)))}{a*v+b*y-y^2*k/2-log(1+sqrt(2pi/k)*exp(a+b^2/(2k)))}
#' The pms parametrization composes of the log odds (p) and the conditional mean (m) and variance (s) of the Gaussian part, i.e. the density is \eqn{1-p} if \eqn{v == 0}, or \eqn{p*dnorm(x, m, sqrt(s))} if \eqn{v == 1}.
#' The two can be related as \code{s = 1 / k}, \code{m = b / k}, \code{p = a + b^2/2/k - log(k/2/pi)/2}. Also see Page 7 of McDavid (2019).
#' @return A vector of 3 numbers.
#' @examples
#' set.seed(1)
#' convert_para(c(0, 0, 1), TRUE) # (-log(2*pi)/2, 0, 1)
#' convert_para(c(-log(2*pi)/2, 0, 1), FALSE) # (0, 0, 1)
#' paras <- c(rnorm(1), rnorm(1), abs(rnorm(1)))
#' convert_para(convert_para(paras, TRUE), FALSE) - paras
#' convert_para(convert_para(paras, FALSE), TRUE) - paras
#'
#' n <- 100
#' v <- stats::rbinom(100, 1, 0.8)
#' y <- stats::rnorm(100) * v
#' a <- paras[1]; b <- paras[2]; k <- paras[3]
#' tmp <- convert_para(paras, FALSE)
#' log_odds <- tmp[1]; mu <- tmp[2]; sigmasq <- tmp[3]
#' prob <- exp(log_odds) / (exp(log_odds)+1)
#' likelihoods <- numeric(100)
#' likelihoods[v == 0] <- 1 - prob
#' likelihoods[v == 1] <- prob * dnorm(y[v == 1], mu, sqrt(sigmasq))
#' mean(log(likelihoods)) - log_dhurdle1d_abk(v, y, a, b, k)
#' @export
convert_para <- function(paras, pms_to_abk){
  if (!pms_to_abk %in% c(FALSE, TRUE)) {stop("pms_to_abk must be TRUE or FALSE.")}
  if (length(paras) != 3) {stop("paras must have length 3.")}
  if (pms_to_abk) {
    log_odds <- paras[1]; mu <- paras[2]; sigmasq <- paras[3]
    k <- 1/sigmasq; b <- mu*k; a <- log_odds + log(k/2/pi)/2 - b^2/2/k
    return (c(a, b, k))
  } else {
    a <- paras[1]; b <- paras[2]; k <- paras[3]
    return (c(a + b^2/2/k - log(k/2/pi)/2, b/k, 1/k))
  }
}


