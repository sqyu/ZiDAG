# library(knitr); library(rmarkdown); library(devtools); library(roxygen2)
# devtools::document(); devtools::build(); devtools::install(); devtools::check()
# tools::package_native_routine_registration_skeleton(".")  # Copy as src/ZiDAG_init.c


#' Gradient of log_dhurdle1d_abk() (with scalar v and y) with respect to a, b, k.
#'
#' Gradient of \code{log_dhurdle1d_abk()} (with scalar \code{v} and \code{y}) with respect to \code{a}, \code{b}, \code{k}.
#'
#' @param V A logical vector, indicating if each entry in \code{Y} is non-zero, i.e. \code{V = (Y != 0)}.
#' @param Y A numerical vector of i.i.d. 1-d Hurdle random variables.
#' @param A A number or a vector of the same length \code{Y}, the \code{a} parameter(s).
#' @param B A number or a vector of the same length \code{Y}, the \code{b} parameter(s).
#' @param k A number, the \code{k} parameter.
#' @details
#' The \code{i}-th row of the returned matrix is the gradient of \code{log_dhurdle1d_abk(V[i], Y[i], A[i], B[i], k)} with respect to \code{A[i]}, \code{B[i]}, \code{k} if \code{A} and \code{B} are vectors, or with respect to \code{A} and \code{B} if they are scalars.
#' @return A matrix with 3 columns and number of rows equal to the length of \code{Y}.
#' @examples
#' if (!requireNamespace("numDeriv", quietly = TRUE))
#'    stop("Please install package \"numDeriv\".")
#' set.seed(1)
#' n <- 100
#' V <- rbinom(n, 1, 0.8)
#' Y <- rnorm(n) * V
#' A <- rnorm(n)
#' B <- rnorm(n)
#' k <- abs(rnorm(1))
#' grad_mat <- grad_a_b_k(V, Y, A, B, k)
#' numer_grad <- t(sapply(1:n,
#'    function(i){numDeriv::grad(function(x){
#'       log_dhurdle1d_abk(V[i], Y[i], x[1], x[2], x[3])}, c(A[i], B[i], k))}))
#' max(abs(grad_mat - numer_grad))
#'
#' a <- rnorm(1)
#' b <- rnorm(1)
#' grad_mat <- grad_a_b_k(V, Y, a, b, k)
#' numer_grad <- t(sapply(1:n,
#'    function(i){numDeriv::grad(function(x){
#'       log_dhurdle1d_abk(V[i], Y[i], x[1], x[2], x[3])}, c(a, b, k))}))
#' max(abs(grad_mat - numer_grad))
#' @export
grad_a_b_k <- function(V, Y, A, B, k){
  if (length(k) != 1) stop("k must be scalar.")
  if (length(V) != length(Y)) stop("V and Y must have the same length.")
  if (length(A) != length(B) || (!length(A) %in% c(1, length(Y))))
    stop("A and B must both be scalar or both have the same length as Y.")
  exp_part <- exp(A+B^2/2/k)
  return (cbind((V-1)+ 1 / (1 + exp_part * sqrt(2*pi / k)),
                Y - B/k + B / (k + exp_part * sqrt(2*pi * k)),
                B^2/2/k^2 + 1/2/k - Y^2/2 -
                  (B^2+k) / (2 * k^(3/2) * (sqrt(k) + exp_part * sqrt(2*pi)))
  ))
}

#' Gradient of log_dhurdle1d_abk() (with scalar v and y) with respect to a, b, k.
#'
#' Gradient of \code{log_dhurdle1d_abk()} (with scalar \code{v} and \code{y}) with respect to \code{a}, \code{b}, \code{k}.
#'
#' @param V A logical vector, indicating if each entry in \code{Y} is non-zero, i.e. \code{V = (Y != 0)}.
#' @param Y A numerical vector of i.i.d. 1-d Hurdle random variables.
#' @param A A number or a vector of the same length \code{Y}, the \code{a} parameter(s). Assumed to be results from \code{sum_A_mat()} run on \code{Vo} and \code{Yo}.
#' @param B A number or a vector of the same length \code{Y}, the \code{b} parameter(s). Assumed to be results from \code{sum_B_mat()} run on \code{Vo} and \code{Yo}.
#' @param Vo A numerical vector of the same dimension as \code{Yo} indicating if each entry in \code{Yo} is non-zero, i.e. \code{Vo = (Yo != 0)}.
#' @param Yo A numerical vector, a sample for the parent nodes (regressors).
#' @param k11 A number, the \code{k} parameter.
#' @param minus_Y A logical, argument as in \code{sum_B_mat()}.
#' @details
#' The derivative of \code{log_dhurdle_vec_abk(V[,1], Y[,1], A, B, k)} with respect to \code{aa}, \code{bb} and \code{k}, where \code{A=sum_A_mat(aa, V[,2:p], Y[,2:p])} and \code{B=sum_B_mat(bb, V[,2:p], Y[,2:p], minus_Y)}.
#' @examples
#' set.seed(1)
#' n <- 100; p <- 10
#' V <- matrix(rbinom(n*p, 1, 0.8), nrow=n, ncol=p)
#' Y <- matrix(rnorm(n*p) * V, nrow=n, ncol=p)
#' aa <- rnorm(2*p-1)
#' bb <- rnorm(2*p-1)
#' k <- abs(rnorm(1))
#' grad_vec <- grad_full_vec(V[,1], Y[,1], sum_A_mat(aa, V[,2:p], Y[,2:p]),
#'    sum_B_mat(bb, V[,2:p], Y[,2:p], minus_Y=TRUE), k, V[,2:p], Y[,2:p], minus_Y=TRUE)
#' numer_grad <- numDeriv::grad(function(x){
#'    log_dhurdle_vec_abk(V[,1], Y[,1], sum_A_mat(x[1:(2*p-1)], V[,2:p], Y[,2:p]),
#'       sum_B_mat(x[(2*p):(4*p-2)], V[,2:p], Y[,2:p], minus_Y=TRUE), x[4*p-1])}, c(aa, bb, k))
#' max(abs(grad_vec - numer_grad))
#'
#' grad_vec <- grad_full_vec(V[,1], Y[,1], sum_A_mat(aa, V[,2:p], Y[,2:p]),
#'    sum_B_mat(bb, V[,2:p], Y[,2:p], minus_Y=FALSE), k, V[,2:p], Y[,2:p], minus_Y=FALSE)
#' numer_grad <- numDeriv::grad(function(x){
#'    log_dhurdle_vec_abk(V[,1], Y[,1], sum_A_mat(x[1:(2*p-1)], V[,2:p], Y[,2:p]),
#'       sum_B_mat(x[(2*p):(4*p-2)], V[,2:p], Y[,2:p], minus_Y=FALSE), x[4*p-1])}, c(aa, bb, k))
#' max(abs(grad_vec - numer_grad))
#' @export
grad_full_vec <- function(V, Y, A, B, k11, Vo, Yo, minus_Y=TRUE){
  grad_tmp <- grad_a_b_k(V, Y, A, B, k11)
  n <- length(V)
  return (c(mean(grad_tmp[,1]), crossprod(grad_tmp[,1], Vo)/n, crossprod(grad_tmp[,1], Yo)/n,
            mean(grad_tmp[,2]), crossprod(grad_tmp[,2], Vo)/n, crossprod(grad_tmp[,2], Yo)/n * (1-2*minus_Y),
            mean(grad_tmp[,3])))
}

#' Fits MLE for 1-d Hurdle model under the abk parametrization.
#'
#' Fits MLE for 1-d Hurdle model under the abk parametrization.
#'
#' @param V A matrix of 0/1s, equal to Y != 0.
#' @param Y A data matrix of the same size as \code{V}.
#' @param left An integer between 1 and \code{ncol(Y)}. The index of the variable to be fit.
#' @param value_only If \code{TRUE}, returns the minimized negative log likelihood only. Defaults to \code{TRUE}.
#' @details
#' Fits an MLE for the 1-d Hurdle model to \code{Y[,left]} under the abk parametrization.
#' @return If \code{value_only == TRUE}, returns the minimized negative log likelihood only.
#' Otherwise, returns
#'   \item{nll}{A number, the minimized negative log likelihood.}
#'   \item{par}{A vector of 3 numbers, the fitted a, b, and k parameters.}
#'   \item{n}{An integer, the sample size.}
#'   \item{effective_df}{Number 3, the effective degree of freedom.}
#' @examples
#' Y <- matrix(rhurdle1d_abk(n=1e4, a=1.2, b=2.3, k=1.9), ncol=1)
#' V <- Y != 0
#' mle1d_abk(V, Y, 1, TRUE)
#' mle1d_abk(V, Y, 1, FALSE)
#' @export
mle1d_abk <- function(V, Y, left, value_only=TRUE){
  n <- nrow(V)
  V <- V[,left]; Y <- Y[,left]
  if (sum(V) == n){
    stop("V contains no 0 value. Stopped.")
  } else if (sum(V) == 0){
    stop("V contains no non-zero value. Stopped")
  } else if (sum(V) == 1){
    stop("V contains only one non-zero value and is not enough for estimation of k. Stopped.")
  }
  phat <- mean(V); muhat <- mean(Y[V]); sigmasqhat <- mean((Y[V]-muhat)^2)
  khat <- 1/sigmasqhat; bhat <- muhat*khat; ahat <- log(phat/(1-phat))+log(khat/2/pi)/2-bhat^2/2/khat
  nll <- -log_dhurdle1d_abk(V,Y,ahat,bhat,khat)
  if (value_only)
    return (nll)
  else{
    pars <- c(ahat, bhat, khat)
    names(pars) <- paste(c("a", "b", "k"), left, sep="")
    return (list("nll"=nll, "par"=pars, "n"=n, "effective_df"=3))
  }
}


#' Fits the Hurdle model assuming linear abk parametrization using \code{stats::optim}.
#'
#' Fits the Hurdle model assuming linear abk parametrization using \code{stats::optim}.
#'
#' @param V A matrix of 0/1s, equal to Y != 0.
#' @param Y A data matrix of the same size as \code{V}.
#' @param left An integer between 1 and \code{ncol(Y)}. The index of the variable to be fit.
#' @param right A vector of integers between 1 and \code{ncol(Y)} different from \code{left}. Indices of the "regressors".
#' @param maxit An integer, the maximum number of integers, argument to \code{stats::optim}.
#' @param tol A positive number, the tolerance passed as the \code{pgtol} argument to \code{stats::optim}.
#' @param runs A positive integer, number of reruns; if larger than \code{1}, \code{stats::optim} will be called multiple times with different initial values, and the fit that maximizes the log likelihood will be used.
#' @param value_only If \code{TRUE}, returns the minimized negative log likelihood only. Defaults to \code{TRUE}.
#' @param report An integer indicating verbosity, argument to \code{stats::optim}. If set to 1 no output is print during estimation.
#' @details
#' Fits the Hurdle model assuming linear abk parametrization, where \code{Y[,left]} conditional on \code{Y[,right]} is a 1-d Hurdle model with respect to the sum of the Lebesgue measure and a point mass at 0 with density
#' \eqn{av+by-y^2k/2-\log(1+\sqrt(2\pi/k)\exp(a+b^2/(2k)))}{a*v+b*y-y^2*k/2-log(1+sqrt(2pi/k)*exp(a+b^2/(2k)))},
#' with \code{a} and \code{b} both linear functions in \code{V[,right]} and \code{Y[,right]}.
#' @return If \code{value_only == TRUE}, returns the minimized negative log likelihood only. Otherwise, returns
#'   \item{nll}{A number, the minimized negative log likelihood.}
#'   \item{par}{A vector of length \code{4*length(right)+3}, the fitted parameters, in the other of: the intercept for the \code{a} (a scalar), linear coefficients on \code{V[,right]} for \code{a}, linear coefficients on \code{Y[,right]} for \code{a}, the intercept for the \code{b} (a scalar), linear coefficients on \code{V[,right]} for \code{b}, linear coefficients on \code{Y[,right]} for \code{b}.}
#'   \item{n}{An integer, the sample size.}
#'   \item{effective_df}{\code{4*length(right)+3}, the effective degree of freedom.}
#' @examples
#' m <- 3; n <- 1000
#' adj_mat <- make_dag(m, "complete")
#' dat <- gen_zero_dat(1, "abk", adj_mat, n, k_mode=1, min_num=10, gen_uniform_degree=1)
#' zi_fit_R(dat$V, dat$Y, 3, 1:2, maxit=1000, runs=2, value_only=TRUE, report=0)
#' zi_fit_R(dat$V, dat$Y, 3, 1:2, maxit=1000, runs=2, value_only=FALSE, report=0)
#' @export
zi_fit_R <- function(V, Y, left, right, maxit=200, tol=1e-8, runs=1, value_only=TRUE, report=10){
  r <- length(right)
  V1 <- V[,left]; Y1 <- Y[,left]
  Vo <- V[,right,drop=FALSE]; Yo <- Y[,right,drop=FALSE]
  nll_optim <- function(x){-log_dhurdle_vec_abk(V1,Y1,sum_A_mat(x[1:(2*r+1)],Vo,Yo),sum_B_mat(x[(2*r+2):(4*r+2)],Vo,Yo),x[4*r+3])}
  grad_optim <- function(x){-grad_full_vec(V1,Y1,sum_A_mat(x[1:(2*r+1)],Vo,Yo),sum_B_mat(x[(2*r+2):(4*r+2)],Vo,Yo),x[4*r+3],Vo,Yo)}
  best_res <- stats::optim(par = c(stats::rnorm(4*r+2),abs(stats::rnorm(1))),
                           fn = nll_optim,
                           gr = grad_optim,
                           method = "L-BFGS-B",
                           lower=c(rep(-Inf, 4*r+2),0.0001),
                           control=list(trace=3*(report!=0),REPORT=max(report,1),maxit=maxit,pgtol=tol))
  if (runs>1){
    for (run in 1:(runs-1)){
      res <- stats::optim(par = c(stats::rnorm(4*r+2),abs(stats::rnorm(1))),
                          fn = nll_optim,
                          gr = grad_optim,
                          method = "L-BFGS-B",
                          lower=c(rep(-Inf, 4*r+2),0.0001),
                          control=list(trace=3*(report!=0),REPORT=max(report,1),maxit=maxit,pgtol=tol))
      if (best_res$value > res$value){
        best_res <- res
      }
    }
  }
  if (value_only)
    return (best_res$value)
  else{
    if (r)
      names(best_res$par) <- c(paste("a",left,sep=""), paste("q",right,sep=""), paste("r",right,sep=""), paste("b",left,sep=""), paste("s",right,sep=""), paste("t",right,sep=""), paste("k",left,sep=""))
    else
      names(best_res$par) <- c(paste(c("a","b","k"),left,sep=""))
    return (list(nll=best_res$value, par=best_res$par, n=nrow(V), effective_df=4*r+3))
  }
}


#' Fits the Hurdle model assuming linear abk parametrization using C code.
#'
#' Fits the Hurdle model assuming linear abk parametrization using C code.
#'
#' @param V A matrix of 0/1s, equal to Y != 0.
#' @param Y A data matrix of the same size as \code{V}.
#' @param left An integer between 1 and \code{ncol(Y)}. The index of the variable to be fit.
#' @param right A vector of integers between 1 and \code{ncol(Y)} different from \code{left}. Indices of the "regressors".
#' @param maxit An integer, the maximum number of integers, argument to \code{stats::optim}.
#' @param runs A positive integer, number of reruns; if larger than \code{1}, \code{stats::optim} will be called multiple times with different initial values, and the fit that maximizes the log likelihood will be used.
#' @param value_only If \code{TRUE}, returns the minimized negative log likelihood only. Defaults to \code{TRUE}.
#' @param verbosity A positive integer, verbosity level. If set to 0 no output is print during estimation.
#' @param step_size A double, the step size of the first trial step. Defaults to 0.1.
#' @param lm_tol A double, accuracy of the line minimization.
#' @param epsabs A double, accuracy of the minimization.
#' @param maxsize A double, the final size of the simplex.
#' @param method A positive integer indicating the method to use. 5 is recommended and default.
#' \describe{
#'   \item{0}{Fletcher-Reeves conjugate gradient}
#'   \item{1}{Polak-Ribiere conjugate gradient}
#'   \item{2}{Vector Broyden-Fletcher-Goldfarb-Shanno method}
#'   \item{3}{Steepest descent algorithm}
#'   \item{4}{Nelder-Mead simplex}
#'   \item{5}{Vector Broyden-Fletcher-Goldfarb-Shanno method version 2}
#'   \item{6}{Simplex algorithm of Nelder and Mead version 2}
#'   \item{7}{Simplex algorithm of Nelder and Mead: random initialization}
#' }
#' @details
#' Fits the Hurdle model assuming linear abk parametrization, where \code{Y[,left]} conditional on \code{Y[,right]} is a 1-d Hurdle model with respect to the sum of the Lebesgue measure and a point mass at 0 with density
#' \eqn{av+by-y^2k/2-\log(1+\sqrt(2\pi/k)\exp(a+b^2/(2k)))}{a*v+b*y-y^2*k/2-log(1+sqrt(2pi/k)*exp(a+b^2/(2k)))},
#' with \code{a} and \code{b} both linear functions in \code{V[,right]} and \code{Y[,right]}.
#' @return If \code{value_only == TRUE}, returns the minimized negative log likelihood only. Otherwise, returns
#'   \item{nll}{A number, the minimized negative log likelihood.}
#'   \item{par}{A vector of length \code{4*length(right)+3}, the fitted parameters, in the other of: the intercept for the \code{a} (a scalar), linear coefficients on \code{V[,right]} for \code{a}, linear coefficients on \code{Y[,right]} for \code{a}, the intercept for the \code{b} (a scalar), linear coefficients on \code{V[,right]} for \code{b}, linear coefficients on \code{Y[,right]} for \code{b}.}
#'   \item{n}{An integer, the sample size.}
#'   \item{effective_df}{\code{4*length(right)+3}, the effective degree of freedom.}
#' @examples
#' m <- 3; n <- 1000
#' adj_mat <- make_dag(m, "complete")
#' dat <- gen_zero_dat(1, "abk", adj_mat, n, k_mode=1, min_num=10, gen_uniform_degree=1)
#' zi_fit_C(dat$V, dat$Y, 3, 1:2, maxit=1000, runs=2, value_only=TRUE, verbosity=0)
#' zi_fit_C(dat$V, dat$Y, 3, 1:2, maxit=1000, runs=2, value_only=FALSE, verbosity=0)
#' @export
#' @useDynLib ZiDAG optim
zi_fit_C <- function(V, Y, left, right, maxit=200, runs=1, value_only=TRUE, verbosity=0, step_size=0.1, lm_tol=0.01, epsabs=1e-10, maxsize=1e-10, method=5){
  r <- length(right)
  par <- c(stats::rnorm(4*r+2),abs(stats::rnorm(1)))
  WY <- cbind(V[,left],Y[,left],V[,right],Y[,right])
  best_res <- .C("optim", sample_size_input=as.integer(nrow(V)), num_parents_input=as.integer(r), value=as.double(0), grad=as.double(numeric(4*r+3)), par=as.double(par), WY=as.double(WY), step_size=as.double(step_size), tol=as.double(lm_tol), maxiter=as.integer(maxit), epsabs=as.double(epsabs), maxsize=as.double(maxsize), method=as.integer(method), verbosity=as.integer(verbosity))
  if (runs>1){
    for (run in 1:(runs-1)){
      par <- c(stats::rnorm(4*r+2),abs(stats::rnorm(1)))
      res <- .C("optim", sample_size_input=as.integer(nrow(V)), num_parents_input=as.integer(r), value=as.double(0), grad=as.double(numeric(4*r+3)), par=as.double(par), WY=as.double(WY), step_size=as.double(step_size), tol=as.double(lm_tol), maxiter=as.integer(maxit), epsabs=as.double(epsabs), maxsize=as.double(maxsize), method=as.integer(method), verbosity=as.integer(verbosity))
      if (best_res$value > res$value)
        best_res <- res
    }
  }
  if (value_only)
    return (best_res$value)
  else{
    if (r)
      names(best_res$par) <- c(paste("a",left,sep=""), paste("q",right,sep=""), paste("r",right,sep=""), paste("b",left,sep=""), paste("s",right,sep=""), paste("t",right,sep=""), paste("k",left,sep=""))
    else
      names(best_res$par) <- c(paste(c("a","b","k"),left,sep=""))
    return (list(nll=best_res$value, par=best_res$par, n=nrow(V), effective_df=4*r+3))
  }
}

# zi_fit_abk2 <- function(V, Y, left, right, use_C=TRUE, runs=1, value_only=TRUE, ...){
#   r <- length(right)
#   V1 <- V[,left]; Y1 <- Y[,left]
#   Vo <- V[,right,drop=FALSE]; Yo <- Y[,right,drop=FALSE]
#   nll_optim <- function(x){-log_dhurdle_vec_abk(V1,Y1,sum_A_mat(x[1:(2*r+1)],Vo,Yo),sum_B_mat(x[(2*r+2):(4*r+2)],Vo,Yo),abs(x[4*r+3]))}
#   grad_optim <- function(x){grad <- -grad_full_vec(V1,Y1,sum_A_mat(x[1:(2*r+1)],Vo,Yo),sum_B_mat(x[(2*r+2):(4*r+2)],Vo,Yo),abs(x[4*r+3]),Vo,Yo); if (x[4*r+3]<0) grad[4*r+3] <- -grad[4*r+3]; grad}
#   xinit <- c(stats::rnorm(4*r+2),-abs(stats::rnorm(1)))
#   best_res <- list(value=Inf)
#   for (run in 1:runs){
#     if (use_C) {
#       nll_grad_optim <- function(x){sumA <- sum_A_mat(x[1:(2*r+1)],Vo,Yo); sumB <- sum_B_mat(x[(2*r+2):(4*r+2)],Vo,Yo);
#         grad <- -grad_full_vec(V1,Y1,sum_A_mat(x[1:(2*r+1)],Vo,Yo),sum_B_mat(x[(2*r+2):(4*r+2)],Vo,Yo),abs(x[4*r+3]),Vo,Yo); if (x[4*r+3]<0) grad[4*r+3] <- -grad[4*r+3];
#         list(f=-log_dhurdle_vec_abk(V1,Y1,sumA,sumB,abs(x[4*r+3])), df=grad)
#       }
#       res <- gsl::multimin(x=xinit, f=nll_optim, df=grad_optim, fdf=nll_grad_optim, ...)
#       if (best_res$value > res$f)
#         best_res <- list(value=res$f, par=res$x)
#     } else {
#       res <- stats::optim(par=xinit, fn=nll_optim, gr=grad_optim,
#                           method="L-BFGS-B", lower=c(rep(-Inf, 4*r+2),0.0001), ...)
#       if (best_res$value > res$value)
#         best_res <- list(value=res$value, par=res$par)
#     }
#   }
#   if (value_only)
#     return (best_res$value)
#   else{
#     best_res$par[4*r+3] <- abs(best_res$par[4*r+3])
#     if (r)
#       names(best_res$par) <- c(paste("a",left,sep=""), paste("q",right,sep=""), paste("r",right,sep=""), paste("b",left,sep=""), paste("s",right,sep=""), paste("t",right,sep=""), paste("k",left,sep=""))
#     else
#       names(best_res$par) <- c(paste(c("a","b","k"),left,sep=""))
#     return (list(nll=best_res$value, par=best_res$par, n=nrow(V), effective_df=4*r+3))
#   }
# }



#' Fits a linear Hurdle conditional model with abk parametrization.
#'
#' Fits a linear Hurdle conditional model with abk parametrization.
#'
#' @param V A matrix of 0/1s, equal to Y != 0.
#' @param Y A data matrix of the same size as \code{V}.
#' @param left An integer between 1 and \code{ncol(Y)}. The index of the variable to be fit.
#' @param right A vector of integers between 1 and \code{ncol(Y)} different from \code{left}. Indices of the "regressors".
#' @param use_C A logical. Function calls \code{zi_fit_C()} if \code{use_C == TRUE}, and otherwise calls \code{zi_fit_R()}.
#' @param tol A positive number. If \code{use_C == FALSE} the tolerance passed as the \code{pgtol} argument to \code{stats::optim} (\code{tol} argument of \code{zi_fit_R()}). If \code{use_C == TRUE}, accuracy of the minimization (\code{epsabs} argument of \code{zi_fit_C()}).
#' @param maxit An integer, the maximum number of integers, argument to \code{stats::optim}.
#' @param runs A positive integer, number of reruns; if larger than \code{1}, \code{stats::optim} will be called multiple times with different initial values, and the fit that maximizes the log likelihood will be used.
#' @param value_only If \code{TRUE}, returns the minimized negative log likelihood only. Defaults to \code{TRUE}.
#' @param report An integer indicating verbosity. If \code{use_C == FALSE}, passed as \code{report} argument to \code{zi_fit_R()}, argument to \code{stats::optim}: if set to 1 no output is print during estimation. If \code{use_C == TRUE}, passed as \code{verbosity} argument to \code{zi_fit_C()}, argument to optim() in C code: a positive integer, verbosity level; if set to 0 no output is print during estimation.
#' @param step_size A double, the step size of the first trial step. Defaults to 0.1.
#' @param lm_tol A double, accuracy of the line minimization.
#' @param maxsize A double, the final size of the simplex.
#' @param method A positive integer indicating the method to use. 5 is recommended and default.
#' @details A linear Hurdle conditional model with abk parametrization for the \code{left} node given those in \code{right} has log density with respect to the sum of the Lebesgue measure and a point mass at 0 equal to (in terms of \code{y})
#' \eqn{av+by-y0^2k/2-\log(1+\sqrt(2\pi/k)\exp(a+b^2/(2k)))}{a*v+b*y-y^2*k/2-log(1+sqrt(2pi/k)*exp(a+b^2/(2k)))},
#' where \code{v = (y != 0)}, \code{k} is constant and \code{a} and \code{b} are linear functions in the values for \code{right} and their indicators.
#' This function fits such a model using \code{Y[,left]}, \code{Y[,right]} and \code{V[,right] = (Y[,right] != 0)}.
#'
#' If \code{right} is empty, fits an unconditional univariate Hurdle using \code{mle1d_abk()}.
#' If \code{use_C == TRUE}, calls \code{zi_fit_C()}.
#' If \code{use_C == FALSE}, calls \code{zi_fit_R()}.
#' @return If \code{value_only == TRUE}, returns the minimized negative log likelihood only. Otherwise, returns
#'   \item{nll}{A number, the minimized negative log likelihood.}
#'   \item{par}{A vector of length \code{4*length(right)+3}, the fitted parameters, in the other of: the intercept for the \code{a} (a scalar), linear coefficients on \code{V[,right]} for \code{a}, linear coefficients on \code{Y[,right]} for \code{a}, the intercept for the \code{b} (a scalar), linear coefficients on \code{V[,right]} for \code{b}, linear coefficients on \code{Y[,right]} for \code{b}.}
#'   \item{n}{An integer, the sample size.}
#'   \item{effective_df}{\code{4*length(right)+3}, the effective degree of freedom.}
#' @examples
#' m <- 3; n <- 1000
#' adj_mat <- make_dag(m, "complete")
#' dat <- gen_zero_dat(1, "abk", adj_mat, n, k_mode=1, min_num=10, gen_uniform_degree=1)
#' zi_fit_abk(dat$V, dat$Y, 3, 1:2, use_C=FALSE, maxit=1000, runs=2, value_only=TRUE, report=0)
#' zi_fit_abk(dat$V, dat$Y, 3, 1:2, use_C=TRUE, maxit=1000, runs=2, value_only=TRUE, report=0)
#' zi_fit_abk(dat$V, dat$Y, 3, 1:2, use_C=FALSE, maxit=1000, runs=2, value_only=FALSE, report=0)
#' zi_fit_abk(dat$V, dat$Y, 3, 1:2, use_C=TRUE, maxit=1000, runs=2, value_only=FALSE, report=0)
#' @export
zi_fit_abk <- function(V, Y, left, right, use_C=TRUE, tol=1e-8, maxit=1000, runs=2, value_only=TRUE, report=0, step_size=0.1, lm_tol=0.01, maxsize=1e-10, method=5){
  V <- V == 1
  if (length(right)==0){
    return (mle1d_abk(V, Y, left, value_only))
  }
  if (use_C){
    return (zi_fit_C(V, Y, left, right, maxit, runs=runs, value_only=value_only, verbosity=report, step_size=step_size, lm_tol=lm_tol, epsabs=tol, maxsize=maxsize, method=method))
  } else {
    return (zi_fit_R(V, Y, left, right, maxit, tol=tol, runs=runs, value_only=value_only, report=report))
  }
}

#' Fits a Hurdle polynomial with either abk or pms parametrization.
#'
#' Fits a Hurdle polynomial with either abk or pms parametrization.
#'
#' @param V A matrix of 0/1s, equal to Y != 0.
#' @param Y A data matrix of the same size as \code{V}.
#' @param left An integer between 1 and \code{ncol(Y)}. The index of the variable to be fit.
#' @param right A vector of integers between 1 and \code{ncol(Y)} different from \code{left}. Indices of the "regressors".
#' @param parametrization A string, either \code{"abk"} (canonical) or \code{"pms"} (moment).
#' @param value_only If \code{TRUE}, returns the minimized negative log likelihood only. Defaults to \code{TRUE}.
#' @param control A list of arguments to be passed to \code{zi_fit_abk()}, \code{zi_fit_pms()} and \code{zi_fit_pms_choose_degree()}.
#' @details
#' For \code{"abk"} parametrization, consult \code{?zi_fit_abk} for \code{control}.
#' For \code{"pms"} parametrization, consult \code{?zi_fit_pms} for model fitting using Hurdle polynomials of an exact degree, or \code{?zi_fit_pms_choose_degree} for automatically determining the degree given a maximum degree allowed.
#'
#' If \code{parametrization == "abk"}, calls \code{zi_fit_abk(V, Y, left, right, value_only)} with other arguments set to those in the \code{control} list.
#' If \code{parametrization == "pms" && is.null(control[["max_uniform_degree"]])}, calls \code{zi_fit_pms(V, Y, left, right, value_only)} with other arguments set to those in the \code{control} list.
#' If \code{parametrization == "pms" && !is.null(control[["max_uniform_degree"]])}, calls \code{zi_fit_pms_choose_degree(V, Y, left, right, value_only)} with other arguments set to those in the \code{control} list.
#' For the list of arguments for \code{zi_fit_abk()}, \code{zi_fit_pms()} and \code{zi_fit_pms_choose_degree()}, please refer to their documentation.
#' @return If \code{value_only == TRUE}, returns the minimized negative log likelihood only. Otherwise, returns
#'   \item{nll}{A number, the minimized negative log likelihood.}
#'   \item{par}{A vector of length \code{4*length(right)+3}, the fitted parameters, in the other of: the intercept for the \code{a} (a scalar), linear coefficients on \code{V[,right]} for \code{a}, linear coefficients on \code{Y[,right]} for \code{a}, the intercept for the \code{b} (a scalar), linear coefficients on \code{V[,right]} for \code{b}, linear coefficients on \code{Y[,right]} for \code{b}.}
#'   \item{n}{An integer, the sample size.}
#'   \item{effective_df}{\code{4*length(right)+3}, the effective degree of freedom.}
#' @examples
#' m <- 3; n <- 1000
#' adj_mat <- make_dag(m, "complete")
#' dat <- gen_zero_dat(1, "abk", adj_mat, n, k_mode=1, min_num=10, gen_uniform_degree=1)
#' zi_fit(dat$V, dat$Y, 3, 1:2, "abk", TRUE, list(use_C=FALSE, maxit=1000, runs=2, report=0))
#' zi_fit(dat$V, dat$Y, 3, 1:2, "abk", TRUE, list(use_C=TRUE, maxit=1000, runs=2, report=0))
#' zi_fit(dat$V, dat$Y, 3, 1:2, "abk", FALSE, list(use_C=FALSE, maxit=1000, runs=2, report=0))
#' zi_fit(dat$V, dat$Y, 3, 1:2, "abk", FALSE, list(use_C=TRUE, maxit=1000, runs=2, report=0))
#'
#' dat <- gen_zero_dat(1, "pms", adj_mat, 1000, k_mode=1, min_num=10, gen_uniform_degree=1)
#' zi_fit(dat$V, dat$Y, 3, 1:2, "pms", value_only=TRUE,
#'     list(p_V_degree=2, p_Y_degree=2, p_Y_V_degree=2, mu_V_degree=2, mu_Y_degree=2))
#' zi_fit(dat$V, dat$Y, 3, 1:2, "pms", value_only=FALSE,
#'     list(p_V_degree=2, p_Y_degree=2, p_Y_V_degree=2, mu_V_degree=2, mu_Y_degree=2))
#'
#' zi_fit(dat$V, dat$Y, 3, 1:2, "pms", value_only=TRUE,
#'     list(max_uniform_degree=2L, print_best_degree = TRUE))
#' zi_fit(dat$V, dat$Y, 3, 1:2, "pms", value_only=FALSE,
#'     list(max_uniform_degree=2L, print_best_degree = TRUE))
#' @export
zi_fit <- function(V, Y, left, right, parametrization="pms", value_only=TRUE, control=list()){
  control[["V"]] <- V
  control[["Y"]] <- Y
  control[["left"]] <- left
  if (length(right) == 0)
    control[["right"]] <- integer(0)
  else
    control[["right"]] <- right
  control[["value_only"]] <- value_only
  if (parametrization == "pms"){
    if (!is.null(control[["max_uniform_degree"]]))
      return (do.call("zi_fit_pms_choose_degree", control))
    else
      return (do.call("zi_fit_pms", control))
  } else if (parametrization == "abk") {
    return (do.call("zi_fit_abk", control))
  } else {
    stop ("parametrization must be either pms or abk.")
  }
}


