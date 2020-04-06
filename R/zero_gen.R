
#' Generates samples from the pms parametrization.
#'
#' Generates samples from the pms parametrization.
#'
#' @param n The number of observations to be drawn.
#' @param log_odds The log odds, the \code{p} parameter.
#' @param mu The conditional mean of the Gaussian part, the \code{m} parameter.
#' @param sigmasq The conditional variance of the Gaussian part, the \code{s} parameter.
#' @details
#' Generates samples from the Hurdle 1d distribution, with probability \eqn{1-p} for \eqn{y == 0}, and density \eqn{p*dnorm(x, m, sqrt(s))} for \eqn{y != 0}.
#' @return A vector of size \code{n}.
#' @examples
#' set.seed(1)
#' log_odds <- stats::rnorm(1)
#' mu <- stats::rnorm(1)
#' sigmasq <- abs(stats::rnorm(1))
#' x <- rhurdle1d_pms(1e5, log_odds, mu, sigmasq)
#' c(mean(x != 0), exp(log_odds) / (exp(log_odds) + 1))
#' hist(x[x != 0], freq=FALSE, breaks=200)
#' curve(dnorm(x, mu, sqrt(sigmasq)), from=min(x), to=max(x), add=TRUE, col="red")
#' @export
rhurdle1d_pms <- function(n, log_odds, mu, sigmasq){
  stats::rbinom(n, 1, exp(log_odds) / (exp(log_odds) + 1)) * stats::rnorm(n, mu, sqrt(sigmasq))
}

#' Generates samples from the abk parametrization.
#'
#' Generates samples from the abk parametrization.
#'
#' @param n The number of observations to be drawn.
#' @param a The \code{a} parameter.
#' @param b The \code{b} parameter.
#' @param k The \code{k} parameter.
#' @details
#' Generates samples from the Hurdle 1d distribution, whose log density with abk parametrization with respect to the sum of the Lebesgue measure and a point mass at 0 is
#' \eqn{av+by-y^2k/2-\log(1+\sqrt(2\pi/k)\exp(a+b^2/(2k)))}{a*v+b*y-y^2*k/2-log(1+sqrt(2pi/k)*exp(a+b^2/(2k)))}
#' where \code{v = (y != 0)}.
#' @return A vector of size \code{n}.
#' @examples
#' set.seed(1)
#' a <- stats::rnorm(1)
#' b <- stats::rnorm(1)
#' k <- abs(stats::rnorm(1))
#' x <- rhurdle1d_abk(1e5, a, b, k)
#' c(mean(x != 0), 1 - 1 / (1 + sqrt(2*pi/k) * exp(a+b^2/(2*k))))
#' hist(x[x != 0], freq=FALSE, breaks=200)
#' curve(exp(b*x-x^2*k/2)/(1+sqrt(2*pi/k)*exp(a+b^2/(2*k))),
#'    from=min(x), to=max(x), add=TRUE, col="red")
#' @export
rhurdle1d_abk <- function(n, a, b, k){
  rhurdle1d_pms(n, a + b^2/2/k - log(k/2/pi)/2, b/k, 1/k)
}

rbb <- Vectorize(stats::rbinom, "prob")
rnn <- Vectorize(stats::rnorm, "mean")

#' Generates one sample for each triplet of (log_odd, mu) with a fixed sigmasq under the pms parametrization.
#'
#' Generates one sample for each pair of (log_odd, mu) with a fixed sigmasq under the pms parametrization.
#'
#' @param log_odds A vector of log odds.
#' @param mus A vector of mus, of equal length as \code{log_odds}.
#' @param sigmasq A positive number, the sigma squared.
#' @details
#' Draws one example from the 1-d Hurdle with parameters \code{log_odds[i]}, \code{mus[i]}, \code{sigmasq}, i.e. a mixture of bernoulli and Gaussian with log odds \code{log_odds[i]} for the bernoulli and conditional mean \code{mus[i]} and variance \code{sigmasq} for the Gaussian part, for each \code{i}.
#' @return A vector of the same size as \code{log_odds}.
#' @examples
#' generate_one_each_logitp_mu_sigmasq(-3:3, 0:6/10, 1)
#' @export
generate_one_each_logitp_mu_sigmasq <- function(log_odds, mus, sigmasq){
  if (length(setdiff(c(length(log_odds),length(mus)), c(1))) > 1)
    stop("log_odds and mus must have the same length or be a scalar.")
  if (sigmasq <= 0) stop ("sigmasq must be positive.")
  rbb(1, 1, exp(log_odds) / (exp(log_odds) + 1)) * rnn(1, mus, sqrt(sigmasq))
}

#' Generates one sample for each triplet of (a, b) with a fixed k under the abk parametrization.
#'
#' Generates one sample for each pair of (a,b) with a fixed k under the abk parametrization.
#'
#' @param A A vector of parameter \code{a}s.
#' @param B A vector of parameter \code{b}s, of equal length as \code{A}.
#' @param k A positive number, the \code{k} parameter.
#' @details
#' Draws one example from the 1-d Hurdle with parameters \code{A[i]}, \code{B[i]}, \code{k} for each \code{i}.
#' @return A vector of the same size as \code{log_odds}.
#' @examples
#' generate_one_each_a_b_k(-3:3, 0:6/10, 1)
#' @export
generate_one_each_a_b_k <- function(A, B, k){
  if (length(setdiff(c(length(A),length(B)), c(1))) > 1)
    stop("A and B must have the same length or be a scalar.")
  if (k <= 0) stop ("k must be positive.")
  generate_one_each_logitp_mu_sigmasq(A+B^2/2/k-log(k/2/pi)/2, B/k, 1/k)
}

#' Returns the adjacency matrix of a chain, complete, or a lattice DAG.
#'
#' Returns the adjacency matrix of a chain, complete, or a lattice DAG.
#'
#' @param m A positive integer, the dimension.
#' @param mode One of \code{"chain"}, \code{"complete"}, \code{"lattice"}.
#' @param nrows A positive integer that divides \code{m}, used for \code{mode == "lattice"} as the number of rows in the graph.
#' @param shuffle A logical, indicates whether the nodes \code{1:m} should be shuffled. Defaults to \code{FALSE}.
#' @param seed A number, the seed passed to \code{set.seed()} if \code{shuffle == TRUE}.
#' @details
#' A chain graph has the form 1 -> 2 -> 3 -> ... -> m.
#' A complete graph has an edge from \code{i} to \code{j} for any \code{1 <= i < j <= m}.
#' A lattice graph first places the nodes in a \code{nrows x (m/nrows)} matrix, row by row. That is, the \code{i}-th node is on the \code{i \%/\% (m/nrows) + 1}-th row and the \code{i \%\% (m/nrows)}-th column (or the \code{(m/nrows)}-th column if \code{i} divides it). Then add an edge between any two adjacent nodes, pointing from the left to the right or from the node above to the node below.
#' If \code{shuffle == TRUE}, the rows/columns of the returning adjacency matrix are simultaneously shuffled.
#' @return The \code{m x m} adjacency matrix of the DAG.
#' @examples
#' make_dag(10, "chain", shuffle=FALSE)
#' make_dag(10, "chain", shuffle=TRUE)
#' make_dag(10, "complete", shuffle=FALSE)
#' make_dag(10, "complete", shuffle=TRUE)
#' make_dag(12, "lattice", nrows=3, shuffle=FALSE)
#' make_dag(12, "lattice", nrows=3, shuffle=TRUE)
#' @export
make_dag <- function(m, mode, nrows=0, shuffle=FALSE, seed=NULL){
  adj_mat <- matrix(0, m, m)
  if (mode == "chain") {
    adj_mat[cbind(c(1:(m-1)),c(2:m))] <- 1
  } else if (mode == "complete") {
    adj_mat[upper.tri(adj_mat)] <- 1
  } else if (mode == "lattice") {
    if (m %% nrows != 0) {stop("Nrows must divide m.")}
    if (nrows == 1 || m == nrows) {stop("Number of rows and columns must be > 1 for lattice")}
    ncols <- m / nrows
    for (r in 1:nrows)
      adj_mat[cbind((r-1)*ncols+1:(ncols-1), (r-1)*ncols+2:ncols)] <- 1
    for (c in 1:ncols)
      adj_mat[cbind(c+(0:(nrows-2))*ncols, c+(1:(nrows-1))*ncols)] <- 1
  } else {
    stop ("Mode not supported.")
  }
  if (shuffle) {
    if (!is.null(seed)) set.seed(seed)
    ord <- sample(m, m)
    adj_mat <- adj_mat[ord,]
    adj_mat <- adj_mat[,ord]
  }
  return (adj_mat)
}

