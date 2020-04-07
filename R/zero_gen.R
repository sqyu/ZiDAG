
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


.adj_to_order <- function(adj_mat){
  ## Adjacency matrix to a topological order
  ## Can be optimized
  m <- ncol(adj_mat)
  in_degs <- colSums(adj_mat)
  ord <- integer(0)
  for (i in 1:m){
    next_node <- which.min(in_degs)
    ord <- c(ord, next_node)
    in_degs[next_node] <- Inf
    children_of_next <- which(adj_mat[next_node,] != 0)
    if (length(children_of_next))
      in_degs[children_of_next] <- in_degs[children_of_next] - 1
  }
  return (ord)
}

#' Generates a sample zero-inflated data matrix.
#'
#' Generates a sample zero-inflated data matrix.
#'
#' @param seed A number, the seed passed to \code{set.seed()}.
#' @param gen_para A string, the parametrization, \code{"abk"} or \code{"pms"}.
#' @param adj_mat A square matrix of 0/1s, the adjacency matrix.
#' @param n A positive integer, the number of samples.
#' @param k_mode One of \code{0}, \code{1}, \code{2}. If \code{k_mode == 1}, the \code{k} (if \code{gen_para == "abk"}) or \code{sigmasq} (if \code{gen_para == "pms"}) parameter in each conditional Hurdle distribution is set to 1. If \code{k_mode == 0}, the \code{k} is set to \code{1/j} (where \code{j} is the topological order of the current node) or \code{sigmasq} is set to \code{j}; if \code{k_mode == 2} they are set in the opposite way.
#' @param min_num A positive integer, the minimum number of 0s/1s in each column. If this is not satisfied, this random generator function will be rerun for up to 1000 times.
#' @param therun A positive integer for internal use. Defaults to 1.
#' @param gen_uniform_degree A positive integer, the maximum degree of the Hurdle polynomials in the conditional Hurdle distributions, defaults to \code{1}. Ignored if \code{parametrization == "abk"}.
#' @details
#' See the simulations section in the paper associated with this package for details.
#' @return A generated data, containing
#'     \item{V}{A matrix of 0/1s of dimension \code{ncol(adj_mat) x n}, equal to Y != 0.}
#'     \item{Y}{A data matrix of dimension \code{ncol(adj_mat) x n}.}
#' @examples
#' x <- gen_zero_dat(1, "abk", make_dag(3, "chain"), 1000)
#' x <- gen_zero_dat(1, "pms", make_dag(3, "chain"), 1000)
#' @export
gen_zero_dat <- function(seed, gen_para, adj_mat, n, k_mode=1, min_num=10, therun=1, gen_uniform_degree=1){
  if (therun == 1001) # If already tried to generate for 1000 times
    stop("Not able to generate good sample after 1000 runs: Some of the columns have fewer than ", min_num, " zeros or fewer than ", min_num, " nonzero values. Stopped.")
  if (therun==1) set.seed(seed) # Only set seed if it is the first trial
  diag(adj_mat) <- 0
  ord <- .adj_to_order(adj_mat)
  adj_mat <- adj_mat[ord,ord]
  m <- dim(adj_mat)[1]
  Y <- matrix(0,n,m)
  V <- matrix(NA,n,m)
  for (j in 1:m){
    jks <- which(adj_mat[,j] != 0)
    if (length(jks) == 0){
      Y[,j] <- rhurdle1d_abk(n,-log(2*pi)/2,0,1) #Equivalent to rhurdle1d_pms(n, 0, 0, 1) #rhurdle1d_abk(n,1,1,1)# rnorm(1,0,1), rnorm(1,0,1), 1)#rgamma(1,2,1))
      V[,j] <- (Y[,j] != 0)
    } else{
      if (gen_para == "abk"){
        if (k_mode == 0) {true_k_j <- 1/j #rgamma(1, 2, 1)#1
        } else if (k_mode == 1) {true_k_j <- 1
        } else if (k_mode == 2) {true_k_j <- j
        } else {stop("Bad k_mode specified.")}
        true_gg_j <- c(0, rep(1, 2 * length(jks))) #c(log(true_k_j/2/pi)/2, rep(g_strength, 2 * length(jks)))
        true_hh_j <- c(0,rep(1 * true_k_j, 2 * length(jks))) ## To stabilize the mean
        gs <- sum_A_mat(true_gg_j, V[,jks], Y[,jks])
        hs <- sum_B_mat(true_hh_j, V[,jks], Y[,jks])
        Y[,j] <- generate_one_each_a_b_k(
          (gs-mean(gs))/stats::sd(gs)-log(2*pi)/2,
          (hs-mean(hs))/stats::sd(hs),
          true_k_j
        )
      } else if (gen_para == "pms"){
        if (k_mode == 0) {true_sigmasq_j <- j  ## Inverse of the abk case
        } else if (k_mode == 1) {true_sigmasq_j <- 1
        } else if (k_mode == 2) {true_sigmasq_j <- 1/j
        } else {stop("Bad k_mode specified.")}
        if (gen_uniform_degree == 1) {
          true_logitp_j <- c(0, rep(1, 2 * length(jks)))
          true_mu_j <- c(0, rep(1, 2 * length(jks)))
          log_odds <- sum_A_mat(true_logitp_j, V[,jks], Y[,jks])
          mus <- sum_B_mat(true_mu_j, V[,jks], Y[,jks])
          Y[,j] <- generate_one_each_logitp_mu_sigmasq(
            (log_odds-mean(log_odds))/stats::sd(log_odds),
            (mus-mean(mus))/stats::sd(mus),
            true_sigmasq_j
          )
        } else {
          rhs_design <- full_design1(
            Vo=V[,jks,drop=FALSE], Yo=Y[,jks,drop=FALSE], right=1:length(jks),
            V_degree=gen_uniform_degree,
            Y_degree=gen_uniform_degree,
            Y_V_degree=gen_uniform_degree)
          true_logitp_j <- c(0, 1-grepl("[\\*\\^]", colnames(rhs_design))*9/10)
          true_mu_j <- c(0, (1 - 2 * grepl("Y", colnames(rhs_design)))*(1-grepl("[\\*\\^]", colnames(rhs_design))*9/10)) ## Conform with the abk case
          log_odds <- sum_A_mat(true_logitp_j, Vo=rhs_design, Yo=NULL)
          mus <- sum_B_mat(true_mu_j, Vo=rhs_design, Yo=NULL, minus_Y=FALSE)
          Y[,j] <- generate_one_each_logitp_mu_sigmasq(
            (log_odds-mean(log_odds))/stats::sd(log_odds),
            (mus-mean(mus))/stats::sd(mus),
            true_sigmasq_j
          )
        }
      } else {stop("Wrong gen_para specified.")}
      V[,j] <- (Y[,j] != 0)
      if (sum(V[,j]) < min_num || sum(1-V[,j]) < min_num)
        return (gen_zero_dat(seed=seed, gen_para=gen_para, adj_mat=adj_mat, n=n, k_mode=k_mode, min_num=min_num, therun=therun+1, gen_uniform_degree=gen_uniform_degree))
    }
  }
  rev_ord <- order(ord)
  V <- V[,rev_ord]
  Y <- Y[,rev_ord]
  return (list("V"=V, "Y"=Y))
}


