#.r.TailProbs and .minD from http://www.statslab.cam.ac.uk/~rds37/papers/r_concave_tail.R
.r.TailProbs <- function(eta, B, r) {
  # TailProbs returns a vector with the tail probability for each \tau = ceil{B*2\eta}/B + 1/B,...,1
  # We return 1 for all \tau = 0, 1/B, ... , ceil{B*2\eta}/B
  MAXa <- 100000
  MINa <- 0.00001
  s <- -1 / r
  etaB <- eta * B
  k_start <- (ceiling(2 * etaB) + 1)
  output <- rep(1, B)
  if (k_start > B)
    return(output)

  a_vec <- rep(MAXa, B)

  Find.a <-
    function(prev_a)
      stats::uniroot(
        Calc.a,
        lower = MINa,
        upper = prev_a,
        tol = .Machine$double.eps ^ 0.75
      )$root
  Calc.a <- function(a) {
    denom <- sum((a + 0:k) ^ (-s))
    num <- sum((0:k) * (a + 0:k) ^ (-s))
    num / denom - etaB
  }

  for (k in k_start:B)
    a_vec[k] <- Find.a(a_vec[k - 1])

  OptimInt <- function(a) {
    num <- (k + 1 - etaB) * sum((a + 0:(t - 1)) ^ (-s))
    denom <- sum((k + 1 - (0:k)) * (a + 0:k) ^ (-s))
    1 - num / denom
  }
  # NB this function makes use of several gloabl variables

  prev_k <- k_start
  for (t in k_start:B) {
    cur_optim <- rep(0, B)
    cur_optim[B] <- OptimInt(a_vec[B])
    if (prev_k <= (B - 1)) {
      for (k in prev_k:(B - 1))
        cur_optim[k] <-
          stats::optimize(
            f = OptimInt,
            lower = a_vec[k + 1],
            upper = a_vec[k],
            maximum  = TRUE
          )$objective
    }
    output[t] <- max(cur_optim)
    prev_k <- which.max(cur_optim)
  }
  return(output)
}

.minD <- function(theta, B, r = c(-1 / 2, -1 / 4)) {
  pmin(c(rep(1, B), .r.TailProbs(theta ^ 2, B, r[1])), .r.TailProbs(theta, 2 *
                                                                    B, r[2]))
}

#' Applies a fit function repeatedly to subsamples drawn as complementary pairs (CPSS).
#'
#' Applies a fit function repeatedly to subsamples drawn as complementary pairs (CPSS).
#'
#' @param n An integer >= 2, the sample size.
#' @param fitting_func A function that takes a vector of indices (a subvector of \code{1:n}) and returns some value that supports the "+" operation (e.g. an estimated adjacency matrix).
#' @param B An integer, number of pairs of subsamples to be drawn. Defaults to \code{50}.
#' @param parallel A logical that indicates whether \code{parallel::mclapply()} should be used for fitting. Otherwise, \code{lapply()} is used.
#' @param num_cores An integer, number of cores to be used if \code{parallel}; will be set to \code{parallel::detectCores() - 1} if \code{NULL}. Defaults to \code{NULL}.
#' @details
#' A complementary pair of sub-indices of \code{1:n} is a pair of disjoint subvectors \code{1:n} each with length \code{floor(n/2)}.
#' The function draws \code{B} such pairs, and apply \code{fitting_func()} to each of these \code{2*B} lists of sub-indices.
#' Finally, the sum of the \code{2*B} returned values/vectors/matrices from \code{fitting_func()} are returned.
#' Assuming \code{fitting_func} takes a list of indices and returns an adjacency matrix, each entry in \code{CPSS_fit(n, fitting_func, B)} returns the number of times that entry equals to \code{1} in the \code{2*B} estimated adjacency matrices.
#'
#' Note: The user is expected to call \code{CPSS_fit()}, followed by \code{CPSS_path()}, and then followed by \code{CPSS_control()}.
#' @return The sum of results from \code{fitting_func()} applied to \code{2*B} lists of sub-indices of \code{1:n}.
#' @examples
#' m <- 3; n <- 200; B <- 5
#' adj_mat <- ZiDAG::make_dag(m, mode = "chain")
#' d <- ZiDAG::gen_zero_dat(seed=1, gen_para="pms", adj_mat=adj_mat, n=n, gen_uniform_degree=1)
#' fitting_func <- function(indices) return (
#'     ZiDAG::ziGDS(V=d$V[indices,], Y=d$Y[indices,], parametrization="pms", verbose=FALSE,
#'     control=list("max_uniform_degree"=1L, "tol"=1e-8, "print_best_degree"=FALSE)))
#' freq_B <- ZiDAG::CPSS_fit(nrow(d$V), fitting_func, B=B, parallel=TRUE, num_cores=2)
#' @export
CPSS_fit <- function(n, fitting_func, B = 50, parallel=TRUE, num_cores=NULL) {
  if (n < 2) stop("The sanple size must be at least 2.")
  # fitting_func: takes a vector of sample indices and returns an adjacency matrix using those samples
  seq_n <- seq_len(n)
  indices_list <-
    lapply(1:B, function(i) {
      sample(seq_n, floor(n / 2))
    })
  if (n %% 2) indices_list <- c(indices_list, lapply(indices_list, function(l) {rest<-setdiff(seq_n, l); rest[-sample(1:length(rest), 1)]}))
  else indices_list <- c(indices_list, lapply(indices_list, function(l) {setdiff(seq_n, l)}))
  one_fit <- function(ii) {
    tmp <- fitting_func(indices_list[[ii]])
    cat(ii, " ")
    tmp
  }
  if (parallel) {
    if (is.null(num_cores)) num_cores <- parallel::detectCores() - 1
    else num_cores <- min(num_cores, parallel::detectCores() - 1)
    res_list <- parallel::mclapply(1:length(indices_list), one_fit, mc.cores = num_cores)
  } else {res_list <- lapply(1:length(indices_list), one_fit)}
  freq_B <- Reduce("+", res_list)
  return (freq_B)
}

#' Generates a path of estimates of adjacency matrix using results from CPSS_fit().
#'
#' Generates a path of estimates of adjacency matrix given the frequencies of each entry in the estimates using subsamples drawn as complementary pairs.
#'
#' @param freq_B A square matrix returned from \code{CPSS_fit()}, the frequency of each entry inn the adjacency matrix in the estimates using the \code{2*B} sets of subsamples.
#' @param B An integer, number of pairs of subsamples to be drawn, should be equal to the \code{B} used when calling \code{CPSS_fit()}. Defaults to \code{50}.
#' @param force_acyclic A logical, whether return estimated adj matrices that correspond to acyclic graphs only.
#' @details
#' Assuming \code{fitting_func} takes a list of indices and returns an adjacency matrix, each entry in \code{freq_B <- CPSS_fit(n, fitting_func, B)} returns the number of times that entry equals to \code{1} in the \code{2*B} estimated adjacency matrices.
#' \code{CPSS_path(freq_B, B)} returns a path of graph estimates by gradually decreasing the frequency threshold.
#' The path is a list containing \code{"graph"} (a list of adjacency matrices), \code{"freq"} (a vector of numbers in [0,1], the frequency threshold), and \code{"actual_fdr"} (a vector of numbers in [0,1], the actual FDR controlled), each corresponding to a final estimate. For each estimate, entries in \code{freq_B} that have frequency above \code{2*B*freq} are set to 1, and 0 otherwise, while \code{actual_fdr} is calculated as in Shah and Samworth (2012).
#' For example, if \code{freq_B} is [0 84 10; 0 0 71; 0 12 0] and \code{force_acyclic == TRUE}, the function returns a path with graphs [0 1 0; 0 0 0; 0 0 0] (\code{freq} 0.84), [0 1 0; 0 0 1; 0 0 0] (\code{freq} 0.71), [0 1 0; 0 0 1; 0 1 0] (\code{freq} 0.12), [0 1 1; 0 0 1; 0 1 0] (\code{freq} 0.10).
#' Fore more details, see Shah and Samworth (2012) and the documentation on \code{CPSS_fit()}.
#'
#' Note: The user is expected to call \code{CPSS_fit()}, followed by \code{CPSS_path()}, and then followed by \code{CPSS_control()}.
#' @return
#' A list containing \code{"graph"} (a list of adjacency matrices), \code{"freq"} (a vector of numbers in [0,1]), and \code{"actual_fdr"} (a vector of numbers in [0,1]), all of the same length; see Details.
#' @examples
#' m <- 3; n <- 100; B <- 10
#' adj_mat <- ZiDAG::make_dag(m, mode = "chain")
#' d <- ZiDAG::gen_zero_dat(seed=1, gen_para="pms", adj_mat=adj_mat, n=n, gen_uniform_degree=1)
#' fitting_func <- function(indices) return (
#'     ZiDAG::ziGDS(V=d$V[indices,], Y=d$Y[indices,], parametrization="pms", verbose=FALSE,
#'     control=list("max_uniform_degree"=1L, "tol"=1e-8, "print_best_degree"=FALSE)))
#' freq_B <- ZiDAG::CPSS_fit(nrow(d$V), fitting_func, B=B, parallel=TRUE, num_cores=2)
#' cpss_path <- CPSS_path(freq_B, B = B)
#' @export
CPSS_path <- function(freq_B, B = 50, force_acyclic=TRUE) {
  theta <- mean(freq_B[diag(nrow(freq_B)) == 0]) / B
  if (theta == 0){
    warning("All estimated graphs in all bootstrapping samples are empty.")
    return (list("graph"=list(freq_B), "freq"=c(0), "actual_fdr"=c(1)))
  }
  mind <- .minD(theta, B)
  counts <- setdiff(sort(unique(freq_B[diag(nrow(freq_B)) == 0]), decreasing = TRUE), 0)
  FDRs <- mind[counts]
  graphs <- list(); actual_fdrs <- c(); freqs <- c()
  for (i in 1:length(counts)) {
    graph <- 1 * (freq_B >= counts[i])
    if (force_acyclic && !ggm::isAcyclic(graph))
      break
    graphs[[length(graphs) + 1]] <- graph
    actual_fdrs <- c(actual_fdrs, mind[counts[i]])
    freqs <- c(freqs, counts[i] / 2 / B)
  }
  if (length(graphs) == 0) {
    warning("Only the empty graph can be produced. This might be because there are multiple edges that appear the most (equal) number of times in the ", 2 * B, " subsamples, and these edges form a cyclic graph.")
    return (list("graph"=list(matrix(0, nrow=nrow(freq_B), ncol=ncol(freq_B))), "freq"=c(max(freq_B)), "actual_fdr"=c(1)))
  }
  return (list("graph"=graphs, "freq"=freqs, "actual_fdr"=actual_fdrs))
}

#' Returns the graph with desired FDR/freq from a path of CPSS estimates from CPSS_path().
#'
#' Returns the graph with desired FDR/freq from a path of CPSS estimates from CPSS_path().
#'
#' @param path A list returned by \code{CPSS_path}
#' @param FDR A number between 0 and 1, the desired FDR (false discovery rate); exactly one of \code{FDR} and \code{freq} must be \code{NULL}.
#' @param freq A number between 0 and 1, the frequency threshold; exactly one of \code{FDR} and \code{freq} must be \code{NULL}.
#' @param B An integer, number of pairs of subsamples to be drawn, should be equal to the \code{B} used when calling \code{CPSS_path()}. Defaults to \code{50}.
#' @details
#' Returns the densest graph along the CPSS path that has FDR <= \code{FDR} or freq >= \code{freq}. If the resulting graph is empty, returns the sparsest non-empty graph.
#' Fore more details, see Shah and Samworth (2012) and the documentation on \code{CPSS_fit()} and \code{CPSS_path()}.
#'
#' Note: The user is expected to call \code{CPSS_fit()}, followed by \code{CPSS_path()}, and then followed by \code{CPSS_control()}.
#' @return A list associated with the selected estimated DAG.
#'    \item{graph}{The estimated adjacency matrix.}
#'    \item{actual_fdr}{The actual FDR controlled.}
#'    \item{freq}{The frequency threshold for the selected graph using \code{2*B} CPSS estimates.}
#' @examples
#' m <- 3; n <- 200; B <- 10
#' adj_mat <- ZiDAG::make_dag(m, mode = "chain")
#' d <- ZiDAG::gen_zero_dat(seed=1, gen_para="pms", adj_mat=adj_mat, n=n, gen_uniform_degree=1)
#' fitting_func <- function(indices) return (
#'     ZiDAG::ziGDS(V=d$V[indices,], Y=d$Y[indices,], parametrization="pms", verbose=FALSE,
#'     control=list("max_uniform_degree"=1L, "tol"=1e-8, "print_best_degree"=FALSE)))
#' freq_B <- ZiDAG::CPSS_fit(nrow(d$V), fitting_func, B=B, parallel=FALSE, num_cores=2)
#' cpss_path <- CPSS_path(freq_B, B = B)
#' CPSS_control(cpss_path, freq=0.5, B=B)
#' CPSS_control(cpss_path, FDR=0.2, B=B)
#' @export
CPSS_control <- function(path, FDR=NULL, freq=NULL, B = 50) {
  if (!xor(is.null(FDR), is.null(freq)))
    stop("One and only one of FDR and freq should be provided. FDR determines the desired FDR to be controlled, freq determines the selection threshold directly.")
  if (length(path) == 0) {
    warning("No DAG available in the path. Returning empty list.")
    return (list())
  }
  if (is.null(FDR)) {
    if (freq < 0 || freq > 1)
      stop("freq must be between 0 and 1.")
    if (freq > path$freq[1]) {
      warning("Requested freq larger than the freq (", path$freq[1], ") for the sparsest non-empty graph. Returning that graph.")
      i <- 1
    } else {i <- max(which(path$freq >= freq))}
  } else {
    if (FDR < 0 || FDR > 1)
      stop("FDR must be between 0 and 1.")
    if (FDR < path$actual_fdr[1]) {
      warning("Requested FDR smaller than the FDR (", path$actual_fdr[1], ") for the sparsest non-empty graph. Returning that graph.")
      i <- 1
    } else {i <- max(which(path$actual_fdr <= FDR))}
  }
  return (list("graph"=path$graph[[i]], "actual_fdr"=path$actual_fdr[i], "freq"=path$freq[i]))
}
