library(bnlearn)
library(parallel)

#r.TailProbs and minD from http://www.statslab.cam.ac.uk/~rds37/papers/r_concave_tail.R
r.TailProbs <- function(eta, B, r) {
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
      uniroot(
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
          optimize(
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

minD <- function(theta, B, r = c(-1 / 2, -1 / 4)) {
  pmin(c(rep(1, B), r.TailProbs(theta ^ 2, B, r[1])), r.TailProbs(theta, 2 *
                                                                    B, r[2]))
}

chooseTau <- function(theta, B, FDR=NULL, freq=NULL) {
  if (!xor(is.null(FDR), is.null(freq)))
    stop("One and only one of FDR and freq should be provided.")
  if (!is.null(FDR) && (FDR < 0 || FDR > 1))
    stop("FDR must be between 0 and 1.")
  if (!is.null(freq) && (freq < 0 || freq > 1))
    stop("freq must be between 0 and 1.")
  mind <- minD(theta, B)
  if (is.null(FDR)) {
    freq <- ceiling(freq * 2 * B)
    return (c(freq, mind[freq]))
  } else {
    if (mind[length(mind)] > FDR) {
      warning("Smallest CPSS bound ", mind[length(mind)],
           " larger than requested FDR ", FDR, ". FDR changed to ", mind[length(mind)], ".")
      FDR <- mind[length(mind)] + 1e-7
    }
    first_index <- min(which(mind <= FDR))
    return (c(first_index, mind[first_index]))
  }
}

moral_graph <- function(adj_mat) {
  rownames(adj_mat) <- colnames(adj_mat) <- paste(1:ncol(adj_mat))
  bn_obj <- bnlearn::empty.graph(rownames(adj_mat))
  amat(bn_obj) <- adj_mat
  mor <- amat(bnlearn::moral(bn_obj))
  rownames(mor) <- colnames(mor) <- NULL
  return (mor)
}

CPSS <- function(n, fitting_func, FDR=NULL, freq=NULL, B = 50, parallel=TRUE) {
  ### TODO!!: If result not DAG, increase threshold until result is DAG
  # fitting_func: takes a vector of sample indices and returns an adjacency matrix using those samples
  if (!xor(is.null(FDR), is.null(freq)))
    stop("One and only one of FDR and freq should be provided. FDR determines the desired FDR to be controlled, freq determines the selection threshold directly.")
  seq_n <- seq_len(n)
  indices_list <-
    lapply(1:B, function(i) {
      sample(seq_n, floor(n / 2))
    })
  indices_list <-
    c(indices_list, lapply(indices_list, function(l) {
      setdiff(seq_n, l)
    }))
  one_fit <- function(ii) {
    tmp <- moral_graph(fitting_func(indices_list[[ii]]))
    tmp[lower.tri(tmp)] <- 0
    cat(ii, " ") # Progress
    tmp
  }
  if (parallel) {
    num_cores <- 1 + (parallel::detectCores() >= 2)
    num_cores <- max(num_cores, parallel::detectCores() - 1)
    res_list <- parallel::mclapply(1:length(indices_list), one_fit, mc.cores = num_cores)
  } else {res_list <- lapply(1:length(indices_list), one_fit)}
  p <- nrow(res_list[[1]])
  freq_B <- Reduce("+", res_list)
  total_freq <- sum(freq_B)
  if (total_freq == 0){
    warning("All estimated graphs in all bootstrapping samples are empty.")
    return (list("graph"=freq_B, "actual_fdr"=1))
  }
  q <- total_freq / length(indices_list)
  theta <- q / (p * (p - 1) / 2)
  min_freq_and_actual_fdr <- chooseTau(theta, B, FDR=FDR, freq=freq)
  min_freq <- min_freq_and_actual_fdr[1]
  actual_fdr <- min_freq_and_actual_fdr[2]
  return (list("graph"=1 * (freq_B >= min_freq), "actual_fdr"=actual_fdr))
}
