########################## BIC scores ##########################

#' Calculates the BIC score from a fit result.
#'
#' Calculates the BIC score from a fit result.
#'
#' @param fit_res A list returned from \code{zi_fit()} or \code{.res_from_model_fits()} with \code{value_only == FALSE}.
#' @details Returns \code{bic_score(fit_res$nll, fit_res$n, fit_res$effective_df)}.
#' @return A number, the BIC score.
#' @examples
#' Y <- matrix(rnorm(500) * rbinom(500, 1, 0.7), nrow=50)
#' res <- zi_fit(Y != 0, Y, 1, 2:5, value_only=FALSE)
#' bic_from_fit(res)
#' bic_score(res$nll, res$n, res$effective_df)
#' @export
bic_from_fit <- function(fit_res){
  return (bic_score(fit_res$nll, fit_res$n, fit_res$effective_df))
}

#' Calculates the BIC score.
#'
#' Calculates the BIC score.
#'
#' @param avgnll Average negative log likelihood over the sample.
#' @param n The sample size.
#' @param numpar The degree of freedom.
#' @details Returns \code{2*avgnll+log(n)/n*numpar}.
#' @return A number, the BIC score.
#' @examples
#' Y <- matrix(rnorm(500) * rbinom(500, 1, 0.7), nrow=50)
#' res <- zi_fit(Y != 0, Y, 1, 2:5, value_only=FALSE)
#' bic_score(res$nll, res$n, res$effective_df)
#' bic_from_fit(res)
#' @export
bic_score <- function(avgnll, n, numpar){
  return (2*avgnll+log(n)/n*numpar)
}


########################## GDS ##########################

.hasPath <- function(v1, v2, in_edges, cur_adj_mat){
  if (length(in_edges[[v2]]) == 0) return (FALSE) # v2 has no parent
  to_search <- c(v2)
  visited <- numeric(ncol(cur_adj_mat))
  visited[v2] <- TRUE
  in_edges[[v2]] <- .rm_elt(in_edges[[v2]], v1) # Do not consider v1->v2 itself; this is safe since in_edges[[v2]] is always nonempty
  cur_adj_mat[v1, v2] <- 0
  while (length(to_search)) {
    v <- to_search[1]
    if (cur_adj_mat[v1, v])
      return (TRUE)
    to_search <- to_search[-1]
    for (vv in in_edges[[v]])
      if (!visited[vv]) {
        to_search <- c(to_search, vv);
        visited[vv] <- TRUE
      }
  }
  return (FALSE)
}

.rm_elt <- function(vec, elt) {return (vec[vec != elt])}
.parent_to_str <- function(parents) {ifelse(length(parents) == 0, "0", paste(parents, collapse=","))}

.add_edge <- function(V, Y, parametrization, in_edges, degrees, cur_adj_mat, nodewise_score, add_dec, maxDegree, maxInDegree, fixedGaps, no_check_DAG, score_cache, control=list()){
  if (!requireNamespace("ggm", quietly = TRUE))
    stop("Please install package \"ggm\".")
  best_in <- best_out <- -1
  best_decrease <- 0
  best_new_score <- Inf
  best_new_right <- NULL
  m <- ncol(V); n <- nrow(V)
  for (ii in 1:m){
    if (degrees[ii] < maxDegree && length(in_edges[[ii]]) < maxInDegree){
      for (oi in setdiff(1:m, c(ii, in_edges[[ii]]))){
        if (degrees[oi] < maxDegree && (is.null(fixedGaps) || !fixedGaps[oi,ii])) {
          if (!(ii %in% in_edges[[oi]]) && (no_check_DAG || !.hasPath(ii, oi, in_edges, cur_adj_mat))){
            cur_adj_mat[oi, ii] <- 1
            new_right <- sort(c(oi, in_edges[[ii]]))
            if (!is.na(add_dec[oi, ii])){
              dec <- add_dec[oi, ii]
              new_score <- nodewise_score[ii] - dec
            } else{
              right_str <- .parent_to_str(new_right)
              if (is.null(score_cache))
                new_score <- bic_from_fit(zi_fit(V, Y, left=ii, right=new_right, parametrization=parametrization, value_only=FALSE, control=control))
              else {
                new_score <- score_cache[[ii]][[right_str]]
                if (is.null(new_score)) {
                  new_score <- bic_from_fit(zi_fit(V, Y, left=ii, right=new_right, parametrization=parametrization, value_only=FALSE, control=control))
                  score_cache[[ii]][[right_str]] <- new_score
                }
              }
              dec <- nodewise_score[ii] - new_score
              add_dec[oi,ii] <- dec
            }
            if (dec > best_decrease){
              best_in <- ii; best_out <- oi;  best_new_score <- new_score
              best_decrease <- dec; best_new_right <- new_right
            }
            cur_adj_mat[oi, ii] <- 0
          }
        }
      }
    }
  }
  return (list("best_in"=best_in, "best_out"=best_out, "best_new_right"=best_new_right, "best_decrease"=best_decrease,
               "best_new_score"=best_new_score, "add_dec"=add_dec, "score_cache"=score_cache))
}

.remove_edge <- function(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, rem_dec, score_cache, control=list()){
  best_in <- best_out <- -1
  best_decrease <- 0
  best_new_score <- Inf
  best_new_right <- NULL
  n <- nrow(V); m <- ncol(V)
  for (ii in 1:m){
    for (oi in in_edges[[ii]]){
      new_right <- .rm_elt(in_edges[[ii]], oi)
      if (!is.na(rem_dec[oi,ii])){
        dec <- rem_dec[oi,ii]
        new_score <- nodewise_score[ii] - dec
      } else {
        right_str <- .parent_to_str(new_right)
        if (is.null(score_cache))
          new_score <- bic_from_fit(zi_fit(V, Y, left=ii, right=new_right, parametrization=parametrization, value_only=FALSE, control=control))
        else {
          new_score <- score_cache[[ii]][[right_str]]
          if (is.null(new_score)) {
            new_score <- bic_from_fit(zi_fit(V, Y, left=ii, right=new_right, parametrization=parametrization, value_only=FALSE, control=control))
            score_cache[[ii]][[right_str]] <- new_score
          }
        }
        dec <- nodewise_score[ii] - new_score
        rem_dec[oi,ii] <- dec
      }
      if (dec > best_decrease){
        best_in <- ii; best_out <- oi;  best_new_score <- new_score
        best_decrease <- dec; best_new_right <- new_right
      }
    }
  }
  return (list("best_in"=best_in, "best_out"=best_out, "best_new_right"=best_new_right, "best_decrease"=best_decrease,
               "best_new_score"=best_new_score, "rem_dec"=rem_dec, "score_cache"=score_cache))
}

.turn_edge <- function(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, add_dec, rem_dec, maxInDegree, no_check_DAG, score_cache, control=list()){
  if (!requireNamespace("ggm", quietly = TRUE))
    stop("Please install package \"ggm\".")
  best_in <- best_out <- -1
  best_decrease <- 0
  best_new_score_in <- best_new_score_out <- Inf
  best_new_right_ii <- best_new_right_oi <- NULL
  m <- ncol(V); n <- nrow(V)
  for (ii in 1:m){
    for (oi in in_edges[[ii]]){
      if (length(in_edges[[oi]]) >= maxInDegree) # Skip if indegree of out node is max
        next
      if (no_check_DAG || !.hasPath(oi, ii, in_edges, cur_adj_mat)){
        cur_adj_mat[ii, oi] <- 1
        cur_adj_mat[oi, ii] <- 0
        new_right_ii <- .rm_elt(in_edges[[ii]], oi)
        if (!is.na(rem_dec[oi,ii])){
          dec_in <- rem_dec[oi,ii]
          new_score_in <- nodewise_score[ii] - dec_in
        } else {
          right_str <- .parent_to_str(new_right_ii)
          if (is.null(score_cache))
            new_score_in <- bic_from_fit(zi_fit(V, Y, left=ii, right=new_right_ii, parametrization=parametrization, value_only=FALSE, control=control))
          else {
            new_score_in <- score_cache[[ii]][[right_str]]
            if (is.null(new_score_in)) {
              new_score_in <- bic_from_fit(zi_fit(V, Y, left=ii, right=new_right_ii, parametrization=parametrization, value_only=FALSE, control=control))
              score_cache[[ii]][[right_str]] <- new_score_in
            }
          }
          dec_in <- nodewise_score[ii] - new_score_in
          rem_dec[oi,ii] <- dec_in
        }
        new_right_oi <- sort(c(ii,in_edges[[oi]]))
        if (!is.na(add_dec[ii,oi])){
          dec_out <- add_dec[ii,oi]
          new_score_out <- nodewise_score[oi] - dec_out
        } else {
          right_str <- .parent_to_str(new_right_oi)
          if (is.null(score_cache))
            new_score_out <- bic_from_fit(zi_fit(V, Y, left=oi, right=new_right_oi, parametrization=parametrization, value_only=FALSE, control=control))
          else {
            new_score_out <- score_cache[[oi]][[right_str]]
            if (is.null(new_score_out)) {
              new_score_out <- bic_from_fit(zi_fit(V, Y, left=oi, right=new_right_oi, parametrization=parametrization, value_only=FALSE, control=control))
              score_cache[[oi]][[right_str]] <- new_score_out
            }
          }
          dec_out <- nodewise_score[oi] - new_score_out
          add_dec[ii,oi] <- dec_out
        }
        if (dec_in+dec_out > best_decrease){
          best_in <- ii; best_out <- oi
          best_new_score_in <- new_score_in
          best_new_score_out <- new_score_out
          best_decrease <- dec_in+dec_out
          best_new_right_ii <- new_right_ii
          best_new_right_oi <- new_right_oi
        }
        cur_adj_mat[oi, ii] <- 1
        cur_adj_mat[ii, oi] <- 0
      }
    }
  }
  return (list("best_in"=best_in, "best_out"=best_out, "best_decrease"=best_decrease,
               "best_new_score_in"=best_new_score_in, "best_new_score_out"=best_new_score_out,
               "add_dec"=add_dec, "rem_dec"=rem_dec, "score_cache"=score_cache,
               "best_new_right_ii"=best_new_right_ii, "best_new_right_oi"=best_new_right_oi))
}

#' Randomly generates DAGs of specified sparsity.
#'
#' Randomly generates DAGs of specified sparsity.
#'
#' @param m An integer, the number of nodes.
#' @param init_spars A vector of numbers between 0 and 1, inclusive. It is assumed that \code{as.integer(m*(m-1)/2*(1-init_spars))} does not contain duplicates.
#' @param init_eachs A positive integer or a vector of positive integers of the same length as \code{init_spars}.
#' @details
#' Randomly generates DAGs of specified sparsity.
#' For each \code{spar} in \code{init_spars}, the function generates a complete graph, then randomly remove \code{as.integer(m*(m-1)/2*(1-spar))} edges. Finally, the ordering of the \code{m} nodes are randomized.
#' For each \code{spar = init_spars[i]} this is repeated \code{init_eachs} (if integer) or \code{init_eachs[i]} (if \code{init_eachs} is a vector of integers) times; duplicate graphs are avoided and at most 200 attempts will be made when generating each new graph.
#' Note that for \code{spar = 0} there will only be one empty graph returned.
#' @return
#' A list of adjacency matrices each of shape \code{m} by \code{m}.
#' @examples
#' rand_init_dags(4, c(0, 0.25, 0.5, 0.75, 1), 5)
#' rand_init_dags(10, c(0, 0.25, 0.5, 0.75, 1), 5)
#' @export
rand_init_dags <- function(m, init_spars, init_eachs) {
  if (length(unique(as.integer(m*(m-1)/2*(1-init_spars)))) != length(init_spars))
    stop("Different init_spars resulted in the same number of edges. Please double check.")
  rand_dag <- function(m, spars, seed) {
    set.seed(seed)
    adj_mat <- ZiDAG::make_dag(m, mode = "complete")
    adj_mat[sample(which(upper.tri(adj_mat)), as.integer(m*(m-1)/2*(1-spars)))] <- 0
    ord <- sample(m)
    return (adj_mat[ord, ord])
  }
  rand_mats <- list()
  if (any(init_eachs %% 1 | init_eachs < 1))
    stop("init_eachs must only contain positive integers.")
  if (length(init_eachs) == 1)
    init_eachs <- rep(init_eachs, length(init_spars))
  else if (length(init_eachs) != length(init_spars))
    stop("init_eachs must be a single integer or a vector of integers of the same length as init_spars.")
  for (spar_i in 1:length(init_spars)) {
    spar <- init_spars[spar_i]
    init_each <- init_eachs[spar_i]
    if (spar == 0)
      rand_mats <- c(rand_mats, list(matrix(0, nrow=m, ncol=m)))
    else {
      if (init_each == 1)
        rand_mats <- c(rand_mats, list(rand_dag(m, spar, spar)))
      else {
        sub_list <- list()
        seed <- 0
        # Generate next matrix of the same sparsity that does not equal to any of the
        # existing ones
        for (i in 1:init_each) {
          this_trial <- 0
          # Give up if already tried 200 times
          while (this_trial <= 200) {
            samp_mat <- rand_dag(m, spar, seed)
            if (i == 1 || !any(sapply(sub_list, function(s){all(s == samp_mat)}))) break
            seed <- seed + 1
            this_trial <- this_trial + 1
          }
          if (this_trial == 200) break
          sub_list <- c(sub_list, list(samp_mat))
        }
        rand_mats <- c(rand_mats, sub_list)
      }
    }
  }
  return (rand_mats)
}

#' Greedy DAG search for DAGs for zero-inflated data based on Hurdle conditionals.
#'
#' Greedy DAG search for DAGs for zero-inflated data based on Hurdle conditionals.
#'
#' @param V A matrix of 0/1s, equal to Y != 0.
#' @param Y A data matrix of the same size as \code{V}.
#' @param parametrization A string, either \code{"abk"} (canonical) or \code{"pms"} (moment).
#' @param fixedGaps A logical square matrix of number of rows/columns equal to number of columns of \code{Y}. Similar to \code{pcalg::gds()}, if its \code{[i,j]} is \code{TRUE}, the result is guaranteed to have no edge between nodes \code{i} and \code{j}.
#' @param maxDegree A positive integer, the maximum degree of a node.
#' @param maxInDegree A positive integer, the maximum in-degree of a node.
#' @param init An list of adjacency matrices or a single adjacency matrix, each of number of rows/columns equal to number of columns of \code{Y}, used as the initial graph. Defaults to \code{NULL} representing the empty graph. If a list is provided, each matrix will be used to initialize the algorithm once and the result with the lowest BIC will be returned.
#' @param verbose A logical, whether to print intermediate steps.
#' @param no_check_DAG A logical. If \code{FALSE}, check on whether the resulting graph is acyclic will not be performed in the forward and edge reversing stages.
#' @param return_BIC A logical. If \code{TRUE}, returns the estimated matrix along with the minimized BIC.
#' @param cache_scores A logical. If \code{TRUE}, nodewise BIC scores will be cached to speed up estimation, at the cost of a memory overhead.
#' @param control A list passed to \code{zi_fit()}. Please consult \code{?zi_fit}.
#' @details
#' Performs greedy DAG search for DAGs for zero-inflated data based on Hurdle conditionals.
#' See Chickering (2003) or \code{?pcalg::gds} for details. However, unlike the Gaussian case in their implementation that returns an equivalence class of DAGs, the DAG estimated by this function is exact.
#' @return If \code{return_BIC} is \code{FALSE}, returns an adjacency matrix; otherwise returns a list with \code{"graph"} element being the estimated matrix and \code{"BIC"} being the estimated BIC.
#' @examples
#' m <- 4; n <- 1000
#' adj_mat <- ZiDAG::make_dag(m, mode = "chain", shuffle=FALSE)
#' d <- ZiDAG::gen_zero_dat(seed=1, gen_para="pms", adj_mat=adj_mat, n=n, gen_uniform_degree=1)
#' est <- ZiDAG::ziGDS(V=d$V, Y=d$Y, parametrization="pms", verbose=FALSE,
#'     control=list("max_uniform_degree"=1L, "tol"=1e-8, "print_best_degree"=FALSE))
#' adj_mat == est
#'
#' d <- ZiDAG::gen_zero_dat(seed=1, gen_para="abk", adj_mat=adj_mat, n=n, gen_uniform_degree=1)
#' est <- ZiDAG::ziGDS(V=d$V, Y=d$Y, parametrization="abk", verbose=FALSE,
#'     control=list(use_C=TRUE, maxit=1000, runs=2, report=0))
#' adj_mat == est
#'
#' est <- ZiDAG::ziGDS(V=d$V, Y=d$Y, parametrization="pms", verbose=FALSE,
#'     control=list("max_uniform_degree"=1L, "tol"=1e-8, "print_best_degree"=FALSE))
#' adj_mat == est
#'
#' est <- ZiDAG::ziGDS(V=d$V, Y=d$Y, parametrization="pms", verbose=FALSE,
#'     control=list("max_uniform_degree"=2L, "tol"=1e-8, "print_best_degree"=FALSE))
#' adj_mat == est
#' @export
ziGDS <- function(V, Y, parametrization, fixedGaps=NULL, maxDegree=Inf, maxInDegree=Inf, init=NULL, verbose=FALSE, no_check_DAG=FALSE, return_BIC=FALSE, cache_scores=TRUE, control=list()){
  if (!requireNamespace("ggm", quietly = TRUE))
    stop("Please install package \"ggm\".")
  if (!all(dim(V)==dim(Y))) stop("V and Y must have the same dimension.")
  m <- ncol(V)
  maxDegree <- min(m-1, maxDegree)
  maxInDegree <- min(m-1, maxInDegree)
  if (maxInDegree > maxDegree)
    stop("maxInDegree must be smaller than maxDegree.")
  if (!is.null(fixedGaps)) {
    if (length(dim(fixedGaps)) != 2 || !all(dim(fixedGaps) == m))
      stop("If provided, fixedGaps must be a ", m, "x", m, " matrix.")
    if (!all(fixedGaps*1 == 0 || fixedGaps*1 == 1))
      stop("If provided, fixedGaps must be a logical or 0/1 matrix.")
    if (sum(abs(fixedGaps - t(fixedGaps))) != 0) {
      stop("fixedGaps must be a symmetric matrix.")
    }
  }
  if (cache_scores)
    score_cache <- lapply(1:m, function(i){hash::hash()})
  else
    score_cache <- NULL
  # If init is not a list of init matrices
  if (!is.list(init))
    init <- list(init)
  est_adj <- NULL; best_BIC <- Inf
  for (init_i in 1:length(init)) {
    if (verbose) message("Initial matrix #", init_i, ":")
    tmp_est <- .ziGDS_help(V, Y, parametrization=parametrization,
                          fixedGaps=fixedGaps, maxDegree=maxDegree, maxInDegree=maxInDegree,
                          init=init[[init_i]], verbose=verbose, no_check_DAG=no_check_DAG,
                          score_cache=score_cache, control=control)
    if (tmp_est$BIC < best_BIC) {
      est_adj <- tmp_est$graph; best_BIC <- tmp_est$BIC
    }
    score_cache <- tmp_est$score_cache
  }
  if (return_BIC)
    return (list("graph"=est_adj, "BIC"=best_BIC))
  return (est_adj)
}

.ziGDS_help <- function(V, Y, parametrization, fixedGaps, maxDegree, maxInDegree, init, verbose, no_check_DAG, score_cache, control){
  m <- ncol(V)
  # If init is NNULL
  if (is.null(init)){
    cur_adj_mat <- matrix(0,m,m)
    in_edges <- lapply(1:m, function(i){integer(0)})
  } else if (is.list(init)){
    if (length(init) != m)
      stop("Initial (in-)edge list must have length ",m,".")
    cur_adj_mat <- matrix(0,m,m)
    for (in_v in 1:m){
      if (any(init[[in_v]] < 1 | init[[in_v]] > m))
        stop("Vertices in edge list must be between 1 and ",m,".")
      if (in_v %in% init[[in_v]]){
        warning("Self edge for", in_v, "ignored.")
        init[[in_v]] <- sort(.rm_elt(init[[in_v]], in_v))
      } else {
        init[[in_v]] <- sort(init[[in_v]])
      }
      cur_adj_mat[init[[in_v]], in_v] <- 1
    }
    in_edges <- init
  } else {
    if (length(dim(init))!=2 || !all(dim(init)==c(m,m)))
      stop("Initial adjacency matrix must have dimension ",m,"*",m,".")
    if (is.logical(init)){init <- init*1}
    else if (!all(unique(c(init)) %in% c(0,1))){
      stop("Initial adjacency matrix must have 0 or 1s only.")
    }
    if (any(diag(init)!=0)){
      if (verbose)
        warning("Self edge(s) ignored.\n")
      diag(init) <- 0
    }
    cur_adj_mat <- init
    in_edges <- lapply(1:m, function(i){which(init[,i]!=0)})
  }
  if (!is.null(fixedGaps) && any(cur_adj_mat[fixedGaps] != 0))
    stop("No edge excluded in fixedGaps should be present in the initial adjacency matrix.")
  if (!no_check_DAG && !ggm::isAcyclic(cur_adj_mat))
    stop("Initial graph must be acyclic.")
  degrees <- rowSums(cur_adj_mat) + colSums(cur_adj_mat)  # only need the degree, not the list
  nodewise_score <- sapply(1:m, function(i){bic_from_fit(zi_fit(V, Y, left=i, right=in_edges[[i]], parametrization=parametrization, value_only=FALSE, control=control))})
  if (!is.null(score_cache))
    for (i in 1:m)
      score_cache[[i]][[.parent_to_str(in_edges[[i]])]] <- nodewise_score[i]
  old_best_sum <- Inf
  best_sum <- sum(nodewise_score)
  iter_count <- 0
  add_dec <- rem_dec <- matrix(NA,m,m)
  while (best_sum < old_best_sum-1e-9){
    old_best_sum <- best_sum
    iter_count <- iter_count + 1
    if (verbose){message("Iteration ", iter_count, ": BIC score ", old_best_sum)}
    while (TRUE){
      res_add <- .add_edge(V, Y, parametrization, in_edges, degrees, cur_adj_mat, nodewise_score, add_dec, maxDegree, maxInDegree, fixedGaps, no_check_DAG, score_cache, control)
      if (res_add$best_in == -1){
        break
      }
      nodewise_score[res_add$best_in] <- res_add$best_new_score
      best_sum <- best_sum - res_add$best_decrease
      add_dec <- res_add$add_dec
      cur_adj_mat[res_add$best_out, res_add$best_in] <- 1
      add_dec[,res_add$best_in] <- NA; rem_dec[,res_add$best_in] <- NA
      add_dec[res_add$best_in, res_add$best_out] <- NA
      in_edges[[res_add$best_in]] <- res_add$best_new_right
      degrees[res_add$best_in] <- degrees[res_add$best_in] + 1
      degrees[res_add$best_out] <- degrees[res_add$best_out] + 1
      score_cache <- res_add$score_cache
      if (verbose) message("Added ", res_add$best_out, "->", res_add$best_in, ", new BIC ", best_sum)
    }
    while (TRUE){
      res_remove <- .remove_edge(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, rem_dec, score_cache, control)
      if (res_remove$best_in == -1){
        break
      }
      nodewise_score[res_remove$best_in] <- res_remove$best_new_score
      best_sum <- best_sum - res_remove$best_decrease
      rem_dec <- res_remove$rem_dec
      cur_adj_mat[res_remove$best_out, res_remove$best_in] <- 0
      add_dec[,res_remove$best_in] <- NA; rem_dec[,res_remove$best_in] <- NA
      in_edges[[res_remove$best_in]] <- res_remove$best_new_right
      degrees[res_remove$best_in] <- degrees[res_remove$best_in] - 1
      degrees[res_remove$best_out] <- degrees[res_remove$best_out] - 1
      score_cache <- res_remove$score_cache
      if (verbose) message("Removed ", res_remove$best_out, "->", res_remove$best_in, ", new BIC ", best_sum)
    }
    while (TRUE){
      res_turn <- .turn_edge(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, add_dec, rem_dec, maxInDegree, no_check_DAG, score_cache, control)
      if (res_turn$best_in == -1){
        break
      }
      nodewise_score[res_turn$best_in] <- res_turn$best_new_score_in
      nodewise_score[res_turn$best_out] <- res_turn$best_new_score_out
      best_sum <- best_sum - res_turn$best_decrease
      add_dec <- res_turn$add_dec; rem_dec <- res_turn$rem_dec
      cur_adj_mat[res_turn$best_out, res_turn$best_in] <- 0
      cur_adj_mat[res_turn$best_in, res_turn$best_out] <- 1
      add_dec[,res_turn$best_in] <- NA; rem_dec[,res_turn$best_in] <- NA
      add_dec[,res_turn$best_out] <- NA; rem_dec[,res_turn$best_out] <- NA
      in_edges[[res_turn$best_in]] <- res_turn$best_new_right_ii
      in_edges[[res_turn$best_out]] <- res_turn$best_new_right_oi
      score_cache <- res_turn$score_cache
      if (verbose) message("Turned ", res_turn$best_out, "->", res_turn$best_in, ", new BIC ", best_sum)
    }
  }
  return (list("graph"=cur_adj_mat, "BIC"=best_sum, "score_cache"=score_cache))
}

########################## SIMY ##########################


.update_progress <- function(current, total) {
  if (current == total){
    message("Complete.")
  }
  if (current > 0 && current < total){
    if (as.integer(current/total*100) != as.integer((current-1)/total*100))
      message(as.integer(current/total*100),"%|")
  }
}

.DecToBin <- function(x) {
  if (x == 0) {return ("0")}
  else if (x < 0) {stop("x needs to be >= 0.")}
  res <- numeric(ceiling(log(x+1)/log(2)))
  for (i in length(res):1){
    res[i] <- x %% 2
    x <- x %/% 2
  }
  return (paste(res[match(1, res):length(res)], collapse=""))
}

.DecToParents <- function(x, node){
  ## Off by one, since lists start with index 1
  ## parents with number > node were -1, so add back +1
  ## node = Inf gives no translation
  ## Example: .DecToParents(14,3)=c(1,4,5) since 13=1101, the digits correspond to nodes c(5,4,2,1)
  x <- x-1
  if (x == 0) {return (integer(0))}
  else if (x < 0) {stop("x needs to be >= 1.")}
  #else if (x >= 2^m) {stop("x must be smaller than 2^m, (x,m)=(", x, ",", m, ") provided.")}
  parents <- integer(0)
  for (i in 1:ceiling(log(x+1)/log(2))){
    if (x %% 2) {
      if (i >= node) {parents <- c(parents, i+1)}
      else {parents <- c(parents, i)}
    }
    x <- x %/% 2
  }
  return (parents)
}

.BinToDec <- function(x) {sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))}
.ParentsToDec <- function(parents, node) {
  ## Off by one, since lists start with index 1
  ## for parents with number > node, subtract 1
  ## node = Inf gives no translation
  if (node %in% parents) {stop("node must not be in parents.")}
  if (length(parents))
    return (sum(2^(sapply(parents,function(i){if(i>node){i-2} else{i-1}})))+1)
  return (1) ## if parents == integer(0)
}

# Could have used bottom-up, but since we will eventually need to call .ParentsToDec
# , it's easier to just iterate over the indices of subsets and call .DecToParents
.getLocalScores <- function(V, Y, parametrization, LS, m, control=list(), pb=NULL, Silander_count=0){
  for (v in 1:m){
    Vout <- .rm_elt(1:m, v)
    for (i in 1:(2^(m-1))){
      parents <- Vout[.DecToParents(i, Inf)]
      res <- zi_fit(V=V, Y=Y, left=v, right=parents, parametrization=parametrization, value_only=FALSE, control=control)
      LS[v, i] <- -bic_from_fit(res)
      if (!is.null(pb)) {
        Silander_count <- Silander_count+1
        .update_progress(Silander_count, pb)
      }
    }
  }
  return (LS)
}

.getBestParents <- function(m, v, LS, pb, Silander_count){
  Vout <- .rm_elt(1:m, v)
  bps <- bss <- numeric(2^(length(Vout)))
  for (csi in 1:(2^(length(Vout)))){ # enumerating subsets of Vout
    # e.g. csi <- 6
    cs <- .DecToParents(csi, v)
    bps[csi] <- csi # e.g. bps[6] <- 19 (c(2,6))
    bss[csi] <- LS[v,csi] # e.g. bss[6] <- LS[3][19] (LS[3][c(2,6)])
    for (cs1_out in cs){ # enumerating LOO-subsets of cs, e.g. cs1_out <- 3
      if (!is.null(pb)){
        Silander_count <- Silander_count+1
        .update_progress(Silander_count, pb)
      }
      if (cs1_out > v)
        cs1 <- csi-2^(cs1_out-2) # index of cs\{cs1_out}
      else
        cs1 <- csi-2^(cs1_out-1)
      if (bss[cs1] > bss[csi]){
        bss[csi] <- bss[cs1]
        bps[csi] <- bps[cs1]
      }
    }
  }
  return (bps)
  # The values of bps are the best parent sets in dec representation
  # The indices of bps are the candidate parent sets in dec representation
}

.getBestSinks <- function(m, bps, LS, pb, Silander_count){
  scores <- numeric(2^m)
  sinks <- rep(-1, 2^m)
  for (wi in 1:(2^m)){ # index as a subset of 1:m
    W <- .DecToParents(wi, Inf)
    for (sink in W){
      upvars <- .rm_elt(W, sink)
      skore <- scores[wi-2^(sink-1)] # upvars
      skore <- skore + LS[sink, bps[sink, .ParentsToDec(upvars, sink)]]
      if (sinks[wi] == -1 || skore > scores[wi]){
        scores[wi] <- skore; sinks[wi] <- sink
      }
      if (!is.null(pb)){
        Silander_count <- Silander_count+1
        .update_progress(Silander_count, pb)
      }
    }
  }
  return (sinks)
  ## indices of sinks correspond to nodes sapply(1:2^m, function(i){.DecToParents(i, Inf)})
}

.Sinks2ord <- function(m, sinks){
  ord <- numeric(m)
  left <- 1:m
  for (i in m:1){
    ord[i] <- sinks[.ParentsToDec(left,Inf)]
    left <- .rm_elt(left, ord[i])
  }
  return (ord)
}

.Ord2net <- function(m, ord, bps){
  predecs <- integer(0)
  adj_mat <- matrix(0, m, m)
  for (i in ord){
    adj_mat[.DecToParents(bps[i, .ParentsToDec(predecs, i)], i), i] <- 1
    predecs <- c(predecs, i)
  }
  return (adj_mat)
}

#' Exhaustive search using BIC for DAGs for zero-inflated data based on Hurdle conditionals.
#'
#' Exhaustive search using BIC for DAGs for zero-inflated data based on Hurdle conditionals.
#'
#' @param V A matrix of 0/1s, equal to Y != 0.
#' @param Y A data matrix of the same size as \code{V}.
#' @param parametrization A string, either \code{"abk"} (canonical) or \code{"pms"} (moment).
#' @param verbose A logical, whether to print intermediate steps.
#' @param control A list passed to \code{zi_fit()}. Please consult \code{?zi_fit}.
#' @details
#' Performs exhaustive DAG search (using BIC) for DAGs for zero-inflated data based on Hurdle conditionals.
#' See Silander and Myllymaki (2006) or \code{?pcalg::simy} for details. However, unlike the Gaussian case in their implementation that returns an equivalence class of DAGs, the DAG estimated by this function is exact.
#' @return An adjacency matrix.
#' @examples
#' m <- 4; n <- 1000
#' adj_mat <- ZiDAG::make_dag(m, mode = "chain", shuffle=TRUE)
#' d <- ZiDAG::gen_zero_dat(seed=1, gen_para="pms", adj_mat=adj_mat, n=n, gen_uniform_degree=1)
#' est <- ZiDAG::ziSIMY(V=d$V, Y=d$Y, parametrization="pms", verbose=FALSE,
#'     control=list("max_uniform_degree"=1L, "tol"=1e-8, "print_best_degree"=FALSE))
#' adj_mat == est
#'
#' d <- ZiDAG::gen_zero_dat(seed=1, gen_para="abk", adj_mat=adj_mat, n=n, gen_uniform_degree=1)
#' est <- ZiDAG::ziSIMY(V=d$V, Y=d$Y, parametrization="abk", verbose=FALSE,
#'     control=list(use_C=TRUE, maxit=1000, runs=2, report=0))
#' adj_mat == est
#'
#' est <- ZiDAG::ziSIMY(V=d$V, Y=d$Y, parametrization="pms", verbose=FALSE,
#'     control=list("max_uniform_degree"=1L, "tol"=1e-8, "print_best_degree"=FALSE))
#' adj_mat == est
#'
#' est <- ZiDAG::ziSIMY(V=d$V, Y=d$Y, parametrization="pms", verbose=FALSE,
#'     control=list("max_uniform_degree"=2L, "tol"=1e-8, "print_best_degree"=FALSE))
#' adj_mat == est
#' @export
ziSIMY <- function(V, Y, parametrization, verbose=FALSE, control=list()){
  n <- nrow(V)
  m <- ncol(V)
  Silander_count <- 0
  if (verbose){
    message("Starting Silander.\nStep 1/5, calculating local scores...")
    pb <- m*2^(m-1)
  }
  else {pb <- NULL}
  LS <- .getLocalScores(V, Y, parametrization, matrix(NA, m, 2^(m-1)), m, control, pb, Silander_count)
  if (verbose){
    message("\nStep 2/5, calculating best parents...")
    Silander_count <- 0
    #pb <- txtProgressBar(min=0, max=m*(m-1)*2^(m-2), style=3)
    pb <- m*(m-1)*2^(m-2)
  } else {pb <- NULL}
  bpss <- t(sapply(1:m, function(i){.getBestParents(m, i, LS, pb, Silander_count)}))
  if (verbose){
    message("\nStep 3/5, calculating best sinks...")
    Silander_count <- 0
    ##pb <- txtProgressBar(min=0, max=m*2^(m-1), style=3)
    pb <- m*2^(m-1)
  } else {pb <- NULL}
  sinks <- .getBestSinks(m, bpss, LS, pb, Silander_count)
  if (verbose)
    message("\nStep 4/5, calculating best order...")
  ord <- .Sinks2ord(m, sinks)
  if (verbose)
    message("\nStep 5/5, calculating best graph...\nSilander done.")
  est_enum <- .Ord2net(m, ord, bpss)
  return (est_enum)
}


########################## Tests if two DAGs are Markov equivalent ##########################


#' Returns if the DAGs corresponding to two adjacency matrices are Markov equivalent.
#'
#' Returns if the DAGs corresponding to two adjacency matrices are Markov equivalent.
#'
#' @param amat1 The first adjacency matrix.
#' @param amat2 The second adjacency matrix.
#' @details
#' First tests if \code{amat1} and \code{amat2} are two 2-d matrices of the same dimensions, then test whether they have the same CPDAG using \code{pcalg::dag2cpdag()}.
#' @return A logical, whether the DAGs corresponding to two adjacency matrices are Markov equivalent.
#' @examples
#' # Complete graphs of the same dimension are always Markov equivalent
#' DAG_eq(make_dag(10, "complete", shuffle=TRUE), make_dag(10, "complete", shuffle=TRUE))
#' # A v structure not Markov equivalent to a chain graph
#' DAG_eq(matrix(c(0,1,0, 0,0,0, 0,1,0), nrow=3, byrow=TRUE), make_dag(3, "chain"))
#' # Reversing a chain graphs gives a Markov equivalent graph
#' DAG_eq(make_dag(10, "chain"), make_dag(10, "chain")[10:1, 10:1])
#' @export
DAG_eq <- function(amat1, amat2){
  if (length(dim(amat1)) != 2 || length(dim(amat2)) != 2 || length(unique(c(dim(amat1),dim(amat2)))) != 1)
    stop("amat1 and amat2 must be square matrices of the same dimension.")
  # Checks if two dags represented by adj matrices are equivalent
  return (all(pcalg::dag2cpdag(amat1) == pcalg::dag2cpdag(amat2)))
}
