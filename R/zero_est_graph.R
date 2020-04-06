#' Calculates the BIC score from a fit result.
#'
#' Calculates the BIC score from a fit result.
#'
#' @param fit_res A list returned from \code{zi_fit()} or \code{.res_from_model_fits()} with \code{value_only == FALSE}.
#' @details Returns \code{bic_score(fit_res$value, fit_res$n, fit_res$effective_df)}.
#' @return A number, the BIC score.
#' @examples
#' Y <- matrix(rnorm(500) * rbinom(500, 1, 0.7), nrow=50)
#' res <- zi_fit(Y != 0, Y, 1, 2:5, value_only=FALSE)
#' bic_from_fit(res)
#' bic_score(res$value, res$n, res$effective_df)
#' @export
bic_from_fit <- function(fit_res){
  return (bic_score(fit_res$value, fit_res$n, fit_res$effective_df))
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
#' bic_score(res$value, res$n, res$effective_df)
#' bic_from_fit(res)
#' @export
bic_score <- function(avgnll, n, numpar){
  return (2*avgnll+log(n)/n*numpar)
}


.add_edge <- function(V, Y, parametrization, in_edges, degrees, cur_adj_mat, nodewise_score, add_dec, maxDegree, maxInDegree, fixedGaps, control=list()){
  if (!requireNamespace("ggm", quietly = TRUE))
    stop("Please install package \"ggm\".")
  best_in <- best_out <- -1
  best_decrease <- 0
  best_new_score <- Inf
  m <- ncol(V); n <- nrow(V)
  for (ii in 1:m){
    if (degrees[ii] < maxDegree && length(in_edges[[ii]]) < maxInDegree){
      for (oi in setdiff(1:m, c(ii, in_edges[[ii]]))){
        if (degrees[oi] < maxDegree && (is.null(fixedGaps) || !fixedGaps[oi,ii])) {
          cur_adj_mat[oi, ii] <- 1
          if (!(ii %in% in_edges[[oi]]) && ggm::isAcyclic(cur_adj_mat)){
            if (!is.na(add_dec[oi, ii])){
              dec <- add_dec[oi, ii]
              new_score <- nodewise_score[ii] - dec
            }
            else{
              new_score <- bic_from_fit(zi_fit(V, Y, left=ii, right=c(oi,in_edges[[ii]]), parametrization=parametrization, value_only=FALSE, control=control))
              dec <- nodewise_score[ii] - new_score
              add_dec[oi,ii] <- dec
            }
            if (dec > best_decrease){
              best_in <- ii; best_out <- oi;  best_new_score <- new_score
              best_decrease <- dec
            }
          }
          cur_adj_mat[oi, ii] <- 0
        }
      }
    }
  }
  return (list("best_in"=best_in, "best_out"=best_out, "best_decrease"=best_decrease,
               "best_new_score"=best_new_score, "add_dec"=add_dec))
}

.remove_edge <- function(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, rem_dec, control=list()){
  best_in <- best_out <- -1
  best_decrease <- 0
  best_new_score <- Inf
  n <- nrow(V); m <- ncol(V)
  for (ii in 1:m){
    for (oi in in_edges[[ii]]){
      if (!is.na(rem_dec[oi,ii])){
        dec <- rem_dec[oi,ii]
        new_score <- nodewise_score[ii] - dec
      } else {
        new_score <- bic_from_fit(zi_fit(V, Y, ii, setdiff(in_edges[[ii]],oi), parametrization=parametrization, value_only=FALSE, control=control))
        dec <- nodewise_score[ii] - new_score
        rem_dec[oi,ii] <- dec
      }
      if (dec > best_decrease){
        best_in <- ii; best_out <- oi;  best_new_score <- new_score
        best_decrease <- dec
      }
    }
  }
  return (list("best_in"=best_in, "best_out"=best_out, "best_decrease"=best_decrease,
               "best_new_score"=best_new_score, "rem_dec"=rem_dec))
}

.turn_edge <- function(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, add_dec, rem_dec, maxInDegree, control=list()){
  if (!requireNamespace("ggm", quietly = TRUE))
    stop("Please install package \"ggm\".")
  best_in <- best_out <- -1
  best_decrease <- 0
  best_new_score_in <- best_new_score_out <- Inf
  m <- ncol(V); n <- nrow(V)
  for (ii in 1:m){
    for (oi in in_edges[[ii]]){
      if (length(in_edges[[oi]]) >= maxInDegree) # Skip if indegree of out node is max
        next
      cur_adj_mat[ii, oi] <- 1
      cur_adj_mat[oi, ii] <- 0
      if (ggm::isAcyclic(cur_adj_mat)){
        if (!is.na(rem_dec[oi,ii])){
          dec_in <- rem_dec[oi,ii]
          new_score_in <- nodewise_score[ii] - dec_in
        } else {
          new_score_in <- bic_from_fit(zi_fit(V, Y, ii, setdiff(in_edges[[ii]],oi), parametrization=parametrization, value_only=FALSE, control=control))
          dec_in <- nodewise_score[ii] - new_score_in
          rem_dec[oi,ii] <- dec_in
        }
        if (!is.na(add_dec[ii,oi])){
          dec_out <- add_dec[ii,oi]
          new_score_out <- nodewise_score[oi] - dec_out
        } else {
          new_score_out <- bic_from_fit(zi_fit(V, Y, oi, c(ii,in_edges[[oi]]), parametrization=parametrization, value_only=FALSE, control=control))
          dec_out <- nodewise_score[oi] - new_score_out
          add_dec[ii,oi] <- dec_out
        }
        if (dec_in+dec_out > best_decrease){
          best_in <- ii; best_out <- oi
          best_new_score_in <- new_score_in
          best_new_score_out <- new_score_out
          best_decrease <- dec_in+dec_out
        }
      }
      cur_adj_mat[oi, ii] <- 1
      cur_adj_mat[ii, oi] <- 0
    }
  }
  return (list("best_in"=best_in, "best_out"=best_out, "best_decrease"=best_decrease,
               "best_new_score_in"=best_new_score_in, "best_new_score_out"=best_new_score_out,
               "add_dec"=add_dec, "rem_dec"=rem_dec))
}










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
