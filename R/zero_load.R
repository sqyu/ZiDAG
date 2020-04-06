# library(knitr); library(rmarkdown); library(devtools); library(roxygen2)
# document(); build(); install(); check()
# tools::package_native_routine_registration_skeleton(".")  # Copy as src/ZiDAG_init.c


#'
#'
#'
#'
#' @param
#' @details
#'
#' @return
#' @examples
#'
#' @export
zi_fit <- function(V, Y, left, right, parametrization="pms", value_only=TRUE, control=list()){
  control[["V"]] <- V
  control[["Y"]] <- Y
  control[["left"]] <- left
  if (length(right) == 0) # Something weird about R
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

#'
#'
#'
#'
#' @param
#' @details
#'
#' @return
#' @examples
#'
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








#'
#'
#'
#'
#' @param
#' @details
#'
#' @return
#' @examples
#'
#' @export
ziGDS <- function(V, Y, parametrization, fixedGaps=NULL, maxDegree=Inf, maxInDegree=Inf, init=NULL, verbose=FALSE, control=list()){
  if (!requireNamespace("ggm", quietly = TRUE))
    stop("Please install package \"ggm\".")
  if (!all(dim(V)==dim(Y))) stop("V and Y must have the same dimension.")
  n <- nrow(V); m <- ncol(V)
  maxDegree <- min(m-1, maxDegree)
  maxInDegree <- min(m-1, maxInDegree)
  if (maxInDegree > maxDegree)
    stop("maxInDegree must be smaller than maxDegree.")
  if (!is.null(fixedGaps)) {
    if (length(dim(fixedGaps)) != 2 || !all(dim(fixedGaps) == m))
      stop("If provided, fixedGaps must be a ", m, "x", m, " matrix.")
    if (!all(fixedGaps*1 == 0 || fixedGaps*1 == 1))
      stop("If provided, fixedGaps must be a boolean or 0/1 matrix.")
    if (sum(abs(fixedGaps - t(fixedGaps))) != 0) {
      stop("fixedGaps must be a symmetric matrix.")
    }
  }
  if (is.null(init)){
    cur_adj_mat <- matrix(0,m,m)
    in_edges <- lapply(1:m, function(i){integer(0)})
  } else if (is.list(init)){
    if (length(init)!=m)
      stop("Initial (in-)edge list must have length ",m,".")
    cur_adj_mat <- matrix(0,m,m)
    for (in_v in 1:m){
      if (any(init[[in_v]] < 1 | init[[in_v]] > m))
        stop("Vertices in edge list must be between 1 and ",m,".")
      if (in_v %in% init[[in_v]]){
        warning("Self edge for", in_v, "ignored.")
        init[[in_v]] <- sort(setdiff(init[[in_v]], c(in_v)))
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
  if (!ggm::isAcyclic(cur_adj_mat))
    stop("Initial graph must be acyclic.")
  degrees <- rowSums(cur_adj_mat) + colSums(cur_adj_mat)  # only need the degree, not the list
  nodewise_score <- sapply(1:m, function(i){bic_from_fit(zi_fit(V, Y, left=i, right=in_edges[[i]], parametrization=parametrization, value_only=FALSE, control=control))})
  old_best_sum <- Inf
  best_sum <- sum(nodewise_score)
  iter_count <- 0
  add_dec <- rem_dec <- matrix(NA,m,m)
  while (best_sum < old_best_sum-1e-9){
    old_best_sum <- best_sum
    iter_count <- iter_count + 1
    if (verbose){message(paste("Iteration ", iter_count, ": BIC score ", old_best_sum, sep=""), "\n")}
    while (TRUE){
      res_add <- .add_edge(V, Y, parametrization, in_edges, degrees, cur_adj_mat, nodewise_score, add_dec, maxDegree, maxInDegree, fixedGaps, control)
      if (res_add$best_in == -1){
        break
      }
      nodewise_score[res_add$best_in] <- res_add$best_new_score
      best_sum <- best_sum - res_add$best_decrease
      add_dec <- res_add$add_dec
      cur_adj_mat[res_add$best_out, res_add$best_in] <- 1
      add_dec[,res_add$best_in] <- NA; rem_dec[,res_add$best_in] <- NA
      add_dec[res_add$best_in, res_add$best_out] <- NA
      in_edges[[res_add$best_in]] <- sort(c(in_edges[[res_add$best_in]], res_add$best_out))
      degrees[res_add$best_in] <- degrees[res_add$best_in] + 1
      degrees[res_add$best_out] <- degrees[res_add$best_out] + 1
      if (verbose) message(paste("Added ", res_add$best_out, "->", res_add$best_in, ", new BIC ", best_sum, sep=""), "\n")
    }
    while (TRUE){
      res_remove <- .remove_edge(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, rem_dec, control)
      if (res_remove$best_in == -1){
        break
      }
      nodewise_score[res_remove$best_in] <- res_remove$best_new_score
      best_sum <- best_sum - res_remove$best_decrease
      rem_dec <- res_remove$rem_dec
      cur_adj_mat[res_remove$best_out, res_remove$best_in] <- 0
      add_dec[,res_remove$best_in] <- NA; rem_dec[,res_remove$best_in] <- NA
      in_edges[[res_remove$best_in]] <- setdiff(in_edges[[res_remove$best_in]], res_remove$best_out)
      degrees[res_add$best_in] <- degrees[res_add$best_in] - 1
      degrees[res_add$best_out] <- degrees[res_add$best_out] - 1
      if (verbose) message(paste("Removed ", res_remove$best_out, "->", res_remove$best_in, ", new BIC ", best_sum, sep=""), "\n")
    }
    while (TRUE){
      res_turn <- .turn_edge(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, add_dec, rem_dec, maxInDegree, control)
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
      in_edges[[res_turn$best_in]] <- setdiff(in_edges[[res_turn$best_in]], res_turn$best_out)
      in_edges[[res_turn$best_out]] <- sort(c(in_edges[[res_turn$best_out]], res_turn$best_in))
      if (verbose) message(paste("Turned ", res_turn$best_out, "->", res_turn$best_in, ", new BIC ", best_sum, sep=""), "\n")
    }
  }
  return (cur_adj_mat)
}
