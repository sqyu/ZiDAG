update_progress <- function(current, total) {
  if (current == total){
    message("Complete.\n")
  }
  if (current > 0 && current < total){
    if (as.integer(current/total*100) != as.integer((current-1)/total*100))
      message(as.integer(current/total*100),"%|", sep="")
  }
}

DecToBin <- function(x) {
  if (x == 0) {return ("0")}
  else if (x < 0) {stop("x needs to be >= 0.")}
  res <- numeric(ceiling(log(x+1)/log(2)))
  for (i in length(res):1){
    res[i] <- x %% 2
    x <- x %/% 2
  }
  return (paste(res[match(1, res):length(res)], collapse=""))
}

DecToParents <- function(x, node){
  ## Off by one, since lists start with index 1
  ## parents with number > node were -1, so add back +1
  ## node = Inf gives no translation
  ## Example: DecToParents(14,3)=c(1,4,5) since 13=1101, the digits correspond to nodes c(5,4,2,1)
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

BinToDec <- function(x) {sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))}
ParentsToDec <- function(parents, node) {
  ## Off by one, since lists start with index 1
  ## for parents with number > node, subtract 1
  ## node = Inf gives no translation
  if (node %in% parents) {stop("node must not be in parents.")}
  if (length(parents))
    return (sum(2^(sapply(parents,function(i){if(i>node){i-2} else{i-1}})))+1)
  return (1) ## if parents == integer(0)
}

# Could have used bottom-up, but since we will eventually need to call ParentsToDec
# , it's easier to just iterate over the indices of subsets and call DecToParents
getLocalScores <- function(V, Y, parametrization, LS, m, control=list(), pb=NULL, Silander_count=0){
  for (v in 1:m){
    Vout <- setdiff(1:m, v)
    for (i in 1:(2^(m-1))){
      parents <- Vout[DecToParents(i, Inf)]
      res <- zi_fit(V=V, Y=Y, left=v, right=parents, parametrization=parametrization, value_only=FALSE, control=control)
      LS[v, i] <- -bic_from_fit(res)
      if (!is.null(pb)) {
        Silander_count <- Silander_count+1
        update_progress(Silander_count, pb)
      }
    }
  }
  return (LS)
}

getBestParents <- function(m, v, LS, pb, Silander_count){
  Vout <- setdiff(1:m, v)
  bps <- bss <- numeric(2^(length(Vout)))
  for (csi in 1:(2^(length(Vout)))){ # enumerating subsets of Vout
    # e.g. csi <- 6
    cs <- DecToParents(csi, v)
    bps[csi] <- csi # e.g. bps[6] <- 19 (c(2,6))
    bss[csi] <- LS[v,csi] # e.g. bss[6] <- LS[3][19] (LS[3][c(2,6)])
    for (cs1_out in cs){ # enumerating LOO-subsets of cs, e.g. cs1_out <- 3
      if (!is.null(pb)){
        Silander_count <- Silander_count+1
        update_progress(Silander_count, pb)
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

getBestSinks <- function(m, bps, LS, pb, Silander_count){
  scores <- numeric(2^m)
  sinks <- rep(-1, 2^m)
  for (wi in 1:(2^m)){ # index as a subset of 1:m
    W <- DecToParents(wi, Inf)
    for (sink in W){
      upvars <- setdiff(W, sink)
      skore <- scores[wi-2^(sink-1)] # upvars
      skore <- skore + LS[sink, bps[sink, ParentsToDec(upvars, sink)]]
      if (sinks[wi] == -1 || skore > scores[wi]){
        scores[wi] <- skore; sinks[wi] <- sink
      }
      if (!is.null(pb)){
        Silander_count <- Silander_count+1
        update_progress(Silander_count, pb)
      }
    }
  }
  return (sinks)
  ## indices of sinks correspond to nodes sapply(1:2^m, function(i){DecToParents(i, Inf)})
}

Sinks2ord <- function(m, sinks){
  ord <- numeric(m)
  left <- 1:m
  for (i in m:1){
    ord[i] <- sinks[ParentsToDec(left,Inf)]
    left <- setdiff(left, ord[i])
  }
  return (ord)
}

Ord2net <- function(m, ord, bps){
  predecs <- integer(0)
  adj_mat <- matrix(0, m, m)
  for (i in ord){
    adj_mat[DecToParents(bps[i, ParentsToDec(predecs, i)], i), i] <- 1
    predecs <- c(predecs, i)
  }
  return (adj_mat)
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
ziSIMY <- function(V, Y, parametrization, verbose=FALSE, control=list()){
  n <- nrow(V)
  m <- ncol(V)
  Silander_count <- 0
  if (verbose){
    message("Starting Silander.\nStep 1/5, calculating local scores...")
    pb <- m*2^(m-1)
  }
  else {pb <- NULL}
  LS <- getLocalScores(V, Y, parametrization, matrix(NA, m, 2^(m-1)), m, control, pb, Silander_count)
  if (verbose){
    message("\nStep 2/5, calculating best parents...")
    Silander_count <- 0
    #pb <- txtProgressBar(min=0, max=m*(m-1)*2^(m-2), style=3)
    pb <- m*(m-1)*2^(m-2)
  } else {pb <- NULL}
  bpss <- t(sapply(1:m, function(i){getBestParents(m, i, LS, pb, Silander_count)}))
  if (verbose){
    message("\nStep 3/5, calculating best sinks...")
    Silander_count <- 0
    ##pb <- txtProgressBar(min=0, max=m*2^(m-1), style=3)
    pb <- m*2^(m-1)
  } else {pb <- NULL}
  sinks <- getBestSinks(m, bpss, LS, pb, Silander_count)
  if (verbose)
    message("\nStep 4/5, calculating best order...")
  ord <- Sinks2ord(m, sinks)
  if (verbose)
    message("\nStep 5/5, calculating best graph...\nSilander done.")
  est_enum <- Ord2net(m, ord, bpss)
  return (est_enum)
}
