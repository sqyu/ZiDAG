update_progress <- function(current, total) {
  if (current == total){
    cat("Complete.\n")
  }
  if (current > 0 && current < total){
    if (as.integer(current/total*100) != as.integer((current-1)/total*100))
      cat(as.integer(current/total*100),"%|", sep="")
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

# SLOW! Top-down, using pseudocode from paper, lots of repetition
getLocalScores0 <- function(V, Y, parametrization, LS, vars, control=list(), pb=NULL){
  for (v in vars){
    LS <- getLocalScores0(V, Y, parametrization, LS, setdiff(vars, v), control, pb)
    #LS[[v]][ParentsToDec(setdiff(vars, v), v)] <- fit(V=V, Y=Y, left=v, right=setdiff(vars, v), parametrization=parametrization, value_only=T, control=control)
    res <- fit(V=V, Y=Y, left=v, right=setdiff(vars, v), parametrization=parametrization, value_only=F, control=control)
    LS[v, ParentsToDec(setdiff(vars, v), v)] <- -bic_from_fit(res)
    if (!is.null(pb)) {
      count <<- count+1;
      #setTxtProgressBar(pb, count)
      update_progress(count, pb)
    }
  }
  return (LS)
}

# Could have used bottom-up, but since we will eventually need to call ParentsToDec
# , it's easier to just iterate over the indices of subsets and call DecToParents
getLocalScores <- function(V, Y, parametrization, LS, m, control=list(), pb=NULL){
  for (v in 1:m){
    Vout <- setdiff(1:m, v)
    for (i in 1:(2^(m-1))){
      parents <- Vout[DecToParents(i, Inf)]
      res <- fit(V=V, Y=Y, left=v, right=parents, parametrization=parametrization, value_only=F, control=control)
      LS[v, i] <- -bic_from_fit(res)
      if (!is.null(pb)) {
        count <<- count+1
        #setTxtProgressBar(pb, count)
        update_progress(count, pb)
      }
    }
  }
  return (LS)
}

## assuming general V, not needed
getBestParents0 <- function(V, v, LS){
  ## e.g. V=c(2,3,5,6), v=3
  if (!v %in% V) {stop("v must be in V.")}
  if (length(V)==1){return (integer(0))}
  Vout <- setdiff(V, v) ## e.g. c(2,5,6)
  bps <- bss <- numeric(2^(length(Vout)))
  for (csi in 1:(2^(length(Vout)))){ # enumerating subsets of Vout
    # e.g. csi <- 6
    csi_par_local <- DecToParents(csi, Inf) # conv csi to positions within Vout, e.g. csi_par_local <- c(1,3)
    cs <- Vout[csi_par_local] # true vertex indices, enumeration over subsets of V, e.g. cs <- c(2,6)
    cs_dec <- ParentsToDec(cs, v) # true vertex indices (dec representation) of cs, e.g. c(2,6) for 3->19
    bps[csi] <- cs_dec # e.g. bps[6] <- 19 (c(2,6))
    bss[csi] <- LS[v,cs_dec] # e.g. bss[6] <- LS[3][19] (LS[3][c(2,6)])
    for (cs1_out in csi_par_local){ # enumerating LOO-subsets of cs, e.g. cs1_out <- 3
      cs1_one_out_i <- csi-2^(cs1_out-1) # index of the new subset of cs within bps and bss, e.g. cs1_one_out_i = 6-2^(3-1)=2, which is index for c(1), 1 as position within Vout, or c(2), 2 as true vertex index
      if (bss[cs1_one_out_i] > bss[csi]){
        bss[csi] <- bss[cs1_one_out_i]
        bps[csi] <- bps[cs1_one_out_i]
      }
    }
  }
  return (bps)
  # The values of bps are in dec rep of true vertex indices
  # The true vertex indices corresponding to the indices of bps can be recovered by
  # sapply(1:2^(length(V)-1), function(csi){setdiff(V, v)[DecToParents(csi, Inf)]})
}

getBestParents <- function(m, v, LS, pb){
  Vout <- setdiff(1:m, v)
  bps <- bss <- numeric(2^(length(Vout)))
  for (csi in 1:(2^(length(Vout)))){ # enumerating subsets of Vout
    # e.g. csi <- 6
    cs <- DecToParents(csi, v)
    bps[csi] <- csi # e.g. bps[6] <- 19 (c(2,6))
    bss[csi] <- LS[v,csi] # e.g. bss[6] <- LS[3][19] (LS[3][c(2,6)])
    for (cs1_out in cs){ # enumerating LOO-subsets of cs, e.g. cs1_out <- 3
      if (!is.null(pb)){
        count <<- count+1
        #setTxtProgressBar(pb, count)
        update_progress(count, pb)
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

# Sanity check
#V <- c(1,2,4);  v <- 1
#tmp <- LS[v, ]
#names(tmp) <- sapply(1:2^(m-1), function(i){paste(DecToParents(i, v), collapse=",")})
# LS
#tmp
# Picked best parents
#sapply(getBestParents0(V, v,LS), function(i){paste(DecToParents(i,v), collapse=",")})
# Within
#sapply(1:2^(length(V)-1), function(csi){paste(setdiff(V, v)[DecToParents(csi, Inf)], collapse=",")})

# Sanity check if getBestParents is getBestParents0 with V=1:m
#sapply(1:m, function(v){all(getBestParents0(1:m, v,LS) == getBestParents(m, v,LS))})



getBestSinks <- function(m, bps, LS, pb){
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
        count <<- count+1
        #setTxtProgressBar(pb, count)
        update_progress(count, pb)
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

count <- 0

Silander <- function(V, Y, parametrization, verbose=F, control=list()){
  n <- nrow(V)
  m <- ncol(V)
  if (verbose){
    cat("Starting Silander.\nStep 1/5, calculating local scores...")
    count <<- 0
    #pb <- txtProgressBar(min=0, max=m*2^(m-1), style=3)
    pb <- m*2^(m-1)
  }
  else {pb <- NULL}
  LS <- getLocalScores(V, Y, parametrization, matrix(NA, m, 2^(m-1)), m, control, pb)
  if (verbose){
    cat("\nStep 2/5, calculating best parents...")
    count <<- 0
    #pb <- txtProgressBar(min=0, max=m*(m-1)*2^(m-2), style=3)
    pb <- m*(m-1)*2^(m-2)
  } else {pb <- NULL}
  bpss <- t(sapply(1:m, function(i){getBestParents(m, i, LS, pb)}))
  if (verbose){
    cat("\nStep 3/5, calculating best sinks...")
    count <<- 0
    ##pb <- txtProgressBar(min=0, max=m*2^(m-1), style=3)
    pb <- m*2^(m-1)
  } else {pb <- NULL}
  sinks <- getBestSinks(m, bpss,LS,pb)
  if (verbose)
    cat("\nStep 4/5, calculating best order...")
  ord <- Sinks2ord(m, sinks)
  if (verbose)
    cat("\nStep 5/5, calculating best graph...\nSilander done.")
  est_enum <- Ord2net(m, ord, bpss)
  return (est_enum)
}
