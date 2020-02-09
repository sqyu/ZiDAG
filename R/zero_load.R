## 310418: 1. Added ..._auto
##         2. Added estimate_build_bhy and estimate_oneside_build_bhy
##         3. estimate and estimate_oneside now return the adj matrices.
## On 20190707: Fixed maxDegree from limiting the indegree to the
##  degree in the adjacency matrix in my_gds. Added implementation for fixedGaps.

#source("http://bioconductor.org/biocLite.R")
#biocLite(c("graph", "RBGL", "Rgraphviz"))
#install.packages("pcalg", dependencies=TRUE)

require(gtools) # permutations
require(sets) # enumeration of DAGs
require(ggm)  # isAcyclic
require(pcalg)
require(sfsmisc) ## For is.whole in validate.vertex inherited from Score

source("Silander.R")
dyn.load("optim.so")
## Change from 300507: (Suppose left=1) g1*2*Vo now becomes g1*Vo, i.e. in 300507 g1 was q/2 but now g1 is q2

logsumexp_vec <- function(xs){ ## Calculates mean(log(1+exp(xs)))
  return (mean(pmax(xs,0) + log(1+exp(-abs(xs)))))
}

logprob1d_ghk <- function(V,Y,g,h,k){ # 1-d hurdle
  # V, Y: vectors
  mean(g*V+h*Y-Y^2*k/2)-logsumexp_vec(g+h^2/2/k-log(k/2/pi)/2)
}

logprob_vec_ghk <- function(V,Y,G,H,k){
  # V, Y, G, H: vectors
  # each entry in (V,Y,G,H) correspond to a different sample with different parameters
  c(crossprod(V,G)+crossprod(Y,H))/length(V)-mean(Y^2)*k/2-logsumexp_vec(G+H^2/2/k-log(k/2/pi)/2)
}

sum_g <- function(gg,vo,yo){
  #gg : c(g11, g1, h1o) ## 1, m-1, m-1
  c(crossprod(gg, c(1,vo,yo)))
}
sum_h <- function(hh,vo,yo){
  #hh : c(h11, ho1, k1) ## 1, m-1, m-1
  c(crossprod(hh, c(1,vo,-yo)))
}
sum_G <- function(gg, VO, YO){
  c(tcrossprod(gg, cbind(1,VO,YO)))
}

sum_H <- function(hh, VO, YO, minus_Y=TRUE){
  if (minus_Y)
    c(tcrossprod(hh, cbind(1, VO, -YO)))
  else
    c(tcrossprod(hh, cbind(1, VO, YO)))
}

hess_g_h_k_block <- function(gg,hh,k11,Vo,Yo){
  ## Please check that the entries with g12 are not off by 2 or 4
  # gg:c(g11, g1, h1o), hh:c(h11, ho1, k1)
  r <- ncol(Vo)
  g <- sum_G(gg, Vo, Yo) # n-vector
  h <- sum_H(hh, Vo, Yo) # n-vector
  exp_part <- exp(g+h^2/2/k11) # n-vector
  common <- exp_part / (sqrt(k11)+exp_part*sqrt(2*pi))^2 # n-vector
  Vo2 <- Vo^2; VoYo <- Vo*Yo; Yo2 <- Yo^2
  dg12_g12 <- -sqrt(k11*2*pi)*colMeans(sweep(Vo2, 1, common, FUN="*")) # (m-1)-vector
  dg12_h12 <- -sqrt(k11*2*pi)*colMeans(sweep(VoYo, 1, common, FUN="*")) # (m-1)-vector
  dg12_h21 <- -sqrt(2*pi/k11)*colMeans(sweep(Vo2*h, 1, common, FUN="*")) # (m-1)-vector
  dg12_k12 <- sqrt(2*pi/k11)*colMeans(sweep(VoYo*h, 1, common, FUN="*")) # (m-1)-vector
  dh12_h12 <- -sqrt(2*pi*k11)*colMeans(sweep(Yo2, 1, common, FUN="*")) # (m-1)-vector
  dh12_h21 <- -dg12_k12 # (m-1)-vector
  dh12_k12 <- sqrt(2*pi/k11)*colMeans(sweep(Yo2*h, 1, common, FUN="*")) # (m-1)-vector
  common2 <- common*(k11)^(-3/2)*(2*exp_part*sqrt(k11)*pi+sqrt(2*pi)*(h^2+k11))
  remove(common)
  dh21_h21 <- -colMeans(sweep(Vo2, 1, common2, FUN="*")) # (m-1)-vector
  dh21_k12 <- colMeans(sweep(VoYo, 1, common2, FUN="*")) # (m-1)-vector
  dk12_k12 <- -colMeans(sweep(Yo2, 1, common2, FUN="*")) # (m-1)-vector
  remove(Vo2, VoYo, Yo2, common2)
  HESS <- matrix(0,4,4*r)
  ind <- 16*(0:(r-1))
  HESS[1+ind] <- dg12_g12
  HESS[2+ind] <- HESS[5+ind] <- dg12_h12
  HESS[3+ind] <- HESS[9+ind] <- dg12_h21
  HESS[4+ind] <- HESS[13+ind] <- dg12_k12
  HESS[6+ind] <- dh12_h12
  HESS[7+ind] <- HESS[10+ind] <- dh12_h21
  HESS[8+ind] <- HESS[14+ind] <- dh12_k12
  HESS[11+ind] <- dh21_h21
  HESS[12+ind] <- HESS[15+ind] <- dh21_k12
  HESS[16+ind] <- dk12_k12
  return (HESS)
}

hess_g_h_k_null <- function(g11,h11,k11,Vo,Yo){
  ## Please check that the entries with g12 are not off by 2 or 4
  r <- ncol(Vo)
  exp_part <- exp(g11+h11^2/2/k11) # n-vector
  common <- exp_part / (sqrt(k11)+exp_part*sqrt(2*pi))^2 # n-vector
  dg12_g12 <- -sqrt(k11*2*pi)*common*colMeans(Vo^2) # (m-1)-vector
  dg12_h12 <- -sqrt(k11*2*pi)*common*colMeans(Vo*Yo) # (m-1)-vector
  dg12_h21 <- dg12_g12/k11*h11# (m-1)-vector
  dg12_k12 <- -h11/k11*dg12_h12# (m-1)-vector
  dh12_h12 <- -sqrt(2*pi*k11)*common*colMeans(Yo^2) # (m-1)-vector
  dh12_h21 <- -dg12_k12 # (m-1)-vector
  dh12_k12 <- -dh12_h12/k11*h11
  common2 <- (h11^2+k11+exp_part*sqrt(k11*2*pi))/(h11*k11)
  dh21_h21 <- common2*dg12_h21 # (m-1)-vector
  dh21_k12 <- common2*dg12_k12 # (m-1)-vector
  dk12_k12 <- -common2*dh12_k12 # (m-1)-vector
  HESS <- matrix(0,4,4*r)
  ind <- 16*(0:(r-1))
  HESS[1+ind] <- dg12_g12
  HESS[2+ind] <- HESS[5+ind] <- dg12_h12
  HESS[3+ind] <- HESS[9+ind] <- dg12_h21
  HESS[4+ind] <- HESS[13+ind] <- dg12_k12
  HESS[6+ind] <- dh12_h12
  HESS[7+ind] <- HESS[10+ind] <- dh12_h21
  HESS[8+ind] <- HESS[14+ind] <- dh12_k12
  HESS[11+ind] <- dh21_h21
  HESS[12+ind] <- HESS[15+ind] <- dh21_k12
  HESS[16+ind] <- dk12_k12
  return (-HESS)
}

grad_g_h_k <- function(v,y,g,h,k){
  exp_part <- exp(g+h^2/2/k)
  return (cbind((v-1)+1/(1+exp_part*sqrt(2*pi/k)),
                y-h/k+h/(k+exp_part*sqrt(2*pi*k)),
                h^2/2/k^2+1/2/k-y^2/2-(h^2+k)/(2*k^(3/2)*(sqrt(k)+exp_part*sqrt(2*pi)))))
}

grad_full_vec <- function(V,Y,G,H,k11,Vo,Yo){
  grad_tmp <- grad_g_h_k(V,Y,G,H,k11)
  n <- length(V)
  return (c(mean(grad_tmp[,1]), crossprod(grad_tmp[,1],Vo)/n, crossprod(grad_tmp[,1],Yo)/n,
            mean(grad_tmp[,2]), crossprod(grad_tmp[,2],Vo)/n, -crossprod(grad_tmp[,2],Yo)/n,
            mean(grad_tmp[,3])))
}

convert_para <- function(paras, lpms_to_ghk){
  if (!lpms_to_ghk %in% c(FALSE, TRUE)) {stop("lpms_to_ghk must be TRUE or FALSE.")}
  if (length(paras) != 3) {stop("paras must have length 3.")}
  if (lpms_to_ghk) {
    lp <- paras[1]; mu <- paras[2]; sigmasq <- paras[3]
    k <- 1/sigmasq; h <- mu*k; g <- lp+log(k/2/pi)/2-h^2/2/k
    return (c(g,h,k))
  } else {
    g <- paras[1]; h <- paras[2]; k <- paras[3]
    return (c(g+h^2/2/k-log(k/2/pi)/2, h/k, 1/k))
  }
}

generate_logitp_mu_sigmasq <- function(n, lp, mu, sigmasq){
  rbinom(n, 1, inv.logit(lp)) * rnorm(n, mu, sqrt(sigmasq))
}
generate_g_h_k <- function(n, g, h, k){
  generate_logitp_mu_sigmasq(n, g+h^2/2/k-log(k/2/pi)/2, h/k, 1/k)
}

rbb <- Vectorize(rbinom, "prob")
rnn <- Vectorize(rnorm, "mean")
generate_one_each_logitp_mu_sigmasq <- function(lps, mus, sigmasq){
  if (length(setdiff(c(length(lps),length(mus)), c(1))) > 1)
    stop("lps and mus must have the same length or be a scalar.")
  rbb(1, 1, inv.logit(lps)) * rnn(1, mus, sqrt(sigmasq))
}
generate_one_each_g_h_k <- function(gs, hs, k){
  if (length(setdiff(c(length(gs),length(hs)), c(1))) > 1)
    stop("gs and hs must have the same length or be a scalar.")
  generate_one_each_logitp_mu_sigmasq(gs+hs^2/2/k-log(k/2/pi)/2, hs/k, 1/k)
}

mle1d_ghk <- function(V,Y,left,value_only=T){
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
  khat <- 1/sigmasqhat; hhat <- muhat*khat; ghat <- log(phat/(1-phat))+log(khat/2/pi)/2-hhat^2/2/khat
  nll <- -logprob1d_ghk(V,Y,ghat,hhat,khat)
  if (value_only)
    return (nll)
  else{
    pars <- c(ghat, hhat, khat)
    names(pars) <- paste(c("g", "h", "k"), left, sep="")
    return (list("value"=nll, "par"=pars, "n"=n, "effective_df"=3))
  }
}

fit_R <- function(V, Y, left, right, maxit=200, tol=1e-8, runs=1, value_only=T, report=10){
  r <- length(right)
  V1 <- V[,left]; Y1 <- Y[,left]
  Vo <- V[,right,drop=F]; Yo <- Y[,right,drop=F]
  #G <- H <- numeric(2*r+1)
  nll_optim <- function(x){-logprob_vec_ghk(V1,Y1,sum_G(x[1:(2*r+1)],Vo,Yo),sum_H(x[(2*r+2):(4*r+2)],Vo,Yo),x[4*r+3])}
  grad_optim <- function(x){-grad_full_vec(V1,Y1,sum_G(x[1:(2*r+1)],Vo,Yo),sum_H(x[(2*r+2):(4*r+2)],Vo,Yo),x[4*r+3],Vo,Yo)}
  #nll_optim <- function(x){G <<- sum_G(x[1:(2*r+1)],Vo,Yo); H <<- sum_H(x[(2*r+2):(4*r+2)],Vo,Yo); -logprob_vec_ghk(V1,Y1,G,H,x[4*r+3])}
  #grad_optim <- function(x){-grad_full_vec(V1,Y1,G,H,x[4*r+3],Vo,Yo)}
  best_res <- optim(par = c(rnorm(4*r+2),abs(rnorm(1))),
                    fn = nll_optim,
                    gr = grad_optim,
                    method = "L-BFGS-B",
                    lower=c(rep(-Inf, 4*r+2),0.0001),
                    control=list(trace=3*(report!=0),REPORT=max(report,1),maxit=maxit,pgtol=tol))
  if (runs>1){
    for (run in 1:(runs-1)){
      res <- optim(par = c(rnorm(4*r+2),abs(rnorm(1))),
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
      names(best_res$par) <- c(paste("g",left,sep=""), paste("q",right,sep=""), paste("r",right,sep=""), paste("h",left,sep=""), paste("s",right,sep=""), paste("t",right,sep=""), paste("k",left,sep=""))
      #c(paste("g",c(left,right),sep=""), paste("h",left,",",right,sep=""), paste("h",left,",",left,sep=""), paste("h",right,",",left,sep=""), paste("k",c(right,left),sep=""))
    else
      names(best_res$par) <- c(paste(c("g","h","k"),left,sep=""))
    best_res$n <- nrow(V)
    best_res$effective_df <- 4*r+3
    return (best_res)
  }
}

fit_C <- function(V, Y, left, right, maxit=200, runs=1, value_only=T, verbosity=0, step_size=0.1, tol=0.01, epsabs=1e-10, maxsize=1e-10, method=5){
  r <- length(right)
  par <- c(rnorm(4*r+2),abs(rnorm(1)))
  WY <- cbind(V[,left],Y[,left],V[,right],Y[,right])
  best_res <- .C("optim", sample_size_input=as.integer(nrow(V)), num_parents_input=as.integer(r), value=as.double(0), grad=as.double(numeric(4*r+3)), par=as.double(par), WY=as.double(WY), step_size=as.double(step_size), tol=as.double(tol), maxiter=as.integer(maxit), epsabs=as.double(epsabs), maxsize=as.double(maxsize), method=as.integer(method), verbosity=as.integer(verbosity))
  if (runs>1){
    for (run in 1:(runs-1)){
      par <- c(rnorm(4*r+2),abs(rnorm(1)))
      res <- .C("optim", sample_size_input=as.integer(nrow(V)), num_parents_input=as.integer(r), value=as.double(0), grad=as.double(numeric(4*r+3)), par=as.double(par), WY=as.double(WY), step_size=as.double(step_size), tol=as.double(tol), maxiter=as.integer(maxit), epsabs=as.double(epsabs), maxsize=as.double(maxsize), method=as.integer(method), verbosity=as.integer(verbosity))
      if (best_res$value > res$value){
        best_res <- res
      }
    }
  }
  if (value_only)
    return (best_res$value)
  else{
    if (r)
      names(best_res$par) <- c(paste("g",left,sep=""), paste("q",right,sep=""), paste("r",right,sep=""), paste("h",left,sep=""), paste("s",right,sep=""), paste("t",right,sep=""), paste("k",left,sep=""))
    #c(paste("g",c(left,right),sep=""), paste("h",left,",",right,sep=""), paste("h",left,",",left,sep=""), paste("h",right,",",left,sep=""), paste("k",c(right,left),sep=""))
    else
      names(best_res$par) <- c(paste(c("g","h","k"),left,sep=""))
    best_res$n <- nrow(V)
    best_res$effective_df <- 4*r+3
    return (best_res)
  }
}

fit_abk <- function(V, Y, left, right, use_C=TRUE, tol=1e-8, maxit=1000, runs=2, value_only=T, report=0){
  V <- V == 1
  if (length(right)==0){
    return (mle1d_ghk(V,Y,left,value_only))
  }
  if (use_C){
    return (fit_C(V, Y, left, right, maxit, runs=runs, value_only=value_only, verbosity=report, step_size=0.1, tol=0.01, epsabs=tol, maxsize=1e-10, method=5))
  } else {
    return (fit_R(V, Y, left, right, maxit, tol=tol, runs=runs, value_only=value_only, report=report))
  }
}


fit <- function(V, Y, left, right, parametrization="pms", value_only=T, control=list()){
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
      return (do.call("fit_pms_choose_degree", control))
    else
      return (do.call("fit_pms", control))
  } else if (parametrization == "abk") {
    return (do.call("fit_abk", control))
  } else {
    stop ("parametrization must be either pms or abk.")
  }
}

KS_stat <- function(x, ydist){ ### !!! ks.test does not have "abs". Why??
  ## sqrt(n)*D_n for KS test
  n <- length(x)
  x <- ydist(sort(x)) - (0:(n - 1))/n
  return (sqrt(n)*max(abs(c(1/n-x, x))))
}

###### Generate DAG
gen_dag <- function(seed, min_per_rank=2, max_per_rank=3, min_ranks=3, max_ranks=3, percent=40, plot=T){
  nodes <- 0
  set.seed(seed)
  if (min_ranks < max_ranks) {ranks <- sample(min_ranks:max_ranks, 1)
  } else {ranks <- min_ranks}
  cat("digraph {\n")
  adj <- NULL
  for (i in 1:ranks){
    if (min_per_rank < max_per_rank) {new_nodes <- sample(min_per_rank:max_per_rank, 1)
    } else {new_nodes <- min_per_rank}
    if (nodes){
      for (j in 1:nodes){
        for (k in 1:new_nodes){
          if (sample(0:99,1) < percent){
            cat(paste(j," -> ",k+nodes, ";\n", sep=""))
            adj <- rbind(adj, c(j, k+nodes))
          }
        }
      }
    }
    nodes <- nodes + new_nodes
  }
  cat("}\n")
  adj_mat <- matrix(0, nodes, nodes)
  adj_mat[adj] <- 1
  diag(adj_mat) <- 0
  if (plot){
    require(igraph)
    gr <- graph_from_adjacency_matrix(adj_mat, mode="directed")
    plot(gr)
    return (list("adj"=adj, "adj_mat"=adj_mat, "gr"=gr, "nodes"=nodes))
  }
  else {
    return (list("adj"=adj, "adj_mat"=adj_mat, "nodes"=nodes))
  }
}

make_dag <- function(m, mode, nrows=0, shuffle=FALSE, seed=NULL){
  adj_mat <- matrix(0,m,m)
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

gen_dag_chain <- function(seed, num_nodes, prop, plot=T){
  # percent: expected percentage of edges
  if (num_nodes <= 1) {return (NULL)}
  set.seed(seed)
  adj <- cbind(1:(num_nodes-1), 2:num_nodes)
  prop <- min((num_nodes-2)/num_nodes,max(0, prop-2/num_nodes)) # subtract the (n-1) edges from the chain
  adj_mat <- matrix(0, num_nodes, num_nodes)
  adj_mat[(1:(num_nodes-1))*(num_nodes+1)] <- 1
  adj <- rbind(which(upper.tri(adj_mat) & adj_mat == 0, arr.ind=T)[sample((num_nodes-1)*(num_nodes-2)/2, prop*num_nodes*(num_nodes-1)/2),], cbind(1:(num_nodes-1), 2:num_nodes))
  adj <- adj[order(adj[,1], adj[,2]), , drop=F]
  adj_mat[adj] <- 1
  if (plot){
    require(igraph)
    gr <- graph_from_adjacency_matrix(adj_mat, mode="directed")
    plot(gr)
    return (list("adj"=adj, "adj_mat"=adj_mat, "gr"=gr, "nodes"=num_nodes))
  }
  else {
    return (list("adj"=adj, "adj_mat"=adj_mat, "nodes"=num_nodes))
  }
}

adj_to_order <- function(adj_mat){
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

gen_dat <- function(seed, gen_para, adj_mat, n, k_mode=1, min_num=10, therun=1, gen_uniform_degree=1){
  if (therun == 1001) # If already tried to generate for 1000 times
    stop("Not able to generate good sample after 1000 runs: Some of the columns have fewer than ", min_num, " zeros or fewer than ", min_num, " nonzero values. Stopped.")
  if (therun==1) set.seed(seed) # Only set seed if it is the first trial
  diag(adj_mat) <- 0
  ord <- adj_to_order(adj_mat)
  adj_mat <- adj_mat[ord,ord]
  m <- dim(adj_mat)[1]
  Y <- matrix(0,n,m)
  V <- matrix(NA,n,m)
  for (j in 1:m){
    jks <- which(adj_mat[,j] != 0)
    if (length(jks) == 0){
      Y[,j] <- generate_g_h_k(n,-log(2*pi)/2,0,1)#Equivalent to generate_logitp_mu_sigmasq(n, 0, 0, 1) #generate_g_h_k(n,1,1,1)# rnorm(1,0,1), rnorm(1,0,1), 1)#rgamma(1,2,1))
      V[,j] <- (Y[,j] != 0)
    } else{
      if (gen_para == "abk"){
        if (k_mode == 0) {true_k_j <- 1/j #rgamma(1, 2, 1)#1
        } else if (k_mode == 1) {true_k_j <- 1
        } else if (k_mode == 2) {true_k_j <- j
        } else {stop("Bad k_mode specified.")}
        true_gg_j <- c(0, rep(1, 2 * length(jks))) #c(log(true_k_j/2/pi)/2, rep(g_strength, 2 * length(jks)))
        true_hh_j <- c(0,rep(1 * true_k_j, 2 * length(jks))) ## To stabilize the mean
        gs <- sum_G(true_gg_j,V[,jks],Y[,jks])
        hs <- sum_H(true_hh_j,V[,jks],Y[,jks])
        Y[,j] <- generate_one_each_g_h_k(
                    (gs-mean(gs))/sd(gs)-log(2*pi)/2, 
                    (hs-mean(hs))/sd(hs), 
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
          lps <- sum_G(true_logitp_j, V[,jks], Y[,jks])
          mus <- sum_H(true_mu_j, V[,jks], Y[,jks])
          Y[,j] <- generate_one_each_logitp_mu_sigmasq(
                    (lps-mean(lps))/sd(lps), 
                    (mus-mean(mus))/sd(mus), 
                    true_sigmasq_j
                  )
        } else {
          rhs_design <- full_design1(
                          Vo=V[,jks,drop=F], Yo=Y[,jks,drop=F], right=1:length(jks), 
                          V_degree=gen_uniform_degree, 
                          Y_degree=gen_uniform_degree, 
                          Y_V_degree=gen_uniform_degree)
          true_logitp_j <- c(0, 1-grepl("[\\*\\^]", colnames(rhs_design))*9/10)
          true_mu_j <- c(0, (1 - 2 * grepl("Y", colnames(rhs_design)))*(1-grepl("[\\*\\^]", colnames(rhs_design))*9/10)) ## Conform with the abk case
          lps <- sum_G(true_logitp_j, VO=rhs_design, YO=NULL)
          mus <- sum_H(true_mu_j, VO=rhs_design, YO=NULL, minus_Y=FALSE)
          Y[,j] <- generate_one_each_logitp_mu_sigmasq(
                    (lps-mean(lps))/sd(lps), 
                    (mus-mean(mus))/sd(mus), 
                    true_sigmasq_j
          )
        }
      } else {stop("Wrong gen_para specified.")}
      V[,j] <- (Y[,j] != 0)
      if (sum(V[,j]) < min_num || sum(1-V[,j]) < min_num)
        return (gen_dat(seed=seed, gen_para=gen_para, adj_mat=adj_mat, n=n, k_mode=k_mode, min_num=min_num, therun=therun+1, gen_uniform_degree=gen_uniform_degree))
    }
  }
  rev_ord <- order(ord)
  V <- V[,rev_ord]
  Y <- Y[,rev_ord]
  return (list("V"=V, "Y"=Y))
}


gen_dat_abk_tri_counterexample <- function(seed, n){
  set.seed(seed)
  m <- 3
  Y <- matrix(0,n,m)
  V <- matrix(NA,n,m)
  Y[,1] <- generate_g_h_k(n,log((1+sqrt(2*exp(1)*pi))/(1+sqrt(2*pi))),0,1)
  V[,1] <- (Y[,1] != 0)
  for (i in 1:n)
    Y[i,2] <- generate_g_h_k(1, sum_g(c(0,0,0),V[i,1],Y[i,1]),
                             sum_h(c(0,1,0),V[i,1],Y[i,1]), 1)
  V[,2] <- (Y[,2] != 0)
  for (i in 1:n)
    Y[i,3] <- generate_g_h_k(1, sum_g(c(0,0,1),V[i,2],Y[i,2]),
                             sum_h(c(0,0,0),V[i,2],Y[i,2]), 1)
  V[,3] <- (Y[,3] != 0)
  if (length(intersect(colSums(V), c(0,n-1,n))))
    warning("Bad sample: Some of the columns have no zeros or have at most 1 nonzero value.")
  return (list("V"=V, "Y"=Y))
}




gen_dat_abk_2d <- function(seed, n, g1,h1,k1,g2,h2,k2,q12,r12,s12,t12){
  set.seed(seed)
  Y <- matrix(0,n,2)
  V <- matrix(NA,n,2)
  Y[,1] <- generate_g_h_k(n,g1,h1,k1)
  V[,1] <- (Y[,1] != 0)
  for (i in 1:n)
    Y[i,2] <- generate_g_h_k(1, sum_g(c(g2,q12,r12),V[i,1],Y[i,1]),
                                 sum_h(c(h2,s12,t12),V[i,1],Y[i,1]), k2)
  V[,2] <- (Y[,2] != 0)
  if (length(intersect(colSums(V), c(0,n-1,n))))
    warning("Bad sample: Some of the columns have no zeros or have at most 1 nonzero value.")
  return (list("V"=V, "Y"=Y))
}

gen_dat_abk_all1 <- function(seed, adj_mat, n){
  set.seed(seed)
  m <- dim(adj_mat)[1]
  Y <- matrix(0,n,m)
  V <- matrix(NA,n,m)
  for (j in 1:m){
    jks <- which(adj_mat[,j] != 0)
    if (length(jks) == 0){
      Y[,j] <- generate_g_h_k(n, 1, 1, 1)#rgamma(1,2,1))
      V[,j] <- (Y[,j] != 0)
    } else{
      true_gg_j <- rep(1,2*length(jks)+1)#rnorm(2*length(jks)+1, 0, 2)
      true_hh_j <- rep(1,2*length(jks)+1)#rnorm(2*length(jks)+1, 0, 2)
      true_k11_j <- 1#rgamma(1, 2, 1)
      for (i in 1:n)
        Y[i,j] <- generate_g_h_k(1, sum_g(true_gg_j,V[i,jks],Y[i,jks]),
                                 sum_h(true_hh_j,V[i,jks],Y[i,jks]), true_k11_j)
      V[,j] <- (Y[,j] != 0)
    }
  }
  return (list("V"=V, "Y"=Y))
}

estimate <- function(dat, cur_set, parametrization, plot=T, verbose=T, control=list()){
  V <- dat$V; Y <- dat$Y
  m <- dim(V)[2]
  left_set <- setdiff(1:m, cur_set)
  while (length(left_set) > 1){
    min_is <- NULL
    which_min_is <- NULL
    for (ii in 1:length(left_set)){
      i <- left_set[ii]
      fi <- fit(V, Y, i, cur_set, parametrization=parametrization, value_only=T, control=control)
      min_i <- Inf
      for (j in left_set){
        if (i == j) next
        if (verbose)
          cat("i=",i,", j=",j, sep="")
        this_j <- fit(V, Y, j, c(i,cur_set), parametrization=parametrization, value_only=T, control=control) + fi
        if (this_j < min_i){
          min_i <- this_j; which_min_i <- j
        }
        if (verbose)
          cat(": ", this_j, "\n", sep="")
      }
      if (verbose)
        cat("****min for ", i, ": ", which_min_i, ", ", min_i, "\n", sep="")
      min_is <- c(min_is, min_i)
      which_min_is <- c(which_min_is, which_min_i)
    }
    for (ii in 1:length(left_set)){
      i <- left_set[ii]
      if (verbose)
        cat("(",i,"->",which_min_is[ii],"): ", min_is[ii], ", ",sep="")
    }
    if (verbose)
      cat("\n")
    winner <- left_set[which.min(min_is)]
    if (verbose)
      cat("Winner:", winner, "\n")
    cur_set <- c(cur_set, winner)
    if (verbose)
      cat("Topological order so far:", cur_set, "\n")
    left_set <- setdiff(left_set, winner)
  }
  cur_set <- c(cur_set,left_set)
  if (verbose)
    cat(paste("n=",n,", estimated topological order: ", paste(cur_set,collapse=" "), sep=""), "\n")
  if (plot && !is.null(dag$gr))
    plot(dag$gr, main = paste("n=",n,", estimated topological order: ", paste(cur_set,collapse=" "), sep=""))
  return (cur_set)
}




estimate_oneside <- function(dat, cur_set, parametrization, plot=T, verbose=T, control=list()){
  V <- dat$V; Y <- dat$Y
  m <- dim(V)[2]
  left_set <- setdiff(1:m, cur_set)
  while (length(left_set) > 1){
    fits <- numeric(length(left_set))
    for (ii in 1:length(left_set)){
      fits[ii] <- fit(V, Y, left_set[ii], cur_set, parametrization=parametrization, value_only=T, control=control)
      if (verbose)
        cat("**** ", left_set[ii], ": ", fits[ii], "\n", sep="")
    }
    if (verbose)
      cat("\n")
    winner <- left_set[which.min(fits)]
    if (verbose)
      cat("Winner:", winner, "\n")
    cur_set <- c(cur_set, winner)
    if (verbose)
      cat("Topological order so far:", cur_set, "\n")
    left_set <- setdiff(left_set, winner)
  }
  cur_set <- c(cur_set,left_set)
  if (verbose)
    cat(paste("n=",n,", estimated topological order: ", paste(cur_set,collapse=" "), sep=""), "\n")
  if (plot && !is.null(dag$gr))
    plot(dag$gr, main = paste("n=",n,", estimated topological order: ", paste(cur_set,collapse=" "), sep=""))
  return (cur_set)
}


estimate_ks_abk <- function(dat, cur_set, plot=T, verbose=T){
  V <- dat$V; Y <- dat$Y
  m <- dim(V)[2]
  left_set <- setdiff(1:m, cur_set)

  while (length(left_set) > 1){
    Ds <- numeric(length(left_set))
    for (ii in 1:length(left_set)){
      if (length(cur_set) == 0){
        resi <- Y[V[,left_set[ii]], left_set[ii]]
        resi <- (resi - mean(resi)) / sd(resi)
      } else {
        resi <- resid(lm(Y[V[,left_set[ii]],left_set[ii]] ~
                            cbind(V[V[,left_set[ii]], cur_set],
                            Y[V[,left_set[ii]], cur_set])))
        resi <- (resi - mean(resi)) / sd(resi)
      }
      Ds[ii] <- KS_stat(resi, function(x){pnorm(x, 0, 1)})
      if (verbose)
        cat("**** ", left_set[ii], ": ", Ds[ii], "\n", sep="")
    }
    if (verbose)
      cat("\n")
    winner <- left_set[which.min(Ds)]
    if (verbose)
      cat("Winner:", winner, "\n")
    cur_set <- c(cur_set, winner)
    if (verbose)
      cat("Topological order so far:", cur_set, "\n")
    left_set <- setdiff(left_set, winner)
  }
  cur_set <- c(cur_set,left_set)
  if (verbose)
    cat(paste("n=",n,", estimated topological order: ", paste(cur_set,collapse=" "), sep=""), "\n")
  if (plot && !is.null(dag$gr))
    plot(dag$gr, main = paste("n=",n,", estimated topological order: ", paste(cur_set,collapse=" "), sep=""))
  return (cur_set)
}

max_ll_complete <- function(dat, top_ord, parametrization, verbose=T, control=list()){
    V <- dat$V; Y <- dat$Y
    m <- dim(V)[2]
    if (!identical(as.integer(sort(top_ord)), 1:m))
      stop("The topological ordering must be a permutation of 1 through ", m, ".")
    nll <- 0
    for (ni in 1:length(top_ord)){
      node <- top_ord[ni]; #parents <- switch((ni==1)+1, top_ord[1:(ni-1)], c())
      if (ni == 1)
        fit_res <- mle1d_ghk(V,Y,node,value_only=!verbose)
      else
        fit_res <- fit(V, Y, node, top_ord[1:(ni-1)], parametrization=parametrization, value_only=!verbose, control=control)
      if (verbose){
        cat(paste(node,"|",paste(top_ord[1:(ni-1)],collapse=","),": ",paste(names(fit_res$par), round(fit_res$par,3), collapse=", ", sep=": "),sep=""), "\n")
        cat(paste("-log(likelihood) =", fit_res$value), "\n")
        nll <- nll + fit_res$value
      } else {nll <- nll+fit_res}
    }
    return (nll)
}

estimate_enumerate_complete <- function(dat, parametrization, verbose=T, control=list()){
  best_ord <- integer(0); best_nll <- Inf
  perms <- permutations(ncol(dat$V), ncol(dat$V))
  for (i in 1:nrow(perms)){
    nll <- max_ll_complete(dat, top_ord=perms[i,], parametrization=parametrization, verbose=verbose, control=control)
    if (verbose)
    cat(paste(paste(perms[i,], collapse=", "),": ", nll, sep=""), "\n")
    if (nll < best_nll){
      best_nll <- nll; best_ord <- perms[i,]
    }
  }
  return (list("order"=best_ord, "value"=best_nll))
}

#### NEEDS OPTIMIZATION!!!
find_adj <- function(m, cur_adj, last_num, perms){
  ## DFS
  #count <<- count+1
  res_set <- set()
  for (i in 1:nrow(perms)){
    res_set <- set_union(res_set, cur_adj[perms[i,], perms[i,]])
  }
  if (last_num < m*(m-1)/2){
    for (l in (last_num+1):(m*(m-1)/2)){
      cur_adj[upper.tri(cur_adj)][l] <- 1
      res_set <- set_union(res_set, find_adj(m, cur_adj, l, perms))
      cur_adj[upper.tri(cur_adj)][l] <- 0
    }
  }
  return (res_set)
}

find_adj2 <- function(m, cur_adj, last_num, perms, res_set){
  ## DFS, so
  new_set <- set()
  if (!set_is_proper_subset(set(cur_adj), res_set)){ # Could compare to upper tri matrices only
    #count <<- count+1
    for (i in 1:nrow(perms)){
      new_set <- set_union(new_set, cur_adj[perms[i,], perms[i,]])
    }
  }
  res_set <- set_union(new_set, res_set)
  if (last_num < m*(m-1)/2){
    for (l in (last_num+1):(m*(m-1)/2)){
      cur_adj[upper.tri(cur_adj)][l] <- 1
      res_set <- set_union(res_set, find_adj2(m, cur_adj, l, perms, res_set))
      cur_adj[upper.tri(cur_adj)][l] <- 0
    }
  }
  return (res_set)
}

find_adj3 <- function(m, cur_adj, last_num, perms, res_set, upper_set){
  ## DFS
  new_set <- new_upper_set <- set()
  if (!set_is_proper_subset(set(cur_adj), upper_set)){ # Could compare to upper tri matrices only
    #count <<- count+1
    for (i in 1:nrow(perms)){
      new_adj <- cur_adj[perms[i,], perms[i,]]
      new_set <- set_union(new_set, new_adj)
      if (sum(new_adj[lower.tri(new_adj)]) == 0)
        new_upper_set <-set_union(new_upper_set, new_adj)
    }
  }
  res_set <- set_union(new_set, res_set)
  upper_set <- set_union(new_upper_set, upper_set)
  if (last_num < m*(m-1)/2){
    for (l in (last_num+1):(m*(m-1)/2)){
      cur_adj[upper.tri(cur_adj)][l] <- 1
      rec_res <- find_adj3(m, cur_adj, l, perms, res_set, upper_set)
      res_set <- set_union(res_set, rec_res$res_set)
      upper_set <- set_union(upper_set, rec_res$upper_set)
      cur_adj[upper.tri(cur_adj)][l] <- 0
    }
  }
  return (list("res_set"=res_set,"upper_set"=upper_set))
}


enumerate_dags <- function(m, method=1){
  # p=5: 3.445176 min for method 1
  #count <- 0
  if (method==1) # Fastest for 3,4, no check for duplicates at all
    return (find_adj(m, matrix(0,m,m), 0, permutations(m, m)))
  if (method==2) # Slowest for 3,4, check for duplicates against all results
    return (find_adj2(m, matrix(0,m,m), 0, permutations(m, m), set()))
  if (method==3) # Medium for 3,4, check for duplicates against all results that are upper tri matrices
    return (find_adj3(m, matrix(0,m,m), 0, permutations(m, m), set(), set())$res_set)
}

bic_score <- function(avgnll, n, numpar){
  return (2*avgnll+log(n)/n*numpar)
}

bic_from_fit <- function(fit_res){
  return (bic_score(fit_res$value, fit_res$n, fit_res$effective_df))
}


estimate_enumerate_all0 <- function(dat, parametrization, verbose=T, all_dags=NULL){
  ## Number of nodes should be limited to 3
  m <- ncol(dat$V)
  if (is.null(all_dags))
    all_dags <- enumerate_dags(m, method=1)
  best_graph <- NULL; best_bic <- Inf
  count <- 0
  for (graph in all_dags){
    count <- count + 1
    bic <- 0
    for (j in 1:m){
      parents <- which(graph[,j]==1)
      if (length(parents)==0) {parents <- integer(0)}
      bic <- bic + bic_from_fit(fit(dat$V, dat$Y, left = j, right = parents, parametrization=parametrization, value_only=F, control=control))
    }
    if (verbose)
      cat(paste("Graph ", count, "/", length(all_dags), ": ", toString(graph),", bic: ", bic, sep=""), "\n")
    if (bic < best_bic){
      best_bic <- bic; best_graph <- graph
    }
  }
  return (list("order"=best_graph, "value"=best_bic))
}

estimate_enumerate_all <- function(dat, parametrization, verbose=T, all_dags=NULL){
  ## Number of nodes should be limited to 3
  m <- ncol(dat$V)
  if (is.null(all_dags))
    all_dags <- enumerate_dags(m, method=1)
  best_graph <- NULL; best_bic <- Inf
  count <- 0
  all_fits <- list()
  for (j in 1:m){
    all_fits[[paste(j,"|",sep="")]] <- bic_from_fit(fit(dat$V, dat$Y, left = j, right = integer(0), parametrization=parametrization, value_only=F, control=control))
    for (r in 1:(m-1)){
      all_combn <- combn((1:m)[-j], r)
      for (ci in 1:ncol(all_combn)){
        all_fits[[paste(j, paste(all_combn[,ci],collapse=","), sep="|")]] <- bic_from_fit(fit(dat$V, dat$Y, left = j, right = all_combn[,ci], parametrization=parametrization, value_only=F, control=control))
      }
    }
  }
  if (verbose)
    cat("All conditionals fitted.\n")
  for (graph in all_dags){
    count <- count + 1
    bic <- 0
    for (j in 1:m)
      bic <- bic + all_fits[[paste(j, paste(which(graph[,j]==1), collapse=","), sep="|")]]
    if (verbose)
      cat(paste("Graph ", count, "/", length(all_dags), ": ", toString(graph),", bic: ", bic, "\n", sep=""))
    if (bic < best_bic){
      best_bic <- bic; best_graph <- graph
    }
  }
  return (list("order"=best_graph, "value"=best_bic))
}

indep_ll_test <- function(x, y, S, suffStat){
  ## elts of suffStat:
  ## V, Y: dat$V, dat$Y
  ## parametrization
  ## conser: if TRUE (reject less->more independent than dependent), pick the min(ll(x|y,S)-ll(x|S), ll(y|x,S)-ll(y|S)), otherwise use the max
  ## control
  ## Returns: p-value for H0: x and y are independent (low p-value: dependent)
  if (length(S) == 0) S <- integer(0)
  xy_diff <- fit(suffStat$V,suffStat$Y,left=x,right=unique(S),parametrization=suffStat$parametrization, value_only=T, control=suffStat$control) -
   fit(suffStat$V,suffStat$Y,left=x,right=unique(c(y,S)), parametrization=suffStat$parametrization, value_only=T, control=suffStat$control)
  yx_diff <- fit(suffStat$V,suffStat$Y,left=y,right=unique(S),parametrization=suffStat$parametrization, value_only=T, control=suffStat$control) -
    fit(suffStat$V,suffStat$Y,left=y,right=unique(c(x,S)), parametrization=suffStat$parametrization, value_only=T, control=suffStat$control)
  if (suffStat$conser)
    return (1-pchisq(2*n*min(xy_diff, yx_diff),4))
  else
    return (1-pchisq(2*n*max(xy_diff, yx_diff),4))
}


build <- function(dat, ord, parametrization, alpha=0.05, verbose=T, control=list()){
  m <- ncol(dat$V); n <- nrow(dat$V)
  if (length(ord) != m) {stop(paste("Ord should have length", m))}
  est_adj <- NULL
  for (ii in 2:m){
    i <- ord[ii]
    if (verbose)
      cat("i = ",i,":\n",sep="")
    for (jj in 1:(ii-1)){
      j <- ord[jj]
      if (2*n*(fit(dat$V, dat$Y, left=i, right=setdiff(ord[1:(ii-1)],j), parametrization=parametrization, value_only=T, control=control) -
        fit(dat$V, dat$Y, left=i, right=ord[1:(ii-1)], parametrization=parametrization, value_only=T, control=control))
        >= qchisq(1-alpha,4)){
        est_adj <- cbind(est_adj, c(j,i))
        if (verbose)
          cat("\t(",j,"->",i,")\n",sep="")
      }
    }
  }
  adj_mat <- matrix(0, m, m)
  adj_mat[t(est_adj)] <- 1
  return (adj_mat)
}

build_bhy <- function(dat, ord, parametrization, alpha=0.05, verbose=T, control=list()){
  m <- ncol(dat$V); n <- nrow(dat$V)
  if (length(ord) != m) {stop(paste("Ord should have length", m))}
  pvals <- numeric(m*(m-1)/2)
  position <- 1
  # Run longer if test statistic negative
  control_max <- control
  control_max[["runs"]] <- max(control[["runs"]], 10) # Okay even if runs not in control
  control_max[["maxit"]] <- max(control[["maxit"]], 10000)
  # If control == control_max, no need to rerun
  control_already_max <- ifelse(length(control)==0, FALSE, identical(control[order(names(control))], control_max[order(names(control_max))]))
  for (ii in 2:m){
    i <- ord[ii]
    if (verbose)
      cat("i = ",i,":\n",sep="")
    full_ll_i <- fit(dat$V,dat$Y,left=i,right=ord[1:(ii-1)],parametrization=parametrization,value_only=T,control=control)
    full_ll_i_max <- NULL
    for (jj in 1:(ii-1)){
      j <- ord[jj]
      test_stat <- 2*n*(fit(dat$V,dat$Y,left=i,right=setdiff(ord[1:(ii-1)],j),parametrization=parametrization,value_only=T,control=control)
                        - full_ll_i) ###- fit(dat$V,dat$Y,left=i,right=ord[1:(ii-1)],parametrization=parametrization,value_only=T,control=control))
      if (test_stat <= 0 && !control_already_max) {#(test_stat >= qchisq(1-alpha, 4) || test_stat <= 0)
        if (verbose)
          cat("Refitting...")
        full_ll_i_max <- fit(dat$V,dat$Y,left=i,right=ord[1:(ii-1)],parametrization=parametrization,value_only=T,control=control_max)
        test_stat <- 2*n*(fit(dat$V,dat$Y,left=i,right=setdiff(ord[1:(ii-1)],j),parametrization=parametrization,value_only=T,control=control_max)
                          - full_ll_i_max)
      }
      pvals[position] <- 1-pchisq(test_stat, 4)
      position <- position+1
    }
  }
  pvals <- p.adjust(p=pvals,method="BY")
  est_adj_mat <- matrix(0,m,m)
  est_adj_mat[upper.tri(est_adj_mat)] <- (pvals <= alpha)
  #est_adj <- apply(which(est_adj_mat == 1, arr.ind=T), c(1,2), function(x){ord[x]})
  #return (est_adj)
  reverse_order <- order(ord)
  est_adj_mat <- est_adj_mat[reverse_order, reverse_order]
  return (est_adj_mat)
}

new_ZeroInfBICScore <- function(d, parametrization, control=list()){
  return(new("ZeroInfBICScore", dataV = d$V, dataY = d$Y,
             parametrization = parametrization,
             targets = list(integer(0)),
             target.index = rep(as.integer(1), nrow(d$Y)),
             nodes = colnames(d$Y),
             control = control))
}

estimate_build_bhy <- function(dat, parametrization, alpha, verbose=F, control=list()){
  ord <- estimate(dat, cur_set=integer(0), parametrization=parametrization, plot=F, verbose=verbose, control=control)
  build_bhy(dat, ord, parametrization=parametrization, alpha=alpha, verbose=verbose, control=control)
}

estimate_oneside_build_bhy <- function(dat, parametrization, alpha, verbose=F, control=list()){
  ord <- estimate_oneside(dat, cur_set=integer(0), parametrization=parametrization, plot=F, verbose=verbose, control=control)
  build_bhy(dat, ord, parametrization=parametrization, alpha=alpha, verbose=verbose, control=control)
}

gds_auto <- function(d, parametrization, maxDegree=Inf, verbose=FALSE, control=list()){
  res <- gds(new_ZeroInfBICScore(d, parametrization, control), fixedGaps=NULL,
             phase=c("forward","backward","turning"),
             iterate=TRUE, maxDegree=maxDegree, verbose=verbose)
  return (as(res$essgraph, "matrix")*1)
}

ges_auto <- function(d, parametrization, maxDegree=Inf, verbose=FALSE, control=list()){
  res <- ges(new_ZeroInfBICScore(d, parametrization, control), fixedGaps=NULL,
             phase=c("forward","backward","turning"),
             iterate=TRUE, maxDegree=maxDegree, verbose=verbose)
  return (as(res$essgraph, "matrix")*1)
}

arges_auto <- function(d, parametrization, maxDegree=Inf, verbose=FALSE, alpha=0.05, conser=F, control=list()){
  m <- ncol(d$V)
  suffStat <- control
  suffStat[["V"]] <- d$V
  suffStat[["Y"]] <- d$Y
  suffStat[["conser"]] <- conser
  antiCIG <- matrix(0, m, m) # non-edges in the CIG
  for (k in 1:(m-1)){
    for (j in (k+1):m){
      antiCIG[k,j] <- antiCIG[j,k] <- (indep_ll_test(x=k, y=j, S=(1:m)[-c(k,j)], parametrization=parametrization, suffStat=suffStat) > alpha)
    }
  }
  res <- ges(new_ZeroInfBICScore(d, parametrization, control), fixedGaps=antiCIG,
             phase=c("forward","backward","turning"),
             iterate=TRUE, maxDegree=maxDegree, verbose=verbose, adapt="vstructure")
  return (as(res$essgraph, "matrix")*1)
}


comp <- function(adj, est_adj, num_nodes){ ## assuming correct topological order
  joint <- unique(rbind(adj, est_adj))
  FP <- nrow(joint) - nrow(adj)
  TP <- nrow(est_adj) - FP
  CP <- nrow(adj)
  CN <- num_nodes*(num_nodes-1)/2 - CP
  return (c("TPR"=TP/CP, "FPR"=FP/CN))
}

comp_fdr <- function(adj, est_adj, num_nodes){ ## assuming correct topological order
  joint <- unique(rbind(adj, est_adj))
  FP <- nrow(joint) - nrow(adj)
  TP <- nrow(est_adj) - FP
  CP <- nrow(adj)
  CN <- num_nodes*(num_nodes-1)/2 - CP
  return (c("TPR"=TP/CP, "FPR"=FP/CN, "FDR"=FP/(TP+FP)))
}


which.group <- function(n, m){
  return (ceiling(n/m))
}
which.group <- Vectorize(which.group)

which.position <- function(n, m){
  if (n%%m == 0) {return (m)}
  else {return (n%%m)}
}
which.position <- Vectorize(which.position)


belong_to_CPDAG <- function(adj_mat, cpdag_res, method=1){
  # Checks if cpdag_res is a CPDAG of adj_mat
  cpdag_amat <- t(as(cpdag_res, "amat"))
  if (length(set(nrow(adj_mat),ncol(adj_mat),nrow(cpdag_amat),ncol(cpdag_amat))) > 1)
    stop("adj_mat and cpdag_res must be square matrices of the same dimensions.")
  m <- nrow(cpdag_amat)
  if (method == 1){
    if (max(cpdag_amat) > 1){
      stop("Max of cpdag_amat is larger than 1. Stopped.")
    }
    for (i in 1:(m-1)){
      for (j in 1:m){
        if (as.logical(adj_mat[i,j]+adj_mat[j,i]) != as.logical(cpdag_amat[i,j]+cpdag_amat[j,i])) # If adj wrong
          return (FALSE)
        if (cpdag_amat[i,j]+cpdag_amat[j,i] == 1){
          if (adj_mat[i,j]!=cpdag_amat[i,j] || adj_mat[j,i]!=cpdag_amat[j,i])
            return (FALSE)
        }
      }
    }
    return (TRUE)
  } else if (method == 2) {
    return (all(dag2cpdag(adj_mat)*1 == cpdag_amat))
  } else {
    stop("method can only be 1 or 2. Stopped.")
  }
}

DAG_eq <- function(amat1, amat2){
  if (length(dim(amat1)) != 2 || length(dim(amat2)) != 2 || length(unique(c(dim(amat1),dim(amat2)))) != 1){
    stop("amat1 and amat2 must be square matrices of the same dimension.")
  }
  # Checks if two dags represented by adj matrices are equivalent
  return (all(dag2cpdag(amat1) == dag2cpdag(amat2)))
}

update_correctness_adj <- function(old_correct_equivalent, est_adj, true_adj, name, verbose=F){
  good <- all(est_adj == true_adj)
  good_eq <- DAG_eq(est_adj, true_adj)
  if (verbose){
    if (!good)
      cat(paste(name, " failed, edges: ", paste(apply(which(est_adj==1, arr.ind=T), 1, function(i){paste(i,collapse="->")}), collapse=", "), sep=""), "\n")
    if (!good_eq)
      cat(paste(name, "also not equivalent to truth."), "\n")
  }
  return (c(good, good_eq)+old_correct_equivalent)
}


##### Score for gds and ges #####
### From https://github.com/cran/pcalg/blob/master/R/AllClasses.R
.tidyTargets <- function(p, targets, target.index = NULL) {
  stopifnot((p <- as.integer(p)) > 0)

  # Check and convert targets
  if (!is.list(targets) || !all(sapply(targets, is.numeric))) {
    stop("Argument 'targets' must be a list of integer vectors.")
  }
  rawTargets <- lapply(targets, function(v) unique(sort(as.integer(v))))
  targets <- unique(rawTargets)
  if (length(targets) < length(rawTargets)) {
    stop("List of targets must be unique.")
  }
  allTargets <- unlist(targets)
  if (length(allTargets) > 0) {
    if (any(is.na(allTargets))) {
      stop("Argument 'targets' must not contain NAs.")
    }
    min.max <- range(allTargets)
    if (min.max[1] <= 0 || min.max[2] > p) {
      stop("Targets are out of range.")
    }
  }

  # Check validity of target index, if provided
  if (!is.null(target.index)) {
    if (!is.numeric(target.index)) {
      stop("Argument 'target.index' must be an integer vector.")
    }
    target.index <- as.integer(target.index)
    min.max <- range(target.index)
    if (min.max[1] <= 0 || min.max[2] > length(targets)) {
      stop("Target index is out of range.")
    }
    # target.index <- match(rawTargets, targets)[target.index]
  }

  # Return value
  if (is.null(target.index)) {
    targets
  } else {
    list(targets = targets, target.index = target.index)
  }
}


setRefClass("ZeroParDAG",
            contains = "ParDAG",
            fields = list(.nodes = "character",
                          .in.edges = "list", .params = "list"),
            methods = list(
              initialize = function(nodes,
                      in.edges=replicate(length(nodes), integer(0), simplify = FALSE),
                      params=lapply(in.edges, function(l) numeric(4*length(l)+3))){
                .nodes <<- nodes
                .in.edges <<- in.edges
                .params <<- params
              },
              node.count = function(){
                return (length(.nodes))
              },
              simulate = function(n, target, int.level){
                stop("simulate() not implemented.")
              },
              edge.count = function(){
                return (sum(sapply(in.edges, length)))
              },
              mle.fit = function(score){
                stop("mle.fit() not implemented.")
              }
            )
)

setRefClass("ZeroInfBICScore",
            contains = "Score",

            fields = list(.nodes = "character", c.fcn = "character",
                          decomp = "logical", .pardag.class = "character",
                          .c.fcn = "character", pp.dat="list"),

            validity = function(object){
              #if (!all(dim(object$pp.dat$dataV) == dim(object$pp.dat$dataY)))
              #  return (FALSE)
              if (ncol(object$pp.dat$data) != 2*object$pp.dat$m)
              return(TRUE)
            },

            methods = list(
              ## overwrites initialize
              initialize = function(dataV,
                                    dataY,
                                    parametrization,
                                    targets = list(integer(0)),
                                    target.index = rep(as.integer(1), nrow(dataY)),
                                    nodes = colnames(dataY),
                                    control = list(),
                                    ...){

                if (is.null(nodes)) {
                  nodes <- as.character(1:ncol(dataV))
                }
                #targetList <- .tidyTargets(ncol(dataV), targets, target.index)
                .nodes <<- nodes
                pp.dat$targets <<- .tidyTargets(length(nodes), targets)

                #if (is.unsorted(targetList$target.index)) {
                #  perm <- order(targetList$target.index)
                #} else {
                #  perm <- seq_along(targetList$target.index)
                #}

                dataV <- as.matrix(dataV)
                dataV <- apply(dataV, 2, as.logical)
                dataY <- as.matrix(dataY)
                m <- pp.dat$m <<- pp.dat$vertex.count <<- ncol(dataV)
                n <- pp.dat$n <<- nrow(dataV)

                ##pp.dat$target.index <<- targetList$target.index[perm]
                ##pp.dat$data <<- cbind(dataV[perm,], dataY[perm,]) ## [V, Y]
                pp.dat$data <<- cbind(dataV, dataY)

                ##A <- !targets2mat(pp.dat$vertex.count, pp.dat$targets, pp.dat$target.index)
                ##pp.dat$non.int <<- lapply(seq_len(ncol(A)), function(i) which(A[, i]))
                # apply() cannot be used since we need a list, not a matrix.
                ##pp.dat$data.count <<- as.integer(colSums(A))
                pp.dat$data.count <<- rep(nrow(dataV), m) ## ??
                pp.dat$total.data.count <<- as.integer(nrow(dataV))

                pp.dat$control <<- control
                pp.dat$parametrization <<- parametrization

                decomp <<- TRUE
                c.fcn <<- "none"
                .pardag.class <<- "ZeroParDAG" ## ???

                pp.dat$local.score <<- function(vertex, parents) local.score(vertex, parents)
                pp.dat$global.score <<- function(edges) global.score(vertex, parents)
                pp.dat$local.fit <<- function(vertex, parents) local.fit(vertex, parents)
                pp.dat$global.fit <<- function(edges) global.fit(vertex, parents)
              }, ## initialize

              local.score = function(vertex, parents, ...){
                validate.vertex(vertex)
                validate.parents(parents)
                res <- fit(pp.dat$data[,1:pp.dat$m], pp.dat$data[,(pp.dat$m+1):(2*pp.dat$m)], left=vertex, right=parents, parametrization=pp.dat$parametrization, value_only=F, control=pp.dat$control)
                return (-bic_from_fit(res))
              },

              local.fit = function(vertex, parents, ...){
                validate.vertex(vertex)
                validate.parents(parents)
                return (fit(pp.dat$data[,1:pp.dat$m], pp.dat$data[,(pp.dat$m+1):(2*pp.dat$m)], left=vertex, right=parents, parametrization=pp.dat$parametrization, value_only=F, control=pp.dat$control)$par)
              }
            )
)


#### Own implementation of gds
add_edge <- function(V, Y, parametrization, in_edges, degrees, cur_adj_mat, nodewise_score, add_dec, maxDegree, maxInDegree, fixedGaps, control=list()){
  best_in <- best_out <- -1
  best_decrease <- 0
  best_new_score <- Inf
  m <- ncol(V); n <- nrow(V)
  for (ii in 1:m){
    if (degrees[ii] < maxDegree && length(in_edges[[ii]]) < maxInDegree){
      for (oi in setdiff(1:m, c(ii, in_edges[[ii]]))){
        if (degrees[oi] < maxDegree && (is.null(fixedGaps) || !fixedGaps[oi,ii])) {
          cur_adj_mat[oi, ii] <- 1
          if (!(ii %in% in_edges[[oi]]) && isAcyclic(cur_adj_mat)){
            if (!is.na(add_dec[oi, ii])){
              dec <- add_dec[oi, ii]
              new_score <- nodewise_score[ii] - dec
            }
            else{
              new_score <- bic_from_fit(fit(V, Y, left=ii, right=c(oi,in_edges[[ii]]), parametrization=parametrization, value_only=F, control=control))
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

remove_edge <- function(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, rem_dec, control=list()){
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
        new_score <- bic_from_fit(fit(V, Y, ii, setdiff(in_edges[[ii]],oi), parametrization=parametrization, value_only=F, control=control))
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

turn_edge <- function(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, add_dec, rem_dec, maxInDegree, control=list()){
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
      if (isAcyclic(cur_adj_mat)){
        if (!is.na(rem_dec[oi,ii])){
          dec_in <- rem_dec[oi,ii]
          new_score_in <- nodewise_score[ii] - dec_in
        } else {
          new_score_in <- bic_from_fit(fit(V, Y, ii, setdiff(in_edges[[ii]],oi), parametrization=parametrization, value_only=F, control=control))
          dec_in <- nodewise_score[ii] - new_score_in
          rem_dec[oi,ii] <- dec_in
        }
        if (!is.na(add_dec[ii,oi])){
          dec_out <- add_dec[ii,oi]
          new_score_out <- nodewise_score[oi] - dec_out
        } else {
          new_score_out <- bic_from_fit(fit(V, Y, oi, c(ii,in_edges[[oi]]), parametrization=parametrization, value_only=F, control=control))
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

my_gds <- function(V, Y, parametrization, fixedGaps=NULL, maxDegree=Inf, maxInDegree=Inf, init=NULL, verbose=F, control=list()){
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
  if (!isAcyclic(cur_adj_mat))
    stop("Initial graph must be acyclic.")
  degrees <- rowSums(cur_adj_mat) + colSums(cur_adj_mat)  # only need the degree, not the list
  nodewise_score <- sapply(1:m, function(i){bic_from_fit(fit(V, Y, left=i, right=in_edges[[i]], parametrization=parametrization, value_only=F, control=control))})
  old_best_sum <- Inf
  best_sum <- sum(nodewise_score)
  iter_count <- 0
  add_dec <- rem_dec <- matrix(NA,m,m)
  while (best_sum < old_best_sum-1e-9){
    old_best_sum <- best_sum
    iter_count <- iter_count + 1
    if (verbose){cat(paste("Iteration ", iter_count, ": BIC score ", old_best_sum, sep=""), "\n")}
    while (TRUE){
      res_add <- add_edge(V, Y, parametrization, in_edges, degrees, cur_adj_mat, nodewise_score, add_dec, maxDegree, maxInDegree, fixedGaps, control)
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
      if (verbose) cat(paste("Added ", res_add$best_out, "->", res_add$best_in, ", new BIC ", best_sum, sep=""), "\n")
    }
    while (TRUE){
      res_remove <- remove_edge(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, rem_dec, control)
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
      if (verbose) cat(paste("Removed ", res_remove$best_out, "->", res_remove$best_in, ", new BIC ", best_sum, sep=""), "\n")
    }
    while (TRUE){
      res_turn <- turn_edge(V, Y, parametrization, in_edges, cur_adj_mat, nodewise_score, add_dec, rem_dec, maxInDegree, control)
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
      if (verbose) cat(paste("Turned ", res_turn$best_out, "->", res_turn$best_in, ", new BIC ", best_sum, sep=""), "\n")
    }
  }
  return (cur_adj_mat)
}
