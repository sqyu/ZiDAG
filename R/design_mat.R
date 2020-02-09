rhs_formula <- function(r, d){
  # "Vo[,1]+...+Vo[,r]+polym(Yo[,1],...,Yo[,r], degree=d)
  #1 old naive: return (list("V1"=paste(paste(paste("Vo[,",1:r,"]",sep=""), collapse=" + "), "+",
  #        paste("polym(", paste(paste("Yo[,",1:r,"]",sep=""), collapse=", "), ", raw=T, degree=", d, ")", sep="")),
  #      "Y1"=paste(paste(paste("Vo[V1,",1:r,"]",sep=""), collapse=" + "), "+",
  #        paste("polym(", paste(paste("Yo[V1,",1:r,"]",sep=""), collapse=", "), ", raw=T, degree=", d, ")", sep=""))
  #))
  #2 old shared design between V and Y, but no interaction on V: return (paste(paste(paste("Vo[,",right,"]",sep=""), collapse=" + "), "+",
  #              paste("polym(", paste(paste("Yo[,",1:r,"]",sep=""), collapse=", "), ", raw=T, degree=", d, ")", sep="")))
  return (paste("0 + polym(", paste(paste("Yo[,",1:r,"]",sep=""), collapse=", "), ", raw=T, degree=", d, ")", sep=""))
}



design_V <- function(Vo, degree, right_names, only=FALSE){
  ## Assumes: degree >= 0. Ensured when calling from fit_pms_fromVY
  ## Returns: interactions of degree <= min(degree, r) if only==FALSE, otherwise = degree if degree <= r else NULL. Also NULL if degree==0.
  if (degree == 0){return(NULL)}
  if (degree < 0){stop("Degree must be >= 0.")}
  r <- ncol(Vo)
  x <- name <- NULL
  if (only && (degree > r))
    return (NULL)
  for (d in ifelse(only, degree, 1):degree){
    if (d == 1){
      x <- cbind(x,Vo)
      name <- c(name, paste("W", right_names, sep=""))
    } else{
      comb <- combn(r,d)
      x <- cbind(x, apply(comb, 2, function(com){apply(Vo[,com,drop=F],1,prod)}))
      name <- c(name, apply(comb, 2, function(com){paste(paste("W",right_names[com],sep=""),collapse="*")}))
    }
  }
  colnames(x) <- name
  return (x)
}

split_sum <- function(n,k){
  ## Assumes: n >= 0, k > 0
  ## Returns: each column represents an ordered partition of n into k non-negative numbers
  if (n < 0 || k <= 0){
    stop("n must be non-negative and k must be positive.")
  }
  if (k > 1){
    return (apply(combn(n+k-1,k-1), 2, function(tmp){c(tmp,n+k)-c(0,tmp)-1}))
  } else {
    return (matrix(n,1,1))
  }
}

#(1,2,"Y") -> "Y1^2"
name_fun <- function(name, deg, prefix=""){
  if(deg==1){
    paste(prefix,name,sep="")
  } else if (deg>1){
    paste(prefix,name,"^",deg,sep="")
  } else {""}
}


design_Y_from_poly <- function(Yo, Y_degree, right_names){
  r <- ncol(Yo)
  if (r != length(right_names)) stop("ncol(Yo) == length(right_names) must be TRUE.")
  if (Y_degree == 1){
    colnames(Yo) <- paste("Y", right_names, sep="")
    return (Yo)
  }
  design_matrix_Y <- model.matrix(formula(paste(" ~", rhs_formula(r, Y_degree))))
  if (ncol(design_matrix_Y) == 1) # If r == 1 and Y_degree == 1, model.matrix does not give the degree after ")"
    colnames(design_matrix_Y)[1] <- paste("raw = T, degree = ", Y_degree, ")1", sep = "")
  #"polym(...)0.1.2.3.4"->c(0,1,2,3,4)
  design_digits <- matrix(sapply(colnames(design_matrix_Y), function(s){as.integer(strsplit(strsplit(s, paste("raw = T, degree = ", Y_degree, ")", sep = ""))[[1]][2], "\\.")[[1]])}), nrow=r)
  #"polym(...)0.1.2.3.4"->c(2,3,4,5)
  design_which_included <- apply(design_digits, 2, function(d){which(d!=0)})
  #"polym(...)0.1.2.3.4"->"Y2*Y3^2*Y4^3*Y4^5"
  design_names <- sapply(1:length(design_which_included), function(di){paste(sapply(design_which_included[[di]],function(dd){name_fun(right_names[dd],design_digits[dd,di],"Y")}), collapse="*")})
  colnames(design_matrix_Y) <- design_names
  return (design_matrix_Y)
}

design_Y <- function(Yo, degree, right_names, only=FALSE){
  ## Assumes: degree >= 0. Ensured when calling from fit_pms_fromVY
  ## Returns: interactions of degree <= degree. NULL if degree == 0.
  ## Slower than model.matrix(formula(paste(" ~", rhs_formula(r, degree))))
  if (degree == 0){return (NULL)}
  if (degree < 0){stop("Degree must be >= 0.")}
  r <- ncol(Yo)
  x <- NULL
  name <- NULL
  if (!only)
    return (design_Y_from_poly(Yo, degree, right_names))
  #for (d in ifelse(only, degree, 1):degree){
  if (degree == 1){
    x <- Yo
    name <- c(name, paste("Y",right_names,sep=""))
  }
  else{
    comb <- split_sum(degree, r)
    x <- cbind(x, apply(comb, 2, function(com){apply(sweep(Yo[,com!=0,drop=FALSE], 2, com[com!=0], "^"), 1, prod)}))
    name <- c(name, apply(comb, 2, function(com){gsub("\\^1", "", paste(paste("Y",right_names[com!=0],"^",com[com!=0],sep=""),collapse="*"))}))
  }
  #}
  colnames(x) <- name
  return (x)
}


quad_design <- function(Vo, Yo, right_names){
  if (length(dim(Vo))!=2 || length(dim(Yo))!=2 || any(dim(Vo)!=dim(Yo)) || ncol(Vo)!=length(right_names))
    stop ("Vo and Yo must be matrices of the same dimension and same number of columns as right_names.")
  r <- ncol(Vo)
  VoYo <- cbind(Vo, Yo)
  colnames(VoYo) <- c(paste("W", right_names, sep=""), paste("Y", right_names, sep=""))
  res <- cbind(VoYo)
  names <- colnames(VoYo)
  res <- cbind(res, do.call("cbind", lapply(seq_len(r-1), function(ri){Vo[,ri]*Vo[,(ri+1):r]})))
  names <- c(names, do.call("c", sapply(seq_len(r-1), function(ri){paste("W",right_names[ri],"*","W",right_names[(ri+1):r],sep="")}, simplify = FALSE)))
  res <- cbind(res, Yo^2)
  names <- c(names, paste("Y", right_names, "^2", sep=""))
  res <- cbind(res, do.call("cbind", lapply(seq_len(r-1), function(ri){Yo[,ri]*Yo[,(ri+1):r]})))
  names <- c(names, do.call("c", sapply(seq_len(r-1), function(ri){paste("Y",right_names[ri],"*","Y",right_names[(ri+1):r],sep="")}, simplify = FALSE)))
  res <- cbind(res, do.call("cbind", lapply(1:r, function(ri){Vo[,ri]*Yo[,-ri]})))
  if (r > 1)
    names <- c(names, c(sapply(1:r, function(ri){paste("W",right_names[ri],"*","Y",right_names[-ri],sep="")})))
  colnames(res) <- names
  return (res)
}

full_design1 <- function(Vo, Yo, right, V_degree=1, Y_degree=1, Y_V_degree=1){
  ## Y_degree >= 0, V_degree >= 0, Y_V_degree can be 0; if Y_degree and V_degree are both non-zero, Y_V_degree can be between 2 and Y_degree+V_degree
  ### Temporarily rules out V1*Y1 by checking names of variables using colnames.
  #rhs_form_Y <- rhs_formula(r, Y_degree)
  #design_matrix_Y <- model.matrix(formula(paste(" ~", rhs_form_Y)))
  #design_matrix_V <- design_V(Vo = Vo, degree=degree, right_names=right)
  r <- ncol(Vo)
  if (r != length(right) || r != ncol(Yo)) stop("ncol(Vo) == ncol(Yo) == length(right) must be TRUE.")
  if (Y_degree==0){
    warning("Y_degree = 0 provided.\n")
  } else if (Y_degree < 0){stop("Negative degree provided for Y. Stopped.")}
  if (V_degree==0){
    warning("V_degree = 0 provided.\n")
  } else if (V_degree < 0){stop("Negative degree provided for V. Stopped.")}
  if (Y_V_degree < 0){stop("Negative degree provided for interactions between Y and V. Stopped.")}
  else if (Y_V_degree > 0) {
    if (Y_V_degree == 1)
      warning("Y_V_degree=1 treated as no interaction between Y and V.")
    if (V_degree == 0 || Y_degree == 0)
      stop("Interaction between V and Y not allowed when V_degree or Y_degree is 0.")
    if (Y_V_degree > V_degree + Y_degree)
      stop("Y_V_degree can be V_degree + Y_degree at most.")
  }
  if (Y_degree == 0 || V_degree == 0 || Y_V_degree <= 1){
    design_matrix_V <- design_V(Vo = Vo, degree=min(r,V_degree), right_names=right, only=FALSE)
    design_matrix_Y <- design_Y(Yo = Yo, degree=Y_degree, right_names=right, only=FALSE)
    design_matrix_Y_V <- NULL
  } else if (V_degree == 2 && Y_degree == 2 && Y_V_degree == 2) {
    return (quad_design(Vo, Yo, right))
  } else {
    ## effective degrees
    V_degree <- min(r, V_degree)
    Y_V_degree <- min(V_degree+Y_degree, Y_V_degree) # In fact, min(Y_V_degree, r+Y_degree, V_degree+Y_degree)
    ## !! To avoid things like W1*Y1, must have Y_V_degree <- min(min(r-1,V_degree)+Y_degree, Y_V_degree)
    design_matrix_V <- lapply(1:min(r,V_degree), function(d){design_V(Vo = Vo, degree=d, right_names=right, only=TRUE)})
    design_matrix_Y <- lapply(1:Y_degree, function(d){design_Y(Yo = Yo, degree=d, right_names=right, only=TRUE)})
    ### Interactions between Y and V, but does not rule out V1*Y1!!!
    design_matrix_Y_V <- do.call("cbind", lapply(2:Y_V_degree, function(td){
      do.call("cbind", lapply(max(td-Y_degree,1):min(td-1,V_degree), function(vd){ ## max(td-Y_degree,1) <= min(td-1,V_degree) is guaranteed
        do.call("cbind", lapply(1:ncol(design_matrix_V[[vd]]), function(vi){
          v_names <- str_to_var(colnames(design_matrix_V[[vd]])[vi])
          good_ys <- sapply(colnames(design_matrix_Y[[td-vd]]), function(s){!length(intersect(str_to_var(s), v_names))})
          tmp <- design_matrix_V[[vd]][,vi] * design_matrix_Y[[td-vd]][,good_ys,drop=FALSE]
          if (length(tmp)){ ## colnames(numeric(0)) will trigger an error
            colnames(tmp) <- paste(colnames(design_matrix_V[[vd]])[vi], colnames(design_matrix_Y[[td-vd]])[good_ys], sep="*")
            return (tmp)
          } else {
            return (NULL)
          }
        })) # Loop thru all possible sums of v_deg and t_deg, then thru v_deg, then thru columns of design_matrix_V[[v_deg]]
      }))
    }))
    design_matrix_V <- do.call("cbind", design_matrix_V)
    design_matrix_Y <- do.call("cbind", design_matrix_Y)
  }
  return (cbind(design_matrix_V, design_matrix_Y, design_matrix_Y_V))
}

full_design2 <- function(Vo, Yo, right, V_degree=1, Y_degree=1, Y_V_degree=1){
  r <- length(right)
  if (r != ncol(Yo) || r != ncol(Vo)) stop("ncol(Vo) == ncol(Yo) == length(right) must be TRUE.")
  if (Y_degree==0){
    warning("Y_degree = 0 provided.\n")
  } else if (Y_degree < 0){stop("Negative degree provided for Y. Stopped.")}
  if (V_degree==0){
    warning("V_degree = 0 provided.\n")
  } else if (V_degree < 0){stop("Negative degree provided for V. Stopped.")}
  if (Y_V_degree < 0){stop("Negative degree provided for interactions between Y and V. Stopped.")}
  else if (Y_V_degree > 0) {
    if (Y_V_degree == 1)
      warning("Y_V_degree=1 treated as no interaction between Y and V.")
    if (V_degree == 0 || Y_degree == 0)
      stop("Interaction between V and Y not allowed when V_degree or Y_degree is 0.")
    if (Y_V_degree > V_degree + Y_degree)
      stop("Y_V_degree can be V_degree + Y_degree at most.")
  }
  design_matrix_Y <- model.matrix(formula(paste(" ~", rhs_formula(r, Y_degree))))
  if (ncol(design_matrix_Y) == 1) # If r == 1 and Y_degree == 1, model.matrix does not give the degree after ")"
    colnames(design_matrix_Y)[1] <- paste("raw = T, degree = ", Y_degree, ")1", sep = "")
  #colnames(design_matrix_Y)[1] <- paste("raw = T, degree = ", Y_degree, ")", paste(numeric(r), collapse="."), sep = "")
  zero_indices <- sapply(1:ncol(Vo), function(i){which(Vo[,i]==0)})

  #"polym(...)0.1.2.3.4"->c(0,1,2,3,4)
  design_digits <- matrix(sapply(colnames(design_matrix_Y), function(s){as.integer(strsplit(strsplit(s, paste("raw = T, degree = ", Y_degree, ")", sep = ""))[[1]][2], "\\.")[[1]])}), nrow=r)
  #"polym(...)0.1.2.3.4"->10
  design_degrees <- colSums(design_digits)
  #"polym(...)0.1.2.3.4"->c(2,3,4,5)
  design_which_included <- apply(design_digits, 2, function(d){which(d!=0)})
  #"polym(...)0.1.2.3.4"->"Y2*Y3^2*Y4^3*Y4^5"
  design_names <- sapply(1:length(design_which_included), function(di){paste(sapply(design_which_included[[di]],function(dd){name_fun(right[dd],design_digits[dd,di],"Y")}), collapse="*")})
  colnames(design_matrix_Y) <- design_names
  design_YV <- Reduce("cbind", lapply(1:length(design_which_included), function(di){
    if (design_degrees[di] >= Y_V_degree || V_degree == 0 || length(design_which_included[[di]]) == r) {NULL}
    else {
      left_overs <- setdiff(1:r, design_which_included[[di]])
      eff_deg <- min(length(left_overs), Y_V_degree-design_degrees[di], V_degree) # effective degrees left
      if (eff_deg == 0) {return (design_which_included[[di]])}
      else {return (Reduce("cbind", lapply(1:eff_deg, function(V_deg){
        if (length(left_overs) == 1) { ## Otherwise utils::combn will treat left_overs as the number of elements to choose from in this case
          tmp <- matrix(design_matrix_Y[,di]*Vo[,left_overs])
          colnames(tmp) <- paste(design_names[di], "*W", right[left_overs], sep="")
          return (tmp)
        } else {
          combns <- utils::combn(left_overs, V_deg)
          which_cols_has_this_index <- list()
          for (ci in 1:ncol(combns))
            for (cc in combns[,ci])
              which_cols_has_this_index[[toString(cc)]] <- c(which_cols_has_this_index[[toString(cc)]],ci)
          tmp <- matrix(rep(design_matrix_Y[,di],each=ncol(combns)), ncol=ncol(combns), byrow=T)
          colnames(tmp) <- rep(design_names[di], ncol(tmp))
          for (i in left_overs){
            tmp[zero_indices[[i]], which_cols_has_this_index[[toString(i)]]] <- 0
            colnames(tmp)[which_cols_has_this_index[[toString(i)]]] <- paste(colnames(tmp)[which_cols_has_this_index[[toString(i)]]],"*W",right[i],sep="")
          }
          if (design_names[di] == "") {colnames(tmp) <- substr(colnames(tmp), 2, sapply(colnames(tmp), nchar))}
          return (tmp)
        }
      })))}
    }
  }))
  if (V_degree == 0){
    design_V <- NULL
  } else if (r == 1) { ## Otherwise utils::combn will treat left_overs as the number of elements to choose from in this case
    design_V <- Vo
  } else {
    design_V <- Reduce("cbind", lapply(1:min(r,V_degree), function(V_deg){
      combns <- utils::combn(r, V_deg)
      which_cols_has_this_index <- list()
      for (ci in 1:ncol(combns))
        for (cc in combns[,ci])
          which_cols_has_this_index[[toString(cc)]] <- c(which_cols_has_this_index[[toString(cc)]],ci)
      tmp <- matrix(1, ncol=ncol(combns), nrow=nrow(Vo), byrow=T)
      colnames(tmp) <- rep("", ncol(tmp))
      for (i in 1:r){
        tmp[zero_indices[[i]], which_cols_has_this_index[[toString(i)]]] <- 0
        colnames(tmp)[which_cols_has_this_index[[toString(i)]]] <- paste(colnames(tmp)[which_cols_has_this_index[[toString(i)]]],"*W",right[i],sep="")
      }
      return (tmp)
    }))
  }
  design_Y_V_YV <- cbind(design_matrix_Y, design_V, design_YV)
  return (design_Y_V_YV)
}


# Test two design matrix constructors are equivalent
#for (r in 1:10){
#  cat(r); right <- (11-r):10
#  for (Y_degree in 1:6){
#    cat(", ", Y_degree, "\n")
#    for (V_degree in 0:6){
#      cat(", , ", V_degree, ": ", switch((V_degree!=0)+1, 0, 0:(Y_degree+V_degree)), "\n")
#      for (Y_V_degree in switch((V_degree!=0)+1, 0, 0:(Y_degree+V_degree))){
#        des1 <- full_design1(d$V[,right,drop=F], d$Y[,right,drop=F], right, V_degree, Y_degree, Y_V_degree)
#        des0 <- full_design2(d$V[,right,drop=F], d$Y[,right,drop=F], right, V_degree, Y_degree, Y_V_degree)
#        dif <- sum(abs(des1[,order(colSums(des1))] - des0[,order(colSums(des0))]))
#        if (abs(dif) > 1e-10){
#          cat("r=",r,", Y_degree=",Y_degree,", V_degree=",V_degree,", Y_V_degree=",Y_V_degree,": ", dif, "\n")
#          stop("")
#        }}}}}




num_par_V <- function(V_degree, r){sapply(1:min(V_degree,r),function(d){choose(r,d)})}
num_par_Y <- function(Y_degree, r){sapply(1:Y_degree,function(d){choose(d+r-1,r-1)})}
num_par_no_Y_V <- function(V_degree, Y_degree, r){(Y_degree+1)/r*choose(Y_degree+r,r-1)+sum(sapply(1:min(V_degree,r),function(i){choose(r,i)}))}
num_par_per_regression <- function(V_degree, Y_degree, Y_V_degree, r, alias=FALSE){
  ### alias: If TRUE, calculate the naive number of parameters
  ### ; if FALSE, exclude terms like W1*Y1 and W1*W2*Y2^2
  if (Y_degree < 0 || V_degree < 0 || Y_V_degree < 0 ||
      Y_V_degree == 1 || (Y_V_degree > 0 && (V_degree == 0 || Y_degree == 0 ||
                                             Y_V_degree > V_degree + Y_degree))){stop("Invalid input.")}
  V_degree <- min(r, V_degree)
  Y_V_degree <- min(min(r-1,V_degree)+Y_degree, Y_V_degree)
  if (Y_V_degree == 0){return (num_par_no_Y_V(V_degree, Y_degree, r))}
  if (alias){
    V_nums <- num_par_V(V_degree, r)
    Y_nums <- num_par_Y(Y_degree, r)
    return(1+sum(V_nums)+sum(Y_nums)+sum(sapply(2:Y_V_degree, function(td){sum(sapply(max(td-Y_degree,1):min(td-1,V_degree), function(vd){V_nums[vd]*Y_nums[td-vd]}))})))
  } else{
    V_nums <- num_par_V(V_degree, r)
    Y_nums <- sapply(1:r, function(ri){num_par_Y(Y_degree, ri)})
    if (is.null(dim(Y_nums))){Y_nums <- matrix(Y_nums, nrow=Y_degree, ncol=r)}
    return(1+sum(V_nums)+sum(Y_nums[,r])+sum(sapply(2:Y_V_degree, function(td){
      sum(sapply(max(td-Y_degree,1):min(td-1,r-1,V_degree), function(vd){V_nums[vd]*Y_nums[td-vd,r-vd]}))})))
  }
}

str_to_var <- function(s){
  ## "W1^4*W2^2*W3" (or V or Y) to c("1","2","3")
  sapply(strsplit(s, "\\*")[[1]], function(ss){strsplit(strsplit(ss, "\\^")[[1]][1], "[VWY]")[[1]][2]})
}



### For producing designs with Y_degree == V_degree == Y_V_degree,
###   separated by degree
flip_1s <- function(vec) {
  # Given a vec with k elements equal to 1, returns the matrix of size
  # length(vec) x 2^k - 1, where each column is the result of flipping
  # the signs of at least one of the 1s
  pos_1s <- which(vec == 1)
  if (length(pos_1s) == 1) {
    vec[pos_1s] <- -1
    return (as.matrix(vec, ncol=1))
  }
  return (do.call("cbind", sapply(1:length(pos_1s), function(num_flips) {
    to_flips <- combn(pos_1s, num_flips)
    apply(to_flips, 2, function(to_flip) {tmp<-vec; tmp[to_flip]<--1; tmp})
  })))
}

split_sum_neg1 <- function(n, k){
  ## Assumes: n >= 0, k > 0
  ## Returns: each column represents an ordered sequence of k integers >= -1 whose abs sum equals to n
  if (n < 0 || k <= 0){
    stop("n must be non-negative and k must be positive.")
  }
  no_neg1 <- split_sum(n, k)
  if (n == 1) # Must have exactly one 1 in each column
    return (cbind(no_neg1, -no_neg1))
  if (n == 0 || k == 1)  # n == 0 split into k integers OR n > 1 split into 1 integer -> no 1 in the result
    return (no_neg1)
  extras <- do.call("cbind", apply(no_neg1, 2, function(vec){
    if (1 %in% vec) flip_1s(vec)
    else NULL
  }))
  return (cbind(no_neg1, extras))
}

one_col <- function(Vo, Yo, vec, right_names){
  pos_vs <- which(vec == -1)
  pos_ys <- which(vec > 0)
  if (length(pos_vs) == 0 && length(pos_ys) == 0)
    stop("vec must not be all zeros.")
  if (length(pos_vs) == 0) {
    x <- rep(1, nrow(Vo))
    name <- c()
  }
  else {
    x <- apply(Vo[, pos_vs, drop=FALSE], 1, prod)
    name <- c(paste("W", right_names[pos_vs], collapse="*", sep=""))
  }
  if (length(pos_ys) > 0) {
    x <- x * apply(sweep(Yo[, pos_ys, drop=FALSE], 2, vec[pos_ys], "^"), 1, prod)
    name <- c(name, gsub("\\^1", "", paste(paste("Y", right_names[pos_ys], sep=""), vec[pos_ys], sep="^", collapse="*")))
    name <- paste(name, collapse="*")
  }
  x <- as.matrix(x, ncol=1)
  colnames(x) <- name
  return (x)
}

full_design_same_exact_degree <- function(Vo, Yo, right_names, exact_degree) {
  if (exact_degree == 0){return (NULL)}
  if (exact_degree < 0){stop("Degree must be >= 0.")}
  r <- ncol(Vo)
  if (r != length(right_names) || r != ncol(Yo)) stop("ncol(Vo) == ncol(Yo) == length(right_names) must be TRUE.")
  x <- NULL
  name <- NULL
  ass <- split_sum_neg1(exact_degree, r)  # -1 corresponds to Vo, >= 1 corresponds to the degree of Yo
  do.call("cbind", lapply(1:ncol(ass), function(i){one_col(Vo, Yo, ass[,i], right_names)}))
}

full_design_same_degree <- function(Vo, Yo, right_names, max_degree) {
  lapply(1:max_degree, function(d){full_design_same_exact_degree(Vo, Yo, right_names, d)})
}

#Test new design matrix constructors equivalent to old when all three degrees are the same
#time_diffs <- c()
#for (r in 1:10){
#  cat(r); right <- (11-r):10
#  for (deg in 1:6){
#    cat(", ", deg, "\n")
#    t1 <- Sys.time()
#    des0 <- full_design1(d$V[,right,drop=F], d$Y[,right,drop=F], right_names[right], deg, deg, deg)
#    t2 <- Sys.time()
#    tmp <- full_design_same_degree(d$V[,right,drop=F], d$Y[,right,drop=F], right_names[right], deg)
#    t3 <- Sys.time()
#    des1 <- do.call("cbind", tmp)
#    # Since both functions produce the same column names
#    dif <- sum(abs(des1[,sort(colnames(des1))] - des0[,sort(colnames(des0))]))
#    if (abs(dif) > 1e-10){
#      cat("r=",r,", deg=",deg,": ", dif, "\n")
#      stop("")
#    }
#    time_diffs <- c(time_diffs, c(t2-t1-(t3-t2)))
#    print(time_diffs)
#}}
