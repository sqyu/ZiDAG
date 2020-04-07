.rhs_formula <- function(r, d){
  # "Vo[,1]+...+Vo[,r]+polym(Yo[,1],...,Yo[,r], degree=d)
  #1 old naive: return (list("V1"=paste(paste(paste("Vo[,",1:r,"]",sep=""), collapse=" + "), "+",
  #        paste("polym(", paste(paste("Yo[,",1:r,"]",sep=""), collapse=", "), ", raw=TRUE, degree=", d, ")", sep="")),
  #      "Y1"=paste(paste(paste("Vo[V1,",1:r,"]",sep=""), collapse=" + "), "+",
  #        paste("polym(", paste(paste("Yo[V1,",1:r,"]",sep=""), collapse=", "), ", raw=TRUE, degree=", d, ")", sep=""))
  #))
  #2 old shared design between V and Y, but no interaction on V: return (paste(paste(paste("Vo[,",right,"]",sep=""), collapse=" + "), "+",
  #              paste("polym(", paste(paste("Yo[,",1:r,"]",sep=""), collapse=", "), ", raw=TRUE, degree=", d, ")", sep="")))
  return (paste("0 + polym(", paste(paste("Yo[,",1:r,"]",sep=""), collapse=", "), ", raw=TRUE, degree=", d, ")", sep=""))
}

.design_V <- function(Vo, degree, right_names, only=FALSE){
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
      name <- c(name, paste("V", right_names, sep=""))
    } else{
      comb <- utils::combn(r,d)
      x <- cbind(x, apply(comb, 2, function(com){apply(Vo[,com,drop=FALSE],1,prod)}))
      name <- c(name, apply(comb, 2, function(com){paste(paste("V",right_names[com],sep=""),collapse="*")}))
    }
  }
  colnames(x) <- name
  return (x)
}

.split_sum <- function(n,k){
  ## Assumes: n >= 0, k > 0
  ## Returns: each column represents an ordered partition of n into k non-negative numbers
  if (n < 0 || k <= 0){
    stop("n must be non-negative and k must be positive.")
  }
  if (k > 1){
    return (apply(utils::combn(n+k-1,k-1), 2, function(tmp){c(tmp,n+k)-c(0,tmp)-1}))
  } else {
    return (matrix(n,1,1))
  }
}

.name_fun <- function(name, deg, prefix=""){#(1,2,"Y") -> "Y1^2"
  if(deg==1){
    paste(prefix,name,sep="")
  } else if (deg>1){
    paste(prefix,name,"^",deg,sep="")
  } else {""}
}

.design_Y_from_poly <- function(Yo, Y_degree, right_names){
  r <- ncol(Yo)
  if (r != length(right_names)) stop("ncol(Yo) == length(right_names) must be TRUE.")
  if (Y_degree == 1){
    colnames(Yo) <- paste("Y", right_names, sep="")
    return (Yo)
  }
  design_matrix_Y <- stats::model.matrix(stats::formula(paste(" ~", .rhs_formula(r, Y_degree))))
  if (ncol(design_matrix_Y) == 1) # If r == 1 and Y_degree == 1, model.matrix does not give the degree after ")"
    colnames(design_matrix_Y)[1] <- paste("raw = TRUE, degree = ", Y_degree, ")1", sep = "")
  #"polym(...)0.1.2.3.4"->c(0,1,2,3,4)
  design_digits <- matrix(sapply(colnames(design_matrix_Y), function(s){as.integer(strsplit(strsplit(s, paste("raw = TRUE, degree = ", Y_degree, ")", sep = ""))[[1]][2], "\\.")[[1]])}), nrow=r)
  #"polym(...)0.1.2.3.4"->c(2,3,4,5)
  design_which_included <- apply(design_digits, 2, function(d){which(d!=0)})
  #"polym(...)0.1.2.3.4"->"Y2*Y3^2*Y4^3*Y4^5"
  design_names <- sapply(1:length(design_which_included), function(di){paste(sapply(design_which_included[[di]],function(dd){.name_fun(right_names[dd],design_digits[dd,di],"Y")}), collapse="*")})
  colnames(design_matrix_Y) <- design_names
  return (design_matrix_Y)
}

.design_Y <- function(Yo, degree, right_names, only=FALSE){
  ## Assumes: degree >= 0. Ensured when calling from fit_pms_fromVY
  ## Returns: interactions of degree <= degree. NULL if degree == 0.
  ## Slower than stats::model.matrix(stats::formula(paste(" ~", .rhs_formula(r, degree))))
  if (degree == 0){return (NULL)}
  if (degree < 0){stop("Degree must be >= 0.")}
  r <- ncol(Yo)
  x <- NULL
  name <- NULL
  if (!only)
    return (.design_Y_from_poly(Yo, degree, right_names))
  #for (d in ifelse(only, degree, 1):degree){
  if (degree == 1){
    x <- Yo
    name <- c(name, paste("Y",right_names,sep=""))
  }
  else{
    comb <- .split_sum(degree, r)
    x <- cbind(x, apply(comb, 2, function(com){apply(sweep(Yo[,com!=0,drop=FALSE], 2, com[com!=0], "^"), 1, prod)}))
    name <- c(name, apply(comb, 2, function(com){gsub("\\^1", "", paste(paste("Y",right_names[com!=0],"^",com[com!=0],sep=""),collapse="*"))}))
  }
  #}
  colnames(x) <- name
  return (x)
}

.quad_design <- function(Vo, Yo, right_names){
  if (length(dim(Vo))!=2 || length(dim(Yo))!=2 || any(dim(Vo)!=dim(Yo)) || ncol(Vo)!=length(right_names))
    stop ("Vo and Yo must be matrices of the same dimension and same number of columns as right_names.")
  r <- ncol(Vo)
  VoYo <- cbind(Vo, Yo)
  colnames(VoYo) <- c(paste("V", right_names, sep=""), paste("Y", right_names, sep=""))
  res <- cbind(VoYo)
  names <- colnames(VoYo)
  res <- cbind(res, do.call("cbind", lapply(seq_len(r-1), function(ri){Vo[,ri]*Vo[,(ri+1):r]})))
  names <- c(names, do.call("c", sapply(seq_len(r-1), function(ri){paste("V",right_names[ri],"*","V",right_names[(ri+1):r],sep="")}, simplify = FALSE)))
  res <- cbind(res, Yo^2)
  names <- c(names, paste("Y", right_names, "^2", sep=""))
  res <- cbind(res, do.call("cbind", lapply(seq_len(r-1), function(ri){Yo[,ri]*Yo[,(ri+1):r]})))
  names <- c(names, do.call("c", sapply(seq_len(r-1), function(ri){paste("Y",right_names[ri],"*","Y",right_names[(ri+1):r],sep="")}, simplify = FALSE)))
  res <- cbind(res, do.call("cbind", lapply(1:r, function(ri){Vo[,ri]*Yo[,-ri]})))
  if (r > 1)
    names <- c(names, c(sapply(1:r, function(ri){paste("V",right_names[ri],"*","Y",right_names[-ri],sep="")})))
  colnames(res) <- names
  return (res)
}

#' Design matrix for pms parametrization given the max degree of the Hurdle polynomial.
#'
#' Design matrix for pms parametrization given the max degree of the Hurdle polynomial.
#'
#' @param Vo A numerical matrix of the same dimension as \code{Yo} indicating if each entry in \code{Yo} is non-zero, i.e. \code{Vo = (Yo != 0)}.
#' @param Yo A numerical matrix, the sample matrix for the parent nodes (regressors).
#' @param right A string vector of length equal to \code{ncol(Vo)}, the names of the parent nodes.
#' @param V_degree A non-negative integer, the degree for the \code{Vo} in the Hurdle polynomial. Defaults to 1.
#' @param Y_degree A non-negative integer, the degree for the \code{Yo} in the Hurdle polynomial. Defaults to 1.
#' @param Y_V_degree A non-negative integer, the degree for interaction between \code{Vo} and \code{Yo} in the Hurdle polynomial. Defaults to 1. If equal to 1, no interaction will be included (since it would be either a pure \code{V} term or a pure \code{Y} term).
#' @details
#' A Hurdle polynomial in parents \code{Yo} is a polynomial in \code{Yo} and their 0/1 indicators \code{Vo}.
#' The \code{V_degree} of a term that is a product of some columns of \code{Vo} only is the number of parents that appears in it. For example, \code{V1 * V2 * V3} has \code{V_degree} equal to 3. Note that \code{V1^p} is equal to \code{V1} for any \code{p >= 1} so it does not make sense to include a power.
#' The \code{Y_degree} of a term that is a product of powers of some columns of \code{Yo} only is the degree of a polynomial in its usual sense. For example, \code{Y1^2 * Y2 * Y3^3} has \code{Y_degree} equal to 2+1+3=6.
#' The \code{Y_V_degree} of a term that involves both some columns of \code{Vo} and some of \code{Yo} is the sum of the \code{V_degree} of the \code{V} part and the \code{Y_degree} of the \code{Y} part. For example, \code{Y1^2 * V2 * Y3^3 * V4 * V5} has \code{Y_V_degree} equal to 2+1+3+1+1=8.
#' The function returns a Hurdle polynomial including all possible terms with \code{V_degree}, \code{Y_degree}, \code{Y_V_degree} less than or equal to those specified.
#' For example, if \code{Vo} and \code{Yo} has two columns and \code{V_degree == 2}, \code{Y_degree == 2}, \code{Y_V_degree == 2}, the function returns a matrix with columns \code{V1}, \code{V2}, \code{V1*V2}, \code{Y1}, \code{Y2}, \code{Y1*Y2}, \code{Y1^2}, \code{Y2^2}, \code{Y1*V2}, \code{Y2*V1}.
#' Note that terms like \code{V1*Y1} are not included as it is equivalent to \code{Y1}.
#' Intercepts (a column with all 1s) are also not included.
#' @return The design matrix with number of rows equal to \code{nrow(Vo)}.
#' @examples
#' Yo <- matrix(rnorm(100) * rbinom(100, 1, 0.7), ncol=2)
#' full_design1(Yo != 0, Yo, c("1", "2"), V_degree=2, Y_degree=2, Y_V_degree=2)
#' full_design1(Yo != 0, Yo, c("1", "2"), V_degree=0, Y_degree=0, Y_V_degree=0)
#' @export
full_design1 <- function(Vo, Yo, right, V_degree=1, Y_degree=1, Y_V_degree=1){
  ## Y_degree >= 0, V_degree >= 0, Y_V_degree can be 0; if Y_degree and V_degree are both non-zero, Y_V_degree can be between 2 and Y_degree+V_degree
  ### Rules out V1*Y1 by checking names of variables using colnames.
  #rhs_form_Y <- .rhs_formula(r, Y_degree)
  #design_matrix_Y <- stats::model.matrix(formula(paste(" ~", rhs_form_Y)))
  #design_matrix_V <- .design_V(Vo = Vo, degree=degree, right_names=right)
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
      message("Y_V_degree=1 treated as no interaction between Y and V.")
    if (V_degree == 0 || Y_degree == 0)
      stop("Interaction between V and Y not allowed when V_degree or Y_degree is 0.")
    if (Y_V_degree > V_degree + Y_degree)
      stop("Y_V_degree can be V_degree + Y_degree at most.")
  }
  if (Y_degree == 0 || V_degree == 0 || Y_V_degree <= 1){
    design_matrix_V <- .design_V(Vo = Vo, degree=min(r,V_degree), right_names=right, only=FALSE)
    design_matrix_Y <- .design_Y(Yo = Yo, degree=Y_degree, right_names=right, only=FALSE)
    design_matrix_Y_V <- NULL
  } else if (V_degree == 2 && Y_degree == 2 && Y_V_degree == 2) {
    return (.quad_design(Vo, Yo, right))
  } else {
    ## effective degrees
    V_degree <- min(r, V_degree)
    Y_V_degree <- min(V_degree+Y_degree, Y_V_degree) # In fact, min(Y_V_degree, r+Y_degree, V_degree+Y_degree)
    ## !! To avoid things like W1*Y1, must have Y_V_degree <- min(min(r-1,V_degree)+Y_degree, Y_V_degree)
    design_matrix_V <- lapply(1:min(r,V_degree), function(d){.design_V(Vo = Vo, degree=d, right_names=right, only=TRUE)})
    design_matrix_Y <- lapply(1:Y_degree, function(d){.design_Y(Yo = Yo, degree=d, right_names=right, only=TRUE)})
    ### Interactions between Y and V, but does not rule out V1*Y1!!!
    design_matrix_Y_V <- do.call("cbind", lapply(2:Y_V_degree, function(td){
      do.call("cbind", lapply(max(td-Y_degree,1):min(td-1,V_degree), function(vd){ ## max(td-Y_degree,1) <= min(td-1,V_degree) is guaranteed
        do.call("cbind", lapply(1:ncol(design_matrix_V[[vd]]), function(vi){
          v_names <- .str_to_var(colnames(design_matrix_V[[vd]])[vi])
          good_ys <- sapply(colnames(design_matrix_Y[[td-vd]]), function(s){!length(intersect(.str_to_var(s), v_names))})
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



.num_par_V <- function(V_degree, r){sapply(1:min(V_degree,r),function(d){choose(r,d)})}
.num_par_Y <- function(Y_degree, r){sapply(1:Y_degree,function(d){choose(d+r-1,r-1)})}
.num_par_no_Y_V <- function(V_degree, Y_degree, r){(Y_degree+1)/r*choose(Y_degree+r,r-1)+sum(sapply(1:min(V_degree,r),function(i){choose(r,i)}))}
.num_par_per_regression <- function(V_degree, Y_degree, Y_V_degree, r, alias=FALSE){
  ### alias: If TRUE, calculate the naive number of parameters
  ### ; if FALSE, exclude terms like W1*Y1 and W1*W2*Y2^2
  if (Y_degree < 0 || V_degree < 0 || Y_V_degree < 0 ||
      Y_V_degree == 1 || (Y_V_degree > 0 && (V_degree == 0 || Y_degree == 0 ||
                                             Y_V_degree > V_degree + Y_degree))){stop("Invalid input.")}
  V_degree <- min(r, V_degree)
  Y_V_degree <- min(min(r-1,V_degree)+Y_degree, Y_V_degree)
  if (Y_V_degree == 0){return (.num_par_no_Y_V(V_degree, Y_degree, r))}
  if (alias){
    V_nums <- .num_par_V(V_degree, r)
    Y_nums <- .num_par_Y(Y_degree, r)
    return(1+sum(V_nums)+sum(Y_nums)+sum(sapply(2:Y_V_degree, function(td){sum(sapply(max(td-Y_degree,1):min(td-1,V_degree), function(vd){V_nums[vd]*Y_nums[td-vd]}))})))
  } else{
    V_nums <- .num_par_V(V_degree, r)
    Y_nums <- sapply(1:r, function(ri){.num_par_Y(Y_degree, ri)})
    if (is.null(dim(Y_nums))){Y_nums <- matrix(Y_nums, nrow=Y_degree, ncol=r)}
    return(1+sum(V_nums)+sum(Y_nums[,r])+sum(sapply(2:Y_V_degree, function(td){
      sum(sapply(max(td-Y_degree,1):min(td-1,r-1,V_degree), function(vd){V_nums[vd]*Y_nums[td-vd,r-vd]}))})))
  }
}

.str_to_var <- function(s){
  ## "W1^4*W2^2*W3" (or V or Y) to c("1","2","3")
  sapply(strsplit(s, "\\*")[[1]], function(ss){strsplit(strsplit(ss, "\\^")[[1]][1], "[VWY]")[[1]][2]})
}


### For producing designs with Y_degree == V_degree == Y_V_degree, separated by degree
.flip_1s <- function(vec) {
  # Given a vec with k elements equal to 1, returns the matrix of size
  # length(vec) x 2^k - 1, where each column is the result of flipping
  # the signs of at least one of the 1s
  pos_1s <- which(vec == 1)
  if (length(pos_1s) == 1) {
    vec[pos_1s] <- -1
    return (as.matrix(vec, ncol=1))
  }
  return (do.call("cbind", sapply(1:length(pos_1s), function(num_flips) {
    to_flips <- utils::combn(pos_1s, num_flips)
    apply(to_flips, 2, function(to_flip) {tmp<-vec; tmp[to_flip]<--1; tmp})
  })))
}

.split_sum_neg1 <- function(n, k){
  ## Assumes: n >= 0, k > 0
  ## Returns: each column represents an ordered sequence of k integers >= -1 whose abs sum equals to n
  if (n < 0 || k <= 0){
    stop("n must be non-negative and k must be positive.")
  }
  no_neg1 <- .split_sum(n, k)
  if (n == 1) # Must have exactly one 1 in each column
    return (cbind(no_neg1, -no_neg1))
  if (n == 0 || k == 1)  # n == 0 split into k integers OR n > 1 split into 1 integer -> no 1 in the result
    return (no_neg1)
  extras <- do.call("cbind", apply(no_neg1, 2, function(vec){
    if (1 %in% vec) .flip_1s(vec)
    else NULL
  }))
  return (cbind(no_neg1, extras))
}

.one_col <- function(Vo, Yo, vec, right_names){
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
    name <- c(paste("V", right_names[pos_vs], collapse="*", sep=""))
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

.full_design_same_exact_degree <- function(Vo, Yo, right_names, exact_degree) {
  if (exact_degree == 0){return (NULL)}
  if (exact_degree < 0){stop("Degree must be >= 0.")}
  r <- ncol(Vo)
  if (r != length(right_names) || r != ncol(Yo)) stop("ncol(Vo) == ncol(Yo) == length(right_names) must be TRUE.")
  x <- NULL
  name <- NULL
  ass <- .split_sum_neg1(exact_degree, r)  # -1 corresponds to Vo, >= 1 corresponds to the degree of Yo
  do.call("cbind", lapply(1:ncol(ass), function(i){.one_col(Vo, Yo, ass[,i], right_names)}))
}

# Returns a list of design matrices
# Each matrix corresponds to an exact degree
# Used for incrementally choosing the best degree using BIC in fit_pms_choose_degree()
.full_design_same_degree <- function(Vo, Yo, right_names, max_degree) {
  lapply(1:max_degree, function(d){.full_design_same_exact_degree(Vo, Yo, right_names, d)})
}

