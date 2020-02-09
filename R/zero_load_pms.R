source("zero_load.R")
source("design_mat.R")
require("glmnet")

generate_p_mu_sigmasq <- function(n, p, mu, sigmasq){
  rbinom(n, 1, p) * rnorm(n, mu, sqrt(sigmasq))
}


logprob1d_pms <- function(V,Y,p,mu,sigmasq){ # 1-d hurdle
  # V, Y: vectors
  V_on <- mean(V)
  return (log(1-p)*(1-V_on)+log(p)*V_on-
            sum((Y[V]-mu)^2)/2/sigmasq/length(V)
          -log(2*pi*sigmasq)/2*V_on)
}


mle1d_pms <- function(V,Y,left,value_only=T){
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
  ll <- -logprob1d_pms(V,Y,phat,muhat,sigmasqhat)
  if (value_only)
    return (ll)
  else{
    pars <- c(phat, muhat, sigmasqhat)
    names(pars) <- paste(c("p", "mu", "sigmasq"), "_V", left, sep="")
    return (list("value"=ll, "par"=pars, "n"=n, "effective_df"=3))
  }
}

make_folds <- function(y, nfits){
  id <- numeric(length(y))
  id[y==0] <- (sample.int(sum(y==0)) %% nfits) + 1
  id[y==1] <- (sample.int(sum(y==1)) %% nfits) + 1
  return (id)
}

fit_lm_sub <- function(y, X, binomial, intercept, penalize, 
                       pen_factors=NULL, tol=1e-8, maxit=100000,
                       seed=NULL, nfits=10, CV_or_BIC=FALSE){
  # nfits: nfolds if CV, nlambda if BIC
  if (!is.null(dim(y))) {stop("y must be a vector, not a matrix.")}
  if (length(dim(X)) != 2) {stop("X must be a 2-d matrix.")}
  if (typeof(X) == "logical")
    X <- 1*X # change to numerical; otherwise glm may fail
  n <- nrow(X); p <- ncol(X)
  if (n != length(y)) {stop("The length of y must be equal to the number of rows in X.")}
  if (!penalize) {
    if (binomial) {
        if (intercept)
          model <- glm(y ~ X, family="binomial", 
                       control=list("epsilon"=tol, "maxit"=maxit))
        else
          model <- glm(y ~ X + 0, family="binomial", 
                       control=list("epsilon"=tol, "maxit"=maxit))
      sigmasq <- NULL
    } else {
      if (intercept)
        model <- lm(y ~ X)
      else
        model <- lm(y ~ X + 0)
      yhat <- predict(model)
      sigmasq <- sum((y-yhat)^2)/nrow(X)
    }
    names(model$coefficients) <- substr(names(model$coefficients), 2, 
                                        sapply(names(model$coefficients), 
                                               nchar))
    model$coefficients <- model$coefficients[!is.na(model$coefficients)]
    names(model$coefficients)[names(model$coefficients) == "Intercept)"] <- "1"
    return (list("coefficients"=model$coefficients, 
                 "logLik"=logLik(model), "sigmasq"=sigmasq, 
                 "df"=length(model$coefficients)+(!is.null(sigmasq))))
  } else { # Penalized regression
    family <- ifelse(binomial, "binomial", "gaussian")
    if (!is.null(pen_factors) && length(pen_factors) != ncol(X))
      stop("In fit_lm_sub: length of pen_factors must be equal to the number of columns in X. Please report the bug.")
    if (intercept) {
      dd <- svd(cbind(1, X))$d
      pen_factors <- c(0, pen_factors) # No penalty on intercept
      pen_factors <- pen_factors / sum(pen_factors) * (1 + ncol(X))
    } else {
      dd <- svd(X)$d
      pen_factors <- pen_factors / sum(pen_factors) * ncol(X)
    }
    if (!is.null(seed)) {set.seed(seed)} # CV folds are random
    if (CV_or_BIC) {
      if (family == "binomial") {
        fraction <- table(y)/length(y)
        weights <- 1 - fraction[as.character(y)] # weights for 0-1
        foldid <- make_folds(y, nfits) # to make folds evenly by taking into account the 0/1s
        thefit.cv <- cv.glmnet(X, y, family=family, 
                               intercept=intercept, nfits=nfits, 
                               penalty.factor=pen_factors,
                               alpha=0, nlambda=100, thresh=tol, 
                               maxit=maxit, standardize.response=FALSE, 
                               foldid = foldid, weights = weights)
      } else  {
        thefit.cv <- cv.glmnet(X, y, family=family, 
                               intercept=intercept, nfits=nfits,
                               penalty.factor=pen_factors,
                               alpha=0, nlambda=100, thresh=tol, 
                               maxit=maxit, standardize.response=FALSE)
      }
      thefit <- glmnet(X, y, lambda=thefit.cv$lambda.min, 
                       family=family, intercept=intercept, 
                       penalty.factor=pen_factors,
                       alpha=0, nlambda=100, thresh=tol, 
                       maxit=maxit, standardize.response=FALSE)
      yhat <- c(predict(thefit, X, 
                        s = thefit.cv$lambda.min, type="response"))
      if (family == "binomial") {
        avgnlls <- c(-mean(y*log(yhat)+(1-y)*log(1-yhat)))
        sigmasqs <- NULL
      } else {
        sigmasqs <- c(mean((y-yhat)^2))
        avgnlls <- c(1/2+log(2*pi*sigmasqs)/2)
      }
      dfs <- c(sum(dd^2/(dd^2 + thefit.cv$lambda.min * pen_factors), na.rm=T) + 
                 (!is.null(sigmasqs)))
      best_li <- 1
    } else {
      thefit <- glmnet(X, y, family=family, intercept=intercept, alpha=0, nlambda=nfits, thresh=tol, maxit=maxit, standardize.response=FALSE)
      yhats <- predict(thefit, X, type="response") # matrix of responses ncol(X) * nlambda
      if (family == "binomial") {
        avgnlls <- -colMeans(y*log(yhats)+(1-y)*log(1-yhats))
        sigmasqs <- NULL
      } else {
        sigmasqs <- colMeans((y-yhats)^2)
        avgnlls <- 1/2+log(2*pi*sigmasqs)/2
      }
      # na.rm = T since the extra regressors might be colinear (e.g. used 
      # specified design matrix for factors without knowing if intercept 
      # is used in the actual fitting) and thus X may have zero eigenvalues
      dfs <- sapply(1:length(thefit$lambda), function(li){
            sum(dd^2/(dd^2 + thefit$lambda[li] * pen_factors), na.rm=T)}) + 
            (!is.null(sigmasqs))
      BICs <- bic_score(avgnlls, n, dfs)
      best_li <- which.min(BICs)
    }
    coefs <- coef(thefit)[, best_li]
    names(coefs)[names(coefs) == "(Intercept)"] <- "1"
    if (!intercept) {coefs <- coefs[names(coefs) != "1"]} # No intercept
    return (list("coefficients"=coefs, "logLik"=-avgnlls[best_li]*n,
                 "sigmasq"=switch(is.null(sigmasqs)+1, sigmasqs[best_li], NULL),
                 "df"=dfs[best_li]))
  }
}

fit_lm <- function(y, X, binomial, intercept=TRUE, tol=1e-8, 
                   maxit=100000, seed=NULL, 
                   penalize_decider=function(X){ncol(X)>=nrow(X)/2}, 
                   pen_factors=NULL, nfits=10, runs=2, 
                   CV_or_BIC=FALSE){
  ## Calls fit_lm_sub for estimation
  penalize <- NULL
  if (is.logical(penalize_decider))
    penalize <- penalize_decider
  else if (is.function(penalize_decider))
    penalize <- penalize_decider(X)
  if (!is.logical(penalize))
    stop("penalize_decider must be a bool or function from X to a bool.")
  if (!penalize) {
    if (binomial) { # If fit naive logistic regression
      res <- tryCatch(
        res <- fit_lm_sub(y, X, binomial=TRUE, intercept=intercept, 
                          penalize=FALSE, pen_factors=NULL,
                          tol=tol, maxit=maxit, CV_or_BIC=CV_or_BIC),
                    warning = function(w){res <- w})
      if ("message" %in% names(res)) # If "fitted probabilities numerically 0 or 1 occurred", fit l2-penalized logistic regression
        penalize <- TRUE
    } else { # If fit naive linear regression
      res <- fit_lm_sub(y, X, binomial=FALSE, intercept=intercept, 
                        penalize=FALSE, pen_factors=NULL,
                        tol=tol, maxit=maxit, CV_or_BIC=CV_or_BIC)
      if (res$sigmasq < min(1e-4,tol*100)) # If fit almost perfect, fit l2-penalized liner regression
        penalize <- TRUE
    }
  }
  if (penalize) { # If fitting penalized regression OR tbe unpenalized regression was bad
    best_res <- fit_lm_sub(y, X, binomial, intercept, penalize=TRUE, 
                           pen_factors=pen_factors, tol=tol, maxit=maxit, 
                           seed=seed, nfits=nfits, CV_or_BIC=CV_or_BIC)
    if (!is.null(seed) || runs==1) # If seed provided, only run once (as the results will be the same for all runs)
      return (best_res)
    for (run in 2:runs){
      res <- fit_lm_sub(y, X, binomial, intercept, penalize=TRUE, 
                        pen_factors=pen_factors, tol=tol, maxit=maxit, 
                        seed=seed, nfits=nfits, CV_or_BIC=CV_or_BIC)
      if (res$logLik > best_res$logLik) # since we MAXIMIZE the log likelihood
        best_res <- res
    }
    return (best_res)
  } else {
    return (res)
  }
}

res_from_model_fits <- function(p_model, mu_model, n, left_name, value_only){
  sigmasq <- mu_model$sigmasq
  best_res <- list()
  best_res$value <- -(p_model$logLik+mu_model$logLik)/n # Negative mean log likelihood
  if (value_only)
    return (best_res$value)
  else{
    best_res$par <- c(p_model$coefficients, mu_model$coefficients, sigmasq)
    names(best_res$par) <- paste(c(paste("p", names(p_model$coefficients), sep="_"), paste("mu", names(mu_model$coefficients), sep="_"), "sigmasq"), "_V", left_name, sep="")
    best_res$n <- n
    best_res$effective_df <- p_model$df+mu_model$df
    return (best_res)
  }
}

fit_pms_input_check <- function(V, Y, extra_regressors,
                                extra_reg_pen_factors, func_name) {
### Checks if inputs are valid and returns updated 
### extra_reg_pen_factors: if NULL, set to 0 of length ncol(extra_regressors);
### else if of length 1, repeated; else program stops if not of the right size.
  if (any(dim(V) != dim(Y)))
    stop("In ", func_name, ": V and Y must have the same dimension.")
  if (!is.null(extra_regressors)) {
    if (!is.matrix(extra_regressors) || nrow(V) != nrow(extra_regressors))
      stop("In ", func_name, ": extra_regressors must be a matrix with the same number of rows as V and Y.")
    if (is.null(extra_reg_pen_factors)) {
      warning("In ", func_name, ": no extra_reg_pen_factors given. Set to 0 for all covariates in extra_regressors.")
      return (rep(0, ncol(extra_regressors)))
    }
    if (any(extra_reg_pen_factors < 0))
      stop("In ", func_name, ": extra_reg_pen_factors must be all non-negative.")
    if (length(extra_reg_pen_factors) == 1)
      return (rep(extra_reg_pen_factors, ncol(extra_regressors)))
    if (length(extra_reg_pen_factors) != ncol(extra_regressors))
      stop("In ", func_name, ": extra_reg_pen_factors must be NULL (in which case it will be set to 0), of length 1, or of length equal to the number of columns in extra_regressors.")
    return (extra_reg_pen_factors)
  }
}

fit_pms <- function(V, Y, left, right, extra_regressors=NULL, 
                    extra_reg_pen_factors=NULL,
                    p_V_degree=1, p_Y_degree=1, p_Y_V_degree=1, 
                    mu_V_degree=1, mu_Y_degree=1, mu_Y_V_degree=1, 
                    value_only=T, tol=1e-8, maxit=100000, seed=NULL, 
                    penalize_decider=function(X){ncol(X)>=nrow(X)/2}, 
                    nfits=10, runs=2, CV_or_BIC=FALSE){
  ## Y_degree >= 0, V_degree >= 0, Y_V_degree can be 0; if Y_degree and V_degree are both non-zero, Y_V_degree can be between 2 and Y_degree+V_degree
  ### Temporarily rules out V1*Y1 by checking names of variables using colnames.
  ## If p degrees and mu degrees match, equivalent to old fit_pms_fromVY(V, Y, left, right, V_degree, Y_degree, Y_V_degree, value_only, seed)
  ## If p_V_degree==p_Y_degree==1 and mu_V_degree==mu_Y_degree==mu_Y_V_degree==2, equivalent to old fit_pms_fromVY_quad(V, Y, left, right, p_quad=FALSE, mu_quad=TRUE, value_only, seed); similarly for other cases
  extra_reg_pen_factors <- fit_pms_input_check(V, Y, extra_regressors, 
                                               extra_reg_pen_factors, "fit_pms")
  r <- length(right)
  if (r == 0){
    return (mle1d_pms(V, Y, left, value_only))
  }
  V1 <- V[, left]; Y1 <- Y[, left]
  Vo <- V[, right, drop=F]; Yo <- Y[, right, drop=F]
  p_design_matrix <- full_design1(Vo, Yo, right, p_V_degree, p_Y_degree, p_Y_V_degree)
  p_pen_factors <- c(rep(1, ncol(p_design_matrix)), extra_reg_pen_factors)
  p_design_matrix <- cbind(p_design_matrix, extra_regressors)
  if (p_V_degree == mu_V_degree && p_Y_degree == mu_Y_degree && p_Y_V_degree == mu_Y_V_degree) {
    mu_design_matrix <- p_design_matrix[V1,,drop=F]
    mu_pen_factors <- p_pen_factors
  }
  else {
    mu_design_matrix <- full_design1(Vo[V1,,drop=F], Yo[V1,,drop=F], right, mu_V_degree, mu_Y_degree, mu_Y_V_degree)
    mu_pen_factors <- c(rep(1, ncol(mu_design_matrix)), extra_reg_pen_factors)
    mu_design_matrix <- cbind(mu_design_matrix, extra_regressors[V1,,drop=F])
  }
  p_model <- fit_lm(y = V1, X = p_design_matrix, binomial = TRUE, 
                    intercept = TRUE, tol=tol, maxit=maxit, 
                    seed=seed, penalize_decider=penalize_decider,
                    pen_factors=p_pen_factors,
                    nfits=nfits, runs=runs, CV_or_BIC=CV_or_BIC)
  mu_model <- fit_lm(y = Y1[V1==TRUE], X = mu_design_matrix, 
                     binomial = FALSE, intercept = TRUE, tol=tol,
                     maxit=maxit, seed=seed, penalize_decider=penalize_decider,
                     pen_factors=mu_pen_factors,
                     nfits=nfits, runs=runs, CV_or_BIC=CV_or_BIC)
  return (res_from_model_fits(p_model, mu_model, nrow(V), left, value_only))
}

fit_pms_choose_degree <- function(V, Y, left, right, max_uniform_degree, 
                                  extra_regressors=NULL, 
                                  extra_reg_pen_factors=NULL,
                                  value_only=T, tol=1e-8, maxit=100000, 
                                  seed=NULL, penalize_decider=function(X){ncol(X)>=nrow(X)/2}, 
                                  nfits=10, runs=2, print_best_degree=FALSE){
  # Forces using BIC for model selection
  # Forces Y_degree == V_degree == Y_V_degree for both p and mu
  # Try Y_degree from 1 to max_uniform_degree and pick the one with the best BIC
  if (typeof(max_uniform_degree) != "integer" || max_uniform_degree < 1)
    stop("max_uniform_degree must be a positive integer.")
  extra_reg_pen_factors <- fit_pms_input_check(V, Y, extra_regressors, extra_reg_pen_factors, "fit_pms_choose_degree")  
  if (max_uniform_degree > 2)
    warning("Using a max_uniform_degree >= 3 may result in very slow execution.")
  r <- length(right)
  if (r == 0){
    return (mle1d_pms(V, Y, left, value_only))
  }
  V1 <- V[,left]; Y1 <- Y[,left]
  Vo <- V[,right,drop=F]; Yo <- Y[,right,drop=F]
  design_matrices_by_degree <- full_design_same_degree(Vo, Yo, right, max_uniform_degree)
  design_matrix <- NULL
  best_bic <- Inf
  best_deg <- 0
  best_res <- list()
  for (deg in 1:max_uniform_degree) {
    design_matrix <- cbind(design_matrix, design_matrices_by_degree[[deg]])
    pen_factors <- c(rep(1, ncol(design_matrix)), extra_reg_pen_factors)
    design_matrix_with_extra <- cbind(design_matrix, extra_regressors)
    p_model <- fit_lm(y = V1, X = design_matrix_with_extra, 
                      binomial = TRUE, intercept = TRUE, 
                      tol=tol, maxit=maxit, seed=seed, 
                      penalize_decider=penalize_decider,
                      pen_factors=pen_factors,
                      nfits=nfits, runs=runs, CV_or_BIC=FALSE)
    mu_model <- fit_lm(y = Y1[V1==TRUE], X = design_matrix_with_extra[V1,,drop=FALSE], 
                       binomial = FALSE, intercept = TRUE, 
                       tol=tol, maxit=maxit, seed=seed, 
                       penalize_decider=penalize_decider, 
                       pen_factors=pen_factors,
                       nfits=nfits, runs=runs, CV_or_BIC=FALSE)
    fit_res <- res_from_model_fits(p_model, mu_model, nrow(V), left, value_only=FALSE)
    res_bic <- bic_from_fit(fit_res)
    if (best_bic > res_bic) {
      best_deg <- deg
      best_bic <- res_bic
      best_res <- fit_res
    }
  }
  if (print_best_degree)
    cat("For left=", left, ", right=", paste(right, collapse=","), ", best degree is ", best_deg, "\n", sep="")
  if (value_only)
    return (best_res$value)
  return (best_res)
}
