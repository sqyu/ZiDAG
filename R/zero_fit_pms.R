.mle1d_pms <- function(V,Y,left,value_only=TRUE){
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
  nll <- -log_dhurdle1d_pms(V,Y,phat,muhat,sigmasqhat)
  if (value_only)
    return (nll)
  else{
    pars <- c(phat, muhat, sigmasqhat)
    names(pars) <- paste(c("p", "mu", "sigmasq"), "_V", left, sep="")
    return (list("nll"=nll, "par"=pars, "n"=n, "effective_df"=3))
  }
}


.make_folds <- function(y, nfits){
  id <- numeric(length(y))
  id[y==0] <- (sample.int(sum(y==0)) %% nfits) + 1
  id[y==1] <- (sample.int(sum(y==1)) %% nfits) + 1
  return (id)
}

# Helper function for zi_fit_lm
# Performs logistic regression if binomial == TRUE and linear regression otherwise
# If penalize == TRUE, use glmnet L2 penalty with cross validation or BIC (depending on CV_or_BIC)
# pen_factors: passed to glmnet (if intercept == FALSE, 0 is assigned for the intercept)
# nfits: number of folds if CV, number of lambdas if BIC
.zi_fit_lm_sub <- function(y, X, binomial, intercept, penalize,
                       pen_factors=NULL, tol=1e-8, maxit=100000,
                       nfits=10, CV_or_BIC=FALSE){
  if (!is.null(dim(y))) {stop("y must be a vector, not a matrix.")}
  if (length(dim(X)) != 2) {stop("X must be a 2-d matrix.")}
  if (typeof(X) == "logical")
    X <- 1*X # change to numerical; otherwise glm may fail
  n <- nrow(X); p <- ncol(X)
  if (n != length(y)) {stop("The length of y must be equal to the number of rows in X.")}
  if (!penalize) {
    if (binomial) {
        if (intercept)
          model <- stats::glm(y ~ X, family="binomial",
                       control=list("epsilon"=tol, "maxit"=maxit))
        else
          model <- stats::glm(y ~ X + 0, family="binomial",
                       control=list("epsilon"=tol, "maxit"=maxit))
      sigmasq <- NULL
    } else {
      if (intercept)
        model <- stats::lm(y ~ X)
      else
        model <- stats::lm(y ~ X + 0)
      yhat <- stats::predict(model)
      sigmasq <- sum((y-yhat)^2)/nrow(X)
    }
    names(model$coefficients) <- substr(names(model$coefficients), 2,
                                        sapply(names(model$coefficients),
                                               nchar))
    model$coefficients <- model$coefficients[!is.na(model$coefficients)]
    names(model$coefficients)[names(model$coefficients) == "Intercept)"] <- "1"
    return (list("coefficients"=model$coefficients,
                 "logLik"=as.numeric(stats::logLik(model)), "sigmasq"=sigmasq,
                 "df"=length(model$coefficients)+(!is.null(sigmasq))))
  } else { # Penalized regression
    family <- ifelse(binomial, "binomial", "gaussian")
    if (is.null(pen_factors))
      pen_factors <- rep(1, ncol(X))
    else if (length(pen_factors) != ncol(X))
      stop("In .zi_fit_lm_sub: length of pen_factors must be equal to the number of columns in X. Please report the bug.")
    if (intercept) {
      pen_factors <- c(0, pen_factors) # No penalty on intercept
      pen_factors <- pen_factors / sum(pen_factors) * (1 + ncol(X))
    } else {
      pen_factors <- pen_factors / sum(pen_factors) * ncol(X)
    }
    if (CV_or_BIC) {
      if (family == "binomial") {
        fraction <- table(y)/length(y)
        weights <- 1 - fraction[as.character(y)] # weights for 0-1
        foldid <- .make_folds(y, nfits) # to make folds evenly by taking into account the 0/1s
        thefit.cv <- glmnet::cv.glmnet(X, y, family=family,
                               intercept=intercept, nfits=nfits,
                               penalty.factor=pen_factors,
                               alpha=0, nlambda=100, thresh=tol,
                               maxit=maxit, standardize.response=FALSE,
                               foldid = foldid, weights = weights)
      } else  {
        thefit.cv <- glmnet::cv.glmnet(X, y, family=family,
                               intercept=intercept, nfits=nfits,
                               penalty.factor=pen_factors,
                               alpha=0, nlambda=100, thresh=tol,
                               maxit=maxit, standardize.response=FALSE)
      }
      thefit <- glmnet::glmnet(X, y, lambda=thefit.cv$lambda.min,
                       family=family, intercept=intercept,
                       penalty.factor=pen_factors,
                       alpha=0, nlambda=100, thresh=tol,
                       maxit=maxit, standardize.response=FALSE)
      yhat <- c(stats::predict(thefit, X, s = thefit.cv$lambda.min, type="response"))
      if (family == "binomial") {
        avgnlls <- c(-mean(y*log(yhat)+(1-y)*log(1-yhat)))
        sigmasqs <- NULL
      } else {
        sigmasqs <- c(mean((y-yhat)^2))
        avgnlls <- c(1/2+log(2*pi*sigmasqs)/2)
      }
      if (intercept) X <- cbind(1, X)
      if (all(pen_factors == 1)) {
        dd <- svd(X)$d
        dfs <- c(sum(dd^2/(dd^2 + thefit.cv$lambda.min), na.rm=TRUE) + (!is.null(sigmasqs)))
      } else { # Expensive
        dfs <- sum(diag(X %*% solve(crossprod(X) + thefit.cv$lambda.min * diag(pen_factors), t(X)))) + (!is.null(sigmasqs))
      }
      best_li <- 1
    } else {
      thefit <- glmnet::glmnet(X, y, family=family, intercept=intercept, alpha=0, nlambda=nfits, thresh=tol, maxit=maxit, standardize.response=FALSE)
      yhats <- stats::predict(thefit, X, type="response") # matrix of responses ncol(X) * nlambda
      if (family == "binomial") {
        avgnlls <- -colMeans(y*log(yhats)+(1-y)*log(1-yhats))
        sigmasqs <- NULL
      } else {
        sigmasqs <- colMeans((y-yhats)^2)
        avgnlls <- 1/2+log(2*pi*sigmasqs)/2
      }
      if (intercept) X <- cbind(1, X)
      if (all(pen_factors == 1)) {
        dd <- svd(X)$d
        dfs <- sapply(thefit$lambda, function(ll){sum(dd^2/(dd^2 + ll))}) + (!is.null(sigmasqs))
      } else { # Expensive
        XX <- crossprod(X)
        dfs <- sapply(thefit$lambda, function(ll){
          sum(diag(X %*% solve(XX + ll * diag(pen_factors), t(X))))}) + (!is.null(sigmasqs))
        remove(XX)
      }
      BICs <- bic_score(avgnlls, n, dfs)
      best_li <- which.min(BICs)
    }
    coefs <- stats::coef(thefit)[, best_li]
    names(coefs)[names(coefs) == "(Intercept)"] <- "1"
    if (!intercept) {coefs <- coefs[names(coefs) != "1"]} # No intercept
    return (list("coefficients"=coefs, "logLik"=as.numeric(-avgnlls[best_li]*n),
                 "sigmasq"=switch(is.null(sigmasqs)+1, sigmasqs[best_li], NULL),
                 "df"=dfs[best_li]))
  }
}

#' Performs unpenalized/ridge logistic/linear regression.
#'
#' Performs unpenalized/ridge logistic/linear regression.
#'
#' @param y The response vector.
#' @param X The design matrix.
#' @param binomial A logical. Fits a logistic regression if \code{TRUE}, otherwise a linear regression.
#' @param intercept A logical, whether to include an intercept in the regression. Defaults to \code{TRUE}.
#' @param tol A number, tolerance. Defaults to \code{1e-8}. Passed to \code{stats::glm()} for penalized logistic regressions, or as the \code{thresh} argument to \code{glmnet::glmnet()} for both logistic and linear regressions if penalized.
#' @param maxit An integer, the maximum number of iterations. Defaults to \code{100000}. Passed to \code{stats::glm()} for penalized logistic regressions, or to \code{glmnet::glmnet()} for both logistic and linear regressions if penalized.
#' @param seed The seed, passed to \code{set.seed()} before model fitting if not \code{NULL}. Defaults to \code{NULL}.
#' @param penalize_decider A logical or a function that takes \code{X} and returns a logical. Defaults to \code{function(X){ncol(X)>=nrow(X)/2}}. Used to decide whether to use penalized l2 (ridge) regression (if \code{TRUE}). Note that if the fits for the unpenalized regression are almost perfect, a penalized regression will be automatically used.
#' @param pen_factors A vector of non-negative numbers, defaults to \code{NULL}. Passed as the \code{penalty.factor} argument to \code{glmnet::glmnet()}. If \code{intercept == TRUE}, a \code{0} will be prepended.
#' @param nfits A positive integer, defaults to \code{10}. Used for penalized regressions, as number of folds if \code{CV_BIC == TRUE} (\code{nfits} argument to \code{glmnet::cv.glmnet()}, with \code{nlambda} set to \code{100}), or the number of lambdas if \code{BIC == FALSE} (as the \code{nlambda} argument to \code{glmnet::glmnet()}).
#' @param runs A positive integer, the number of reruns. The fit with the maximum likelihood will be returned. Defaults to \code{2}.
#' @param CV_or_BIC A logical, whether to use cross validation or BIC to choose the model from the path of penalized fits. Defaults to \code{FALSE}.
#' @details See description of each argument.
#' @return A list of following elements for the fitted model:
#'     \item{coefficients}{The coefficients on columns of \code{X}. The coefficient on the intercept is prepended if \code{intercept == TRUE}.}
#'     \item{logLik}{The log likelihood.}
#'     \item{sigmasq}{The estimated variance (sigma squared). \code{NULL} if \code{binomial == FALSE}.}
#'     \item{df}{The effective degree of freedom.}
#' @examples
#' m <- 5; n <- 100;
#' adj_mat <- ZiDAG::make_dag(m, mode = "chain")
#' d <- ZiDAG::gen_zero_dat(seed=1, gen_para="pms", adj_mat=adj_mat, n=n, gen_uniform_degree=1)
#' X <- full_design1(d$V[,2:m], d$Y[,2:m], right=paste(2:m),
#'          V_degree = 2, Y_degree = 2, Y_V_degree = 2)
#' zi_fit_lm(d$V[,1], X, TRUE)
#' zi_fit_lm(d$V[,1], X, penalize_decider=TRUE, TRUE, nfits=5, runs=1)
#' zi_fit_lm(d$Y[,1], X, FALSE)
#' zi_fit_lm(d$Y[,1], X, FALSE, penalize_decider=function(X){nrow(X)<100*ncol(X)}, nfits=5, runs=1)
#' @export
zi_fit_lm <- function(y, X, binomial, intercept=TRUE, tol=1e-8, maxit=100000, seed=NULL,
                   penalize_decider=function(X){ncol(X)>=nrow(X)/2},
                   pen_factors=NULL, nfits=10, runs=2, CV_or_BIC=FALSE){
  ## Calls .zi_fit_lm_sub for estimation
  penalize <- NULL
  if (is.logical(penalize_decider))
    penalize <- penalize_decider
  else if (is.function(penalize_decider))
    penalize <- penalize_decider(X)
  if (!is.logical(penalize))
    stop("penalize_decider must be a logical or function that takes X and returns a logical.")
  if (!is.null(seed)) {set.seed(seed)}
  if (!penalize) {
    if (binomial) { # If fit naive logistic regression
      res <- tryCatch(
        res <- .zi_fit_lm_sub(y, X, binomial=TRUE, intercept=intercept,
                          penalize=FALSE, pen_factors=NULL,
                          tol=tol, maxit=maxit, CV_or_BIC=CV_or_BIC),
                    warning = function(w){res <- w})
      if (!"message" %in% names(res)) # Otherwise, warnings like "fitted probabilities numerically 0 or 1 occurred" occurred, so fit l2-penalized logistic regression
        return (res)
    } else { # If fit naive linear regression
      res <- .zi_fit_lm_sub(y, X, binomial=FALSE, intercept=intercept,
                        penalize=FALSE, pen_factors=NULL,
                        tol=tol, maxit=maxit, CV_or_BIC=CV_or_BIC)
      if (res$sigmasq >= min(1e-4, tol*100)) # Otherwise, fit almost perfect, so fit l2-penalized liner regression
        return (res)
    }
  }
  best_res <- list(logLik = -Inf)
  for (run in seq_len(runs)) {
    # If fit penalized regression (including cases where the unpenalized regression was bad)
    res <- .zi_fit_lm_sub(y, X, binomial, intercept, penalize=TRUE,
                          pen_factors=pen_factors, tol=tol, maxit=maxit,
                          nfits=nfits, CV_or_BIC=CV_or_BIC)
    if (res$logLik > best_res$logLik) # since we MAXIMIZE the log likelihood
      best_res <- res
  }
  return (best_res)
}

# Combine results from the model for p and the model for mu and return
# the sum of their likelihoods, as well as the parameters if not value_only
.res_from_model_fits <- function(p_model, mu_model, n, left_name, value_only){
  sigmasq <- mu_model$sigmasq
  best_res <- list()
  best_res$nll <- -(p_model$logLik+mu_model$logLik)/n # Negative mean log likelihood
  if (value_only)
    return (best_res$nll)
  else{
    best_res$par <- c(p_model$coefficients, mu_model$coefficients, sigmasq)
    names(best_res$par) <- paste(c(paste("p", names(p_model$coefficients), sep="_"), paste("mu", names(mu_model$coefficients), sep="_"), "sigmasq"), "_V", left_name, sep="")
    best_res$n <- n
    best_res$effective_df <- p_model$df+mu_model$df
    return (best_res)
  }
}

.fit_pms_input_check <- function(V, Y, extra_regressors,
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

#' Fits a Hurdle conditional model with pms parametrization of specified degree.
#'
#' Fits a Hurdle conditional model with pms parametrization of specified degree.
#'
#' @param V A matrix of 0/1s, equal to Y != 0.
#' @param Y A data matrix of the same size as \code{V}.
#' @param left An integer between 1 and \code{ncol(Y)}. The index of the variable to be fit.
#' @param right A vector of integers between 1 and \code{ncol(Y)} different from \code{left}. Indices of the "regressors".
#' @param extra_regressors A matrix with the same number of rows as \code{V} and \code{Y}, extra regressors to be included in both regressions (conditional log odds/conditional mean). Defaults to \code{NULL}.
#' @param extra_reg_pen_factors A vector of non-negative numbers, defaults to \code{NULL}. Penalty factors for \code{extra_regressors}. If the main design matrix has \code{d} columns, \code{c(rep(1, d), extra_reg_pen_factors)} will be passed as the \code{penalty.factor} argument to \code{glmnet::glmnet()}. If \code{intercept == TRUE}, a \code{0} will also be prepended.
#' @param p_V_degree A non-negative integer, the degree for the \code{Vo} in the Hurdle polynomial for the conditional log odds. Defaults to 1.
#' @param p_Y_degree A non-negative integer, the degree for the \code{Yo} in the Hurdle polynomial for the conditional log odds. Defaults to 1.
#' @param p_Y_V_degree A non-negative integer, the degree for interaction between \code{Vo} and \code{Yo} in the Hurdle polynomial for the conditional log odds. Defaults to 1. If equal to 1, no interaction will be included (since it would be either a pure \code{V} term or a pure \code{Y} term).
#' @param mu_V_degree A non-negative integer, the degree for the \code{Vo} in the Hurdle polynomial for the conditional mean. Defaults to 1.
#' @param mu_Y_degree A non-negative integer, the degree for the \code{Yo} in the Hurdle polynomial for the conditional mean. Defaults to 1.
#' @param mu_Y_V_degree A non-negative integer, the degree for interaction between \code{Vo} and \code{Yo} in the Hurdle polynomial for the conditional mean. Defaults to 1. If equal to 1, no interaction will be included (since it would be either a pure \code{V} term or a pure \code{Y} term).
#' @param value_only If \code{TRUE}, returns the minimized negative log likelihood only. Defaults to \code{TRUE}.
#' @param tol A number, tolerance. Defaults to \code{1e-8}. Passed to \code{stats::glm()} for penalized logistic regressions, or as the \code{thresh} argument to \code{glmnet::glmnet()} for both logistic and linear regressions if penalized.
#' @param maxit An integer, the maximum number of iterations. Defaults to \code{100000}. Passed to \code{stats::glm()} for penalized logistic regressions, or to \code{glmnet::glmnet()} for both logistic and linear regressions if penalized.
#' @param seed A number, the random seed passed to \code{zi_fit_lm()} for both regressions (conditional log odds/conditional mean).
#' @param penalize_decider A logical or a function that takes a design matrix and returns a logical. Defaults to \code{function(X){ncol(X)>=nrow(X)/2}}. Used to decide whether to use penalized l2 (ridge) regression (if \code{TRUE}) when fitting each conditional distribution. Note that for either regression (conditional log odds/conditional mean), if the fits for unpenalized regressions are almost perfect, penalized regressions will be automatically used.
#' @param nfits A positive integer, defaults to \code{10}. Used for penalized regressions, as number of folds if \code{CV_BIC == TRUE} (\code{nfits} argument to \code{glmnet::cv.glmnet()}, with \code{nlambda} set to \code{100}), or the number of lambdas if \code{BIC == FALSE} (as the \code{nlambda} argument to \code{glmnet::glmnet()}).
#' @param runs A positive integer, the number of reruns. The fit with the maximum likelihood will be returned. Defaults to \code{2}.
#' @details
#' A Hurdle conditional model with pms parametrization for the \code{left} node given those in \code{right} has log density with respect to the sum of the Lebesgue measure and a point mass at 0 equal to (in terms of \code{y})
#' \eqn{\log(1-p)}{log(1-p)} if \code{y == 0}, or \eqn{\log(p)-(y-mu)^2/2/sigmasq}{log(p)-(y-mu)^2/2/sigmasq} otherwise. That is, it is a mixture of a binomial with probability of success \code{p} and a Gaussian with conditional mean \code{mu} and conditional variance \code{sigmasq}.
#' Here \code{sigmasq} is assumed constant, and parameters \code{log(p/(1-p))} and \code{mu} are Hurdle polynomials, i.e. polynomials in the values for \code{right} and their indicators.
#' This function thus fits such a model using \code{Y[,left]}, \code{Y[,right]} and \code{V[,right] = (Y[,right] != 0)}, using a logistic for the log odds \code{log(p/(1-p))} and a linear regression for \code{mu}.
#'
#' Writing \code{Yo <- Y[,right]}, a Hurdle polynomial in parents \code{Yo} is a polynomial in \code{Yo} and their 0/1 indicators \code{Vo}.
#' The \code{V_degree} of a term that is a product of some columns of \code{Vo} only is the number of parents that appears in it. For example, \code{V1 * V2 * V3} has \code{V_degree} equal to 3. Note that \code{V1^p} is equal to \code{V1} for any \code{p >= 1} so it does not make sense to include a power.
#' The \code{Y_degree} of a term that is a product of powers of some columns of \code{Yo} only is the degree of a polynomial in its usual sense. For example, \code{Y1^2 * Y2 * Y3^3} has \code{Y_degree} equal to 2+1+3=6.
#' The \code{Y_V_degree} of a term that involves both some columns of \code{Vo} and some of \code{Yo} is the sum of the \code{V_degree} of the \code{V} part and the \code{Y_degree} of the \code{Y} part. For example, \code{Y1^2 * V2 * Y3^3 * V4 * V5} has \code{Y_V_degree} equal to 2+1+3+1+1=8.
#' The design matrix thus includes all possible terms with \code{V_degree}, \code{Y_degree}, \code{Y_V_degree} less than or equal to those specified.
#' For example, if \code{Vo} and \code{Yo} has two columns and \code{V_degree == 2}, \code{Y_degree == 2}, \code{Y_V_degree == 2}, the design matrix has columns \code{V1}, \code{V2}, \code{V1*V2}, \code{Y1}, \code{Y2}, \code{Y1*Y2}, \code{Y1^2}, \code{Y2^2}, \code{Y1*V2}, \code{Y2*V1}. Note that terms like \code{V1*Y1} are not included as it is equivalent to \code{Y1}.
#' Parameters \code{p_V_degree}, \code{p_Y_degree}, \code{p_Y_V_degree}, \code{mu_V_degree}, \code{mu_Y_degree}, and \code{mu_Y_V_degree} specify these degrees for the regressions for the log odds \code{log(p/(1-p))} and the conditional mean \code{mu}, respectively.
#'
#' For automatically choosing a uniform degree <= a specified maximum degree, please use \code{zi_fit_pms_choose_degree()}.
#' @return If \code{value_only == TRUE}, returns the minimized negative log likelihood only. Otherwise, returns
#'   \item{nll}{A number, the minimized negative log likelihood.}
#'   \item{par}{A vector of length \code{4*length(right)+3}, the fitted parameters, in the other of: the intercept for the \code{a} (a scalar), linear coefficients on \code{V[,right]} for \code{a}, linear coefficients on \code{Y[,right]} for \code{a}, the intercept for the \code{b} (a scalar), linear coefficients on \code{V[,right]} for \code{b}, linear coefficients on \code{Y[,right]} for \code{b}.}
#'   \item{n}{An integer, the sample size.}
#'   \item{effective_df}{\code{4*length(right)+3}, the effective degree of freedom.}
#' @examples
#' m <- 3; n <- 1000
#' adj_mat <- make_dag(m, "complete")
#' dat <- gen_zero_dat(1, "pms", adj_mat, n, k_mode=1, min_num=10, gen_uniform_degree=1)
#' extra_regressors <- matrix(rnorm(n * 4), nrow=n)
#' extra_reg_pen_factors <- c(1, 2, 3, 4) / sum(c(1, 2, 3, 4))
#' zi_fit_pms(dat$V, dat$Y, 3, 1:2, extra_regressors=extra_regressors,
#'     extra_reg_pen_factors=extra_reg_pen_factors, p_V_degree=2, p_Y_degree=2,
#'     p_Y_V_degree=2, mu_V_degree=2, mu_Y_degree=2, mu_Y_V_degree=2, value_only=TRUE)
#' zi_fit_pms(dat$V, dat$Y, 3, 1:2, extra_regressors=extra_regressors,
#'     extra_reg_pen_factors=extra_reg_pen_factors, p_V_degree=2, p_Y_degree=2,
#'     p_Y_V_degree=2, mu_V_degree=2, mu_Y_degree=2, mu_Y_V_degree=2, value_only=FALSE)
#' @export
zi_fit_pms <- function(V, Y, left, right, extra_regressors=NULL,
                    extra_reg_pen_factors=NULL, p_V_degree=1, p_Y_degree=1,
                    p_Y_V_degree=1, mu_V_degree=1, mu_Y_degree=1, mu_Y_V_degree=1,
                    value_only=TRUE, tol=1e-8, maxit=100000, seed=NULL,
                    penalize_decider=function(X){ncol(X)>=nrow(X)/2},
                    nfits=10, runs=2){
  ## Y_degree >= 0, V_degree >= 0, Y_V_degree can be 0; if Y_degree and V_degree are both non-zero, Y_V_degree can be between 2 and Y_degree+V_degree
  ### Temporarily rules out V1*Y1 by checking names of variables using colnames.
  ## If p degrees and mu degrees match, equivalent to old fit_pms_fromVY(V, Y, left, right, V_degree, Y_degree, Y_V_degree, value_only, seed)
  ## If p_V_degree==p_Y_degree==1 and mu_V_degree==mu_Y_degree==mu_Y_V_degree==2, equivalent to old fit_pms_fromVY_quad(V, Y, left, right, p_quad=FALSE, mu_quad=TRUE, value_only, seed); similarly for other cases
  extra_reg_pen_factors <- .fit_pms_input_check(V, Y, extra_regressors,
                                               extra_reg_pen_factors, "zi_fit_pms")
  r <- length(right)
  if (r == 0){
    return (.mle1d_pms(V, Y, left, value_only))
  }
  V1 <- V[, left]; Y1 <- Y[, left]
  Vo <- V[, right, drop=FALSE]; Yo <- Y[, right, drop=FALSE]
  p_design_matrix <- full_design1(Vo, Yo, right, p_V_degree, p_Y_degree, p_Y_V_degree)
  p_pen_factors <- c(rep(1, ncol(p_design_matrix)), extra_reg_pen_factors)
  p_design_matrix <- cbind(p_design_matrix, extra_regressors)
  if (p_V_degree == mu_V_degree && p_Y_degree == mu_Y_degree && p_Y_V_degree == mu_Y_V_degree) {
    mu_design_matrix <- p_design_matrix[V1,,drop=FALSE]
    mu_pen_factors <- p_pen_factors
  }
  else {
    mu_design_matrix <- full_design1(Vo[V1,,drop=FALSE], Yo[V1,,drop=FALSE], right, mu_V_degree, mu_Y_degree, mu_Y_V_degree)
    mu_pen_factors <- c(rep(1, ncol(mu_design_matrix)), extra_reg_pen_factors)
    mu_design_matrix <- cbind(mu_design_matrix, extra_regressors[V1,,drop=FALSE])
  }
  p_model <- zi_fit_lm(y = V1, X = p_design_matrix, binomial = TRUE,
                    intercept = TRUE, tol=tol, maxit=maxit,
                    seed=seed, penalize_decider=penalize_decider,
                    pen_factors=p_pen_factors,
                    nfits=nfits, runs=runs, CV_or_BIC=FALSE)
  mu_model <- zi_fit_lm(y = Y1[V1==TRUE], X = mu_design_matrix,
                     binomial = FALSE, intercept = TRUE, tol=tol,
                     maxit=maxit, seed=seed, penalize_decider=penalize_decider,
                     pen_factors=mu_pen_factors,
                     nfits=nfits, runs=runs, CV_or_BIC=FALSE)
  return (.res_from_model_fits(p_model, mu_model, nrow(V), left, value_only))
}

#' Fits and chooses a Hurdle conditional model with pms parametrization of degree <= a maximum degree.
#'
#' Fits and chooses a Hurdle conditional model with pms parametrization of degree <= a maximum degree.
#'
#' @param V A matrix of 0/1s, equal to Y != 0.
#' @param Y A data matrix of the same size as \code{V}.
#' @param left An integer between 1 and \code{ncol(Y)}. The index of the variable to be fit.
#' @param right A vector of integers between 1 and \code{ncol(Y)} different from \code{left}. Indices of the "regressors".
#' @param max_uniform_degree A positive integer, the maximum degree for the Hurdle polynomials.
#' @param extra_regressors A matrix with the same number of rows as \code{V} and \code{Y}, extra regressors to be included in both regressions (conditional log odds/conditional mean). Defaults to \code{NULL}.
#' @param extra_reg_pen_factors A vector of non-negative numbers, defaults to \code{NULL}. Penalty factors for \code{extra_regressors}. If the main design matrix has \code{d} columns, \code{c(rep(1, d), extra_reg_pen_factors)} will be passed as the \code{penalty.factor} argument to \code{glmnet::glmnet()}. If \code{intercept == TRUE}, a \code{0} will also be prepended.
#' @param value_only If \code{TRUE}, returns the minimized negative log likelihood only. Defaults to \code{TRUE}.
#' @param tol A number, tolerance. Defaults to \code{1e-8}. Passed to \code{stats::glm()} for penalized logistic regressions, or as the \code{thresh} argument to \code{glmnet::glmnet()} for both logistic and linear regressions if penalized.
#' @param maxit An integer, the maximum number of iterations. Defaults to \code{100000}. Passed to \code{stats::glm()} for penalized logistic regressions, or to \code{glmnet::glmnet()} for both logistic and linear regressions if penalized.
#' @param seed A number, the random seed passed to \code{zi_fit_lm()} for both regressions (conditional log odds/conditional mean).
#' @param penalize_decider A logical or a function that takes a design matrix and returns a logical. Defaults to \code{function(X){ncol(X)>=nrow(X)/2}}. Used to decide whether to use penalized l2 (ridge) regression (if \code{TRUE}) when fitting each conditional distribution. Note that for either regression (conditional log odds/conditional mean), if the fits for unpenalized regressions are almost perfect, penalized regressions will be automatically used.
#' @param nfits A positive integer, defaults to \code{10}. Used for penalized regressions, as number of folds if \code{CV_BIC == TRUE} (\code{nfits} argument to \code{glmnet::cv.glmnet()}, with \code{nlambda} set to \code{100}), or the number of lambdas if \code{BIC == FALSE} (as the \code{nlambda} argument to \code{glmnet::glmnet()}).
#' @param runs A positive integer, the number of reruns. The fit with the maximum likelihood will be returned. Defaults to \code{2}.
#' @param print_best_degree A logical, whether to print the degree (1, ..., \code{max_uniform_degree}) that minimizes the BIC.
#' @details
#' A Hurdle conditional model with pms parametrization for the \code{left} node given those in \code{right} has log density with respect to the sum of the Lebesgue measure and a point mass at 0 equal to (in terms of \code{y})
#' \eqn{\log(1-p)}{log(1-p)} if \code{y == 0}, or \eqn{\log(p)-(y-mu)^2/2/sigmasq}{log(p)-(y-mu)^2/2/sigmasq} otherwise. That is, it is a mixture of a binomial with probability of success \code{p} and a Gaussian with conditional mean \code{mu} and conditional variance \code{sigmasq}.
#' Here \code{sigmasq} is assumed constant, and parameters \code{log(p/(1-p))} and \code{mu} are Hurdle polynomials, i.e. polynomials in the values for \code{right} and their indicators.
#' This function thus fits such a model using \code{Y[,left]}, \code{Y[,right]} and \code{V[,right] = (Y[,right] != 0)}, using a logistic for the log odds \code{log(p/(1-p))} and a linear regression for \code{mu}.
#'
#' Writing \code{Yo <- Y[,right]}, a Hurdle polynomial in parents \code{Yo} is a polynomial in \code{Yo} and their 0/1 indicators \code{Vo}.
#' The degree of a term in a Hurdle polynomial is the number of \code{V} terms plus the sum of the degrees of the \code{Y} terms. For example, \code{Y1^2 * V2 * Y3^3 * V4 * V5} has degree equal to 2+1+3+1+1=8.
#' Given a degree, the design matrix thus includes all possible terms with degree less than or equal to the specified degree.
#' For example, if \code{Vo} and \code{Yo} has two columns and if we choose degree \code{2}, the design matrix has columns \code{V1}, \code{V2}, \code{V1*V2}, \code{Y1}, \code{Y2}, \code{Y1*Y2}, \code{Y1^2}, \code{Y2^2}, \code{Y1*V2}, \code{Y2*V1}. Note that terms like \code{V1*Y1} are not included as it is equivalent to \code{Y1}.
#'
#' This function fits models using Hurdle polynomials with degrees 1, 2, ..., \code{max_uniform_degree}, and automatically chooses the degree that minimizes the BIC.
#' It is equivalent to calling \code{zi_fit_pms()} with all degree arguments equal to \code{d}, with \code{d} in 1, ..., \code{max_uniform_degree}, and returning the one with the smallest BIC.
#' @return If \code{value_only == TRUE}, returns the minimized negative log likelihood only. Otherwise, returns
#'   \item{nll}{A number, the minimized negative log likelihood.}
#'   \item{par}{A vector of length \code{4*length(right)+3}, the fitted parameters, in the other of: the intercept for the \code{a} (a scalar), linear coefficients on \code{V[,right]} for \code{a}, linear coefficients on \code{Y[,right]} for \code{a}, the intercept for the \code{b} (a scalar), linear coefficients on \code{V[,right]} for \code{b}, linear coefficients on \code{Y[,right]} for \code{b}.}
#'   \item{n}{An integer, the sample size.}
#'   \item{effective_df}{\code{4*length(right)+3}, the effective degree of freedom.}
#' @examples
#' m <- 3; n <- 1000
#' adj_mat <- make_dag(m, "complete")
#' dat <- gen_zero_dat(1, "pms", adj_mat, n, k_mode=1, min_num=10, gen_uniform_degree=1)
#' extra_regressors <- matrix(rnorm(n * 4), nrow=n)
#' extra_reg_pen_factors <- c(1, 2, 3, 4) / sum(c(1, 2, 3, 4))
#' zi_fit_pms_choose_degree(dat$V, dat$Y, 3, 1:2, max_uniform_degree=2L,
#'     extra_regressors=extra_regressors, extra_reg_pen_factors=extra_reg_pen_factors,
#'     value_only=TRUE, print_best_degree=TRUE)
#' zi_fit_pms_choose_degree(dat$V, dat$Y, 3, 1:2, max_uniform_degree=2L,
#'     extra_regressors=extra_regressors, extra_reg_pen_factors=extra_reg_pen_factors,
#'     value_only=FALSE, print_best_degree=TRUE)
#' @export
zi_fit_pms_choose_degree <- function(V, Y, left, right, max_uniform_degree,
                                  extra_regressors=NULL,
                                  extra_reg_pen_factors=NULL,
                                  value_only=TRUE, tol=1e-8, maxit=100000,
                                  seed=NULL, penalize_decider=function(X){ncol(X)>=nrow(X)/2},
                                  nfits=10, runs=2, print_best_degree=FALSE){
  # Forces using BIC for model selection
  # Forces Y_degree == V_degree == Y_V_degree for both p and mu
  # Try Y_degree from 1 to max_uniform_degree and pick the one with the best BIC
  if (typeof(max_uniform_degree) != "integer" || max_uniform_degree < 1)
    stop("max_uniform_degree must be a positive integer.")
  extra_reg_pen_factors <- .fit_pms_input_check(V, Y, extra_regressors, extra_reg_pen_factors, "zi_fit_pms_choose_degree")
  if (max_uniform_degree > 2)
    warning("Using a max_uniform_degree >= 3 may result in very slow execution.")
  r <- length(right)
  if (r == 0){
    return (.mle1d_pms(V, Y, left, value_only))
  }
  V1 <- V[,left]; Y1 <- Y[,left]
  Vo <- V[,right,drop=FALSE]; Yo <- Y[,right,drop=FALSE]
  design_matrices_by_degree <- .full_design_same_degree(Vo, Yo, right, max_uniform_degree)
  design_matrix <- NULL
  best_bic <- Inf
  best_deg <- 0
  best_res <- list()
  for (deg in 1:max_uniform_degree) {
    design_matrix <- cbind(design_matrix, design_matrices_by_degree[[deg]])
    pen_factors <- c(rep(1, ncol(design_matrix)), extra_reg_pen_factors)
    design_matrix_with_extra <- cbind(design_matrix, extra_regressors)
    p_model <- zi_fit_lm(y = V1, X = design_matrix_with_extra,
                      binomial = TRUE, intercept = TRUE,
                      tol=tol, maxit=maxit, seed=seed,
                      penalize_decider=penalize_decider,
                      pen_factors=pen_factors,
                      nfits=nfits, runs=runs, CV_or_BIC=FALSE) # Use BIC since it is used for choosing the degree anyways
    mu_model <- zi_fit_lm(y = Y1[V1==TRUE], X = design_matrix_with_extra[V1,,drop=FALSE],
                       binomial = FALSE, intercept = TRUE,
                       tol=tol, maxit=maxit, seed=seed,
                       penalize_decider=penalize_decider,
                       pen_factors=pen_factors,
                       nfits=nfits, runs=runs, CV_or_BIC=FALSE)
    fit_res <- .res_from_model_fits(p_model, mu_model, nrow(V), left, value_only=FALSE)
    res_bic <- bic_from_fit(fit_res)
    if (best_bic > res_bic) {
      best_deg <- deg
      best_bic <- res_bic
      best_res <- fit_res
    }
  }
  if (print_best_degree)
    message("For left=", left, ", right=", paste(right, collapse=","), ", best degree is ", best_deg, "\n", sep="")
  if (value_only)
    return (best_res$nll)
  return (best_res)
}
