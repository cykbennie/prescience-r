#' @title Approximate Best Subset Maximum Binary Prediction Rule (PRESCIENCE)
#'
#' @description
#' \code{select} employs the Approximate Best Subset Maximum Binary Prediction Rule
#' (PRESCIENCE) on the given variable selection problem.
#'
#' @param formula an object of class formula, relating the dependent variable to the grouping variable.
#' @param data a data frame containing the variables in the model.
#' @param nfoc integer. The number of focus variable(s) (excluding the intercept).
#' @param q integer. The cardinality constraint for the covariate selection.
#' @param bound numeric. The maximum absolute value of the bounds for all variables.
#' @param beta0 integer. The coefficient taking value either 1 or -1 to normalize the scale for the first focus variable. Default = 1.
#' @param warmstart logical. If \code{TRUE}, use the warm start strategy.
#' @param tau the tuning parameter for enlarging the estimated bounds. Default = 1.5.
#' @param mio integer. 1 for MIO method 1 and 2 for method 2 in the paper. Default = 1.
#' @param tlim time limit (in seconds) specified for the MIO solver. Default = 86400.
#' @import stats
#' @import slam
#' @import gurobi
#' @export
#' @return a list with 7 elements:
#' \item{tolerance}{tolerance level}
#' \item{status}{optimization status}
#' \item{score}{Gurobi score}
#' \item{gap}{the MIO optimization gap value in case of early termination (0 if optimal solution is found within the time limit)}
#' \item{rtime}{time used by the MIO solver in the estimation procedure}
#' \item{ncount}{the number of Branch-and-bound nodes}
#' \item{bhat}{maximum score estimates for coefficients}
#' @author Yankang (Bennie) Chen <yankang.chen@@columbia.edu>
#' @references \emph{Best Subset Binary Prediction} by Le-Yu Chen and Sokbae Lee (2018).
#' \url{https://arxiv.org/abs/1610.02738}
#' @examples
#' results <- select(auto ~ dcost + cars + dovtt + divtt,
#' data = transportation, nfoc = 1, q = 1, bound = 10)
#'
#' summary(results)
#' coef(results)


select <- function(formula, data, nfoc, q, bound, beta0 = 1, warmstart = TRUE, tau = 1.5, mio = 1, tlim = 86400){
  # Error checking
  if (missing(formula) || class(formula) != "formula")
    stop("'formula' is missing or incorrect.")

  requireNamespace("gurobi")
  gc()

  df <- model.frame(formula, data)
  varnames <- colnames(df)
  varnames <- varnames[2:length(varnames)]
  varnames <- c(varnames[1:nfoc], '(Intercept)', varnames[(nfoc+1):length(varnames)])

  Y_tr <- as.matrix(df[1])
  temp <- df[2:ncol(df)]
  n_tr <- nrow(Y_tr)
  x_std <- scale(temp)
  x_foc <- cbind(as.matrix(x_std[,1:nfoc]), matrix(rep(1, n_tr), ncol = 1))  # [focus variables, intercept]
  x_aux <- x_std[,(nfoc+1):ncol(x_std)]  # [auxiliary variables]

  remove(df, temp, x_std)
  gc()

  d <- ncol(x_aux)
  bnd <- cbind(matrix(rep(-bound, nfoc+d), ncol = 1), matrix(rep(bound, nfoc+d), ncol = 1))  # set the initial parameter bounds
  tol <- floor(sqrt(log(n_tr) * n_tr) / 2)  # set the tolerance level value

  # Results
  results <- list()
  results$tolerance <- tol/n_tr

  if (warmstart) { # warm start MIO
    max_score <- warm_start_max_score(Y_tr, x_foc, x_aux, beta0, q, tlim, tol, bnd, mio, tau)
    results$status <- max_score$status
    bhat <- c(1, max_score$bhat)
    results$score <- max_score$score/n_tr
    results$gap <- max_score$gap/n_tr
    results$rtime <- max_score$rtime
    results$ncount <- max_score$ncount
  } else { # cold start MIO
    max_score <- max_score_constr_fn(Y_tr, x_foc, x_aux, beta0, q, tlim, tol, bnd, mio)
    results$status <- max_score$status
    bhat <- c(1, max_score$bhat)
    results$score <- max_score$score/n_tr
    results$gap <- max_score$gap/n_tr
    results$rtime <- max_score$rtime
    results$ncount <- max_score$ncount
  }

  bhat <- as.data.frame(bhat)
  rownames(bhat) <- varnames
  colnames(bhat) <- "Estimate"
  results$bhat <- bhat

  class(results) <- c("select", "list")
  return(results)
}


miobnd_fn <- function(x, beta0, bnd){

    n <- nrow(x)
    k <- ncol(x) - 1

    # Build model
    model <- list()
    model$modelsense <- "max"
    model$sense <- ">"
    model$lb <- bnd[,1]
    model$ub <- bnd[,2]

    # Set parameters
    params <- list()
    tol <- 1e-6
    params$outputflag <- 0
    params$OptimalityTol <- tol
    params$FeasibilityTol <- tol
    params$IntFeasTol <- tol

    v <- rep(0, 2)
    value <- rep(0, n)

    for (i in 1:n) {

        alpha <- beta0 * x[i,1]

        model$obj <- x[i,2:(k+1)]
        model$objcon <- alpha

        model$A <- as.simple_triplet_matrix(matrix(x[i,2:(k+1)], nrow = 1))
        model$rhs <- -alpha
        try(result <- gurobi(model, params))
        try(v[1] <- result$objval, silent = TRUE)

        model$obj <- -x[i,2:(k+1)]
        model$objcon <- -alpha

        model$A <- as.simple_triplet_matrix(matrix(-x[i,2:(k+1)], nrow = 1))
        model$rhs <- alpha
        try(result <- gurobi(model, params))
        try(v[2] <- result$objval, silent = TRUE)

        value[i] <- max(v)
    }

    return(value)
}


get_bnd <- function(y, x, beta0, bnd){

    k <- ncol(x) - 1
    n <- k+1
    logit <- matrix(coefficients(glm.fit(x, y, family = binomial(link = "logit"), intercept = FALSE)), ncol = 1)
    p_hat <- 1/(1+exp(-x %*% logit))

    constr <- matrix(rep(p_hat-0.5, n), ncol = n) * x
    bound <- bnd

    remove(logit, p_hat)
    gc()

    # Build model
    model <- list()
    model$sense <- ">"
    model$A <- as.simple_triplet_matrix(constr[,2:n])
    model$rhs <- -constr[,1]*beta0

    model$lb = bound[,1]
    model$ub = bound[,2]

    # Set parameters
    params <- list()
    tol <- 1e-6
    params$outputflag <- 0
    params$OptimalityTol <- tol
    params$FeasibilityTol <- tol
    params$IntFeasTol <- tol

    for (i in 1:k) {

        objcoef <- rep(0, k)
        objcoef[i] <- 1
        model$obj <- objcoef
        model$modelsense <- "min"

        try(result <- gurobi(model, params))
        try(bound[i,1] <- result$objval, silent = TRUE)

        model$modelsense <- "max"
        model$lb <- bound[,1]

        try(result <- gurobi(model, params))
        try(bound[i,2] <- result$objval, silent = TRUE)

        model$ub = bound[,2]
    }

    return(bound)
}


max_score_constr_fn <- function(y, x_foc, x_aux, beta0, q, tlim, abgap, bnd, mio){

    N <- nrow(y)
    k <- ncol(x_foc) - 1
    d <- ncol(x_aux)

    bhat <- matrix(rep(0, k+d), ncol = 1)  # (k+d) vector
    gap <- 0
    rtime <- 0

    miobnd <- miobnd_fn(cbind(x_foc, x_aux), beta0, bnd)  # (N) vector

    # Build model
    model <- list()
    model$sense <- "<"
    model$modelsense <- "max"
    model$lb <- c(rep(0, N), bnd[,1], rep(0, d))  # (N+k+2*d) vector
    model$ub <- c(rep(1, N), bnd[,2], rep(1, d))  # (N+k+2*d) vector
    model$vtype <- c(rep("B", N), rep("C", k+d), rep("B", d))  # (N+k+2*d) vector

    # Set parameters
    params <- list()
    tol <- 1e-6
    params$outputflag <- 0
    params$OptimalityTol <- tol
    params$FeasibilityTol <- tol
    params$IntFeasTol <- tol

    if (tlim > 0) params$TimeLimit <- tlim
    if (abgap > 0) params$MIPGapAbs <- abgap

    ztemp1 <- matrix(rep(0, N*d), ncol = d)  # N*d matrix
    ztemp2 <- matrix(rep(0, N*(2*d+1)), ncol = N)  # (2*d+1)*N matrix
    htemp <- rbind(diag(d), -diag(d), matrix(rep(0, d), ncol = d))  # (2*d+1)*d matrix
    etemp <- rbind(-diag(bnd[(k+1):(k+d), 2]), diag(bnd[(k+1):(k+d), 1]), matrix(rep(1, d), ncol = d))  # (2*d+1)*d matrix
    mtemp1 <- cbind(ztemp2, matrix(rep(0, k*(2*d+1)), ncol = k), htemp, etemp)  # (2*d+1)*(N+k+2*d) matrix
    mtemp2 <- c(rep(0, 2*d), q)  # (2*d+1) vector

    if (mio == 1) {
        model$obj <- c(2*y-1, rep(0, k+2*d))  # Method 1 formulation
        model$objcon <- sum(1-y)
        miobnd_bar <- miobnd+tol  # (N) vector
        mtemp3 <- cbind(diag(-miobnd_bar), x_foc[,2:(k+1)], x_aux, ztemp1)  # N*(N+k+2*d) matrix
        mtemp4 <- cbind(diag(miobnd), -x_foc[,2:(k+1)], -x_aux, ztemp1)  # N*(N+k+2*d) matrix
        model$A <- as.simple_triplet_matrix(rbind(mtemp4, mtemp3, mtemp1))  # (2*N+2*d+1)*(N+k+2*d) matrix
        model$rhs <- c(miobnd*(1-tol)+beta0*x_foc[,1], -tol*miobnd_bar-beta0*x_foc[,1], mtemp2)  # (2*N+2*d+1) vector
    } else {
        model$obj <- c(rep(1, N), rep(0, k+2*d))  # Method 2 formulation
        temp2 <- 1-2*y  # (N) vector
        temp3 <- matrix(rep(temp2, k), ncol = k)  # N*k matrix
        temp4 <- matrix(rep(temp2, d), ncol = d)  # N*d matrix
        model$A <- as.simple_triplet_matrix(rbind(cbind(diag(miobnd), temp3*x_foc[,2:(k+1)], temp4*x_aux, ztemp1), mtemp1))
        # (N+2*d+1)*(N+k+2*d) sparse matrix
        model$rhs <- c(miobnd*(1-tol) - beta0*temp2*x_foc[,1], mtemp2)  # (N+2*d+1) vector
    }

    remove(ztemp1, ztemp2, htemp, etemp, mtemp1, mtemp2)
    gc()

    results <- list()
    try(result <- gurobi(model, params))
    try(results$status <- result$status, silent = TRUE)
    try(results$bhat <- result$x[(N+1):(N+k+d)], silent = TRUE)
    try(results$score <- result$objval, silent = TRUE)
    try(results$gap <- result$objbound - results$score, silent = TRUE)
    try(results$rtime <- result$runtime, silent = TRUE)
    try(results$ncount <- result$nodecount, silent = TRUE)

    return(results)
}


warm_start_max_score <- function(y, x_foc, x_aux, beta0, q, tlim, tol, bnd, mio, tau){

    bnd_h <- get_bnd(y, cbind(x_foc, x_aux), beta0, bnd)  # (k-1+d)*2 matrix
    bnd_abs <- tau * apply(abs(bnd_h), 1, max)  # (k-1+d) vector
    bnd0 <- cbind(apply(cbind(-bnd_abs, bnd[,1]), 1, max), apply(cbind(bnd_abs, bnd[,2]), 1, min))  # (k-1+d)*2 matrix
    # this is the refind bound used for warm start MIO

    max_score <- max_score_constr_fn(y, x_foc, x_aux, beta0, q, tlim, tol, bnd0, mio)
    results <- list()
    results$status <- max_score$status
    results$bhat <- max_score$bhat
    results$score <- max_score$score
    results$gap <- max_score$gap
    results$rtime <- max_score$rtime
    results$ncount <- max_score$ncount

    return(results)
}
