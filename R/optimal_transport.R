#' solve the Optimal transport problem
#' @description Optimally transport according to wasserstein cost
#' C_{ij} = d(QX_i, Y_j)^2.  this is the usual optimal transport problem
#' except with an allowance of a parameter Q for an orthogonal matrix.
#' Give two matrices of dimensions n x d and m x y respectively
#' @param X an n x d matrix of points where each row is a point
#' @param Y similar, except possibly a different number of points
#' @param Q An orthogonal matrix
#' @param lambda the parameter to penalize in sinkhorn divergence
#' @param eps the tolerance
#' @return P, the n x m matrix of assignments
#' @export
optimal_transport <- function(X,Y, Q= NULL,lambda = .1,eps = .01) {
  if(is.null(Q)) {
    d <- dim(X)[2]
    Q <- diag(1,d,d)
  }
  X <- X%*%Q

  n <- dim(X)[1]
  m <- dim(Y)[1]

  tmp1 <- X%*%t(Y)
  tmp2 <- outer(rep(1, n), rowSums(Y^2))
  tmp3 <- outer(rowSums(X^2), rep(1,m))
  C <- tmp2 - 2*tmp1 + tmp3
  #C <- as.matrix(rect.dist(X,Y))^2


  r <- 1/n
  c <- 1/m
  P <- exp(-lambda * C)
  u <- rep(0,n)
  while (max(abs(u - rowSums(P))) > eps) {
    u = rowSums(P)
    P <- r*P / u
    v <- colSums(P)
    P <- c*t( t(P)/v)
  }

  return(P)

}

#' Procrustes
#' @description Function to find the optimal Q.  Classic Procrustes, but with the allowance
#' of an assignment matrix.
#' @param X the vectors X
#' @param Y the vectors Y
#' @param Pi the matrix of transports
#' @param p the dimension of the positive component
#' @param q the dimension of the negative component
#' @import Matrix
#' @export
#' @return the procrustes solution
procrustes <- function(X,Y, Pi = NULL,p = dim(X)[2],q=0) {
  if (is.null(Pi)) {
    stop("need to run optimal transport first")
  }
  if(q == 0) {
    vals <- svd(t(X)%*% Pi %*% Y)
    toReturn <- vals$u %*% t(vals$v)
  } else {
    Xp <- X[,c(1:p)]
    Yp <- Y[,c(1:p)]
    Xq <- X[,c((p+1):(p+q))]
    Yq <- Y[,c((p+1):(p+q))]
    valsp <- svd(t(Xp)%*% Pi %*% Yp)
    valsp <- valsp$u %*% t(valsp$v)
    valsq <- svd(t(Xq)%*% Pi %*% Yq)
    valsq <- valsq$u %*% t(valsq$v)
    toReturn <- Matrix::bdiag(valsp,valsq)
  }


  return(toReturn)

}

#' Iterate the optimal transport problem with penalization for different lambda.
#' @description Function to iterate optimal transport based on sinkhorn divergence, for a fixed
#' penalization parameter
#' @param X the n x d matrix of vectors
#' @param Y the m x d matrix of vectors
#' @param Q optional, an initialization point
#' @param lambda the penalization parameter
#' @param eps tolerance for computing sinkhorn divergence
#' @param numReps when to stop
#' @param eps_OT tolerance for the individual optimal transport problem
#' @param p the dimension of the positive component
#' @param q the dimension of the negative component
#' @return a list of the final Pi and Q
#' @export
#' @examples
#' library(rstiefel)
#' set.seed(2019)
#' X <- matrix(rnorm(1000,1,.2),ncol= 4)
#' Y <- rbind(X,X)
#' W <- rustiefel(4,4)
#' Y <- Y %*% W
#' test <- iterative_optimal_transport(X,Y,numReps = 1000,lambda = .0001)
#' norm(test$Q - W,"2")
#' X <- matrix(rnorm(5000,.2,.02),ncol= 5)
#' Y <- rbind(X,X)
#' W <- rustiefel(5,5)
#' Y <- Y %*% W
#' Y <- matrix(rnorm(200,.7),ncol =5)
#' test2 <- iterative_optimal_transport(X,Y,numReps = 1000,lambda = .0001)
#' norm(test2$Q - W,"2")
iterative_optimal_transport <-function(X,Y, Q = NULL,
                                        lambda = .01,eps = .01,
                                        numReps =1000,eps_OT = .01

                                        ,p = dim(X)[2],q=0) {
  if(is.null(Q)) {
    d <- dim(X)[2]
    Q <- diag(1,d,d)
  }
  Pi <- optimal_transport(X,Y,Q,eps = eps_OT,lambda=lambda)
  Q <- procrustes(X,Y,Pi,p,q)
  c <- norm(X %*% Q - Pi%*% Y,"F")
  i <- 1
  while ( i < numReps) {
    if( c > eps ) {

      Pi <- optimal_transport(X,Y,Q,eps = eps_OT,lambda = lambda)
      Q <- procrustes(X,Y,Pi,p,q)
      c <- norm(X %*% Q - Pi%*% Y,"F")
      i <- i+1
    } else {
      break
    }

  }
  toReturn <- list(Pi,Q,c)
  names(toReturn) <- c("Pi", "Q","obj.value")

  return(toReturn)

}

#' Match datasets
#' @description Function to iterate over a decreasing sequence of lambdas, the penalization parameters.
#' If lambda is big, the function is more concave, so we iterate, starting
#' from Lambda = .5, and decreasing each time by alpha.
#' This method takes longer, but is more likely to converge to the true solution,
#' since we start from a more concave problem and iteratively solve it by setting
#' lambda = alpha * lambda, for alpha in (0,1).
#' @param X an n x d dataset of vectors
#' @param Y an m x d dataset of vectors
#' @param Q an initial guess
#' @param lambda_init the initial value of lambda for penalization
#' @param lambda_final For termination
#' @param p the dimension of the positive component
#' @param q the dimension of the negative component
#' @param alpha the parameter for which lambda is multiplied by
#' @param eps the tolerance for the iterative optimal transport problem
#' @param numReps the number of reps for each subiteration
#' @param eps_OT the tolerance for the individual optimal transport problem
#' @return a list of the final orthogonal matrix and the assignment matrix
#' @export
#' @import rstiefel
#' @examples
#' library(rstiefel)
#'set.seed(2030)
#'X <- matrix(rnorm(100,1,.2),ncol= 4)
#' Y <- rbind(X,X)
#' W <- rustiefel(4,4)
#' Y <- Y %*% W
#' test <- match_support(X,Y,numReps = 10,Q = diag(c(1,-1,-1,1)))
#' #cheating a bit by starting with Q with the diagonals equal to
#' # the sign of the true matrix
#' # others have pointed out that initializing Q at all 2^d sign matrices ( diagonal matrices
#'  # whose entries are plus or minus one) might have better global
#'  # convergence
#'set.seed(2010)
#' X <- matrix(rnorm(100,.2,.02),ncol= 5)
#' Y <- rbind(X,X)
#' W <- rustiefel(5,5)
#' Y <- Y %*% W
#' test2 <- match_support(X,Y)
#' norm(test2$Q - W)
#'
match_support <- function(X,Y,
                           Q = NULL,lambda_init = .5, lambda_final = .01,alpha = .95,
                           eps = .01,numReps =100,eps_OT = .01,p = dim(X)[2],q=0) {

  lambda <- lambda_init
  if (is.null(Q)) {
    d <- dim(X)[2]
    Q <- diag(1,d,d)
  }
  while(lambda > lambda_final) {
    Q <- iterative_optimal_transport(X,Y,Q,lambda = lambda,eps = eps,numReps = numReps,eps_OT = eps_OT)
    Pi <- Q$Pi
    c <- Q$`obj.value`
    Q <- Q$Q
    lambda <- alpha*lambda
  }

  toReturn <- list(Pi,Q,c)
  #toReturn <- Q
  names(toReturn) <- c("Pi", "Q","obj.value")
  return(toReturn)

}

#' minimize over all possible orthogonal matrices
#'
#'#' Match datasets
#' @description find the orthogonal matrix that minimizes match_support starting
#' from all 2^d sign matrix initializations
#' @param X an n x d dataset of vectors
#' @param Y an m x d dataset of vectors
#' @param Q an initial guess
#' @param lambda_init the initial value of lambda for penalization
#' @param lambda_final For termination
#' @param p the dimension of the positive component
#' @param q the dimension of the negative component
#' @param alpha the parameter for which lambda is multiplied by
#' @param eps the tolerance for the iterative optimal transport problem
#' @param numReps the number of reps for each subiteration
#' @param eps_OT the tolerance for the individual optimal transport problem
#' @param costType one of "kernel", "obj.value", or "both", to determine the
#' minimum cost
#' @return a list of the final orthogonal matrix and the assignment matrix
#' @export
#' @import rstiefel
match_support_min <-function(X,Y,Q = NULL,lambda_init = .5, lambda_final = .01,
                             alpha = .95,eps = .01,numReps =100,eps_OT = .01,
                             p = dim(X)[2],q=0,costType = "kernel") {

  d <- dim(X)[2]

  #get all the sign matrices:
  ds <- list()
  for ( i in 1:d) {
    ds[[i]] <- c(-1,1)
  }
  signs <- expand.grid(ds)
  costs_kernel <- rep(0,d)
  costs_obj <- rep(0,d)
  get_matched <- list()

  #iterate over all sign matrices using match_support():
  for ( i in c(1:nrow(signs))) {
    print(paste("trying sign",i," of ", 2^d))
    currentsign <- diag(signs[i,])
    get_matched[[i]] <- match_support(X = X,Y = Y,lambda_init = lambda_init,
                                        lambda_fina= lambda_final,eps = eps,
                                        numReps = numReps,eps_OT = eps_OT,
                                        alpha = alpha,Q = currentsign,
                                        p=p,q=q)
    costs_obj[i] <-  get_matched[[i]]$obj.value
    costs_kernel[i] <- kernel.stat(Xhat%*% get_matched[[i]]$Q,Yhat)
  }

  #find the smallest cost value and return the corresponding result
  if (costType == "kernel") {
    minval <- which.min(costs_kernel)
    final <- get_matched[[minval]]
  } else {
    minval <- which.min(costs_obj)
    final <- get_matched[[minval]]
  }

  return(final)

}

#' minimize over all possible orthogonal matrices
#'
#'#' Match datasets
#' @description find the orthogonal matrix according to the initialization via
#' the median sign flip heuristic
#' @param X an n x d dataset of vectors
#' @param Y an m x d dataset of vectors
#' @param Q an initial guess
#' @param lambda_init the initial value of lambda for penalization
#' @param lambda_final For termination
#' @param p the dimension of the positive component
#' @param q the dimension of the negative component
#' @param alpha the parameter for which lambda is multiplied by
#' @param eps the tolerance for the iterative optimal transport problem
#' @param numReps the number of reps for each subiteration
#' @param eps_OT the tolerance for the individual optimal transport problem
#' @param costType one of "kernel", "obj.value", or "both", to determine the
#' minimum cost
#' @return a list of the final orthogonal matrix and the assignment matrix
#' @export
#' @import rstiefel
match_support_med <-function(X,Y,Q = NULL,lambda_init = .5, lambda_final = .01,
                             alpha = .95,eps = .01,numReps =100,eps_OT = .01,
                             p = dim(X)[2],q=0,costType = "kernel") {

  d <- dim(X)[2]

  meds1 <- rep(0,d)

  for ( dim in c(1:d)) {
    meds1[dim] <- sign(median(X[,dim]))
  }

  meds2 <- rep(0,d)
  for ( dim in c(1:d)) {
    meds2[i] <- sign(median(Y[,dim]))
  }

  Q_median <- diag(meds1*meds2)

  #get the median:
  get_matched <- match_support(X = X,Y = Y,lambda_init = lambda_init,
                                      lambda_fina= lambda_final,eps = eps,
                                      numReps = numReps,eps_OT = eps_OT,
                                      alpha = alpha,Q = Q_median,
                                      p=p,q=q)

  return(get_matched)
}

