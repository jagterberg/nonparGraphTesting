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
#' @import pdist
#' @return P, the n x m matrix of assignments
#' @export
optimal_transport <- function(X,Y, Q= NULL,lambda = .1,eps = .01) {
  if(is.null(Q)) {
    d <- dim(X)[2]
    Q <- diag(1,d,d)
  }
  X <- X%*%Q
  C <- as.matrix(pdist::pdist(X,Y))^2

  n <- dim(X)[1]
  m <- dim(Y)[1]
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
#' @export
#' @return the procrustes solution
procrustes <- function(X,Y, Pi = NULL) {
  if (is.null(Pi)) {
    stop("need to run optimal transport first")
  }
  vals <- svd(t(X)%*% Pi %*% Y)
  toReturn <- vals$u %*% t(vals$v)
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
#' norm(test$`Orthogonal Matrix` - W,"2")
#' X <- matrix(rnorm(5000,.2,.02),ncol= 5)
#' Y <- rbind(X,X)
#' W <- rustiefel(5,5)
#' Y <- Y %*% W
#' Y <- matrix(rnorm(200,.7),ncol =5)
#' test2 <- iterative_optimal_transport(X,Y,numReps = 1000,lambda = .0001)
#' norm(test2$`Orthogonal Matrix` - W,"2")
iterative_optimal_transport <-function(X,Y, Q = NULL,lambda = .01,eps = .01,numReps =1000) {
  if(is.null(Q)) {
    d <- dim(X)[2]
    Q <- diag(1,d,d)
  }
  Pi <- optimal_transport(X,Y,Q)
  Q <- procrustes(X,Y,Pi)
  c <- norm(X %*% Q - Pi%*% Y,"F")
  i <- 1
  while ( i < numReps) {
    if( c > eps ) {

      Pi <- optimal_transport(X,Y,Q)
      Q <- procrustes(X,Y,Pi)
      c <- norm(X %*% Q - Pi%*% Y)
      i <- i+1
    } else {
      break
    }

  }
  toReturn <- list(Pi,Q,c)
  names(toReturn) <- c("Pi", "Orthogonal Matrix","obj.value")

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
#' @param alpha the parameter for which lambda is multiplied by
#' @param eps the tolerance for the optimal transport problem
#' @param numReps the number of reps for each subiteration
#' @return a list of the final orthogonal matrix and the assignment matrix
#' @export
#' @import pdist
#' @import rstiefel
#' @examples
#' library(rstiefel)
#'set.seed(2019)
#'X <- matrix(rnorm(100,1,.2),ncol= 4)
#' Y <- rbind(X,X)
#'
#'
#' W <- rustiefel(4,4)
#' Y <- Y %*% W
#' test <- match_support(X,Y)
#'
#'
#' # others have pointed out that initializing Q at all 2^d sign matrices ( diagonal matrices
#'  # whose entries are plus or minus one) might have better global
#'  # convergence
#'
#' X <- matrix(rnorm(100,.2,.02),ncol= 5)
#' Y <- rbind(X,X)
#' W <- rustiefel(5,5)
#' Y <- Y %*% W
#' test2 <- match_support(X,Y)
#'
match_support <- function(X,Y,
                  Q = NULL,lambda_init = .5, lambda_final = .01,alpha = .95,
                  eps = .01,numReps =100) {

  lambda <- lambda_init
  if (is.null(Q)) {
    d <- dim(X)[2]
    Q <- diag(1,d,d)
  }
  while(lambda > lambda_final) {
    Q <- iterative_optimal_transport(X,Y,Q,lambda = lambda,eps = eps,numReps = numReps)
    Pi <- Q$Pi
    c <- Q$`obj.value`
    Q <- Q$`Orthogonal Matrix`
    lambda <- alpha*lambda
  }

  toReturn <- list(Pi,Q,c)
  #toReturn <- Q
  names(toReturn) <- c("Pi", "Q","obj.value")
  return(toReturn)

}





