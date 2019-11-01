#'Match Datasets probabilistically
#'@description Matches the maps using probabilistic entropy maximization
#'@param X an n x d matrix of vectors
#'@param Y an n x d matrix of vectors
#'@param sigma a tuning parameter
#'@param numReps the number of iterations.
#'@return The final orthogonal matrix
#'@export
#'@examples
#'library(rstiefel)
#' set.seed(2019)
#' #generate a bunch of normals
#' X <- matrix(rnorm(1000,1,.2),ncol= 4)
#' #implicit assignment from X to Y
#' Y <- rbind(X,X)
#' #generate a random 4 x 4 orthogonal matrix
#' W <- rustiefel(4,4)
#'
#' #hit Y by the orthogonal matrix
#' Y <- Y %*% W
#'
#' #Solve the problem using both max entropy and optimal transport
#' test <- match_orthogonal(X,Y,numReps = 500)
#' test2 <- solve_optimal_transport(X,Y,numReps = 50, alpha =.5)
#' norm(test$`Orthogonal Matrix` - W,"2")
#'
#' X <- matrix(rnorm(1000,.2,.02),ncol= 5)
#' Y <- rbind(X,X)
#' W <- rustiefel(5,5)
#' Y <- Y %*% W
#' test <- match_orthogonal(X,Y,numReps = 100)
#'
#' test2 <- solve_optimal_transport(X,Y)
#' norm(test2$`Orthogonal Matrix` - W,"2")
#'
#' set.seed(2018)
#' X <- matrix(rnorm(1800,1,.1),ncol= 9)
#' Y <- rbind(X,X)
#' W <- rustiefel(9,9)
#' Y <- Y %*% W
#' test <- match_orthogonal(X,Y,numReps = 50)
#' test2 <- solve_optimal_transport(X,Y,lambda_init = 4,numReps = 200,Q=W)
#' norm(test2$`Orthogonal Matrix` - W,"2")
#'
#' D <- 3
#' N <- 50
#' M <- 60
#' X <- matrix(rnorm(N*D),N,D)
#' X <- t(t(X)/sqrt(colSums(X^2)))
#' Y <- matrix(rnorm(M*D),M,D)
#' Y = t(t(Y)/sqrt(colSums(Y^2)))
#' X = abs(X)
#' Y = -abs(Y)
#' sigma = 0.1
#' niter = 50
#' test <- match_orthogonal(X,Y,sigma=.1,numReps = 50)
#' test2 <- solve_optimal_transport(X,Y,lambda_init = 1)

match_support_entropy <- function(X,Y,sigma = .1,numReps=100) {
  d = dim(X)[2]
  A <- diag(1,d,d)

  for (it in c(1:numReps)) {
    X <- X%*%A
    C <- as.matrix(pdist::pdist(X,Y))^2
    K <- exp(-C/(2/sigma^2))
    P <- t(K / rowSums(K))
    Q <-  t(t(K) / colSums(K))
    M <-  t(Y) %*% (P+ t(Q)) %*% X
    svdVal <- svd(M)

    A <- svdVal$v %*% t(svdVal$u)
  }

  return(A)

}


