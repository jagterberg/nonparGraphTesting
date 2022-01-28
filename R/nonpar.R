#' Nonparametric Random Graph Testing
#' Wrapper function to take in two graphs and run the nonpar
#' test with alignment; however it may be more useful to use nonpar.test
#' directly after alignment
#' @param A An adjacency matrix
#' @param B another adjacency matrix
#' @param d the (known) dimension
#' @param nsims  the number of permutations to run
#' @param Q_init initial guess for orthogonal matrix
#' @param lambda_init initial value of lambda for regularized sinkhorn problem
#' @param lambda_final final value of lambda
#' @param eta decay rate of the lambdas
#' @param eps tolerance
#' @param numReps reps in the iterative sinkhorn divergence problem
#' @param alpha the type one error
#' @param match If true, run a optimal-transport procrustes
#' @export
nonpar <- function(A,B,d = 5,nsims =100,Q_init = NULL,lambda_init = .5
                   ,lambda_final = .01,eta=.5,eps=.01,numReps = 20,match = FALSE,alpha=.05) {
  Xhat <- ase(A,d)
  Yhat <- ase(B,d)
  if(match) { #todo: add more initializations
    matched <- match_support(Xhat, Yhat,Q = Q_init
                             ,lambda_init = lambda_init
                             ,lambda_final = lambda_final
                             ,alpha=eta
                             ,eps=eps
                             ,numReps=numReps)
    Ynew <- Yhat %*% matched$Q
  } else {
    Ynew <- Yhat
  }

  result <- nonpar.test(Xhat,Ynew,nsims,alpha=alpha)
  return(result)
}




#' actually run the test
#' @param Xhat the embedded first graph
#' @param Yhat the embedded second graph
#' @param nsims  the number of permutations to do
#' @param alpha the rejection type 1 error
#' @export
nonpar.test <- function(Xhat,Yhat,nsims = 100,alpha = .05) {


  dist.mat <- get_dist_matrix(Xhat,Yhat)
  i2 <- setdiff( c(1:(nrow(Xhat) + nrow(Yhat))), c(1:nrow(Xhat)) )
  U <- kernel.stat(Xhat,Yhat,dist=dist.mat,i1 = c(1:nrow(Xhat)),i2 = i2)
  #U <- kernel.stat(Xhat, Ynew)
  testresult <- run_perm_test(U,nsims,Xhat,Yhat,dist.mat = dist.mat,alpha=alpha)
  return(testresult)
}

#' Compute the Adjacency Spectral embedding
#' @param A a symmetric matrix
#' @param d the known dimension
#' @import irlba
#' @return adjacency spectral embedding
#' @export
ase <- function(A,d ) {
  A_svd <- irlba(A,d)
  Xhat <- A_svd$u %*% diag(A_svd$d)^(1/2)
  return(Xhat)
}


#' Function to get pairwise distances.
#' @param Z1 the n x d dataset
#' @param Z2 the m x d dataset
#' @param sigma the kernel sigma
#' @import stats
#' @return the kernel matrix
#' @export
get_dist_matrix <- function(Z1,Z2,sigma = .5) {
  #new_dat <- rbind(Z1,Z2)
  #D1 <- exp(-(as.matrix(stats::dist(new_dat))^2)/(2*sigma^2))
  m <- nrow(Z2)
  n <- nrow(Z1)
  D1 <- exp(-(as.matrix(stats::dist(Z1))^2)/(2*sigma^2))
  D2 <- exp(-(as.matrix(stats::dist(Z2))^2)/(2*sigma^2))
  D3 <- exp(-rect.dist(Z1,Z2)/(2*sigma^2))
  i1 <- c(1:nrow(Z1))
  i2 <- setdiff(c(1:(m+n)),c(1:nrow(Z1)))
  D <- matrix(0,n+m,n+m)
  D[i1,i1] <- D1
  D[i2,i2] <- D2
  D[i1,i2] <- D3
  D[i2,i1] <- t(D3)
  return(D)
}

#' Function to get the distance.  Original code by Youngser Park.
#' @param X n x d dataset
#' @param Y m x dataset
#' @return pairwise kernel n x m matrix
#' @export
rect.dist <- function(X,Y){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  m <- nrow(Y)
  tmp1 <- X%*%t(Y)
  tmp2 <- outer(rep(1, n), rowSums(Y^2))
  tmp3 <- outer(rowSums(X^2), rep(1,m))
  D <- tmp2 - 2*tmp1 + tmp3
  #D <- exp(-D/(2*(.5^2)))
  return(D)
}

#' Function to generate the U-statistic given the latent positions X and Y
#' and an optional choice of sigma.
#' Original code by Youngser Park.
#' @param X n x d datset
#' @param Y m x d dataset
#' @param sigma for the kernel
#' @param dist optional, the precalculated kernel matrix
#' @param i1 sets of indices of X to run the permutation test
#' @param i2 sets of indices of X to run the permutation test
#' @return the kernel statistic evaluated at the specified indices
#' @export
kernel.stat <- function(X,Y,sigma=0.5,dist = NULL,i1=c(1:nrow(X)),
                        i2=c((nrow(X) + 1):(nrow(X)+nrow(Y)))){

  n <- nrow(X)
  m <- nrow(Y)

  if (is.null(dist)) {
    tmpXX <- sum(exp(-(as.matrix(stats::dist(X))^2)/(2*sigma^2))) - n
    tmpYY <- sum(exp(-(as.matrix(stats::dist(Y))^2)/(2*sigma^2))) - m
    tmpXY <- sum(exp(-(rect.dist(X,Y))/(2*sigma^2)))

    tmp <- tmpXX/(n*(n-1)) + tmpYY/(m*(m-1)) - 2*tmpXY/(m*n)

    return((m+n)*tmp)
  } else {
    tmpXX <- sum(dist[i1,i1]) - n
    tmpYY <-  sum(dist[i2,i2]) - m
    tmpXY <- sum(dist[i1,i2])
    tmp<- tmpXX /(n*(n-1)) +tmpYY/(m*(m-1)) - 2*tmpXY/(m*n)
    return((m+n)*tmp)
  }


}


#' Helper function for a normalizing constant
#' @param X n x d dataset
get_s <- function(X) {
  return((1/nrow(X)^(.5))*norm(X,"F"))
}

#' Runs the permutation test given the value of the U-statistic
#' the number of repetitions, and the estimated latent positions
#' Code written in Fall of 2018, so may be buggy
#' @param U value of null U statistic
#' @param nsims number of simulations
#' @param X estimates for latent positions 1
#' @param Y estimates for latent positions 2
#' @param dist.mat Optional, a distance matrix precalculated already so as not to recalculate.
#' If null, the function is recursively run to calculate it
#' @param alpha the rejection level
#' @return either an estimated p-value under the null or the vector of permuted test statistics
#' @export
run_perm_test <- function(U,nsims,X,Y,dist.mat = NULL,alpha = .05) {
  toReturn <- rep(-1.0,nsims)
  m <- nrow(X)
  n <- nrow(Y)



  for (i in 1:nsims) {
    #cat(i," out of ",nsims,"\r")
    indices_1 <- sample(c(1:(m + n)),size=n,replace = FALSE)
    indices_2 <- setdiff( c(1:(m+n)), indices_1 )

    if(is.null(dist.mat)) {
      dist.mat <- get_dist_matrix(X,Y)
    }

    Uhat <- kernel.stat(X=X,Y=Y,i1=indices_1,i2=indices_2,dist=dist.mat)
    #if (!p.val) {
    toReturn[i] <- Uhat
    #} else {
    #  if (Uhat > U) {
    #    toReturn[i] <- 1.0
    #  } else {
    #    toReturn[i] <- 0.0
    #  }
    #}


  }

  crit_val <- sort(toReturn)[floor((1-alpha)*length(toReturn))]
  if (U > crit_val) {
    reject <- "reject the null"
  } else {
    reject <- "do not reject the null"
  }

  pval <- sum(toReturn > U)/nsims

  final_return <- list(toReturn,reject,crit_val,U,pval)
  names(final_return) <- c("permutation_results","reject","critical_value","test statistic","estimated p-value")
  return(final_return)

}










