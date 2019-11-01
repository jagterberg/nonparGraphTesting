#' Nonparametric Random Graph Testing
#' @param A An adjacency matrix
#' @param B another adjacency matrix
#' @param d the (known) dimension
#' @export
nonpar <- function(A,B,d = 5) {
  Xhat <- ase(A,d)
  Yhat <- ase(B,d)
  result <- nonpar.test(Xhat,Yhat)
}

#' actually run the test
nonpar.test <- function(Xhat,Yhat) {

}

#' Compute the Adjacency Spectral embedding
#' @param A a symmetric matrix
#' @param d the known dimension
#' @export
ase <- function(A,d) {


}
