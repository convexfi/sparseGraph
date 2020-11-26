library(spectralGraphTopology)

#' @title Learns sparse Laplacian matrix of a connected graph
#'
#' Learns a connected graph via non-convex, sparse promoting regularization
#' functions such as MCP, SCAD, and re-weighted l1-norm.
#'
#' @param S a pxp sample covariance/correlation matrix, where p is the number
#'        of nodes of the graph
#' @param w0 initial estimate for the weight vector the graph or a string
#'        selecting an appropriate method. Available methods are: "qp": finds w0
#'        that minimizes ||ginv(S) - L(w0)||_F, w0 >= 0; "naive": takes w0 as
#'        the negative of the off-diagonal elements of the pseudo inverse,
#'        setting to 0 any elements s.t. w0 < 0
#' @param alpha hyperparameter to control the level of sparsiness of the
#'        estimated graph
#' @param sparsity_type type of non-convex sparsity regularization. Available
#'        methods are: "mcp", "scad", "re-l1", and "none"
#' @param eps hyperparameter for the re-weighted l1-norm
#' @param eta learning rate
#' @param backtrack whether to update the learning rate using backtrack line
#'        search
#' @param maxiter maximum number of iterations
#' @param reltol relative tolerance on the Frobenius norm of the estimated
#'        Laplacian matrix as a stopping criteria
#' @param verbose whether or not to show a progress bar displaying the
#'        iterations
#' @return A list containing possibly the following elements:
#' \item{\code{laplacian}}{the estimated Laplacian Matrix}
#' \item{\code{adjacency}}{the estimated Adjacency Matrix}
#' \item{\code{maxiter}}{number of iterations taken to converge}
#' \item{\code{convergence}}{boolean flag to indicate whether or not the optimization converged}
#' \item{\code{elapsed_time}}{elapsed time recorded at every iteration}
#' @author Ze Vinicius, Jiaxi Ying, and Daniel Palomar
#' @export
learn_laplacian_pgd_connected <- function(S, w0 = "naive", alpha = 0, sparsity_type = "none",
                                          eps = 1e-4, gamma = 2.001, q = 1, backtrack = TRUE,
                                          maxiter = 10000, reltol = 1e-5, verbose = TRUE) {
  # number of nodes
  p <- nrow(S)
  Sinv <- MASS::ginv(S)
  # initial estimate
  w <- spectralGraphTopology:::w_init(w0, Sinv) + 1e-4
  J <- matrix(1, p, p) / p
  Lw <- L(w)
  H <- alpha * (diag(p) - p * J)
  # use safe initial learning rate
  eta <- 1 / (2*p)
  Lw <- L(w)
  K <- S + H
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = maxiter, clear = FALSE, width = 80)
  time_seq <- c(0)
  start_time <- proc.time()[3]
  for (i in 1:maxiter) {
    w_prev <- w
    tryCatch(
      {gradient <- spectralGraphTopology:::Lstar(K - spectralGraphTopology:::inv_sympd(Lw + J))},
        error = function(err) {
          results <- list(laplacian = L(w_prev), adjacency = A(w_prev), maxiter = i,
                     convergence = FALSE, elapsed_time = time_seq)
          return(results)
      }
    )
    if (backtrack) {
      fun <- connected_pgd.obj(Lw = Lw, J = J, K = K)
      while(1) {
        wi <- w - eta * gradient
        wi[wi < 0] <- 0
        # compute the objective function at the updated value of w
        fun_t <- connected_pgd.obj(Lw = L(wi), J = J, K = K)
        # check whether the previous value of the objective function is
        # smaller than the current one
        if (fun < fun_t - sum(gradient * (wi - w)) - (.5/eta)*norm(wi - w, '2')^2) {
          eta <- .5 * eta
        } else {
          eta <- 2 * eta
          break
        }
      }
    } else {
        wi <- w - eta * gradient
        wi[wi < 0] <- 0
    }
    if (verbose)
      pb$tick()
    has_converged <- (norm(L(wi) - Lw, 'F') / norm(Lw, 'F') < reltol) && (i > 1)
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    if (has_converged)
      break
    w <- wi
    Lw <- L(w)
    if (sparsity_type == "mcp") {
      H <- -(alpha + Lw / gamma) * (Lw >= -alpha*gamma)
      diag(H) <- rep(0, p)
      K <- S + H
    } else if (sparsity_type == "re-l1") {
      K <- S + H / ((-Lw + eps) ^ q)
    } else if (sparsity_type == "scad") {
      H <- -alpha * (Lw >= - alpha)
      H <- H + (-gamma * alpha - Lw) / (gamma - 1) * (Lw > -gamma*alpha) * (Lw < -alpha)
      diag(H) <- rep(0, p)
      K <- S + H
    }
  }
  results <- list(laplacian = L(wi), adjacency = A(wi), maxiter = i,
                  convergence = has_converged, elapsed_time = time_seq)
  return(results)
}


connected_pgd.obj <- function(Lw, J, K) {
  eigvals <- eigen(Lw + J, symmetric = TRUE, only.values = TRUE)$values
  eigvals <- eigvals[eigvals > 1e-8]
  return(sum(diag(Lw %*% K)) - sum(log(eigvals)))
  #chol_factor <- chol(Lw + J)
  #return(sum(diag(Lw %*% K)) - 2*sum(log(diag(chol_factor))))
}
