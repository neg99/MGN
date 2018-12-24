# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

eval_basis_compensated <- function(N, glrr) {
    # cat(sprintf("eval_basis_compensated\n"))
    N <- as.integer(N)
    glrr <- as.vector(glrr)
    dim(glrr) <- length(glrr)
    r <- length(glrr) - 1

    result <- .Call("eval_basis_compensatedC", N, glrr, PACKAGE = "rhlra2")
    # cat(sprintf("~eval_basis_compensated\n"))
    list(basis = matrix(result[1 : (N*r)], nrow = N),
         basis_fourier = matrix(result[(N*r + 1) : (2*N*r)], nrow = N),
         unitroots = result[(2*N*r + 1) : (N*(2 * r + 1))],
         A_f = result[(N*(2 * r + 1) + 1) : (N*(2 * r + 2))],
         alpha = Re(result[N*(2 * r + 2) + 1]),
         qrinvmat = matrix(result[(N*(2 * r + 2) + 2) : (N*(2 * r + 2) + r * r + 1)], nrow = r),
         glrr = glrr)
}

eval_basis <- function(N, glrr) {
    N <- as.integer(N)
    glrr <- as.vector(glrr)
    dim(glrr) <- length(glrr)
    r <- length(glrr) - 1

    result <- .Call("eval_basisC", N, glrr, PACKAGE = "rhlra2")

    list(basis = matrix(result[1 : (N*r)], nrow = N),
         basis_fourier = matrix(result[(N*r + 1) : (2*N*r)], nrow = N),
         unitroots = result[(2*N*r + 1) : (N*(2 * r + 1))],
         A_f = result[(N*(2 * r + 1) + 1) : (N*(2 * r + 2))],
         alpha = Re(result[N*(2 * r + 2) + 1]),
         glrr = glrr)
}

eval_pseudograd <- function(N, glrr, signal, tau, basis_obj) {
    # cat(sprintf("eval_pseudograd\n"))
    N <- as.integer(N)
    tau <- as.integer(tau)
    glrr <- as.numeric(glrr)
    dim(glrr) <- length(glrr)

    signal <- as.numeric(signal)
    dim(signal) <- length(signal)

    r <- length(glrr) - 1

    result <- .Call("eval_pseudogradC", N, glrr, basis_obj$A_f, basis_obj$alpha, tau, signal, PACKAGE = "rhlra2")
    # cat(sprintf("~eval_pseudograd\n"))
    matrix(result, nrow = N)
}

eval_pseudograd_compensated <- eval_pseudograd

eval_tangent_basis_compensated <- function(N, glrr) {
    N <- as.integer(N)
    glrr <- as.vector(glrr)
    dim(glrr) <- length(glrr)

    result <- .Call("eval_tangent_basis_compensatedC", N, glrr, PACKAGE = "rhlra2")

    matrix(result, nrow = N)
}

eval_tangent_basis <- function(N, glrr) {
    N <- as.integer(N)
    glrr <- as.vector(glrr)
    dim(glrr) <- length(glrr)

    result <- .Call("eval_tangent_basisC", N, glrr, PACKAGE = "rhlra2")

    matrix(result, nrow = N)
}

glue_series_lists <- function(series_list) {
    .Call("glue_series_lists", series_list, PACKAGE = "rhlra2")
}

generic_mul <- function(v, left_chol_mat,
                        right_chol_mat, series_fft) {
    .Call("generic_mulC", v, left_chol_mat, right_chol_mat, series_fft, PACKAGE = "rhlra2")
}

generic_diag_one_triple <- function(u, v, left_chol_mat,
                        right_chol_mat) {
    .Call("generic_diag_one_tripleC", u, v, left_chol_mat, right_chol_mat, PACKAGE = "rhlra2")
}
