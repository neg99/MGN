library(svd)
library(fftw)
library(Matrix)
library(rhlra2)

sapply_ns <- function(...) sapply(..., simplify = FALSE)

get_rev_row_form <- function(x) {
    trip <- as(x, "dgTMatrix")

    trip_i <- trip@i
    trip_j <- trip@j
    trip_data <- trip@x

    trip_diagonal <- (trip_i - trip_j) + 1
    trip_pos <- trip_i + 1
    trip <- matrix(0, ncol = max(trip_diagonal),
                        nrow = max(trip_pos))
    trip[cbind(trip_pos, trip_diagonal)] <- trip_data
    trip
}

inv_ac_diags <- function(N, coefs) {
    p <- length(coefs)
    convv <- c(1, -coefs)

    get_small_vec <- function(i) {
        cumsum(convv[(i + 1):length(convv)] * convv[1:(length(convv) - i)])
    }

    small_list <- lapply(0:p, get_small_vec)

    get_big_vec <- function(small_vec) {
        diag_num <- p + 1 - length(small_vec)
        sapply(1:(N - diag_num), function(i)
            small_vec[min(i, N - diag_num - i + 1, length(small_vec))])
    }

    lapply(small_list, get_big_vec)
}

band_mat_from_diags <- function(diags) {
    bandSparse(length(diags[[1]]), k = c(0:(length(diags) - 1)),
               diag = diags, symm=TRUE)
}

get_matrix_weight_matrix <- function(left_diags, right_diags) {
    L <- length(left_diags[[1]])
    K <- length(right_diags[[1]])
    N <- L + K - 1
    pl <- length(left_diags) - 1
    pr <- length(right_diags) - 1

    # i - pl - pr - 1

    answer <- lapply(1:(2 * pl + 2 * pr + 1),
                     function(i) numeric(N - abs(i - pl - pr - 1)))

    outer(-pl:pl, -pr:pr, Vectorize(function(i, j) {
        padding <- numeric((abs(i) + abs(j) - abs(i + j)) / 2)
        answer[[i + j + pl + pr + 1]] <<- answer[[i + j + pl + pr + 1]] +
            c(padding,
              ssa_convolve(left_diags[[abs(i) + 1]], right_diags[[abs(j) + 1]]),
              padding)
        0
    }))

    answer <- answer[(pl + pr + 1) : (2 * pl + 2 * pr + 1)]
    bandSparse(N, k = c(0:(pl + pr)), diag = answer, symm=TRUE)
}


ssa_convolve <- function(u, v) {
    p <- planFFT(length(u) + length(v) - 1)
    l_fft <- FFT(c(u, numeric(length(v) - 1)), plan = p)
    r_fft <- FFT(c(v, numeric(length(u) - 1)), plan = p)
    Re(IFFT(l_fft * r_fft, plan = p))
}

trmat_indices <- function(N, L) {
    K <- N - L + 1
    i <- rep(1:L, K)
    j <- rep(1:K, each = L)
    matrix(i + j - 1, nrow = L)
}

traj_matrix <- function(series, L) {
    outer(1:L, 1:(length(series) - L + 1),
          function(i, j) series[i + j - 1])
}

neg_fourier_short_mat <- function(N, indices) {
    full_indx <- outer(0:(N-1), indices)
    matrix(complex(argument = (full_indx %% N)*2*pi/N), nrow = N)
}

get_comp_space_by_v <- function(N, v, v_2 = FALSE, compensated = TRUE,
                                return_pack = FALSE) {
    if (!v_2) {
        result <- NULL
        if (compensated) {
            result <- eval_basis_compensated(N, v)
        } else {
            result <- eval_basis(N, v)
        }

        if (return_pack) {
            return(result)
        } else {
            return(result$basis)
        }
    }
    if (v_2) {
        result <- NULL
        if (compensated) {
            result <- eval_tangent_basis_compensated(N, v)
        } else {
            result <- eval_tangent_basis(N, v)
        }

        if (return_pack) {
            stop("Unsupported")
        }

        return(result)
    }
}

get_pseudograd <- function(N, vspace_pack, signal, j) {
    eval_pseudograd(N, vspace_pack$glrr, signal, j, vspace_pack)
}

get_space_by_v <- function(N, v, v_2 = FALSE) {
    pre_answer <- get_comp_space_by_v(N, v, v_2 = v_2)
    if (v_2) {
        v <- ssa_convolve(v, v)
    }
    r <- length(v) - 1
    j <- which.max(abs(v))
    begin_indices <- numeric(0)
    if (j > 1) begin_indices <- 1:(j - 1)
    end_indices <- numeric(0)
    if (j < r + 1) end_indices <- (N - r + j):N
    a_coefs <- solve(pre_answer[c(begin_indices, end_indices), ])
    Re(pre_answer %*% a_coefs)
}

get_v_from_nonlrf_series <- function(series, r) {
    L <- r + 1
    if (!is.list(series)) {
        N <- length(series)
        trmat <- traj_matrix(series[!is.na(series)], L)
    } else {
        Ns <- sapply(series, length)
        trmat <- do.call(cbind,
                         sapply_ns(series, function(x) traj_matrix(x[!is.na(x)], L)))
    }

    svd(trmat)$u[, L]
}

weighted_project_onto_vspace <- function(series, vspace, weights) {
    real_part <- as.matrix(weights %*% Re(vspace))
    imag_part <- as.matrix(weights %*% Im(vspace))
    wpvs <- real_part + 1i * imag_part
    as.numeric(Re(vspace %*% solve(t(Conj(vspace)) %*% wpvs,
                                   t(Conj(vspace)) %*% as.numeric(weights %*% series))))
}

weighted_project_rotate_basis <- function(vspace, chol_weights) {
    real_part <- as.matrix(chol_weights %*% Re(vspace))
    imag_part <- as.matrix(chol_weights %*% Im(vspace))
    vspace_rot <- real_part + 1i * imag_part
    qr(vspace_rot)
}

weighted_project_onto_vspace_qr <- function(series, qrobj, vspace, chol_weights) {
    series_rot <- as.numeric(chol_weights %*% series)

    coef <- qr.coef(qrobj, series_rot)

    Re(vspace %*% coef)
}

weighted_project_onto_vspace_coef <- function(series, vspace, chol_weights, tol = 1e-14) {
    if (!is.list(series)) {
        vspace_rot <- as.matrix(chol_weights %*% vspace)

        series_rot <- as.numeric(chol_weights %*% series)
    } else {
        vspace_rot <- do.call(rbind, sapply_ns(seq_along(series),
                    function(i) as.matrix(chol_weights[[i]] %*% vspace[[i]])))

        series_rot <- glue_series_lists(sapply_ns(seq_along(series),
                    function(i) as.numeric(chol_weights[[i]] %*% series[[i]])))
    }

    svdobj <- svd(vspace_rot)

    svdobj$d[svdobj$d < svdobj$d[1] * tol] <- Inf
    
    svdobj$v %*% ((t(svdobj$u) %*% series_rot) / svdobj$d)
}


find_step <- function(signal, series, v, vspace_pack, weights_chol, debug = FALSE, ...) {
    r <- length(v) - 1
    answer <- NA

    j <- which.max(abs(v))

    if (!is.list(signal)) {
        noise <- series - signal
        N <- length(series)

        K <- N - r

        pseudograd <- Re(get_pseudograd(N, vspace_pack, signal, j))
        pseudograd_minus <- sapply(1:r, function(i)
            weighted_project_onto_vspace_qr(Re(pseudograd[, i]), vspace_pack$qrobj,
                                            vspace_pack$basis, weights_chol))
        pseudograd <- pseudograd - pseudograd_minus
    } else {
        noise <- mapply("-", series, signal, SIMPLIFY = FALSE)
        Ns <- sapply(series, length)
        Ks <- Ns - r
        pseudograd <- sapply_ns(seq_along(series), function(k) {
            pseudograd_cur <- Re(get_pseudograd(Ns[[k]], vspace_pack[[k]], signal[[k]], j))
            pseudograd_minus <- sapply(1:r, function(i)
                weighted_project_onto_vspace_qr(Re(pseudograd_cur[, i]), vspace_pack[[k]]$qrobj,
                                                vspace_pack[[k]]$basis, weights_chol[[k]]))
            pseudograd_cur - pseudograd_minus
        })
    }

    used_coefs <- weighted_project_onto_vspace_coef(noise, pseudograd, weights_chol)

    ans <- numeric(r+1)

    ans[-j] <- used_coefs

    ans
}
