library(orthopolynom)
set.seed(15)

vspace_pack <- NULL
mats <- NULL
cur_v_global <- NULL

# NAIVE IMPLEMENTATION: absolutely unstable
compl_matrices <- function(N, v, inv_weights_chol) {
    r <- length(v) - 1
    K <- N - r
    qop <- bandSparse(N, K, 0:-r, matrix(rep(v, each = K), ncol = r + 1))
    if (is.null(inv_weights_chol)) {
        return(list(s = qop))
    }
    pregam <- inv_weights_chol %*% qop
    gam <- Cholesky(t(pregam) %*% pregam) #TODO: convolve by hand instead of multiplication
    list(s = qop, sts = gam)
}

project_onto_a_vp <- function(series, v, weights_chol, inv_weights_chol, keep_mats = FALSE,
                              compensated = NULL) {
    if (!keep_mats) {
        mats <<- compl_matrices(length(series), v, inv_weights_chol)
    }
    
    series <- as.numeric(series)
    series - t(inv_weights_chol) %*% 
        (inv_weights_chol %*% (as.numeric(mats$s %*% solve(mats$sts, system = "A",
                as.numeric(series %*% mats$s)))))
}

project_onto_a_mgn <- function(series, v, weights_chol, inv_weights_chol, 
                               keep_mats = NULL, compensated = TRUE) {
    vspace_pack <<- get_comp_space_by_v(length(series), v, return_pack = TRUE, 
                                        compensated = compensated)
    vspace_pack$qrobj <<- weighted_project_rotate_basis(vspace_pack$basis, weights_chol)
    weighted_project_onto_vspace_qr(series, vspace_pack$qrobj, vspace_pack$basis, weights_chol)
}


find_step_vp <- function(signal, series, v, debug = FALSE,
                         project_onto = project_onto_a_vp, 
                         weights_chol = Diagonal(length(series)), 
                         inv_weights_chol = Diagonal(length(series)),
                         svd_test = FALSE,
                         ...) {
    signal <- as.numeric(signal)
    series <- as.numeric(series)
    
    noise <- series - signal
    N <- length(series)
    r <- length(v) - 1
    K <- N - r

    answer <- NA
    
    if (is.null(mats)) {
        mats <- compl_matrices(N, v, inv_weights_chol)
    }
    
    make_grad <- function(i) {
        unit_i <- numeric(r + 1)
        unit_i[i] <- 1
        # cmat_unit_i <- compl_matrices(N, unit_i, NULL)
        h <- -as.numeric(t(inv_weights_chol) %*% (inv_weights_chol %*% 
                (mats$s %*% solve(mats$sts, system = "A", as.numeric(signal[i:(i+K-1)])))))
        
        prehath <- solve(mats$sts, system = "A",
                         as.numeric(series %*% mats$s))
        prehath_long <- numeric(N)
        prehath_long[i:(i+K-1)] <- prehath
        
        prehath2 <- as.numeric(t(inv_weights_chol) %*% (inv_weights_chol %*% prehath_long))
        
        hath <- -project_onto(prehath2, v, weights_chol, inv_weights_chol, TRUE)
        as.numeric(h + hath)
        # as.numeric(h)
    }

    do_business <- function() {
        j <- which.max(abs(v))

        ans <- numeric(r+1)

        indices <- (1:(r+1))[-j]

        grad <- sapply(indices, make_grad)

        if (svd_test) {
            return(singular_values_test(grad, 2 * r))
        }

        ans[-j] <- qr.coef(qr(as.matrix(weights_chol %*% grad)), 
                           as.numeric(weights_chol %*% noise))

        ans
    }

    answer <- do_business()
    # tryCatch({answer <- do_business()}, warning = function(x) {print(x)}, error = function(x) {print(x)})
    answer
}

# perform HSLRA

# series -- initial point

# v_init -- initial GLRR

# it -- limit amount of iterations

# objective -- used for comparison for problems with known
# local minimum to make plots

# opt_method = TRUE -- return approximation instead of
# technical information

# compensated = FALSE -- do not use compensated Horner scheme
# when calculating projections in GLRR (enabled by default)

# step_search -- use MGN or VP algorithm for calculating step

# project_onto -- use MGN or VP algorithm for calculating
# projections onto Z(a)

# weights -- weights matrix of type Matrix

# subit, criterion_split -- parameters for backtracking

run_hlra <- function(series, v_init, it, objective, subit = 50,
                          step_search = "mgn",
                          opt_method = FALSE, project_onto = project_onto_a_mgn,
                          criterion_split = 0, weights = Diagonal(length(series)),
                          compensated = TRUE, svd_test = FALSE) {
    N <- length(series)
    
    weights_chol <- NULL
    inv_weights_chol <- NULL
    
    if (any(diag(weights) == 0)) {
        mask <- which(diag(weights) == 0)
        weights[cbind(mask, mask)] <- 1
        weights_chol <- chol(weights)
        inv_weights_chol <- chol(solve(weights))
        weights[cbind(mask, mask)] <- 0
        weights_chol[cbind(mask, mask)] <- 0
        inv_weights_chol[cbind(mask, mask)] <- 0
    } else {
        weights_chol <- chol(weights)
        inv_weights_chol <- chol(solve(weights))
    }
    
    signal <- NULL
    mats <<- NULL
    vspace_pack <<- NULL
    signal <- project_onto(series, v_init, weights_chol, inv_weights_chol, 
                           compensated = compensated)
    
    dists <- rep(NA, it)
    steps <- rep(NA, it)
    v <- v_init
    prev_step <- 1

    closer <- function(signal, prev_signal, step, prev_step) {
        coef <- sqrt(sum(as.numeric(weights_chol %*% (prev_signal - signal))^2)) / sqrt(sum(as.numeric(weights_chol %*% prev_signal)^2))
        if (is.na(coef)) {
            return(FALSE)
        }
        if (coef > criterion_split) {
            return(sum(as.numeric(weights_chol %*%(prev_signal - signal))*as.numeric(weights_chol %*%(signal + prev_signal - 2*series))) > 0)
        }
        else {
            return(sum(step^2) < sum(prev_step^2))
        }

    }



    for (i in 1:it) {
        step <- NULL
        if (step_search == "vp") {
            step <- find_step_vp(signal, series, v, project_onto = project_onto, 
                                 weights_chol = weights_chol,
                                 inv_weights_chol = inv_weights_chol)
        } else {
            step <- find_step(signal, series, v, vspace_pack, weights_chol)
        }

        steps[i] <- sqrt(sum(step^2))

        improve_flag = FALSE
        dstep <- step

        for (j in 1:subit) {
            new_v <- v + dstep
            dstep <- dstep / 2
            new_signal <- project_onto(series, new_v, weights_chol, inv_weights_chol, 
                                       compensated = compensated)

            if (closer(new_signal, signal, step, prev_step)) {
                improve_flag <- TRUE
                signal <- new_signal
                v <- new_v
                cur_v_global <<- v
                prev_step <- step
                break
            }
        }

        if (improve_flag) {
            ru <- length(v)
            # print(svd(traj_matrix(signal, ru))$d[ru])
            
            dists[i] <- sqrt(sum((signal - series)^2)) - sqrt(sum((series - objective)^2))
        } else {
            break
        }
    }

    if (svd_test) {
        test_info <- NULL

        if (step_search == "vp") {
            test_info <- find_step_vp(signal, series, v, project_onto = project_onto, 
                                 weights_chol = weights_chol,
                                 inv_weights_chol = inv_weights_chol, svd_test = TRUE)
        } else {
            test_info <- find_step(signal, series, v, vspace_pack, weights_chol, svd_test = TRUE)
        }

        return(list(signal = signal, test_info = test_info))
    }

    if (opt_method) {
        return(signal)
    }
    dists
}

vp_demo <- function() {
    N <- 1000
    
    signal <- seq(-1, 1, length.out = N)^2
    norm_const <- sqrt(sum(signal^2))
    signal <- signal / norm_const
    
    net <- seq(from = -1, to = 1, length.out = N)
    
    pre_noise <- abs(net) / sqrt(sum((abs(net))^2))
    polys <- legendre.polynomials(5)
    tang_space_basis <- cbind(as.function(polys[[1]])(net) / sqrt(N),
                              as.function(polys[[2]])(net) / sqrt(N),
                              as.function(polys[[3]])(net) / sqrt(N),
                              as.function(polys[[4]])(net) / sqrt(N),
                              as.function(polys[[5]])(net) / sqrt(N),
                              as.function(polys[[6]])(net) / sqrt(N))
    
    weights <- Diagonal(N)
    # weights <- band_mat_from_diags(inv_ac_diags(N, c(0.9)))
    
    weights_chol <- chol(weights)
    
    noise <<- pre_noise -
        weighted_project_onto_vspace_qr(pre_noise, 
                                        weighted_project_rotate_basis(tang_space_basis, weights_chol), tang_space_basis, weights_chol)
    
    series <- signal + noise
    
    matplot(1:N, cbind(signal, series), type = "l")
    
    v <- c(1, -3, 3, -1) + c(1, 1, 1, 1) * 1e-6
    
    it <- 40
    
    data <- cbind(run_hlra(series, v, it, signal, step_search = "vp",
                                project_onto = project_onto_a_vp, weights = weights),
                  run_hlra(series, v, it, signal, step_search = "vp", weights = weights),
                  run_hlra(series, v, it, signal, step_search = "mgn", weights = weights,
                                compensated = FALSE),
                  run_hlra(series, v, it, signal, step_search = "mgn", weights = weights,
                                compensated = TRUE))
    symbs <- c(0,1,2,5)
    symbs_all <- c(0,1,2,5,15,16,17,18)
    colors <- c("black", "red", "green3", "blue")
    
    symb <- matrix((data > 0) * 4, nrow = it)
    
    matplot(1:it, abs(data),
            log = "y", type = "b", pch = " ", col = colors, xlab = "Iteration", ylab = "Distance to series")
    points(rep(1:it, 4), as.numeric(abs(data)), pch = symbs_all[symb + rep(1:4, each = it)], 
           col = rep(colors, each = it))
    
    legend("topright", c("VPGN", "S-VPGN", "MGN", "S-MGN"), pch = symbs, col = colors)
}

