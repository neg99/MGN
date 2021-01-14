source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")

# library(Rslra)
library(orthopolynom)

r <- 3
coefs <- numeric(0)
sigma <- 0.01
max_step <- 1

distance <- function(series, signal) {
    sqrt(sum((series - signal)^2))
}

distance_diff <- function(estimate, objective, input) {
    as.numeric(sqrt(sum((estimate - input)^2)) - sqrt(sum((input - objective)^2)))
}

eval_methods <- function(N, error = NULL) {
    print(N)
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
    weights_chol <- Diagonal(N)

    series <- signal + (pre_noise -
       weighted_project_onto_vspace_qr(pre_noise,
        weighted_project_rotate_basis(tang_space_basis, weights_chol),
        tang_space_basis, weights_chol))

    # matplot(1:N, cbind(signal, series), type = "l")

    ideal_weights <- weights
    if (is.null(error)) {
        error <- runif(4, -1, 1)
    }
    v <- c(1, -3, 3, -1) + 1e-6 * error

    vpgn_approx <- NULL
    v_vp <- NULL

    # tryCatch({kostya_approx <- slra(series, list(m=r+1), r, opt = list(disp='off',
    #                                                                    Rini=Rini),
    #                                compute.Rh=TRUE, compute.ph=TRUE)$ph},
    #          warning = function(x) {print(x)}, error = function(x) {print(x)})
    tryCatch({vpgn_approx <- run_hlra(series, v, 100, 0, opt_method = TRUE,
                                             step_search = "vp",
                                             project_onto = project_onto_a_vp)
             v_vp <- cur_v_global},
             warning = function(x) {print(x)}, error = function(x) {print(x)})


    s_vpgn_approx <- NULL
    v_vp_comp <- NULL
    tryCatch({s_vpgn_approx <- run_hlra(series, v, 100, 0, opt_method = TRUE,
                                             step_search = "vp",
                                             compensated = FALSE)
             v_vp_comp <- cur_v_global},
             warning = function(x) {print(x)}, error = function(x) {print(x)})

    s_vpgn_h_approx <- NULL
    v_vp_comp <- NULL
    tryCatch({s_vpgn_h_approx <- run_hlra(series, v, 100, 0, opt_method = TRUE,
                                             step_search = "vp",
                                             compensated = TRUE)
    v_vp_comp <- cur_v_global},
    warning = function(x) {print(x)}, error = function(x) {print(x)})


    # list(nikita = distance(my_approx, series), kostya = distance(kostya_approx, series),
    # true = distance(signal, series))

    c(vpgn = distance(vpgn_approx, signal),
         s_vpgn = distance(s_vpgn_approx, signal),
         s_vpgn_h = distance(s_vpgn_h_approx, signal),
         vpgn_dd = distance_diff(vpgn_approx, signal, series),
         s_vpgn_dd = distance_diff(s_vpgn_approx, signal, series),
         s_vpgn_h_dd = distance_diff(s_vpgn_h_approx, signal, series))
}

set.seed(15)

M <- 15
net <- as.integer(exp(seq(log(20), log(5000), length.out = M)))

# M <- 30
# net <- as.integer(exp(seq(log(20), log(50000), length.out = M)))

It <- 15

result_sum <- NULL

for (k in 0:It) {
    print(paste("Iteration:", k))
    # result <- sapply(net, eval_methods)
    result <- sapply(net, function(i) {
        eval_methods(i, 2 * c(bitwShiftR(k, 3)%%2, bitwShiftR(k, 2)%%2, bitwShiftR(k, 1)%%2, k %% 2) - 1)
    })
    if (is.null(result_sum)) {
        result_sum <- result
    } else {
        result_sum <- result_sum + result
    }
}

result <- result_sum / It

result <- sapply(net, eval_methods)

data <- data.frame(t(result))

# pdf("kostya_comp.pdf", width = 2.1, height = 2.1, pointsize = 4)
# par(mar = c(4.6,3.9,1.2,1.2))

colors <- c("black", "red", "green3", "blue")

matplot(net, cbind(data[, c(1:3)]),
        log = "xy", type = "l", xlab = "N",
        ylab = "Distance to solution",
        lty=c(1,2,1,2), col = colors)


legend("topleft", c("VPGN", "S-VPGN", "S-VPGN-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

# dev.off()

# pdf("dist_diff.pdf", width = 4.2, height = 2.1, pointsize = 3)
# par(mar = c(4.6,3.9,1.2,1.2))

disp_data <- cbind(as.numeric(data[, 4]), as.numeric(data[, 5]),
                   as.numeric(data[, 6]))

symbs <- c(0,1,2,5)
symbs_all <- c(0,1,2,5,15,16,17,18)

symb <- (disp_data > 0) * 3

matplot(net, abs(disp_data), log = "xy", type = "b", xlab = "N",
        pch = " ",
        ylab = "Abs. value of Distance Difference", lty = c(1,2,1,2), col = colors)
points(rep(net, 4), as.numeric(abs(disp_data)), pch = symbs_all[symb + rep(1:3, each = M)],
       col = rep(colors, each = M))


legend("topleft", c("VPGN", "S-VPGN", "S-VPGN-H"), bg = rgb(1, 1, 1, 0.85), col = colors,
       pch = symbs)
# dev.off()

stop()

# eval_hessian <- function(N, h = 1e-6) {
#     print(N)
#     signal <- seq(-1, 1, length.out = N)^2
#     norm_const <- sqrt(sum(signal^2))
#     signal <- signal / norm_const
#
#     net <- seq(from = -1, to = 1, length.out = N)
#
#     pre_noise <- abs(net) / sqrt(sum((abs(net))^2))
#     polys <- legendre.polynomials(5)
#     tang_space_basis <- cbind(as.function(polys[[1]])(net) / sqrt(N),
#                               as.function(polys[[2]])(net) / sqrt(N),
#                               as.function(polys[[3]])(net) / sqrt(N),
#                               as.function(polys[[4]])(net) / sqrt(N),
#                               as.function(polys[[5]])(net) / sqrt(N),
#                               as.function(polys[[6]])(net) / sqrt(N))
#
#     weights <- Diagonal(N)
#     weights_chol <- chol(weights)
#
#     series <- signal + (pre_noise -
#                             weighted_project_onto_vspace_qr(pre_noise,
#                             weighted_project_rotate_basis(tang_space_basis, weights_chol),
#                             tang_space_basis, weights_chol))
#
#     v <- c(1, -3, 3, -1)
#     hes_fun <- function(x) {
#         nx <- numeric(4)
#         nx[c(1, 3, 4)] <- x
#         v <- c(1, -3, 3, -1) + nx
#         space <- get_comp_space_by_v(N, v)
#
#         sum((series - weighted_project_onto_vspace_qr(series,
#                             weighted_project_rotate_basis(space, weights_chol), space, weights_chol))^2)
#     }
#
#     hessian(hes_fun, numeric(3), h = h)
# }

# eigen(eval_hessian(20, h = 1e-7))
# eigen(eval_hessian(50, h = 1e-7))
# eigen(eval_hessian(100, h = 1e-7))


compare_projectors <- function(N, debug = FALSE) {
    print(N)
    signal <- seq(-1, 1, length.out = N)^2
    norm_const <- sqrt(sum(signal^2))
    signal <- signal / norm_const

    net <- seq(from = -1, to = 1, length.out = N)

    pre_noise <- abs(net) / sqrt(sum((abs(net))^2))
    polys <- legendre.polynomials(2)
    space_basis <- cbind(as.function(polys[[1]])(net) / sqrt(N),
                              as.function(polys[[2]])(net) / sqrt(N),
                              as.function(polys[[3]])(net) / sqrt(N))

    weights <- Diagonal(N)
    weights_chol <- Diagonal(N)

    series <- signal + (pre_noise -
                            weighted_project_onto_vspace_qr(pre_noise,
                                                            weighted_project_rotate_basis(space_basis, weights_chol),
                                                            space_basis, weights_chol))

    # matplot(1:N, cbind(signal, series), type = "l")

    ideal_weights <- weights
    v <- c(1, -3, 3, -1)


    mgn_proj <- project_onto_a_mgn(series, v, weights_chol, weights_chol, compensated = FALSE)

    mgn_proj_comp <- project_onto_a_mgn(series, v, weights_chol, weights_chol, compensated = TRUE)

    vp_proj <- project_onto_a_vp(series, v, weights_chol, weights_chol)


    c(mgn = distance(mgn_proj, signal),
         mgn_comp = distance(mgn_proj_comp, signal),
         vp = distance(vp_proj, signal))
}

# M <- 5
# net <- as.integer(exp(seq(log(20), log(10000), length.out = M)))

M <- 30
net <- as.integer(exp(seq(log(20), log(50000), length.out = M)))

result <- sapply(net, compare_projectors)

data <- data.frame(t(result))

# pdf("kostya_comp.pdf", width = 2.1, height = 2.1, pointsize = 4)
# par(mar = c(4.6,3.9,1.2,1.2))

colors <- c("black", "red", "green3", "blue")

matplot(net, cbind(data[, c(1:3)]),
        log = "xy", type = "l", xlab = "N",
        ylab = "Distance to real projection",
        lty=c(1,2,1), col = colors)


legend("topleft", c("S-VPGN", "S-VPGN-H", "VPGN"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

# dev.off()
