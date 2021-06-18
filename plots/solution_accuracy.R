source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")

# library(Rslra)
library(orthopolynom)
set.seed(15)

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

svt_projection <- function(signal) {
    singular_values_test(as.vector(signal), r)
}

eval_disperancy <- function(N, debug = FALSE) {
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

    matplot(1:N, cbind(signal, series), type = "l")

    ideal_weights <- weights
    v <- c(1, -3, 3, -1) + 1e-6 * c(1, 1, 1, 1)


    mgn_data <- run_hlra(series, v, 100, 0, svd_test = TRUE,
                                compensated = FALSE)
    mgn_approx <- mgn_data$signal

    v_mgn <- cur_v_global

    s_mgn_data <- run_hlra(series, v, 100, 0, svd_test = TRUE,
                                compensated = TRUE)
    s_mgn_approx <- s_mgn_data$signal

    v_s_mgn <- cur_v_global


    # v_basis <- get_comp_space_by_v(N, v)
    # v_basis <- cbind(as.function(polys[[1]])(net) / sqrt(N),
    #                           as.function(polys[[2]])(net) / sqrt(N),
    #                           as.function(polys[[3]])(net) / sqrt(N))

    Rini <- t(v)

    kostya_approx <- NULL
    kostya_approx_data <- NULL
    v_vp <- NULL

    # tryCatch({kostya_approx <- slra(series, list(m=r+1), r, opt = list(disp='off',
    #                                                                    Rini=Rini),
    #                                compute.Rh=TRUE, compute.ph=TRUE)$ph},
    #          warning = function(x) {print(x)}, error = function(x) {print(x)})
    tryCatch({kostya_approx_data <- run_hlra(series, v, 100, 0, svd_test = TRUE,
                                             step_search = "vp",
                                             project_onto = project_onto_a_vp)
             kostya_approx <- kostya_approx_data$signal
             v_vp <- cur_v_global},
             warning = function(x) {print(x)}, error = function(x) {print(x)})


    kostya_approx_comp <- NULL
    kostya_approx_comp_data <- NULL
    v_vp_comp <- NULL
    tryCatch({kostya_approx_comp_data <- run_hlra(series, v, 100, 0, svd_test = TRUE,
                                             step_search = "vp")
             kostya_approx_comp = kostya_approx_comp_data$signal
             v_vp_comp <- cur_v_global},
             warning = function(x) {print(x)}, error = function(x) {print(x)})

    print(v_vp_comp)

    # list(nikita = distance(my_approx, series), kostya = distance(kostya_approx, series),
    # true = distance(signal, series))

    c(mgn = distance(mgn_approx, signal),
         mgn_comp = distance(s_mgn_approx, signal),
         vp = distance(kostya_approx, signal),
         vp_comp = distance(kostya_approx_comp, signal),

         mgn_svtp = svt_projection(mgn_approx),
         mgn_comp_svtp = svt_projection(s_mgn_approx),
         vp_svtp = svt_projection(kostya_approx),
         vp_comp_svtp = svt_projection(kostya_approx_comp),

         mgn_dd = distance_diff(mgn_approx, signal, series),
         mgn_comp_dd = distance_diff(s_mgn_approx, signal, series),
         vp_dd = distance_diff(kostya_approx, signal, series),
         vp_comp_dd = distance_diff(kostya_approx_comp, signal, series),

         mgn_svtg = mgn_data$test_info$singular_part,
         mgn_comp_svtg = s_mgn_data$test_info$singular_part,
         vp_svtg = kostya_approx_data$test_info$singular_part,
         vp_comp_svtg = kostya_approx_comp_data$test_info$singular_part,

         mgn_cng = mgn_data$test_info$cond_number_part,
         mgn_comp_cng = s_mgn_data$test_info$cond_number_part,
         vp_cng = kostya_approx_data$test_info$cond_number_part,
         vp_comp_cng = kostya_approx_comp_data$test_info$cond_number_part)
}

M <- 30
net <- as.integer(exp(seq(log(20), log(50000), length.out = M)))

# M <- 5
# net <- as.integer(exp(seq(log(20), log(1000), length.out = M)))


result <- sapply(net, eval_disperancy)

data <- data.frame(t(result))

# pdf("kostya_comp.pdf", width = 2.1, height = 2.1, pointsize = 4)
# par(mar = c(4.6,3.9,1.2,1.2))

colors <- c("black", "red", "green3", "blue")

matplot(net, cbind(data[, c(1:4)]),
        log = "xy", type = "l", xlab = "N",
        ylab = "Distance to solution",
        lty=c(1,2,1,2), col = colors)


legend("topleft", c("MGN", "S-MGN", "VPGN", "S-VPGN-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

# dev.off()

alive_amount <- 1:length(net)

if (any(is.na(data[, 7]))) {
    alive_amount <- 1:(which.max(is.na(data[, 7])) - 1)
}


# pdf("projection_svt.pdf", width = 3.2, height = 2.1, pointsize = 3)
# par(mar = c(4.6,3.9,1.2,1.2))

matplot(net[alive_amount], data[alive_amount, c(7:8, 5:6)], log = "xy", type = "l", xlab = "N",
        ylab = "Proportion of residual singular values, projection", lty = c(1,2,1,2), col = colors)

legend("topleft", c("VPGN",
                    "S-VPGN-H", "MGN", "MGN-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

# dev.off()

# pdf("gradient_svt.pdf", width = 3.2, height = 2.1, pointsize = 3)
# par(mar = c(4.6,3.9,1.2,1.2))


matplot(net[alive_amount], data[alive_amount, c(15:16, 13:14)], log = "xy", type = "l", xlab = "N",
        ylab = "Proportion of residual singular values, gradient", lty = c(1,2,1,2), col = colors)

legend("topleft", c("VPGN",
                    "S-VPGN-H", "MGN", "MGN-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

# dev.off()

# pdf("gradient_cn.pdf", width = 3.2, height = 2.1, pointsize = 3)
# par(mar = c(4.6,3.9,1.2,1.2))

matplot(net[alive_amount], data[alive_amount, c(19:20, 17:18)], log = "xy", type = "l", xlab = "N",
        ylab = "Condition number, gradient", lty = c(1,2,1,2), col = colors)

legend("topleft", c("VPGN",
                    "S-VPGN-H", "MGN", "MGN-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

# dev.off()


# pdf("gradient_overall.pdf", width = 3.2, height = 2.1, pointsize = 3)
# par(mar = c(4.6,3.9,1.2,1.2))

matplot(net[alive_amount], cbind( as.matrix(data[alive_amount, c(15:16, 13:14)]) * as.matrix(data[alive_amount, c(19:20, 17:18)]), 1), log = "xy", type = "l", xlab = "N",
        ylab = "Summary conditioning, gradient", lty = c(1,2,1,2,3), col = colors)

legend("topleft", c("VPGN",
                    "S-VPGN-H", "MGN", "MGN-H", "1"),
       lty = c(1,2,1,2,3), bg = rgb(1, 1, 1, 0.85), col = colors)

# dev.off()

# pdf("dist_diff.pdf", width = 4.2, height = 2.1, pointsize = 3)
# par(mar = c(4.6,3.9,1.2,1.2))
# 
# disp_data <- cbind(as.numeric(data[, 9]), as.numeric(data[, 10]),
#                    as.numeric(data[, 11]), as.numeric(data[, 12]))
# 
# symbs <- c(0,1,2,5)
# symbs_all <- c(0,1,2,5,15,16,17,18)
# 
# symb <- (disp_data > 0) * 4
# 
# matplot(net, abs(disp_data), log = "xy", type = "b", xlab = "N",
#         pch = " ",
#         ylab = "Abs. value of Distance Difference", lty = c(1,2,1,2), col = colors)
# points(rep(net, 4), as.numeric(abs(disp_data)), pch = symbs_all[symb + rep(1:4, each = M)],
#        col = rep(colors, each = M))
# 
# 
# legend("topleft", c("MGN", "S-MGN", "VPGN",
#                     "S-VPGN-H"), bg = rgb(1, 1, 1, 0.85), col = colors,
#        pch = symbs)
# dev.off()

stop()

eval_hessian <- function(N, h = 1e-6) {
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
    weights_chol <- chol(weights)

    series <- signal + (pre_noise -
                            weighted_project_onto_vspace_qr(pre_noise,
                            weighted_project_rotate_basis(tang_space_basis, weights_chol),
                            tang_space_basis, weights_chol))

    v <- c(1, -3, 3, -1)
    hes_fun <- function(x) {
        nx <- numeric(4)
        nx[c(1, 3, 4)] <- x
        v <- c(1, -3, 3, -1) + nx
        space <- get_comp_space_by_v(N, v)

        sum((series - weighted_project_onto_vspace_qr(series,
                            weighted_project_rotate_basis(space, weights_chol), space, weights_chol))^2)
    }

    hessian(hes_fun, numeric(3), h = h)
}

# eigen(eval_hessian(20, h = 1e-7))
# eigen(eval_hessian(50, h = 1e-7))
# eigen(eval_hessian(100, h = 1e-7))
