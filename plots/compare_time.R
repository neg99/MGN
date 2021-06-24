source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")

compare_time <- function(series, v_init, it, objective, subit = 16,
                         step_search = "mgn",
                         opt_method = FALSE, project_onto = project_onto_a_mgn,
                         weights = Diagonal(length(series)),
                         compensated = TRUE, vp_drop = Inf, cycles = 10) {
  N <- length(series)
  print(N)
  print(step_search)
  
  weights_chol <- chol(weights)
  
  if (step_search == "vp" && N >= vp_drop) {
    return(NA)
  }
  
  inv_weights_chol <- NULL
  
  if (step_search == "vp") {
    inv_weights_chol <- chol(solve(weights))
  }
  
  signal <- NULL
  mats <<- NULL
  vspace_pack <<- NULL
  
  v <- v_init
  
  fun_to_eval <- function() {
    signal <- project_onto(series, v_init, weights_chol, inv_weights_chol,
                           compensated = compensated)
    
    if (step_search == "vp") {
      step <- find_step_vp(signal, series, v, project_onto = project_onto,
                           weights_chol = weights_chol,
                           inv_weights_chol = inv_weights_chol)
    } else {
      step <- find_step(signal, series, v, vspace_pack, weights_chol)
    }
  }
  
  times <- sapply(1:cycles, function(x) {
    print(x)
    system.time(fun_to_eval())["elapsed"]
  })
  
  mean(times, trim = .2)
}

eval_for_one_N <- function(N, weights, vp_drop = Inf) {
  
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
  
  # weights <- Diagonal(N)
  # weights <- band_mat_from_diags(inv_ac_diags(N, c(0.9)))
  
  weights_chol <- chol(weights)
  
  noise <<- pre_noise -
    weighted_project_onto_vspace_qr(pre_noise,
                                    weighted_project_rotate_basis(tang_space_basis, weights_chol), tang_space_basis, weights_chol)
  
  series <- signal + noise
  
  v <- c(1, -3, 3, -1) + c(1, 1, 1, 1) * 1e-6
  
  
  c(compare_time(series, v, it, signal, step_search = "vp",
                 project_onto = project_onto_a_vp, weights = weights,
                 vp_drop = vp_drop),
    compare_time(series, v, it, signal, step_search = "vp", weights = weights,
                 vp_drop = vp_drop),
    compare_time(series, v, it, signal, step_search = "mgn", weights = weights,
                 compensated = FALSE),
    compare_time(series, v, it, signal, step_search = "mgn", weights = weights))
  
}

colors <- c("black", "red", "green3", "blue")

M <- 15
net <- as.integer(exp(seq(log(100), log(100000), length.out = M)))

white_noise_data <- sapply(net,
                           function(N) {
                             weights <- Diagonal(N)
                             eval_for_one_N(N, weights)
                           })

options(scipen=5)

matplot(net, t(white_noise_data), log = "xy", type = "l",
        xlab = "N", ylab = "Time, sec", col = colors,
        lty = c(1,2,1,2))

legend("topleft", c("VPGN", "S-VPGN-H", "MGN", "MGN-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

white_noise_data_moved <- t(white_noise_data)
white_noise_data_moved <- sapply(1:4, function(i) white_noise_data_moved[, i]/white_noise_data_moved[1, i])

pdf("time_diagonal.pdf", width = 2.1, height = 2.1, pointsize = 4)
par(mar = c(4.6,3.9,1.2,1.2))

matplot(net, white_noise_data_moved, log = "xy", type = "l",
        xlab = "N", ylab = "Times slower", col = colors,
        lty = c(1,2,1,2))

legend("topleft", c("VPGN", "S-VPGN-H", "MGN", "MGN-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

dev.off()

M <- 15
net <- as.integer(exp(seq(log(100), log(100000), length.out = M)))

red_noise_data <- sapply(net,
                         function(N) {
                           weights <- band_mat_from_diags(inv_ac_diags(N, c(0.9)))
                           eval_for_one_N(N, weights, vp_drop = 2000)
                         })



par(mar = c(4.6,3.9,1.2,1.2))

matplot(net, t(red_noise_data), log = "xy", type = "l",
        xlab = "N", ylab = "Time, sec", col = colors,
        lty = c(1,2,1,2))

legend("topright", c("VPGN", "S-VPGN-H", "MGN", "MGN-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

red_noise_data_moved <- t(red_noise_data)
red_noise_data_moved <- sapply(1:4, function(i) red_noise_data_moved[, i]/red_noise_data_moved[1, i])

pdf("time_3diagonal.pdf", width = 2.1, height = 2.1, pointsize = 4)
par(mar = c(4.6,3.9,1.2,1.2))

matplot(net, red_noise_data_moved, log = "xy", type = "l",
        xlab = "N", ylab = "Times slower", col = colors,
        lty = c(1,2,1,2))

legend("topright", c("VPGN", "S-VPGN-H", "MGN", "MGN-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)
dev.off()
