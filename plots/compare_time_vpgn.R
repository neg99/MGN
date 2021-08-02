source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")

library(orthopolynom)

r <- 3
coefs <- numeric(0)
sigma <- 0.01

compare_time <- function(N, I = 10, B = 10, debug = FALSE) {
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
  
  
  mgn_proj <- function() for (j in 1:B) project_onto_a_mgn(series, v, weights_chol, weights_chol, compensated = FALSE)
  
  mgn_proj_comp <- function() for (j in 1:B) project_onto_a_mgn(series, v, weights_chol, weights_chol, compensated = TRUE)
  
  vp_proj <- function() for (j in 1:B) project_onto_a_vp(series, v, weights_chol, weights_chol)
  
  time_result <- function(fun_to_eval) {
    times <- sapply(1:I, function(x) {
      system.time(fun_to_eval())["elapsed"]/B
    })
    
    mean(times, trim = .2)
  }
  
  c(mgn = time_result(mgn_proj),
    mgn_comp = time_result(mgn_proj_comp),
    vp = time_result(vp_proj))
}

# M <- 5
# net <- as.integer(exp(seq(log(200), log(2000), length.out = M)))

M <- 15
net <- as.integer(exp(seq(log(100), log(100000), length.out = M)))

result_proj <- sapply(net, compare_time)

data_proj <- data.frame(t(result_proj))

save.image("compare_time_vpgn.RData")
# load("compare_time_vpgn.RData")


# pdf("fig_time.pdf", width = 3.1, height = 2.1, pointsize = 4)
# par(mar = c(4.6,3.9,1.2,1.2))

colors <- c("black", "red", "green3", "blue")

matplot(net, cbind(data_proj[, c(3,1,2)]),
        log = "xy", type = "l", xlab = "N",
        ylab = "Average time, sec.",
        lty=c(1,2,1), col = colors)


legend("topleft", c("VP", "S-VP", "S-VP-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

# dev.off()

result_moved <- sapply(1:3, function(i) result_proj[i, ]/result_proj[i, 1])

# pdf("fig_time_moved.pdf", width = 3.1, height = 2.1, pointsize = 4)
# par(mar = c(4.6,3.9,1.2,1.2))

colors <- c("black", "red", "green3", "blue")

matplot(net, cbind(result_moved[, c(3,1,2)]),
        log = "xy", type = "l", xlab = "N",
        ylab = "Times slower",
        lty=c(1,2,1), col = colors)


legend("topleft", c("VP", "S-VP", "S-VP-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

# dev.off()
