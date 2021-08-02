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

result_proj <- sapply(net, compare_projectors)

data_proj <- data.frame(t(result_proj))

pdf("fig_proj.pdf", width = 3.1, height = 2.1, pointsize = 4)
par(mar = c(4.6,3.9,1.2,1.2))

colors <- c("black", "red", "green3", "blue")

library(latex2exp)
matplot(net, cbind(data_proj[, c(3,1,2)]),
        log = "xy", type = "l", xlab = "N",
        ylab = TeX('$\\rho_{proj}$'),#"Distance to real projection",
        lty=c(1,2,1), col = colors)


legend("topleft", c("VP", "S-VP", "S-VP-H"),
       lty = c(1,2,1,2), bg = rgb(1, 1, 1, 0.85), col = colors)

dev.off()
