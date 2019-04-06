source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")

set.seed(15)

N <- 50

part_1 = .9^(1:N) * cos(pi * 1:N / 5)
part_2 = .2 * 1.05^(1:N) * cos(pi * 1:N / 12 + pi/4)

signal <- part_1 + part_2
r <- 4
noise <- rnorm(N)
series <- signal + 0.2 * sqrt(sum(signal ^ 2)) *  noise / sqrt(sum(noise ^ 2))


v_init <- svd(traj_matrix(series, r + 1))$u[, r + 1]
weights <- Diagonal(N)

answer <- run_hlra(series = series, 
                   v_init = v_init,
                   it = 100,
                   objective = NULL,
                   opt_method = TRUE,
                   weights = weights)

pdf("model_wo_gap.pdf", width = 4.2, height = 2.1, pointsize = 4)
par(mar = c(4.6,3.9,1.2,1.2))

matplot(1:N, cbind(series, answer, signal), pch = c(1, 26, 26),
        xlab = "Index", ylab = "Value")
lines(1:N, answer, col = "red", lty = 1)
lines(1:N, signal, col = "blue", lty = 2)
legend("top", c("Series", "Approximation", "Signal"), col=c("black", "red", "blue"),
       pch = c(1, 26, 26), lty = c(0, 1, 2))

dev.off()

# construct time series with gaps
series <- signal + 0.2 * sqrt(sum(signal ^ 2)) *  noise / sqrt(sum(noise ^ 2))
series[10:19] <- NA
series[35:39] <- NA
mask <- 1 - is.na(series)

# fill series gaps with mean
series_fill <- series
series_fill[!mask] <- mean(series, na.rm = TRUE)

weights <- Diagonal(N, mask)
# make appropriate weights matrix

v_init <- svd(traj_matrix(series_fill, r + 1))$u[, r + 1]
# take initial GLRR from SVD of trajectory matrix of series

answer <- run_hlra(series = series_fill, 
              v_init = v_init,
              it = 100,
              objective = NULL,
              opt_method = TRUE,
              weights = weights)

# perform HSLRA

# series -- initial point

# v_init -- initial GLRR

# it -- limit amount of iterations

# objective -- used for comparison for problems with known
# local minimum to make plots, not needed here

# opt_method = TRUE -- return approximation instead of
# technical information

# weights -- weights matrix of type Matrix

series_fill_to_plot <- signal
series_fill_to_plot[!!mask] <- NA

pdf("model_gap.pdf", width = 4.2, height = 2.1, pointsize = 4)
par(mar = c(4.6,3.9,1.2,1.2))

matplot(1:N, cbind(series, answer, signal, series_fill_to_plot), pch = c(1, 26, 26, 8),
        xlab = "Index", ylab = "Value", col = 1)
lines(1:N, answer, col = "red", lty = 1)
lines(1:N, signal, col = "blue", lty = 2)
legend("top", c("Series", "Approximation", "Signal", "Missing value"), 
       col=c("black", "red", "blue", "black"),
       pch = c(1, 26, 26, 8), lty = c(0, 1, 2, 0))

dev.off()