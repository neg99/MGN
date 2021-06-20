source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")

set.seed(15)

N <- 50

part_1 = .9^(1:N) * cos(pi * 1:N / 5)
part_2 = .2 * 1.05^(1:N) * cos(pi * 1:N / 12 + pi/4)

signal <- part_1 + part_2
r <- 4
ar <- .9
magnitude_wogap <- 0.15
magnitude_gap <- 0.15
skip_indices <- c(10:19, 35:39)


noise <- as.numeric(arima.sim(n = N, list(ar = ar), sd = 1))
series <- signal + magnitude_wogap * sqrt(sum(signal ^ 2)) *  noise / sqrt(sum(noise ^ 2))


v_init <- svd(traj_matrix(series, r + 1))$u[, r + 1]
weights <- bandSparse(length(series), length(series), 0:length(ar), inv_ac_diags(length(series), ar), symmetric = TRUE)

answer <- run_hlra(series = series, 
                   v_init = v_init,
                   it = 100,
                   objective = NULL,
                   opt_method = TRUE,
                   weights = weights)

# pdf("model_wo_gap_rn.pdf", width = 4.2, height = 2.1, pointsize = 4)
# par(mar = c(4.6,3.9,1.2,1.2))

matplot(1:N, cbind(series, answer, signal), pch = c(1, 26, 26),
        xlab = "Index", ylab = "Value")
lines(1:N, answer, col = "red", lty = 1)
lines(1:N, signal, col = "blue", lty = 2)
legend("top", c("Series", "Approximation", "Signal"), col=c("black", "red", "blue"),
       pch = c(1, 26, 26), lty = c(0, 1, 2))

# dev.off()

# construct time series with gaps
series <- signal + magnitude_gap * sqrt(sum(signal ^ 2)) *  noise / sqrt(sum(noise ^ 2))
series[skip_indices] <- NA
mask <- 1 - is.na(series)

weights_gap <- weights
weights_gap[skip_indices, ] <- 0
weights_gap[, skip_indices] <- 0

# fill series gaps with mean
series_fill <- series
series_fill[!mask] <- mean(series, na.rm = TRUE)

# v_init <- svd(traj_matrix(series_fill, r + 1))$u[, r + 1]
# take initial GLRR from SVD of trajectory matrix of series

answer <- run_hlra(series = series_fill, 
              v_init = v_init,
              it = 100,
              objective = NULL,
              opt_method = TRUE,
              weights = weights_gap)

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

# pdf("model_gap_rn.pdf", width = 4.2, height = 2.1, pointsize = 4)
# par(mar = c(4.6,3.9,1.2,1.2))

matplot(1:N, cbind(series, answer, signal, series_fill_to_plot), pch = c(1, 26, 26, 8),
        xlab = "Index", ylab = "Value", col = 1)
lines(1:N, answer, col = "red", lty = 1)
lines(1:N, signal, col = "blue", lty = 2)
legend("top", c("Series", "Approximation", "Signal", "Missing value"), 
       col=c("black", "red", "blue", "black"),
       pch = c(1, 26, 26, 8), lty = c(0, 1, 2, 0))

# dev.off()