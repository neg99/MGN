source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")

set.seed(15)

N <- 50
signal <- 5 * sin(2 * pi * 1:N / 12)
series <- signal + runif(N)
# construct time series

r <- 2

weights <- band_mat_from_diags(inv_ac_diags(N, c(0.9)))

# make 3-diagonal inverted autocovariation matrix

v_init <- svd(traj_matrix(series, r + 1))$u[, r + 1]

# take initial GLRR from SVD of trajectory matrix of series

answer <- run_hlra(series = series, 
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

matplot(1:N, cbind(answer, signal, series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "S-MGN solution")


answer_mgn <- run_hlra(series = series, 
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        compensated = FALSE,
                        weights = weights)

matplot(1:N, cbind(answer, signal, series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "MGN solution")


answer_vpgn <- run_hlra(series = series, 
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        step_search = "vp",
                        project_onto = project_onto_a_vp,
                        weights = weights)

matplot(1:N, cbind(answer_vpgn, signal, series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "VPGN solution")

answer_svpgn <- run_hlra(series = series, 
                             v_init = v_init,
                             it = 100,
                             objective = NULL,
                             opt_method = TRUE,
                             step_search = "vp",
                             weights = weights)

matplot(1:N, cbind(answer_svpgn, signal, series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "S-VPGN solution")
