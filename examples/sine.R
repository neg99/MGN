source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")

set.seed(15)

N <- 50
signal <- 10 * sin(2 * pi * 1:N / 6)
series <- signal + rnorm(N)
r <- 3


v_init <- svd(traj_matrix(signal, r + 1))$u[r + 1, ]

answer <- compare_steps(series = series, 
              v_init = v_init,
              it = 100,
              objective = NULL,
              opt_method = TRUE)

matplot(1:N, cbind(answer, signal, series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "S-MGN solution")


answer_mgn <- compare_steps(series = series, 
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        compensated = FALSE)

matplot(1:N, cbind(answer, signal, series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "MGN solution")


answer_vpgn <- compare_steps(series = series, 
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        step_search = "vp",
                        project_onto = project_onto_a_vp)

matplot(1:N, cbind(answer_vpgn, signal, series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "VPGN solution")

answer_svpgn <- compare_steps(series = series, 
                             v_init = v_init,
                             it = 100,
                             objective = NULL,
                             opt_method = TRUE,
                             step_search = "vp")

matplot(1:N, cbind(answer_svpgn, signal, series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "S-VPGN solution")
