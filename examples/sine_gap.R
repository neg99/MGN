source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")

set.seed(15)

N <- 50
Ng <- 10

signal <- 10 * sin(2 * pi * 1:N / 6)
series <- c(signal, numeric(Ng)) + c(rnorm(N), numeric(Ng))
r <- 3

mask <- c(rep(1, N), rep(0, Ng))
weights <- Diagonal(N+Ng, mask)

v_init <- svd(traj_matrix(signal, r + 1))$u[r + 1, ]

answer <- compare_steps(series = series, 
              v_init = v_init,
              it = 100,
              objective = NULL,
              opt_method = TRUE,
              weights = weights)

matplot(1:(N + Ng), cbind(answer, c(signal, numeric(Ng)), series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "S-MGN solution")


answer_mgn <- compare_steps(series = series, 
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        compensated = FALSE,
                        weights = weights)

matplot(1:(N + Ng), cbind(answer, c(signal, numeric(Ng)), series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "MGN solution")
