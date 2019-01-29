source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")

set.seed(15)

N <- 10
Ng <- 1

signal <- c(numeric(N-1), 1)
series <- c(signal, numeric(Ng)) + c(numeric(N), numeric(Ng))
# construct time series with gap at its end
r <- 1

mask <- c(rep(1, N), rep(0, Ng))
weights <- Diagonal(N+Ng, mask)
# make appropriate weights matrix

v_init <- c(.5, -.5)
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

matplot(1:(N + Ng), cbind(answer, c(signal, numeric(Ng)), series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "S-MGN solution")


answer_mgn <- run_hlra(series = series, 
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        compensated = FALSE,
                        weights = weights)

matplot(1:(N + Ng), cbind(answer, c(signal, numeric(Ng)), series), type = "l",
        lty = 1, col = c("black", "red", "blue"),
        xlab = "index", ylab = "value", main = "MGN solution")


print(answer)
