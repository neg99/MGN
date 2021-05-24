source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")
library(rhlra)
library(Rssa)
data(co2)

rank <- 12

series <- as.numeric(co2)
L <- 234
result <- hlra_cadzow(series, rank, L = L, left_diags=inv_ac_diags(L, .9), right_diags=inv_ac_diags(L+1, numeric(0)))
F <- result$signal
plot(F, type="l")

print(svd(traj_matrix(F, rank + 1))$d)
print(svd(traj_matrix(series, rank + 1))$d)


v_init <- svd(traj_matrix(F, rank + 1))$u[, rank+1]

# v_init <- c(0.1767979, -0.3472488, 0.3255328, -0.2953839, 0.2646262, -0.2422574, 0.2341790, -0.2421206, 0.2642842, -0.2947319, 0.3243734, -0.3454727, 0.1774180)

print(v_init)
plot(v_init)
plot(abs(v_init), log="y")

Fssa <- ssa(F, L=rank+1, svd.method = "svd")
print(parestimate(Fssa, groups = list(1:rank)))

weights <- bandSparse(length(series), length(series), 0:1, inv_ac_diags(length(series), .9), symmetric = TRUE)
skip_indices <- 50:250
series_wo_gap <- series
series[skip_indices] <- 0
weights[skip_indices, ] <- 0
weights[, skip_indices] <- 0

weights_for_vpgn <- weights
weights_for_vpgn[cbind(skip_indices, skip_indices)] <- 1e-10

answer_mgn <- run_hlra(series = series,
                       v_init = v_init,
                       it = 100,
                       objective = NULL,
                       opt_method = TRUE,
                       weights = weights,
                       compensated = FALSE)


#v_init <- c(1, numeric(11), -1)

answer_mgn_h <- run_hlra(series = series,
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        weights = weights,
                        compensated = TRUE)

mgn_final_point <- cur_v_global

#v_init <- c(1, numeric(11), -1)

answer_vpgn <- run_hlra(series = series,
                       v_init = v_init,
                       it = 100,
                       objective = NULL,
                       opt_method = TRUE,
                       compensated = FALSE,
                       weights = weights_for_vpgn,
                       step_search = "vp",
                       project_onto = project_onto_a_vp)

#v_init <- c(1, numeric(11), -1)

answer_vpgn_h <- run_hlra(series = series,
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        compensated = TRUE,
                        weights = weights_for_vpgn,
                        step_search = "vp")

vpgn_final_point <- cur_v_global

series[skip_indices] <- NA
series_gap <- series_wo_gap
series_gap[-skip_indices] <- NA

matplot(1:length(series), cbind(series, answer_mgn_h, answer_vpgn_h), type = "l",
        lty = 1:3, col = c("blue", "black", "red"),
        xlab = "index", ylab = "value")
legend("topleft", c("Series", "MGN-H", "VPGN-H"), col=c("blue", "black", "red"), lty=1:3)

matplot(1:length(series), cbind(series, answer_mgn, answer_vpgn), type = "l",
        lty = 1:3, col = c("blue", "black", "red"),
        xlab = "index", ylab = "value")
legend("topleft", c("Series", "MGN", "VPGN"), col=c("blue", "black", "red"), lty=1:3)

matplot(1:length(series), cbind(series, answer_mgn_h, answer_vpgn_h, series_gap), type = "l",
        lty = 1:3, col = c("blue", "black", "red", "orange"),
        xlab = "index", ylab = "value")
legend("topleft", c("Series", "MGN-H", "VPGN-H", "True value"), col=c("blue", "black", "red", "orange"), lty=1:3)

print(mean((answer_mgn[skip_indices] - series_wo_gap[skip_indices])^2))
print(mean((answer_mgn_h[skip_indices] - series_wo_gap[skip_indices])^2))
print(mean((answer_vpgn[skip_indices] - series_wo_gap[skip_indices])^2))
print(mean((answer_vpgn_h[skip_indices] - series_wo_gap[skip_indices])^2))

eigens_count <- 1:20
eigens_mgn <- svd(traj_matrix(answer_mgn, length(series) %/% 2))$d[eigens_count]
eigens_mgn_h <- svd(traj_matrix(answer_mgn_h, length(series) %/% 2))$d[eigens_count]
eigens_vpgn <- svd(traj_matrix(answer_vpgn, length(series) %/% 2))$d[eigens_count]
eigens_vpgn_h <- svd(traj_matrix(answer_vpgn_h, length(series) %/% 2))$d[eigens_count]

matplot(eigens_count, cbind(eigens_mgn, eigens_mgn_h, eigens_vpgn, eigens_vpgn_h), log="y", col=c("black", "black", "red", "red"),
        pch=1:4, xlab="Eigen number", ylab="Eigenvalue")
legend("topright", c("MGN", "MGN-H", "VPGN", "VPGN-H"), pch=1:4, col=c("black", "black", "red", "red"))
