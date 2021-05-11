source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")
library(rhlra)
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
                       weights = weights,
                       step_search = "vp",
                       project_onto = project_onto_a_vp)

#v_init <- c(1, numeric(11), -1)

answer_vpgn_h <- run_hlra(series = series,
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        compensated = TRUE,
                        weights = weights,
                        step_search = "vp")

vpgn_final_point <- cur_v_global

matplot(1:length(series), cbind(series, answer_mgn_h, answer_vpgn_h), type = "l",
        lty = 1:3, col = c("blue", "black", "red"),
        xlab = "index", ylab = "value")
legend("topleft", c("Series", "MGN-H", "VPGN-H"), col=c("blue", "black", "red"), lty=1:3)

matplot(1:length(series), cbind(series, answer_mgn, answer_vpgn), type = "l",
        lty = 1:3, col = c("blue", "black", "red"),
        xlab = "index", ylab = "value")
legend("topleft", c("Series", "MGN", "VPGN"), col=c("blue", "black", "red"), lty=1:3)

print(mean((answer_mgn - series)^2))
print(mean((answer_mgn_h - series)^2))
print(mean((answer_vpgn - series)^2))
print(mean((answer_vpgn_h - series)^2))

eigens_count <- 1:20
eigens_mgn <- svd(traj_matrix(answer_mgn, length(series) %/% 2))$d[eigens_count]
eigens_mgn_h <- svd(traj_matrix(answer_mgn_h, length(series) %/% 2))$d[eigens_count]
eigens_vpgn <- svd(traj_matrix(answer_vpgn, length(series) %/% 2))$d[eigens_count]
eigens_vpgn_h <- svd(traj_matrix(answer_vpgn_h, length(series) %/% 2))$d[eigens_count]

matplot(eigens_count, cbind(eigens_mgn, eigens_mgn_h, eigens_vpgn, eigens_vpgn_h), log="y", col=c("black", "black", "red", "red"),
        pch=1:4, xlab="Eigen number", ylab="Eigenvalue")
legend("topright", c("MGN", "MGN-H", "VPGN", "VPGN-H"), pch=1:4, col=c("black", "black", "red", "red"))

answer_mgn_new <- run_hlra(series = series,
                             v_init = vpgn_final_point, #point after VPGN-H
                             it = 100,
                             objective = NULL,
                             opt_method = TRUE,
                             compensated = FALSE,
                             weights = weights)

answer_mgn_h_new <- run_hlra(series = series,
                       v_init = vpgn_final_point, #point after VPGN-H
                       it = 100,
                       objective = NULL,
                       opt_method = TRUE,
                       compensated = TRUE,
                       weights = weights)


mean((answer_mgn_new - series)^2)
mean((answer_mgn_h_new - series)^2)

print(mean((answer_mgn - answer_mgn_h)^2))
print(mean((answer_vpgn - answer_vpgn_h)^2))
print(mean((answer_mgn_h - answer_vpgn_h)^2))

s.mgn <- ssa(as.vector(answer_mgn_h), L=rank+1, svd.method = "svd")
print(parestimate(s.mgn, groups = list(1:rank)))
s.mgn$sigma
rf.mgn <- rforecast(s.mgn, groups = list(1:rank), len =120, only.new = FALSE)
plot(rf.mgn, type = "l")

s.vpgn <- ssa(as.vector(answer_vpgn_h), L=rank+1, svd.method = "svd")
print(parestimate(s.vpgn, groups = list(1:rank)))
s.vpgn$sigma
rf.vpgn <- rforecast(s.vpgn, groups = list(1:rank), len =120, only.new = FALSE)
plot(rf.vpgn, type = "l")
plot(rf.mgn, type = "l")
lines(rf.vpgn, type = "l", col = "red")

matplot(1:(rank + 1), cbind(v_init/-v_init[rank+1], mgn_final_point/-mgn_final_point[rank+1], vpgn_final_point/-mgn_final_point[rank+1]),
        pch = 1:3)

print(v_init)
print(mgn_final_point)
print(vpgn_final_point)
