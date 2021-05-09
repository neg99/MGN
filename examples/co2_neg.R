source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")
library(Rssa)
data(co2)

rank <- 12

series <- as.numeric(co2)
#series <- as.numeric(co2[1:246])
s.init <- ssa(series)
# Now make rank 3 approximation using the Cadzow iterations
F <- cadzow(s.init, rank = rank, tol = 1e-10, maxiter = 10)
s <- ssa(F, L = rank+1)
lr <- lrr(s, groups = list(1:rank))
#plot(s, type = "vectors")
#series <- as.numeric(co2[1:246])
# plot(series)

#v_init <- c(1, numeric(rank - 1), -1)
v_init <- c(lr, -1)

answer_mgn <- run_hlra(series = series,
                       v_init = v_init,
                       it = 100,
                       objective = NULL,
                       opt_method = TRUE,
                       compensated = FALSE)


#v_init <- c(1, numeric(11), -1)

answer_mgn_h <- run_hlra(series = series,
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        compensated = TRUE)

#v_init <- c(1, numeric(11), -1)

answer_vpgn <- run_hlra(series = series,
                       v_init = v_init,
                       it = 100,
                       objective = NULL,
                       opt_method = TRUE,
                       compensated = FALSE,
                       step_search = "vp",
                       project_onto = project_onto_a_vp)

#v_init <- c(1, numeric(11), -1)

answer_vpgn_h <- run_hlra(series = series,
                        v_init = v_init,
                        it = 100,
                        objective = NULL,
                        opt_method = TRUE,
                        compensated = TRUE,
                        step_search = "vp")

vpgn_final_point <- cur_v_global

matplot(1:length(series), cbind(series, answer_mgn_h, answer_vpgn_h), type = "l",
        lty = 1:3, col = c("blue", "black", "red"),
        xlab = "index", ylab = "value")
legend("topleft", c("Series", "MGN-H", "VPGN-H"), col=c("blue", "black", "red"), lty=1:3)

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
                             compensated = FALSE)

answer_mgn_h_new <- run_hlra(series = series,
                       v_init = vpgn_final_point, #point after VPGN-H
                       it = 100,
                       objective = NULL,
                       opt_method = TRUE,
                       compensated = TRUE)


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

