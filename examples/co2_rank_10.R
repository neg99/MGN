source("../auxiliary/1_mgn.R")
source("../auxiliary/2_vp.R")
set.seed(1)
library(Rssa)
library(plotrix)
data(co2)

rank <- 10

series <- as.numeric(co2)
# Now make rank ... approximation using the Cadzow iterations
F <- cadzow(s.init, rank = rank, tol = 0, maxiter = 500)
plot(F, type="l")

print(svd(traj_matrix(F, rank + 1))$d)
print(svd(traj_matrix(series, rank + 1))$d)


v_init <- svd(traj_matrix(F, rank + 1))$u[, rank+1]
print(v_init)

roots <- polyroot(v_init)
plot(Re(roots), Im(roots), xlab = "Re", ylab = "Im")
draw.circle(0, 0, 1, lty = 2)


Fssa <- ssa(F, L=rank+1, svd.method = "svd")
print(parestimate(Fssa, groups = list(1:rank)))

pdf("co2.pdf", width = 2.1, height = 2.1, pointsize = 4)
par(mar = c(4.6,3.9,1.2,1.2))
plot(co2)
dev.off()

pdf("co2_roots.pdf", width = 2.1, height = 2.1, pointsize = 4)
par(mar = c(4.6,3.9,1.2,1.2))
plot(Re(roots), Im(roots), xlab = "Re", ylab = "Im", pch = 4, col = "red")
draw.circle(0, 0, 1, lty = 2)
dev.off()


stop()

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

mgn_final_point <- cur_v_global

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

pdf("co2_approx.pdf", width = 2.1, height = 2.1, pointsize = 4)
par(mar = c(4.6,3.9,1.2,1.2))
matplot(1:length(series), cbind(series, answer_mgn, answer_vpgn), type = "l",
        lty = 1:3, col = c("blue", "black", "red"),
        xlab = "index", ylab = "value")
legend("topleft", c("Series", "MGN", "VPGN"), col=c("blue", "black", "red"), lty=1:3)
dev.off()

print(mean((answer_mgn - series)^2))
print(mean((answer_mgn_h - series)^2))
print(mean((answer_vpgn - series)^2))
print(mean((answer_vpgn_h - series)^2))

eigens_count <- 1:15
eigens_mgn <- svd(traj_matrix(answer_mgn, length(series) %/% 2))$d[eigens_count]
eigens_mgn_h <- svd(traj_matrix(answer_mgn_h, length(series) %/% 2))$d[eigens_count]
eigens_vpgn <- svd(traj_matrix(answer_vpgn, length(series) %/% 2))$d[eigens_count]
eigens_vpgn_h <- svd(traj_matrix(answer_vpgn_h, length(series) %/% 2))$d[eigens_count]

matplot(eigens_count, cbind(eigens_mgn, eigens_mgn_h, eigens_vpgn, eigens_vpgn_h), log="y", col=c("black", "black", "red", "red"),
        pch=1:4, xlab="No. of singular number", ylab="Singular value")
legend("topright", c("MGN", "MGN-H", "VPGN", "S-VPGN-H"), pch=1:4, col=c("black", "black", "red", "red"))

pdf("co2_singular.pdf", width = 2.1, height = 2.1, pointsize = 4)
par(mar = c(4.6,3.9,1.2,1.2))

matplot(eigens_count, cbind(eigens_mgn, eigens_mgn_h, eigens_vpgn, eigens_vpgn_h), log="y", col=c("black", "black", "red", "red"),
        pch=1:4, xlab="No. of singular number", ylab="Singular value")
legend("bottomleft", c("MGN", "MGN-H", "VPGN", "S-VPGN-H"), pch=1:4, col=c("black", "black", "red", "red"))

dev.off()

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

