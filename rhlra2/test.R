N <- 100
v <- c(1, -2, 1)

ordinary_basis_obj <- eval_basis_compensated(N, v)
matplot(1:N, ordinary_basis_obj$basis)

tan_basis <- eval_tangent_basis_compensated(N, v)
matplot(1:N, tan_basis)



series <- seq(-3, 5, length.out = N)
matplot(1:N, eval_pseudograd_compensated(N, v, series, 2, basisobj))

tan_basis <- get_comp_space_by_v(N = N, v = v, basis_use_c = FALSE, v_2 = TRUE)
matplot(1:N, tan_basis)
