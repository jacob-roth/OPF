
j = 2
for i in 1:ny
    println(Gamma[i,j] - Gamma_FD[sigma[i],j])
end
println(argmax(Gamma[:,j] - Gamma_FD[sigma,j]))

[norm(Gamma[:,i] - Gamma_FD[:,i]) for i = 1:nx]

sigma = sortperm(data[:f_RGL_idx])
[norm(Gamma[:,i] - Gamma_FD[sigma,i]) for i = 1:nx]





### scratch ###
bus_idx, gen_idx = RGL_idx(opfdata)
bus_RGL = [bus_idx[:R]; bus_idx[:G]; bus_idx[:L]]
bus_RGL_inv = sortperm(bus_RGL)
G_RGL = G[bus_RGL, bus_RGL]
B_RGL = B[bus_RGL, bus_RGL]
vars, pars, eqns = jac_idx(opfdata)

## test PF equations consistency
Vtilde, Itilde, Stilde = PF_real(sxbar, G_RGL, B_RGL, pars)
Vtilde_, Itilde_, Stilde_ = PF(sxbar, Y, idx)
@test norm(real.(Vtilde_[bus_RGL]) - Vtilde.R) <= tol; @test norm(real.(Vtilde_) - Vtilde.R[bus_RGL_inv]) <= tol
@test norm(imag.(Vtilde_[bus_RGL]) - Vtilde.I) <= tol; @test norm(imag.(Vtilde_) - Vtilde.I[bus_RGL_inv]) <= tol
@test norm(real.(Itilde_[bus_RGL]) - Itilde.R) <= tol; @test norm(real.(Itilde_) - Itilde.R[bus_RGL_inv]) <= tol
@test norm(imag.(Itilde_[bus_RGL]) - Itilde.I) <= tol; @test norm(imag.(Itilde_) - Itilde.I[bus_RGL_inv]) <= tol
@test norm(real.(Stilde_[bus_RGL]) - Stilde.R) <= tol; @test norm(real.(Stilde_) - Stilde.R[bus_RGL_inv]) <= tol
@test norm(imag.(Stilde_[bus_RGL]) - Stilde.I) <= tol; @test norm(imag.(Stilde_) - Stilde.I[bus_RGL_inv]) <= tol
