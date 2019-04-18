@testset "jacobian" begin
Y = computeAdmittanceMatrix(opfdata)
z_idx = OPF.om_z_idx(opfdata)
J_RGL_idx = OPF.om_jac_RGL_idx(opfdata)

data = Dict()
data[:Y]       = Y
data[:opfdata] = opfdata
data[:z_idx]   = z_idx
data[:model]   = sm_eval

## numerical (`opfmodel_z` order)
J_num = OPF.jac_z(sm_zbar, data, :num)[1:(2nbus), :]
## algebraic entrywise  (`opfmodel_z` order)
J_alg_ew = OPF.jac_z(sm_zbar, data, :ew)
## algebraic vectorized  (`opfmodel_z` order)
J_alg_vec = OPF.jac_z(sm_zbar, data, :vec)
## test
@test norm(J_num - J_alg_ew) <= tol
@test norm(J_num - J_alg_vec) <= tol
@test norm(J_alg_ew - J_alg_vec) <= tol
end # testset