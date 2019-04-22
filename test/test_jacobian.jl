@testset "jacobian" begin
Y = computeAdmittanceMatrix(opfdata)
xtilde = true
m_idx = OPF.model_idx(opfdata, xtilde)
data = Dict()
data[:Y]             = Y
data[:m_idx]         = m_idx
data[:opfdata]       = opfdata
data[:BusIdx]        = opfdata.BusIdx
data[:BusGenerators] = opfdata.BusGenerators
data[:model]         = sm_eval

## numerical
J_num = OPF.jac_z(sm_zbar, data, :num)[1:(2nbus), :]
## algebraic entrywise
J_alg_ew,_,_ = OPF.jac_z(sm_zbar, data, :ew)
## algebraic vectorized
J_alg_vec,_,_ = OPF.jac_z(sm_zbar, data, :vec)
## test
@test norm(J_num - J_alg_ew) <= tol
@test norm(J_num - J_alg_vec) <= tol
@test norm(J_alg_ew - J_alg_vec) <= tol

## compare to MatPower
@testset "matpower" begin
if case == "case30"
    dS_dVm, dS_dVa = OPF.dStilde_dVtilde(z, Y, m_idx)
    ## MatPower
    dS_dVa_ = matread("cases/SVa.mat")
    dS_dVm_ = matread("cases/SVm.mat")
    ## test
    @test norm(dS_dVm - dS_dVm_["SVm"]) <= tol
    @test norm(dS_dVa - dS_dVa_["SVa"]) <= tol
end
end # testset

end # testset