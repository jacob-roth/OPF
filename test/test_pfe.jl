@testset "pfe" begin
for full in [true, false]
    Y = computeAdmittanceMatrix(opfdata)
    z_idx = OPF.om_z_idx(opfdata)
    J_RGL_idx = OPF.om_jac_RGL_idx(opfdata)
    x_RGL_idx = OPF.om_x_RGL_idx(J_RGL_idx);       nx = length(x_RGL_idx)
    y_RGL_idx = OPF.om_y_RGL_idx(J_RGL_idx, full); ny = length(y_RGL_idx)
    f_RGL_idx = OPF.om_f_RGL_idx(J_RGL_idx, full); nf = length(f_RGL_idx)
    y_RGL = sm_zbar[y_RGL_idx]
    Pnet, Qnet = OPF.PQnet(sm, opfdata)
    ny = length(y_RGL)

    data = Dict()
    data[:Y]          = Y
    data[:opfdata]    = opfdata
    data[:z_idx]      = z_idx
    data[:model]      = sm_eval
    data[:y_RGL_idx]  = y_RGL_idx
    data[:f_RGL_idx]  = f_RGL_idx
    data[:opfmodel_z] = sm_zbar
    data[:Pnet]       = Pnet
    data[:Qnet]       = Qnet

    ## vectorized algebraic
    Vtilde, Itilde, Stilde = OPF.PF(y_RGL, data)  ## Vtilde, Itilde, Stilde = OPF.PF(sm_zbar, Y, z_idx)
    ## test PF
    F_P = real.(Stilde) - data[:Pnet]
    F_Q = imag.(Stilde) - data[:Qnet]
    @test norm(F_P) <= tol*1000
    @test norm(F_Q) <= tol*1000

    F = fill(NaN, ny)
    OPF.PFE_RGL!(F, y_RGL, data)
    @test norm(F) <= tol*1000

    ## compare to MatPower
    if case == "case30"
        dS_dVm, dS_dVa = OPF.dStilde_dVtilde(Vtilde, Y)
        ## MatPower
        dS_dVa_ = matread("cases/SVa.mat")
        dS_dVm_ = matread("cases/SVm.mat")
        ## test
        @test norm(dS_dVm - dS_dVm_["SVm"]) <= tol
        @test norm(dS_dVa - dS_dVa_["SVa"]) <= tol
    end
end
end # pfe testset