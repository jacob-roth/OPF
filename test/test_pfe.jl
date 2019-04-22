@testset "pfe" begin
for xtilde in [true, false]
    Y = computeAdmittanceMatrix(opfdata)
    Pnet, Qnet = OPF.PQnet(sm, opfdata)
    m_idx = OPF.model_idx(opfdata, xtilde)
    nbus = length(opfdata.buses)
    ngen = length(opfdata.generators)
    nload = sum(opfdata.buses.bustype .== 1)
    nx = length(m_idx[:x])
    ny = length(m_idx[:y])
    nd = length(m_idx[:d])

    @testset "vectorized" begin
    ## vectorized algebraic
    rP, rQ = OPF.PFE_vec(sm_zbar, Y, Pnet, Qnet, m_idx)
    @test norm(rP) <= tol*1000
    @test norm(rQ) <= tol*1000
    end # testset

    @testset "entrywise" begin
    ## entrywise algebraic
    rP, rQ = OPF.PFE(sm_zbar, Y, Pnet, Qnet, m_idx)
    @test norm(rP) <= tol*1000
    @test norm(rQ) <= tol*1000
    ## entrywise algebraic in place
    F = fill(NaN, nx)
    OPF.PFE!(F, sm_zbar, Y, Pnet, Qnet, m_idx)
    @test norm(F) <= tol*1000
    end # testset
end
end # pfe testset