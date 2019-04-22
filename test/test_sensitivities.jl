@testset "sensitivities" begin
## parameters
epsilon = 1e-3
xtilde = false

## setup
Y = computeAdmittanceMatrix(opfdata)
m_idx = OPF.model_idx(opfdata, xtilde)
Pnet, Qnet = OPF.PQnet(sm, opfdata)

## data
data = Dict()
data[:Y]             = Y
data[:m_idx]         = m_idx
data[:opfdata]       = opfdata
data[:Pnet]          = Pnet
data[:Qnet]          = Qnet
data[:BusIdx]        = opfdata.BusIdx
data[:BusGenerators] = opfdata.BusGenerators

## algebraic
Γ_alg = OPF.get_Gamma(sm_zbar, data, :ew, epsilon)
## finite difference
Γ_fd  = OPF.get_Gamma(sm_zbar, data, :fd, epsilon)

## test
@test norm(Γ_alg - Γ_fd) <= epsilon
for i = 1:length(m_idx[:x])
    @test norm(Γ_alg[:,i] - Γ_fd[:,i]) <= epsilon
end
end # testset