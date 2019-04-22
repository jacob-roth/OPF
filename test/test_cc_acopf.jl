@testset "cc-acopf" begin
## data
data = Dict()
data[:Sigma_d]           = Matrix(Diagonal(ones(2nbus)))
data[:Va_min]            = -pi * ones(nbus)
data[:Va_max]            =  pi * ones(nbus)

## options
options = Dict()
options[:lossless]       = false
options[:current_rating] = false
options[:epsilon_Vm]     = 0.0001
options[:epsilon_Va]     = 0.05
options[:epsilon_Qg]     = 0.05
options[:gamma]          = 1.0
options[:relax_Gamma]    = false
options[:print_level]    = 5

cm = OPF.cc_acopf_model(opfdata, options, data)
cm = OPF.cc_acopf_solve(cm, opfdata)

@test norm(getvalue(cm.m[:Gamma]) - (getvalue(cm.m[:dF_dx]) \ -getvalue(cm.m[:dF_dy]))) <= tol
end # testset

# getvalue(cm.m[:zeta])
# norm(getvalue(cm.m[:zeta]))

# cm_eval = setup(cm.m)
# m_idx = OPF.model_idx(opfdata)
# z_idx = OPF.om_z_idx(opfdata)
# J, JJ, dF = OPF.jac_z_alg(sm_zbar, Y, opfdata.BusIdx, opfdata.BusGenerators, z_idx, m_idx, false)
# Gamma = dF[:dF_dx] \ -dF[:dF_dy]
#
# err = 0
# for i = 1:nx
#     for j = 1:ny
#         global err += sum(dF[:dF_dx][i, k] * Gamma[k, j] for k = 1:nx) + dF[:dF_dy][i, j]
#     end
# end