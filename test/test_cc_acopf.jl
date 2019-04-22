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
options[:epsilon_Vm]     = 0.01
options[:epsilon_Va]     = 0.05
options[:epsilon_Qg]     = 0.05
options[:gamma]          = 10.0
options[:relax_Gamma]    = false
options[:print_level]    = 5

cm = OPF.cc_acopf_model(opfdata, options, data)
cm = OPF.cc_acopf_solve(cm, opfdata)

@test norm(getvalue(cm.m[:Gamma]) - (getvalue(cm.m[:dF_dx]) \ -getvalue(cm.m[:dF_dy]))) <= tol
end # testset
