@testset "cc-acopf" begin
for Gamma_type = [:y, :d]
    for xtilde = [true, false]
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
        options[:gamma]          = 1.0
        options[:relax_Gamma]    = false
        options[:Gamma_type]     = Gamma_type
        options[:xtilde]         = xtilde
        options[:print_level]    = 5

        cm = OPF.cc_acopf_model(opfdata, options, data)
        cm = OPF.cc_acopf_solve(cm, opfdata, options)

        if Gamma_type == :y
            @test norm(getvalue(cm.m[:Gamma]) - (getvalue(cm.m[:dF_dx]) \ -getvalue(cm.m[:dF_dy]))) <= tol
        else
            @test norm(getvalue(cm.m[:Gamma]) - (getvalue(cm.m[:dF_dx]) \ -getvalue(cm.m[:dF_dd]))) <= tol
        end
    end
end
end # testset
