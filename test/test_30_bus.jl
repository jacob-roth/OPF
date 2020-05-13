const path = pwd() * "/cases/"
const tol = 1e-9
const plotting = false

import Pkg; Pkg.activate(".."); Pkg.instantiate()
include("../src/OPF.jl")
using DelimitedFiles


options = DefaultOptions()
options[:constr_limit_scale] = 1.44
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
options[:shed_load]      = true
options[:ratelimit]      = 1.0
options[:print_level]    = 1
casedata = load_case("case30", path; other=true)
opfdata  = casedata.opf

for ratelimit in [Inf, 1e-0, 1e-1, 1e-2, 1e-3, 1e-4]
    options[:ratelimit] = ratelimit
    opfmodel = acopf_model(opfdata, options)
    opfmodel_exitrates = acopf_solve_exitrates(opfmodel, casedata, options)
end
