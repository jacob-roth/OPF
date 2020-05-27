## -----------------------------------------------------------------------------
## constants
## -----------------------------------------------------------------------------
const cases_path = pwd() * "/cases/"
const tol = 1e-9

## -----------------------------------------------------------------------------
## environment
## -----------------------------------------------------------------------------
import Pkg
Pkg.activate("..")
Pkg.instantiate()
using Test
using MPCCases, Printf
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra
using NLsolve
include("../src/OPF.jl") # stand-in for "using OPF"; remember to comment out the `module` calls

## -----------------------------------------------------------------------------
## test cases
## -----------------------------------------------------------------------------

case_name = "case1354pegase"
case_path = cases_path * "1354_id/"
# case_path = cases_path * "1354_structjump/"
# case_path = cases_path * "1354_pglib/"

options = DefaultOptions()
options[:current_rating] = false
options[:lossless]       = false
options[:remove_Bshunt]  = false
options[:remove_tap]     = false
options[:print_level]    = 5

options[:current_rating] = false
opfdata = load_case(case_name, case_path; other=false);
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodeldata = get_opfmodeldata(opfdata,options)
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)

options[:current_rating] = true
opfdata = load_case(case_name, case_path; other=false);
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodeldata = get_opfmodeldata(opfdata,options)
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)

options[:current_rating] = true
options[:lossless] = true
opfdata = load_case(case_name, case_path; other=false);
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodeldata = get_opfmodeldata(opfdata,options)
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)

options[:current_rating] = true
options[:lossless] = true
options[:remove_Bshunt] = true
options[:remove_tap] = true
opfdata = load_case(case_name, case_path; other=false);
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodeldata = get_opfmodeldata(opfdata,options)
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)











opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)
point_power = get_point(opfmodel_acopf)


point_pm = Dict()
point_pm[:Vm] = Vm_pm = vec(readdlm("/Users/jakeroth/Desktop/opensource/PowerModels.jl/Vm.csv"))
point_pm[:Va] = Va_pm = vec(readdlm("/Users/jakeroth/Desktop/opensource/PowerModels.jl/Va.csv"))
point_pm[:Pg] = Pg_pm = vec(readdlm("/Users/jakeroth/Desktop/opensource/PowerModels.jl/Pg.csv"))
point_pm[:Qg] = Qg_pm = vec(readdlm("/Users/jakeroth/Desktop/opensource/PowerModels.jl/Qg.csv"))
point_pm[:curr_fr_2] = curr_fr_2_pm = vec(readdlm("/Users/jakeroth/Desktop/opensource/PowerModels.jl/curr_fr_2.csv"))
point_pm[:curr_to_2] = curr_to_2_pm = vec(readdlm("/Users/jakeroth/Desktop/opensource/PowerModels.jl/curr_to_2.csv"))
point_pm[:Gbus] = G_pm = readdlm("/Users/jakeroth/Desktop/opensource/PowerModels.jl/Gbus.csv")
point_pm[:Bbus] = B_pm = readdlm("/Users/jakeroth/Desktop/opensource/PowerModels.jl/Bbus.csv")
Ybus = G_pm - im*B_pm
Ybus - Matrix(opfmodeldata[:Y])

check_feasibility(point_pm, opfdata, options)
curr2_Yl = get_flowmag2s(Vm_pm, Va_pm, opfmodeldata[:Y], opfdata, options, :Yl).flowmag2
curr2_Yft = get_flowmag2s(Vm_pm, Va_pm, opfmodeldata[:Y], opfdata, options, :Yft).flowmag2
idx_fr = findall(point_pm[:curr_fr_2]-curr2 .>= 1e-7)
idx_to = findall(point_pm[:curr_to_2]-curr2 .>= 1e-7)

Yls = zeros(Complex,length(opfdata.lines))
Yfts = zeros(Complex,length(opfdata.lines))
for i in eachindex(opfdata.lines)
    line = opfdata.lines[i]
    Yls[i] = line.r / (line.r^2 + line.x^2) - im * (line.x / (line.r^2 + line.x^2))
    Yfts[i] = opfmodeldata[:Y][line.to,line.from]
end
Yls2 = abs2.(Yls)
Yfts2 = abs2.(Yfts)

i=1988

lines = opfdata.lines
Yl = lines.r ./ (lines.r.^2 .+ lines.x.^2) .- im .* (lines.x ./ (lines.r.^2 .+ lines.x.^2))

options[:current_rating] = true
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata, point_power)
options[:feasibility]=true
opfmodel_acopf = acopf_solve(opfmodel, opfdata)
point_curr = get_point(opfmodel_acopf)

c = get_flowmag2s(point_curr[:Vm], point_curr[:Va], opfmodeldata[:Y], opfdata, options, :Yl).flowmag2




point_power = get_point(opfmodel_acopf)
opfmodeldata = get_opfmodeldata(opfdata, options)
aps = get_apparentpower(point_power[:Vm], point_power[:Va], opfmodeldata)
cur = get_flowmag2s(point_power[:Vm], point_power[:Va], opfmodeldata[:Y], opfdata, options).flowmag2
flo = get_PM_flow(point_power[:Vm], point_power[:Va], opfmodeldata)


options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = false
options[:remove_tap]     = false
opfdata = load_case(case_name, case_path; other=false);
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)
point_curr  = get_point(opfmodel_acopf)

opfdata = load_case(case_name, case_path; other=false);
opfdata.lines.rateA .*= 2.5
opfdata.lines.rateA .*= 0
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)


opfdata = load_case(case_name, case_path; other=false);
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)

function adjust_ratings!(opfdata, options, feas_tol::Float64=1e-6, scale::Float64=1.5)
    opfmodel = acopf_model(opfdata, options)
    opfmodel_acopf = acopf_solve(opfmodel, opfdata)

    point = get_point(opfmodel_acopf)
    PG = point[:Pg]
    QG = point[:Qg]
    VM = point[:Vm]
    VA = point[:Va]
    Y  = computeAdmittanceMatrix(opfdata, options)
    flowmag2s = get_flowmag2s(VM, VA, Y, opfdata, options)
    flowmax   = (opfdata.lines.rateA ./ opfdata.baseMVA).^2
    flowmax_adj = deepcopy(flowmax)
    flowmax_adj[flowmax_adj .== 0] .= Inf
    flows = flowmag2s.flowmag2 .<= flowmax_adj .+ feas_tol
    offend = findall(flows .!= 1)
    flowmax[offend] .= flowmag2s.flowmag2[offend] .* scale
    opfdata.lines.rateA .= sqrt.(flowmax) .* opfdata.baseMVA
end
opfdata = load_case(case_name, case_path; other=false);
adjust_ratings!(opfdata, options)
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)