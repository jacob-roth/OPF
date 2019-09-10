## -----------------------------------------------------------------------------
## constants
## -----------------------------------------------------------------------------
const path = pwd() * "/cases/"
const tol = 1e-9
const plotting = false

## -----------------------------------------------------------------------------
## environment
## -----------------------------------------------------------------------------
import Pkg
Pkg.activate(dirname(dirname(dirname(path))))
Pkg.instantiate()
using Test
using MPCCases, Printf, MAT
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra, Distributions
using NLsolve
# using OPF
include("../src/OPF.jl")

################################################################################
case    = "case118"
opfdata = load_case(case, path, other=false);
nbus    = length(opfdata.buses);
ngen    = length(opfdata.generators);
nload   = sum(opfdata.buses.bustype .== 1);
options = DefaultOptions();
options[:current_rating] = true
adjustments = DefaultAdjustments()

# opfdata0 = deepcopy(opfdata);
opfdata  = deepcopy(opfdata0);
update_ratings_max!(opfdata, options)
max_ratings = deepcopy(opfdata.lines.rateA)

## update loading
update_loadings!(opfdata, options, 0.75, DefaultAdjPF())

## get contingencies
contingencies = get_contingencies(opfdata, options)
c = Dict()
for k in collect(keys(contingencies))
    c[k] = contingencies[k]
end

c_opfdata = OPFData[]
for k in collect(keys(contingencies))
    c = contingencies[k]
    rl = remove_line!(opfdata, c.asset.id);
    push!(c_opfdata, deepcopy(opfdata))
    reinstate_line!(opfdata, rl.id, rl);
end

## solve scacopf
scm = scacopf_model(opfdata, options, adjustments, c)
@time scm = scacopf_solve(scm, opfdata, options, c)

## get dispatch point
x = MathProgBase.getsolution(scm.m.internalModel)
dp = Dict()
dp[:Pg] = deepcopy(x[[x.col for x in getindex(scm.m, :Pg)]])
dp[:Qg] = deepcopy(x[[x.col for x in getindex(scm.m, :Qg)]])
dp[:Vm] = deepcopy(x[[x.col for x in getindex(scm.m, :Vm)]])
dp[:Va] = deepcopy(x[[x.col for x in getindex(scm.m, :Va)]])

## get operating points & flow limits
ops          = get_operating_points(scm, c)
flows        = [get_flowmag2s(op, c_data, options) for (op, c_data) in zip(ops, c_opfdata)]
base_ratings = [get_ratings(x.flowmag2, opfdata.baseMVA) for x in flows]
base_ratings = cat(base_ratings..., dims=2)
adj_ratings  = map(maximum, eachrow(base_ratings))

opfdata.lines.rateA .= adj_ratings .* 1.25
scm_new = scacopf_model(opfdata, options, adjustments, c)
@time scm_new = scacopf_solve(scm, opfdata, options, c)


n1_opfdata = get_n1_limits("case118", path)