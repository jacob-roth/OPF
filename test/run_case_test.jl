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

# case_name = "case118"
# case_path = cases_path
case_name = "case1354pegase"
# case_path = cases_path * "1354_id/"
# case_path = cases_path * "1354_structjump/"
case_path = cases_path * "1354_pglib/"
# case_name = "case9241pegase"
# case_path = cases_path * "9241_id/"
# case_path = cases_path * "1354_structjump/"
# case_path = cases_path * "1354_pglib/"

## infeasible (~491 iterations)
options = DefaultOptions()
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = false
options[:remove_tap]     = false
options[:print_level]    = 5
options[:shed_load]      = false

opfdata = load_case(case_name, case_path; other=false);
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodeldata = get_opfmodeldata(opfdata,options)
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)

## feasible (~18 iterations)
options = DefaultOptions()
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = false
options[:remove_tap]     = true
options[:print_level]    = 5
options[:shed_load]      = false

opfdata = load_case(case_name, case_path; other=false);
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodeldata = get_opfmodeldata(opfdata,options)
opfmodel = acopf_model(opfdata, options)
opfmodel_acopf = acopf_solve(opfmodel, opfdata)

## 118bus n-1
options[:current_rating] = true
options[:ramp_pct] = 0.01
opfdata = load_case(case_name, case_path; other=false);
casedata = CaseData(opfdata, get_phys(opfdata.buses, Dg=1.0, Mg=1.0, Dl=1.0, Dv=1.0));
opfmodeldata = get_opfmodeldata(opfdata,options)
contingencies = get_all_contingencies(opfdata, options)
scopfmodel = scacopf_model(opfdata, options, DefaultAdjustments(), contingencies)
scopfmodel_acopf = scacopf_solve(scopfmodel, opfdata, options, contingencies)
