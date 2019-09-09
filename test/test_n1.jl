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
set_acopf_n1_limits!(opfdata, options, adjustments)