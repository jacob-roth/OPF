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
using MPCCases, Printf
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra, Distributions
using NLsolve
include("../src/OPF.jl") # stand-in for "using OPF"

## -----------------------------------------------------------------------------
## test cases
## -----------------------------------------------------------------------------

include("test_exitrates.jl")
