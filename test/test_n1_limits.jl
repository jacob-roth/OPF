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
include("../src/OPF.jl") # stand-in for "using OPF"

## -----------------------------------------------------------------------------
## test
## -----------------------------------------------------------------------------
options = DefaultOptions()
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true
options[:print_level]    = 1

adjustments = DefaultAdjustments()

case_path     = "/Users/jakeroth/Desktop/opensource/OPF/test/cases/";
case_name     = "case30";
casedata      = load_case(case_name, case_path);
opfdata       = casedata.opf;
contingencies = get_all_contingencies(opfdata, options)
# ctg2          = get_nonislanding_contingencies(opfdata, options)
# c13           = Dict(13 => contingencies[13])
# c16           = Dict(16 => contingencies[16])
# c34           = Dict(34 => contingencies[34])

##
## test n-1
##

max_iter = 10
viol_scale = 1.0
nonviol_scale = 1.0
set_n1_limits!(opfdata, options, adjustments, max_iter, viol_scale, nonviol_scale)
