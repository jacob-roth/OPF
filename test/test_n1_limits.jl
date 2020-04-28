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
using DelimitedFiles

## -----------------------------------------------------------------------------
## test
## -----------------------------------------------------------------------------
options = DefaultOptions()
options[:current_rating] = true
options[:lossless]       = true
# options[:remove_Bshunt]  = true
options[:remove_Bshunt]  = false
options[:remove_tap]     = true
options[:print_level]    = 5

adjustments = DefaultAdjustments()

case_path     = "/Users/jakeroth/Desktop/opensource/OPF/test/cases/";
# case_name     = "case30";
# case_name     = "case118";
case_name     = "case300";
phys = hcat(collect(1:300), collect(ones(300)), collect(ones(300)), collect(ones(300)))
writedlm("cases/case300.phys", phys)
casedata      = load_case(case_name, case_path);
opfdata       = casedata.opf;
contingencies = get_all_contingencies(opfdata, options)

options[:ctg_Vm_adj] = 0.5
options[:sc_total_obj] = true

c = Dict()
for i in 1:100
    c[i] = contingencies[i]
end
opfdata.lines.rateA .= 0.0
point_0, M = get_scacopf_point(opfdata, options, adjustments, c)

c = Dict()
for i in 101:200
    c[i] = contingencies[i]
end
opfmodeldata = get_opfmodeldata(opfdata, options)
opfdata.lines.rateA .= 0.0
point_0, M = get_scacopf_point(opfdata, options, adjustments, c)

c = Dict()
for i in 201:300
    c[i] = contingencies[i]
end
opfmodeldata = get_opfmodeldata(opfdata, options)
opfdata.lines.rateA .= 0.0
point_0, M = get_scacopf_point(opfdata, options, adjustments, c)

c = Dict()
for i in 301:400
    c[i] = contingencies[i]
end
opfmodeldata = get_opfmodeldata(opfdata, options)
opfdata.lines.rateA .= 0.0
point_0, M = get_scacopf_point(opfdata, options, adjustments, c)

c = Dict()
for i in 401:411
    c[i] = contingencies[i]
end
opfmodeldata = get_opfmodeldata(opfdata, options)
opfdata.lines.rateA .= 0.0
point_0, M = get_scacopf_point(opfdata, options, adjustments, c)

# c = Dict()
# for i in 101:150
#     c[i] = contingencies[i]
# end
# opfmodeldata = get_opfmodeldata(opfdata, options)
# opfdata.lines.rateA .= 0.0
# point_0, M = get_scacopf_point(opfdata, options, adjustments, c)
#
# c = Dict()
# for i in 101:125
#     c[i] = contingencies[i]
# end
# opfmodeldata = get_opfmodeldata(opfdata, options)
# opfdata.lines.rateA .= 0.0
# point_0, M = get_scacopf_point(opfdata, options, adjustments, c)
#
# c = Dict()
# for i in 111:115
#     c[i] = contingencies[i]
# end
# opfmodeldata = get_opfmodeldata(opfdata, options)
# opfdata.lines.rateA .= 0.0
# point_0, M = get_scacopf_point(opfdata, options, adjustments, c)
#
# options[:ramup]      = 100.0
# options[:ramdn]      = 0.0
# options[:ctg_Vm_adj] = 0.25
# c = Dict()
# for i in 114:114
#     c[i] = contingencies[i]
# end
# opfmodeldata = get_opfmodeldata(opfdata, options)
# opfdata.lines.rateA .= 0.0
# point_0, M = get_scacopf_point(opfdata, options, adjustments, c)
#
# opfd = deepcopy(opfdata)
# removed_line = remove_line!(opfd, 114)
# c_opfmodeldata = get_opfmodeldata(opfd, options)
# islanding_buses = get_islanding_buses(c_opfmodeldata, options)
#
# pc = get_contingency_point(M, 114)
#
# ctg2          = get_nonislanding_contingencies(opfdata, options)
# c13           = Dict(13 => contingencies[13])
# c16           = Dict(16 => contingencies[16])
# c34           = Dict(34 => contingencies[34])

##
## test n-1
##

max_iter = 10
viol_scale = 1.2
nonviol_scale = 1.05
set_n1_limits!(opfdata, options, adjustments, c, max_iter, viol_scale, nonviol_scale)