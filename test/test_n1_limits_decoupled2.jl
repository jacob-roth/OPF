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
include("../src/OPF_decoupled.jl") # stand-in for "using OPF"

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

lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
nbus = length(buses); nline = length(lines); ngen = length(generators)

viol_scale   = 1.0
nonviol_scale = 1.0

function update_limits!(opfdata, ratings, viol_scale::Float64=1.005, nonviol_scale::Float64=1.0)
    viol_idx = ratings .>= opfdata.lines.rateA
    opfdata.lines.rateA[viol_idx]   .= ratings[viol_idx]   .* viol_scale
    opfdata.lines.rateA[.!viol_idx] .= max.(ratings[.!viol_idx] .* nonviol_scale, opfdata.lines.rateA[.!viol_idx])
    println("highest limit on line $(argmax(opfdata.lines.rateA))")
end

fig, ax = subplots()

opfdata.lines.rateA .= 0.0
point_0, M = get_scacopf_point(opfdata, options, adjustments, contingencies)
@assert(M.status == :Optimal)
flowmag2s_0 = get_flowmag2s(point_0, opfdata, options)
ratings_0   = get_ratings(flowmag2s_0.flowmag2, opfdata.baseMVA)
ratings     = deepcopy(ratings_0)
opfdata.lines.rateA .= max.(ratings, 1.0)
ax[:plot](opfdata.lines.rateA, alpha=1.0, color="black")
for c_id in keys(contingencies)
    point_c     = get_contingency_point(M, c_id)
    flowmag2s_c = get_flowmag2s(point_c, opfdata, options)
    ratings_c   = get_ratings(flowmag2s_c.flowmag2, opfdata.baseMVA)
    update_limits!(opfdata, ratings_c, viol_scale, nonviol_scale)
    ax[:plot](opfdata.lines.rateA, alpha=1.0/c_id, color="black")
end

solved = false; iter = 0
while solved != true && iter <= 10

    point_0, M  = get_scacopf_point(opfdata, options, adjustments, contingencies)
    if M.status == :Optimal; break; end

    flowmag2s_0 = get_flowmag2s(point_0, opfdata, options)
    ratings_0   = get_ratings(flowmag2s_0.flowmag2, opfdata.baseMVA)
    ratings     = deepcopy(ratings_0)
    update_limits!(opfdata, ratings, viol_scale, nonviol_scale)

    for c_id in keys(contingencies)
        opfd = deepcopy(opfdata)
        removed_line = remove_line!(opfd, c_id)
        point_c      = get_contingency_point(M, c_id)
        flowmag2s_c  = get_flowmag2s(point_c, opfd, options).flowmag2; splice!(flowmag2s_c, c_id:c_id-1, 0.0)
        ratings_c    = get_ratings(flowmag2s_c, opfd.baseMVA)
        update_limits!(opfdata, ratings_c, viol_scale, nonviol_scale)
        reinstate_line!(opfd, c_id, removed_line)
    end

    println("rateA = ", opfdata.lines.rateA)
    global iter += 1
end

opfdata.lines.rateA
maximum(opfdata.lines.rateA)
minimum(opfdata.lines.rateA)