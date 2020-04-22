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

# set ratings for no line outages
opfdata.lines.rateA .= 0.0
point_0, M = get_scacopf_point(opfdata, options, adjustments, contingencies)
@assert(M.status == :Optimal)
flowmag2s_0 = get_flowmag2s(point_0, opfdata, options)
ratings_0   = get_ratings(flowmag2s_0.flowmag2, opfdata.baseMVA)
ratings     = deepcopy(ratings_0)
ratings     = max.(ratings, 1.0)
for c_id in keys(contingencies)
    point_c     = get_contingency_point(M, c_id)
    flowmag2s_c = get_flowmag2s(point_c, opfdata, options)
    ratings_c   = get_ratings(flowmag2s_c.flowmag2, opfdata.baseMVA)
    opfdata.lines.rateA = max.(max.(ratings, ratings_0), opfdata.lines.rateA)
    # global ratings = max.(ratings, ratings_0)
end
opfdata.lines.rateA .= ratings

solved = false; iter = 0
while solved != true && iter <= 10

    point_0, M  = get_scacopf_point(opfdata, options, adjustments, contingencies)
    if M.status == :Optimal; break; end
    flowmag2s_0 = get_flowmag2s(point_0, opfdata, options)
    ratings_0   = get_ratings(flowmag2s_0.flowmag2, opfdata.baseMVA)
    ratings     = deepcopy(ratings_0)
    viol_idx    = ratings .> opfdata.lines.rateA
    opfdata.lines.rateA[viol_idx]   .*= 1.5
    opfdata.lines.rateA .= max.(ratings .* 1.5, opfdata.lines.rateA)
    # opfdata.lines.rateA[viol_idx]   .*= 1.05
    # opfdata.lines.rateA[.!viol_idx] .*= 1.01
    # opfdata.lines.rateA .= max.(ratings .* 1.05, opfdata.lines.rateA)

    for c_id in keys(contingencies)
        opfd = deepcopy(opfdata)
        removed_line = remove_line!(opfd, c_id)
        point_c      = get_contingency_point(M, c_id)
        flowmag2s_c  = get_flowmag2s(point_c, opfd, options).flowmag2
        splice!(flowmag2s_c, c_id:c_id-1, 0.0)
        ratings_c    = get_ratings(flowmag2s_c, opfd.baseMVA)
        viol_idx     = ratings_c .> opfdata.lines.rateA
        opfdata.lines.rateA[viol_idx]   .*= 1.5
        opfdata.lines.rateA .= max.(ratings .* 1.5, opfdata.lines.rateA)
        # opfdata.lines.rateA[viol_idx]   .*= 1.05
        # opfdata.lines.rateA[.!viol_idx] .*= 1.01
        ratings      = max.(ratings_c, ratings_0)
        # opfdata.lines.rateA .= max.(ratings .* 1.05, opfdata.lines.rateA)
        reinstate_line!(opfd, c_id, removed_line)
    end

    # for c_id in keys(contingencies)
    #     point_c     = get_contingency_point(M, c_id)
    #     flowmag2s_c = get_flowmag2s(point_c, opfdata, options)
    #     ratings_c   = get_ratings(flowmag2s_c.flowmag2, opfdata.baseMVA)
    #     viol_idx    = ratings_c .> opfdata.lines.rateA
    #     println(findall(viol_idx))
    #     ratings_c[viol_idx] .*= 1.25
    #     println(ratings_c[viol_idx])
    #     ratings     = max.(ratings_c, ratings_0)
    # end

    opfdata.lines.rateA .= max.(ratings, opfdata.lines.rateA)
    println("rateA = ", opfdata.lines.rateA)
    global iter += 1
end

opfdata.lines.rateA


function diagnose_infeasibility()