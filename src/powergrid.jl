using DelimitedFiles
using PowerDynamics

powergrid_path = "/"
json_file = "ieee118_v4.json"
powergrid = read_powergrid(powergrid_path * json_file, Json)

# Check branch data was entered correctly
branch_idx = [(l.from, l.to) for l in powergrid.lines]
admittance = [(l.Y.re, l.Y.im) for l in powergrid.lines]

casefile_path = "/"
branch_file = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.branch")
branch_line_idx = Int.(branch_file[:,1:2])
deduped_branches = Set([branch_line_idx[row_idx,1] < branch_line_idx[row_idx,2] ? (branch_line_idx[row_idx,1], branch_line_idx[row_idx,2]) : (branch_line_idx[row_idx,2], branch_line_idx[row_idx,1]) 
                        for row_idx in 1:size(branch_line_idx,1)])
unpacked_deduped_branches = [(branch[1], branch[2]) for branch in deduped_branches]

@assert length(branch_idx) == length(unpacked_deduped_branches)
@assert sort(branch_idx) == sort(unpacked_deduped_branches)

admittance_file = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.Y", '\t', ComplexF64)
Y_vec = [(real(admittance_file[branch_idx[line_idx][1], branch_idx[line_idx][2]]), imag(admittance_file[branch_idx[line_idx][1], branch_idx[line_idx][2]])) 
         for line_idx in 1:length(branch_idx)]

@assert admittance == Y_vec

# Check bus data was entered correctly
bus_file = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.bus")
bus_type_vals = bus_file[:,2]
bus_P_vals = bus_file[:,3] / 100
bus_Q_vals = bus_file[:,4] / 100

gen_file = readdlm(casefile_path * "pglib/pglib_opf_case118_ieee_lowdamp.gen")
gen_ids = Int.(gen_file[:,1])
gen_P_vals = gen_file[:,2] / 100

M = 0.0531
Ω = 2 * π * 50
H = M * Ω / 2
D = 0.05

nodes = powergrid.nodes
for idx in 1:length(nodes)
    try 
        node = nodes[idx]
        if isa(node, SwingEq)
            gen_idx = findfirst(isequal(idx), gen_ids)
            @assert bus_type_vals[idx] == 2
            @assert abs(node.P - (gen_P_vals[gen_idx] - bus_P_vals[idx])) <= eps()
            @assert node.H == H
            @assert node.Ω == Ω
            @assert node.D == D
        elseif isa(node, PQAlgebraic)
            @assert bus_type_vals[idx] == 1
            @assert node.P == -bus_P_vals[idx]
            @assert node.Q == -bus_Q_vals[idx]
        end
    catch e
        if isa(e, AssertionError)
            println("idx:", idx)
        end
    end
end
