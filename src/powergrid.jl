using DelimitedFiles
using PowerDynamics

powergrid_path = "/"
json_file = "ieee118_v6.json"
powergrid = read_powergrid(powergrid_path * json_file, Json)

# Check branch data was entered correctly
lines_arr = collect(values(powergrid.lines))
line_types = Set([typeof(l) for l in lines_arr])
branch_idx = [(parse(Int, l.from[4:end]), parse(Int, l.to[4:end])) for l in lines_arr]
if (length(line_types) == 1) & (first(line_types) == PiModelLine)
    admittance = [(l.y.re, l.y.im) for l in lines_arr]
    shunt = [(l.y_shunt_km, l.y_shunt_mk) for l in lines_arr]
elseif (length(line_types) == 1) & (first(line_types) == StaticLine)
    admittance = [(l.Y.re, l.Y.im) for l in lines_arr]
end

casefile_path = "/"
branch_file = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.branch")
branch_line_idx = Int.(branch_file[:,1:2])
deduped_branches = Set([branch_line_idx[row_idx,1] < branch_line_idx[row_idx,2] ? (branch_line_idx[row_idx,1], branch_line_idx[row_idx,2]) : (branch_line_idx[row_idx,2], branch_line_idx[row_idx,1]) 
                        for row_idx in 1:size(branch_line_idx,1)])
unpacked_deduped_branches = [(branch[1], branch[2]) for branch in deduped_branches]

@assert length(branch_idx) == length(unpacked_deduped_branches)
@assert sort(branch_idx) == sort(unpacked_deduped_branches)

Y_arr = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.Y", '\t', ComplexF64)
Y_vec = [(real(Y_arr[line_pairs[1], line_pairs[2]]), -imag(Y_arr[line_pairs[1], line_pairs[2]])) for line_pairs in branch_idx]
@assert admittance == Y_vec

if (length(line_types) == 1) & (first(line_types) == PiModelLine)
    shunt_vec = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.YshI") 
    shunt_pairs = [(complex(0,shunt_vec[line_pairs[1]]), complex(0,shunt_vec[line_pairs[2]])) for line_pairs in branch_idx]
    @assert shunt == shunt_pairs
end

# Check bus data was entered correctly
operatingdata_path = "/"
nodes = powergrid.nodes
num_nodes = length(nodes)
Pnet = reshape(readdlm(operatingdata_path * "Pnet.csv"), num_nodes)
Qnet = reshape(readdlm(operatingdata_path * "Qnet.csv"), num_nodes)

bus_file = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.bus")
bus_type_vals = bus_file[:,2]

M = 0.0531
Ω = 2 * π * 50
H = M * Ω / 2
D = 0.05

for idx in 1:length(nodes)
    try 
        node = nodes[idx]
        if isa(node, SwingEq)
            @assert bus_type_vals[idx] == 2
            @assert node.P == Pnet[idx]
            @assert node.H == H
            @assert node.Ω == Ω
            @assert node.D == D
        elseif isa(node, PQAlgebraic)
            @assert bus_type_vals[idx] == 1
            @assert node.P == Pnet[idx]
            @assert node.Q == Qnet[idx]
        end
    catch e
        if isa(e, AssertionError)
            println("idx:", idx)
        end
    end
end
