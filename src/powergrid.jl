using DelimitedFiles
using PowerDynamics

powergrid_path = "/"
powergrid = read_powergrid(powergrid_path * "ieee118_v3.json", Json)

# Check branch data was entered correctly
branch_idx = [(l.from, l.to) for l in powergrid.lines]
admittance = [(l.Y.re, l.Y.im) for l in powergrid.lines]

casefile_path = "/"
branch_file = readdlm(casefile_path * ".branch")
branch_line_idx = Int.(branch_file[:,1:2])
deduped_branches = Set([branch_line_idx[row_idx,1] < branch_line_idx[row_idx,2] ? (branch_line_idx[row_idx,1], branch_line_idx[row_idx,2]) : (branch_line_idx[row_idx,2], branch_line_idx[row_idx,1]) 
                        for row_idx in 1:size(branch_line_idx,1)])
unpacked_deduped_branches = [(branch[1], branch[2]) for branch in deduped_branches]

@assert length(branch_idx) == length(unpacked_deduped_branches)
@assert sort(branch_idx) == sort(unpacked_deduped_branches)

admittance_file = readdlm(casefile_path * ".Y", '\t', ComplexF64)
Y_vec = [(real(admittance_file[branch_idx[line_idx][1], branch_idx[line_idx][2]]), imag(admittance_file[branch_idx[line_idx][1], branch_idx[line_idx][2]])) 
         for line_idx in 1:length(branch_idx)]

@assert admittance == Y_vec

# Check bus data was entered correctly
bus_file = readdlm(casefile_path * ".bus")
bus_type_vals = bus_file[:,2]
bus_P_vals = bus_file[:,3] / 100
bus_Q_vals = bus_file[:,4] / 100

M = 0.0531
Ω = 2 * π * 50
H = M * Ω / 2
D = 0.05

nodes = powergrid.nodes
for idx in 1:length(nodes)
    try 
        node = nodes[idx]
        if isa(node, SwingEq)
            @assert bus_type_vals[idx] == 2
            @assert node.P == -bus_P_vals[idx]
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

# Get counts of bus types
nodes = powergrid.nodes
bus_array = collect(values(nodes))
type_dict = Dict{Type{<:AbstractNode}, Int}()
type_dict[SwingEq] = 0
type_dict[PQAlgebraic] = 0
type_dict[SlackAlgebraic] = 0

for bus in bus_array
    bus_type = typeof(bus)
    type_dict[bus_type] += 1
end

# symbolsof(pq) == [:u_r, :u_i]
# symbolsof(swing) == [:u_r, :u_i, :ω]
# symbolsof(slack) == [:u_r, :u_i]


# Import operatingpoint from FPACOPF into PowerDynamics framework
operatingdata_path = "/"
Va = reshape(readdlm(operatingdata_path * "Va.csv"), 118)
Vm = reshape(readdlm(operatingdata_path * "Vm.csv"), 118)

u_r = Vm .* cos.(Va)
u_i = Vm .* sin.(Vm)
ω = 0

num_gens = type_dict[SwingEq]
num_loads = type_dict[PQAlgebraic]
num_slack = type_dict[SlackAlgebraic]
num_vars = (num_gens * 3) + (num_loads * 2) + (num_slack * 2) 

fpacopf_op = Array{Float64,1}(undef, num_vars)

vec_idx = 1
for node_idx in 1:length(nodes)
    node = nodes[node_idx]
    if isa(node, SwingEq)
        fpacopf_op[vec_idx] = u_r[node_idx]
        fpacopf_op[vec_idx+1] = u_i[node_idx]
        fpacopf_op[vec_idx+2] = ω
        vec_idx += 3
    else
        fpacopf_op[vec_idx] = u_r[node_idx]
        fpacopf_op[vec_idx+1] = u_i[node_idx]
        vec_idx += 2
    end
end

# FPACOPF is not a root it seems
# du = similar(fpacopf_op)
# rpg = rhs(powergrid)
# rpg(du, fpacopf_op, nothing, 0.0)
# if maximum(abs.(du)) > 1e-4
#     @warn "The operationpoint search did not converge in a fixed point!"
# end

# Operationpoint from PowerDynamics
operationpoint = find_operationpoint(powergrid)

# Try using fpacopf_op as ic_guess
fpacopf_operationpoint = find_operationpoint(powergrid, fpacopf_op)

# Or we could just create the operationpoint using fpacopf_op as the .vec of State
fpacopf_operationpoint = State(powergrid, fpacopf_op)

# Try simulate a line fault
solution2 = simulate(LineFault(from=4,to=5), powergrid, fpacopf_operationpoint, timespan = (0.0,5.0))
# plot2 = create_plot(solution2)
# savefig(plot2, "ieee14-minimal-line-fault.pdf")

