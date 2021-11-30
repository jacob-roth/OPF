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

