function get_pg_branches(powergrid::PowerGrid)
    pg_lines = collect(values(powergrid.lines))
    pg_branches = [(parse(Int, l.from[4:end]), parse(Int, l.to[4:end])) for l in pg_lines]
    return pg_branches
end

function get_line_types(powergrid::PowerGrid)
    pg_lines = collect(values(powergrid.lines))
    line_types = Set([typeof(l) for l in pg_lines])
    return line_types
end

function get_pg_admittance_and_shunt(line_types::Set)
    pg_lines = collect(values(powergrid.lines))
    if (length(line_types) == 1) & (first(line_types) == PiModelLine)
        pg_admittance = [(l.y.re, l.y.im) for l in pg_lines]
        pg_shunt = [(l.y_shunt_km, l.y_shunt_mk) for l in pg_lines]
        return pg_admittance, pg_shunt
    elseif (length(line_types) == 1) & (first(line_types) == StaticLine)
        pg_admittance = [(l.Y.re, l.Y.im) for l in pg_lines]
        return pg_admittance
    end
end

# Check branch indices were entered correctly
function check_branch_idx(casefile_path::String, powergrid::PowerGrid)
    branch_file = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.branch")
    branches = Int.(branch_file[:,1:2])
    deduped_branches = Set([(minimum(branches[row_idx,:]), maximum(branches[row_idx,:])) for row_idx in 1:size(branches,1)])
    deduped_branches_vec = collect(deduped_branches)

    pg_branches = get_pg_branches(powergrid)
    @assert length(pg_branches) == length(deduped_branches_vec)
    @assert sort(pg_branches) == sort(deduped_branches_vec)
end

# Check branch admittance and shunts were entered correctly
function check_admittance_and_shunts(casefile_path::String, powergrid::PowerGrid)
    pg_branches = get_pg_branches(powergrid)
    Y_file = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.Y", '\t', ComplexF64)
    Y_vec = [(real(Y_file[line_pairs[1], line_pairs[2]]), -imag(Y_file[line_pairs[1], line_pairs[2]])) for line_pairs in pg_branches]

    line_types = get_line_types(powergrid)
    if (length(line_types) == 1) & (first(line_types) == PiModelLine)
        pg_admittance, pg_shunt = get_pg_admittance_and_shunt(line_types)
        y_shunt_file = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.y_shunt_arr", '\t', ComplexF64) 
        y_shunt_vec = [(y_shunt_file[line_pairs[1], line_pairs[2]], y_shunt_file[line_pairs[2], line_pairs[1]]) for line_pairs in pg_branches]
        @assert pg_admittance == Y_vec
        @assert pg_shunt == y_shunt_vec
    elseif (length(line_types) == 1) & (first(line_types) == StaticLine)
        pg_admittance = get_pg_admittance_and_shunt(line_types)
        @assert pg_admittance == Y_vec
    end
end

function get_bus_types(casefile_path::String)
    bus_file = readdlm(casefile_path * "mpc_lowdamp_pgliblimits.bus")
    bus_types = bus_file[:,2]
    return bus_types
end

function get_Pnet_and_Qnet(operatingdata_path::String, powergrid::PowerGrid)
    nodes = powergrid.nodes
    num_nodes = length(nodes)
    Pnet = reshape(readdlm(operatingdata_path * "Pnet.csv"), num_nodes)
    Qnet = reshape(readdlm(operatingdata_path * "Qnet.csv"), num_nodes)
    return Pnet, Qnet
end

# Check bus data was entered correctly
function check_bus_data(casefile_path::String, operatingdata_path::String, powergrid::PowerGrid; 
                        Ω::Real, H::Real, D::Real)
    bus_types = get_bus_types(casefile_path)
    Pnet, Qnet = get_Pnet_and_Qnet(operatingdata_path, powergrid)
    nodes = powergrid.nodes

    for idx in 1:length(nodes)
        try 
            node = nodes[idx]
            if isa(node, SwingEq)
                @assert bus_types[idx] == 2
                @assert node.P == Pnet[idx]
                @assert node.H == H
                @assert node.Ω == Ω
                @assert node.D == D
            elseif isa(node, PQAlgebraic)
                @assert bus_types[idx] == 1
                @assert node.P == Pnet[idx]
                @assert node.Q == Qnet[idx]
            end
        catch e
            if isa(e, AssertionError)
                println("idx:", idx)
            end
        end
    end
end
