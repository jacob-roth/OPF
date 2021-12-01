# Get counts of bus types
# symbolsof(pq) == [:u_r, :u_i]
# symbolsof(swing) == [:u_r, :u_i, :ω]
# symbolsof(slack) == [:u_r, :u_i]
function get_type_dict(powergrid::PowerGrid)
    nodes = powergrid.nodes
    bus_array = collect(values(nodes))
    bus_types = Set(typeof.(bus_array))
    type_dict = Dict{Type{<:AbstractNode}, Int}()
    for bus_type in bus_types
        type_dict[bus_type] = 0
    end
    for bus in bus_array
        bus_type = typeof(bus)
        type_dict[bus_type] += 1
    end
    return type_dict
end

# Get number of variables for the operatingpoint
function get_num_vars(type_dict::Dict{Type{<:AbstractNode}, Int})
    num_gens = type_dict[SwingEq]
    num_loads = type_dict[PQAlgebraic]
    num_slack = type_dict[SlackAlgebraic]
    num_vars = (num_gens * 3) + (num_loads * 2) + (num_slack * 2) 
    return num_vars
end

# Import operatingpoint from FPACOPF into PowerDynamics framework
function import_own_operatingpoint(powergrid::PowerGrid, operatingdata_path::String)
    nodes = powergrid.nodes
    num_buses = length(nodes)
    type_dict = get_type_dict(powergrid)
    num_vars = get_num_vars(type_dict)

    Va = reshape(readdlm(operatingdata_path * "Va.csv"), num_buses)
    Vm = reshape(readdlm(operatingdata_path * "Vm.csv"), num_buses)

    u_r = Vm .* cos.(Va)
    u_i = Vm .* sin.(Vm)
    ω = 0

    own_op = Array{Float64,1}(undef, num_vars)
    vec_idx = 1
    for node_idx = 1:length(nodes)
        node = nodes["bus$node_idx"]
        if isa(node, SwingEq)
            own_op[vec_idx] = u_r[node_idx]
            own_op[vec_idx+1] = u_i[node_idx]
            own_op[vec_idx+2] = ω
            vec_idx += 3
        else
            own_op[vec_idx] = u_r[node_idx]
            own_op[vec_idx+1] = u_i[node_idx]
            vec_idx += 2
        end
    end
    return own_op
end

# Import operatingpoint from FPACOPF into PowerDynamics framework
function import_own_operatingpoint(powergrid::PowerGrid, optimal_values::Dict)
    nodes = powergrid.nodes
    type_dict = get_type_dict(powergrid)
    num_vars = get_num_vars(type_dict)

    Va = optimal_values[:Va]
    Vm = optimal_values[:Vm]

    u_r = Vm .* cos.(Va)
    u_i = Vm .* sin.(Vm)
    ω = 0

    own_op = Array{Float64,1}(undef, num_vars)
    vec_idx = 1
    for node_idx = 1:length(nodes)
        node = nodes["bus$node_idx"]
        if isa(node, SwingEq)
            own_op[vec_idx] = u_r[node_idx]
            own_op[vec_idx+1] = u_i[node_idx]
            own_op[vec_idx+2] = ω
            vec_idx += 3
        else
            own_op[vec_idx] = u_r[node_idx]
            own_op[vec_idx+1] = u_i[node_idx]
            vec_idx += 2
        end
    end
    return own_op
end
