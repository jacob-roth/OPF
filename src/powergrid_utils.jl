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

function get_cumulative_num_vars(powergrid::PowerGrid)
    nodes = powergrid.nodes
    bus_array = collect(values(nodes))
    bus_types = typeof.(bus_array)

    cumulative_num_vars = Int[]
    var_idx = 1
    for bus_type in bus_types
        push!(cumulative_num_vars, var_idx)
        if bus_type == SwingEq
            var_idx += 3
        else
            var_idx += 2
        end
    end
    return cumulative_num_vars
end

# Import operatingpoint from FPACOPF into PowerDynamics framework
function import_own_operatingpoint(powergrid::PowerGrid, optimal_values::Dict)
    nodes = powergrid.nodes
    type_dict = get_type_dict(powergrid)
    num_vars = get_num_vars(type_dict)

    Va = optimal_values[:Va]
    Vm = optimal_values[:Vm]

    u_r = Vm .* cos.(Va)
    u_i = Vm .* sin.(Va)
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

function import_own_operatingpoint(powergrid::PowerGrid, operatingdata_path::String)
    nodes = powergrid.nodes
    num_buses = length(nodes)

    optimal_values = Dict()
    optimal_values[:Va] = reshape(readdlm(operatingdata_path * "Va.csv"), num_buses)
    optimal_values[:Vm] = reshape(readdlm(operatingdata_path * "Vm.csv"), num_buses)

    return import_own_operatingpoint(powergrid, optimal_values)
end

function err_from_root(powergrid::PowerGrid, operatingpoint_vec::AbstractArray)
    du = similar(operatingpoint_vec)
    rpg = rhs(powergrid)
    rpg(du, operatingpoint_vec, nothing, 0.0)
    return du
end

function err_from_root(operatingpoint::State)
    grid = operatingpoint.grid
    vec = operatingpoint.vec
    return err_from_root(grid, vec)
end

function import_initial_contingencies(ic_file_path::String)
    ic_ID_arr = readdlm(ic_file_path)
    num_ics = size(ic_ID_arr,1)
    ic_ID_vec = [(Int(ic_ID_arr[ic_idx,1]), Int(ic_ID_arr[ic_idx,2])) for ic_idx in 1:num_ics]
    return ic_ID_vec
end
