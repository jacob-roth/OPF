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
    v = Vm .* exp.(im .* Va)

    # u_r = Vm .* cos.(Va)
    # u_i = Vm .* sin.(Va)
    u_r = real.(v)
    u_i = imag.(v)
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

function vecidx2busvar(seek_vec_idx,operationpoint::State)
    powergrid = operationpoint.grid
    point = operationpoint.vec
    nodes = powergrid.nodes
    type_dict = get_type_dict(powergrid)
    num_vars = get_num_vars(type_dict)
    vec_idx = 0
    for node_idx = 1:length(nodes)
        node = nodes["bus$node_idx"]
        if isa(node, SwingEq)
            if seek_vec_idx == vec_idx+1
                return "bus$node_idx" * "_ur"
            elseif seek_vec_idx == vec_idx+2
                return "bus$node_idx" * "_ui"
            elseif seek_vec_idx == vec_idx+3
                return "bus$node_idx" * "_omega"
            end
            vec_idx += 3
        else
            if seek_vec_idx == vec_idx+1
                return "bus$node_idx" * "_ur"
            elseif seek_vec_idx == vec_idx+2
                return "bus$node_idx" * "_ui"
            end
            vec_idx += 2
        end
    end
end


function get_components(operationpoint::State)
    powergrid = operationpoint.grid
    point = operationpoint.vec
    nodes = powergrid.nodes
    type_dict = get_type_dict(powergrid)
    num_vars = get_num_vars(type_dict)

    # nswing = sum(typeof(nodes[k]) == SwingEq for k in keys(nodes))
    # npq = length(nodes) - ngen
    nbus = length(nodes)
    omega = zeros(nbus)
    ur = zeros(nbus)
    ui = zeros(nbus)
    vec_idx = 0
    for node_idx = 1:length(nodes)
        node = nodes["bus$node_idx"]
        if isa(node, SwingEq)
            ur[node_idx] = point[vec_idx+1]
            ui[node_idx] = point[vec_idx+2]
            omega[node_idx] = point[vec_idx+3]
            vec_idx += 3
        else
            ur[node_idx] = point[vec_idx+1]
            ui[node_idx] = point[vec_idx+2]
            vec_idx += 2
        end
    end
    return omega, ur, ui
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

function simulate_initial_contingencies(ic_ID_vec::AbstractArray{NTuple{N,Int}}, 
                                        timespan::Tuple{Real, Real}, tspan_fault::Tuple{Real, Real}, 
                                        powergrid::PowerGrid, operationpoint::State,
                                        timespan_cap::Real=20.0, resolution::Real=1e-3) where {N}
    results_dict = Dict()
    for ic_ID in ic_ID_vec
        branches = collect("branch" .* string.(ic_ID))
        fault = LineFailures(line_names=branches, tspan_fault=tspan_fault)
        fault_solution = simulate(fault, powergrid, operationpoint, timespan);
        generator_indices = findall(bus -> isa(bus, SwingEq), collect(values(powergrid.nodes)))
        generator_names = "bus" .* string.(generator_indices)
        
        results_dict[(ic_ID,:v)] = fault_solution(timespan[1]:resolution:timespan_cap, :, :v)[generator_indices,:]
        results_dict[(ic_ID,:p)] = fault_solution(timespan[1]:resolution:timespan_cap, :, :p)[generator_indices,:]
        results_dict[(ic_ID,:q)] = fault_solution(timespan[1]:resolution:timespan_cap, :, :q)[generator_indices,:]
        results_dict[(ic_ID,:ω)] = fault_solution(timespan[1]:resolution:timespan_cap, generator_names, :ω)
    end
    return results_dict
end

function get_abs_delta(results_dict::Dict)
    abs_delta_dict = Dict()
    for (ic_ID, variable) in keys(results_dict)
        results_arr = results_dict[(ic_ID, variable)]
        abs_delta_arr = abs.((results_arr[:, 2:end] - results_arr[:, 1:end-1]) ./ results_arr[:, 1:end-1])
        abs_delta_dict[(ic_ID, variable)] = abs_delta_arr
    end
    return abs_delta_dict
end

function findall_big_delta(abs_delta_dict::Dict, tol::Real=0.01)
    all_big_delta_dict = Dict()
    for (ic_ID, variable) in keys(abs_delta_dict)
        abs_delta_arr = abs_delta_dict[(ic_ID, variable)]
        num_rows = size(abs_delta_arr,1)
        all_big_delta_idx = Array{Any}(undef, num_rows)
        for row_idx in 1:num_rows
            row = abs_delta_arr[row_idx,:]
            all_big_delta_idx[row_idx] = findall(delta -> delta >= tol, row)
        end
        all_big_delta_dict[(ic_ID, variable)] = all_big_delta_idx
    end
    return all_big_delta_dict
end
