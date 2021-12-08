using DelimitedFiles
using PowerDynamics

powergrid_path = "/"
json_file = "ieee118_v9.json"
powergrid = read_powergrid(powergrid_path * json_file, Json)

# Check branch data
casefile_path = "/"
check_branch_idx(casefile_path, powergrid)
nobus_Bshunt = true
check_admittance_and_shunts(casefile_path, powergrid, nobus_Bshunt)

M = 0.0531
Ω = 2 * π * 50
H = M * Ω / 2
D = 0.05
Γ = 1e-6

# Check bus data
operatingdata_path = "/"
check_bus_data(casefile_path, operatingdata_path, powergrid, Ω, H, D, Γ)

import_path = "/"
imported_vec = import_own_operatingpoint(powergrid, import_path)
imported_operationpoint = State(powergrid, imported_vec)

timespan = (0.0, 101.0)
tspan_fault = (1.0, 100.0) 
timespan_cap = 50.0
resolution = 1e-3
tol = 0.01

ic_file_path = "/ic_IDs.txt"
ic_ID_vec = import_initial_contingencies(ic_file_path)

results_dict = simulate_initial_contingencies(ic_ID_vec[1:2], timespan, tspan_fault, powergrid, operationpoint, timespan_cap, resolution)
abs_delta_dict = get_abs_delta(results_dict)
all_big_delta_dict = findall_big_delta(abs_delta_dict, tol)
