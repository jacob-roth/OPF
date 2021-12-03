using DelimitedFiles
using PowerDynamics

powergrid_path = "/"
json_file = "ieee118_v6.json"
powergrid = read_powergrid(powergrid_path * json_file, Json)

# Check branch data
casefile_path = "/"
check_branch_idx(casefile_path, powergrid)
check_admittance_and_shunts(casefile_path, powergrid)

M = 0.0531
Ω = 2 * π * 50
H = M * Ω / 2
D = 0.05

# Check bus data
operatingdata_path = "/"
check_bus_data(casefile_path, operatingdata_path, powergrid, Ω=Ω, H=H, D=D)
