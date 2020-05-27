## -----------------------------------------------------------------------------
## environment
## -----------------------------------------------------------------------------
using PowerModels
using Ipopt
using DelimitedFiles

## -----------------------------------------------------------------------------
## test cases
## -----------------------------------------------------------------------------
datapath = "/Users/jakeroth/Desktop/opensource/pglib-opf/pglib_opf_case1354_pegase.m"

## -----------------------------------------------------------------------------
## solve
## -----------------------------------------------------------------------------
network_data = PowerModels.parse_file(datapath);
pm_base = instantiate_model(network_data, ACPPowerModel, PowerModels.build_opf);
pm_curr = instantiate_model(network_data, IVRPowerModel, PowerModels.build_opf_iv);
result_base = optimize_model!(pm_base, optimizer=with_optimizer(Ipopt.Optimizer))
result_curr = optimize_model!(pm_curr, optimizer=with_optimizer(Ipopt.Optimizer))

# constraint_current_limit  # https://lanl-ansi.github.io/PowerModels.jl/stable/constraints/#PowerModels.constraint_current_limit

## is current model best solved as Vreal, Vimag not polar??
r = result_curr;

bus = r["solution"]["bus"];
bus_id_k = collect(keys(bus))
bus_id = parse.(Int64, bus_id_k)
σ = sortperm(bus_id)
bus_id_k = bus_id_k[σ]
bus_id = bus_id[σ]
bus_vr = [bus[k]["vr"] for k in bus_id_k]
bus_vi = [bus[k]["vi"] for k in bus_id_k]
bus_vm = abs.(bus_vr .+ im .* bus_vi)
bus_va = angle.(bus_vr .+ im .* bus_vi)

bus_map = Dict()
for i in eachindex(σ)
    v = i
    k = bus_id_k[i]
    bus_map[k] = v
    bus_map[v] = k
end

gen_id_k = collect(keys(network_data["gen"]))
gen_id = parse.(Int64, gen_id_k)
σ_g = sortperm(gen_id)
gen_id_k = gen_id_k[σ_g]
gen_id = gen_id[σ_g]

gen = r["solution"]["gen"];
gen_pg = [gen[id_k]["pg"] for id_k in gen_id_k]
gen_qg = [gen[id_k]["qg"] for id_k in gen_id_k]

# gen_map = Dict()
# for i in eachindex(gen_id_k)
#     k = gen_id_k[i]
#     v = string(network_data["gen"][k]["gen_bus"])
#     gen_map[k] = v
# end

# gen_pg = zeros(length(bus)) #[gen[k]["pg"] for k in gen_id_k]
# gen_qg = zeros(length(bus)) #[gen[k]["qg"] for k in gen_id_k]
# for id_k in gen_id_k
#     gen_pg[bus_map[gen_map[id_k]]] = gen[id_k]["pg"]
#     gen_qg[bus_map[gen_map[id_k]]] = gen[id_k]["qg"]
# end

branch = r["solution"]["branch"]
branch_id_k = collect(keys(branch))
branch_id = parse.(Int64, branch_id_k)
σ_branch = sortperm(branch_id)
branch_id_k = branch_id_k[σ_branch]
branch_id = branch_id[σ_branch]
branch_rateA = [network_data["branch"][id_k]["rate_a"] for id_k in branch_id_k]
branch_cr_fr = [branch[id_k]["cr_fr"] for id_k in branch_id_k]
branch_cr_to = [branch[id_k]["cr_to"] for id_k in branch_id_k]
branch_ci_fr = [branch[id_k]["ci_fr"] for id_k in branch_id_k]
branch_ci_to = [branch[id_k]["ci_to"] for id_k in branch_id_k]
branch_csi_fr = [branch[id_k]["csi_fr"] for id_k in branch_id_k]
branch_csr_fr = [branch[id_k]["csr_fr"] for id_k in branch_id_k]
# branch_csi_to = [branch[id_k]["csi_to"] for id_k in branch_id_k]
# branch_csr_to = [branch[id_k]["csr_to"] for id_k in branch_id_k]

# branch_c_fr_2 = (abs.(branch_cr_fr) + abs.(branch_cr_fr)).^2
# branch_c_to_2 = (abs.(branch_cr_to) + abs.(branch_ci_to)).^2
# branch_cs_fr_2 = (abs.(branch_csr_fr) + abs.(branch_csi_fr)).^2
branch_c_fr_2 = (abs2.(branch_cr_fr) + abs2.(branch_cr_fr))
branch_c_to_2 = (abs2.(branch_cr_to) + abs2.(branch_ci_to))
branch_cs_fr_2 = (abs2.(branch_csr_fr) + abs2.(branch_csi_fr))
branch_c_fr_2 .<= branch_rateA.^2


## Vm, Va, Pg, Qg
point = Dict()
point[:Vm] = bus_vm
point[:Va] = bus_va
point[:Pg] = gen_pg
point[:Qg] = gen_qg

Vm = point[:Vm]
Va = point[:Va]
Pg = point[:Pg]
Qg = point[:Qg]

writedlm("Vm.csv", Vm)
writedlm("Va.csv", Va)
writedlm("Pg.csv", Pg)
writedlm("Qg.csv", Qg)

curr_fr_2 = zeros(length(branch_id_k))
curr_to_2 = zeros(length(branch_id_k))
for i in eachindex(branch_id_k)
    id_k = branch_id_k[i]
    line = network_data["branch"][id_k]
    fr = bus_map[string(network_data["branch"][id_k]["f_bus"])]
    to = bus_map[string(network_data["branch"][id_k]["t_bus"])]
    curr2[i] = Vm[fr]^2 + Vm[to]^2


    Yl = line["br_r"] / (line["br_r"]^2 + line["br_x"]^2) - im * (line["br_x"] / (line["br_r"]^2 + line["br_x"]^2))
    tl = network_data["branch"][id_k]["tap"]
    ts = network_data["branch"][id_k]["shift"]
    Tl = tl * cos(ts) + im * tl * sin(ts)
    Vf = Vm[fr] * exp(im * Va[fr])
    Vt = Vm[to] * exp(im * Va[to])
    curr_fr_2[i] = abs2((Yl + im * 0.0/2) * Vf/abs2(Tl) - Yl * Vt/conj(Tl))
    curr_to_2[i] = abs2((Yl + im * 0.0/2) * Vt - Yl * Vf/Tl)
end

writedlm("curr_fr_2.csv", curr_fr_2)
writedlm("curr_to_2.csv", curr_to_2)
Y = Matrix(calc_admittance_matrix(network_data).matrix)
writedlm("Gbus.csv", real.(Y))
writedlm("Bbus.csv", imag.(Y))