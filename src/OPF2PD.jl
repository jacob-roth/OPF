const physDefault = Dict()
const physDefault[:Ω] = 2 * π * 50
const physDefault[:H] = 0.0531 * physDefault[:Ω] / 2 # M*Ω/2
const physDefault[:Dg] = 0.05
const physDefault[:Dv] = 0.005
# opf2pd("/Users/jakeroth/git/OPF/test/test.json", optimal_values, opfmodeldata, phys)
function opf2pd(fileout::String, optimal_values::Dict, opfmodeldata::Dict, line_type::String, 
                phys::Dict = physDefault)
  # setup
  @assert line_type ∈ Set(["StaticLine", "PiModelLine"])
  nbus = length(optimal_values[:Vm])
  nline = length(opfmodeldata[:lines])
  out = Dict()
  out["nodes"] = []
  out["lines"] = []
  out["version"] = "1"

  # buses
  for i = 1:nbus
    # common
    B = Dict()
    B["name"] = "bus" * string(i)
    B["params"] = Dict()
    B["params"]["Y_n"] = 0.0 # not sure why but https://github.com/jacob-roth/PowerDynamics.jl/blob/main/examples/ieee14bus/grid.json has this for all nodes

    # slack
    if opfmodeldata[:buses][i].bustype == 3
      V = optimal_values[:Vm][i]
      theta = optimal_values[:Va][i]
      B["params"]["U"] = real(V * exp(im * theta))
      B["type"] = "SlackAlgebraic"

    # generator
    elseif opfmodeldata[:buses][i].bustype == 2
      B["params"]["P"] = +1optimal_values[:Pnet][i]
      # B["params"]["Q"] = optimal_values[:Qnet][i]
      # B["params"]["V"] = optimal_values[:Vm][i]
      B["params"]["H"] = phys[:H]
      B["params"]["D"] = phys[:Dg]
      B["params"]["Ω"] = phys[:Ω]
      B["type"] = "SwingEq"

    # load
    elseif opfmodeldata[:buses][i].bustype == 1
      B["params"]["P"] = +1optimal_values[:Pnet][i]
      # B["params"]["V"] = optimal_values[:Vm][i]
      B["params"]["Q"] = +1optimal_values[:Qnet][i]
      B["type"] = "PQAlgebraic"
    end

    # add to `nodes` array
    push!(out["nodes"], B)
  end

  # lines
  from_to_pairs = Set{Tuple{Int, Int}}()
  for l = 1:nline
    f = opfmodeldata[:lines][l].from
    t = opfmodeldata[:lines][l].to
    ff = min(f,t) # need to have from.id <= to.id for PD lightgraphs dependency
    tt = max(f,t) # need to have from.id <= to.id for PD lightgraphs dependency

    if (ff,tt) ∈ from_to_pairs
      continue
    else
      new_pair = (ff,tt)
      push!(from_to_pairs, new_pair)
    end

    L = Dict()
    L["name"] = "branch" * string(l)
    L["params"] = Dict()
    L["params"]["from"] = "bus" * string(ff)
    L["params"]["to"] = "bus" * string(tt)
    L["type"] = line_type
    if line_type == "PiModelLine"
      L["params"]["y"] = Dict()
      L["params"]["y"]["re"] = 0.0
      L["params"]["y"]["im"] = -opfmodeldata[:Y][ff, tt]
      L["params"]["y_shunt_km"] = complex(0.0, opfmodeldata[:YshI][ff])
      L["params"]["y_shunt_mk"] = complex(0.0, opfmodeldata[:YshI][tt])
      # L["params"]["t_shunt_km"] =
      # L["params"]["t_shunt_mk"] = 
    else
      L["params"]["Y"] = Dict()
      L["params"]["Y"]["re"] = 0.0
      L["params"]["Y"]["im"] = -opfmodeldata[:Y][ff, tt]
    end

    # add to `nodes` array
    push!(out["lines"], L)
  end

  # write
  open(fileout, "w") do f
    JSON.print(f, out, 1)
  end
  # return JSON.print(out, 1)
end

function opf2pd(fileout::String, operatingdata_path::String, 
                opfmodeldata::Dict, line_type::String, 
                phys::Dict=physDefault)
  Y = opfmodeldata[:Y]
  if isa(Y, AbstractArray{<:Complex})
    opfmodeldata[:Y] = imag.(Y)
  end

  num_buses = opfmodeldata[:nbus]
  optimal_values = Dict()
  optimal_values[:Vm] = reshape(readdlm(operatingdata_path * "Vm.csv"), num_buses)
  optimal_values[:Va] = reshape(readdlm(operatingdata_path * "Va.csv"), num_buses)
  optimal_values[:Pnet] = reshape(readdlm(operatingdata_path * "Pnet.csv"), num_buses)
  optimal_values[:Qnet] = reshape(readdlm(operatingdata_path * "Qnet.csv"), num_buses)
  return opf2pd(fileout, optimal_values, opfmodeldata, line_type, phys)
end

function opf2pd(fileout::String, casedata_path::String, operatingdata_path::String, 
                opfmodeldata::Dict, line_type::String, 
                phys::Dict=physDefault)
  Y = opfmodeldata[:Y]
  if isa(Y, AbstractArray{<:Complex})
    opfmodeldata[:Y] = imag.(Y)
  end

  num_buses = opfmodeldata[:nbus]
  optimal_values = Dict()
  optimal_values[:Vm] = reshape(readdlm(operatingdata_path * "Vm.csv"), num_buses)
  optimal_values[:Va] = reshape(readdlm(operatingdata_path * "Va.csv"), num_buses)

  bus_file = readdlm(casedata_path * "mpc_lowdamp_pgliblimits.bus")
  gen_file = readdlm(casedata_path * "pglib_opf_case118_ieee_lowdamp.gen")
  baseMVA = 100.0

  bus_type = Int.(bus_file[:,2])
  Pd = bus_file[:,3] / baseMVA
  Qd = bus_file[:,4] / baseMVA
  gen_id = Int.(gen_file[:,1])
  Pg = gen_file[:,2] / baseMVA
  Qg = gen_file[:,3] / baseMVA

  Pnet = -Pd
  Qnet = -Qd
  for idx in 1:length(gen_id)
    if bus_type[gen_id[idx]] == 2
      Pnet[gen_id[idx]] += Pg[idx]
      Qnet[gen_id[idx]] += Qg[idx]
    end
  end
  optimal_values[:Pnet] = Pnet
  optimal_values[:Qnet] = Qnet
  return opf2pd(fileout, optimal_values, opfmodeldata, line_type, phys)
end
