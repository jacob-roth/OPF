const physDefault = Dict()
const physDefault[:Ω] = 2 * π * 50
const physDefault[:H] = 0.0531 * physDefault[:Ω] / 2 # M*Ω/2
const physDefault[:Dg] = 0.05
const physDefault[:Dv] = 0.005
# opf2pd("/Users/jakeroth/git/OPF/test/test.json", optimal_values, opfmodeldata, phys)
function opf2pd(fileout::String, optimal_values::Dict, opfmodeldata::Dict, phys::Dict = physDefault)
  # setup
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
      B["params"]["U"] = V * exp(im * theta)
      B["type"] = "SlackAlgebraic"

    # generator
    elseif opfmodeldata[:buses][i].bustype == 2
      B["params"]["P"] = optimal_values[:Pnet][i]
      # B["params"]["Q"] = optimal_values[:Qnet][i]
      # B["params"]["V"] = optimal_values[:Vm][i]
      B["params"]["H"] = phys[:H]
      B["params"]["D"] = phys[:Dg]
      B["params"]["Ω"] = phys[:Ω]
      B["type"] = "SwingEq"

    # load
    elseif opfmodeldata[:buses][i].bustype == 1
      B["params"]["P"] = optimal_values[:Pnet][i]
      # B["params"]["V"] = optimal_values[:Vm][i]
      B["params"]["Q"] = optimal_values[:Qnet][i]
      B["type"] = "PQAlgebraic"
    end

    # add to `nodes` array
    push!(out["nodes"], B)
  end

  # lines
  for l = 1:nline
    f = opfmodeldata[:lines][l].from
    t = opfmodeldata[:lines][l].to
    ff = min(f,t) # need to have from.id <= to.id for PD lightgraphs dependency
    tt = max(f,t) # need to have from.id <= to.id for PD lightgraphs dependency
    L = Dict()
    L["name"] = "branch" * string(l)
    L["params"] = Dict()
    L["params"]["from"] = "bus" * string(ff)
    L["params"]["to"] = "bus" * string(tt)
    L["params"]["Y"] = Dict()
    L["params"]["Y"]["re"] = 0.0
    L["params"]["Y"]["im"] = opfmodeldata[:Y][ff, tt] # or should it be [f,t]?
    L["type"] = "StaticLine"
    # L["params"]["y"]["re"] = 0.0
    # L["params"]["y"]["im"] = opfmodeldata[:Y][f, t]
    # L["params"]["y_shunt_mk"] = 
    # L["params"]["y_shunt_km"] = 
    # L["params"]["t_shunt_mk"] = 
    # L["params"]["t_shunt_km"] = 
    # L["type"] = "PiModelLine"

    # add to `nodes` array
    push!(out["lines"], L)
  end

  # write
  open(fileout, "w") do f
    JSON.print(f, out, 1)
  end
  # return JSON.print(out, 1)
end
