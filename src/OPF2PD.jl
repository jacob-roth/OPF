const physDefault = Dict()
const physDefault[:Ω] = 2 * π * 50
const physDefault[:H] = 0.0531 * physDefault[:Ω] / 2 # M*Ω/2
const physDefault[:Dg] = 0.05
const physDefault[:Dv] = 0.005
const physDefault[:Γ] = 1e-6
# opf2pd("/Users/jakeroth/git/OPF/test/test.json", optimal_values, opfmodeldata, phys)
function opf2pd(fileout::String, optimal_values::Dict, opfmodeldata::Dict, line_type::String, bus_type::String,
                nobus_Bshunt::Bool=false, phys::Dict = physDefault)
  # setup
  @assert line_type ∈ Set(["StaticLine", "PiModelLine", "PiModelLineTapratio"])
  @assert bus_type ∈ Set(["SwingEq", "SwingEqLVS"])
  nbus = length(optimal_values[:Vm])
  nline = length(opfmodeldata[:lines])
  lines = opfmodeldata[:lines]
  buses = opfmodeldata[:buses]
  baseMVA = opfmodeldata[:baseMVA]
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
      if bus_type == "SwingEq"
        B["params"]["P"] = +1optimal_values[:Pnet][i]
        # B["params"]["Q"] = optimal_values[:Qnet][i]
        # B["params"]["V"] = optimal_values[:Vm][i]
        B["params"]["H"] = phys[:H]
        B["params"]["D"] = phys[:Dg]
        B["params"]["Ω"] = phys[:Ω]
        B["type"] = "SwingEq"
      elseif bus_type == "SwingEqLVS"
        B["params"]["P"] = +1optimal_values[:Pnet][i]
        B["params"]["H"] = phys[:H]
        B["params"]["D"] = phys[:Dg]
        B["params"]["Ω"] = phys[:Ω]
        B["params"]["Γ"] = phys[:Γ]
        B["params"]["V"] = optimal_values[:Vm][i] * exp(im * optimal_values[:Va][i])
        B["type"] = "SwingEqLVS"
      end

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

  ## best so far
  # lines
  if haskey(optimal_values, :ic_IDs)
    ic_IDs = optimal_values[:ic_IDs]
    println("Initial contingencies: ", ic_IDs)
  else
    ic_IDs = Set()
  end

  if haskey(optimal_values, :IDs) 
    failed_line = optimal_values[:IDs]
    println("Failed line: ", failed_line)
  else
    failed_line = Set()
  end

  from_to_pairs = Set{Tuple{Int, Int}}()
  for l = setdiff(1:nline, failed_line, ic_IDs)
    f = opfmodeldata[:lines][l].from
    t = opfmodeldata[:lines][l].to
    ff = min(f,t) # need to have from.id <= to.id for PD lightgraphs dependency
    tt = max(f,t) # need to have from.id <= to.id for PD lightgraphs dependency
    # good_order = (t == tt) * (f == ff)

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
      # y = Ys = 1/(lines[l].x*im)
      # y_shunt_km = im*(lines[l].b/2 + (buses[ff].Bs / baseMVA))
      # y_shunt_mk = im*(lines[l].b/2 + (buses[tt].Bs / baseMVA))
      from_bus = (f == ff) ? f : t # are these correct or do they need to depend on whether t,f are swapped
      to_bus = (t == tt) ? t : f # are these correct or do they need to depend on whether t,f are swapped
    
      line_b = sum(lines[(lines.from .== from_bus) .& (lines.to .== to_bus)].b) # duplicate line shunts
      # line_b = 1.0 / sum(1.0 ./ lines.b[(lines.from .== from_bus) .& (lines.to .== to_bus)]) # duplicate line shunts
      # if isinf(line_b)
      #   line_b=0
      # end
      from_Bs = nobus_Bshunt ? 0.0 : first(buses[buses.bus_i .== from_bus].Bs) # don't sum bus shunts
      to_Bs = nobus_Bshunt ? 0.0 : first(buses[buses.bus_i .== to_bus].Bs) # don't sum bus shunts
      # from_Bs = sum(buses[buses.bus_i .== from_bus].Bs)
      # to_Bs = sum(buses[buses.bus_i .== to_bus].Bs)
      y_shunt_km = im*(line_b / 2 + (from_Bs / baseMVA))
      y_shunt_mk = im*(line_b / 2 + (to_Bs / baseMVA))
      

      L["params"]["y"] = Dict()
      L["params"]["y"]["re"] = 0.0
      # L["params"]["y"]["im"] = imag(Ys)
      L["params"]["y_shunt_km"] = +1y_shunt_km # -1 +1 -1 +1
      L["params"]["y_shunt_mk"] = +1y_shunt_mk # +1 +1 -1 -1
      # L["params"]["y"] = Dict()
      # L["params"]["y"]["re"] = 0.0
      L["params"]["y"]["im"] = -opfmodeldata[:Y][ff, tt]
      # L["params"]["y_shunt_km"] = complex(0.0, opfmodeldata[:YshI][ff])
      # L["params"]["y_shunt_mk"] = complex(0.0, opfmodeldata[:YshI][tt])
    elseif line_type == "PiModelLineTapratio"
      L["params"]["y"] = Dict()
      L["params"]["y"]["re"] = 0.0
      L["params"]["y"]["im"] = -opfmodeldata[:Y][ff, tt]
      L["params"]["y_shunt_km"] = complex(0.0, opfmodeldata[:YshI][ff])
      L["params"]["y_shunt_mk"] = complex(0.0, opfmodeldata[:YshI][tt])
      tap = (opfmodeldata[:lines][l].ratio == 0) ? (1.0) : (opfmodeldata[:lines][l].ratio)
      tap *= exp(opfmodeldata[:lines][l].angle * pi/180 * im)
      L["params"]["t_shunt_km"] = tap
      L["params"]["t_shunt_mk"] = tap
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
                opfmodeldata::Dict, line_type::String, bus_type::String,
                nobus_Bshunt::Bool=false, phys::Dict=physDefault)
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
  return opf2pd(fileout, optimal_values, opfmodeldata, line_type, bus_type, nobus_Bshunt, phys)
end



  # # lines
  # for l = 1:nline
  #   f = opfmodeldata[:lines][l].from
  #   t = opfmodeldata[:lines][l].to
  #   ff = min(f,t) # need to have from.id <= to.id for PD lightgraphs dependency
  #   tt = max(f,t) # need to have from.id <= to.id for PD lightgraphs dependency
  #   good_order = (ff==f)*(tt==t)
    
  #   L = Dict()
  #   L["name"] = "branch" * string(l)
  #   L["params"] = Dict()
  #   L["params"]["from"] = "bus" * string(ff)
  #   L["params"]["to"] = "bus" * string(tt)
  #   L["type"] = line_type
  #   if line_type == "PiModelLine"
  #     y = Ys = 1/(lines[l].x*im)
  #     y_shunt_km = im*(lines[l].b/2 + (buses[ff].Bs / baseMVA))
  #     y_shunt_mk = im*(lines[l].b/2 + (buses[tt].Bs / baseMVA))
  #     from_bus = f
  #     to_bus = t
    
  #     line_b = sum(lines[(lines.from .== from_bus) .& (lines.to .== to_bus)].b) # duplicate line shunts
  #     from_Bs = first(buses[buses.bus_i .== from_bus].Bs) # don't sum bus shunts
  #     to_Bs = first(buses[buses.bus_i .== to_bus].Bs) # don't sum bus shunts
  #     # from_Bs = sum(buses[buses.bus_i .== from_bus].Bs)
  #     # to_Bs = sum(buses[buses.bus_i .== to_bus].Bs)
  #     y_shunt_km = im*(line_b / 2 + (from_Bs / baseMVA))
  #     y_shunt_mk = im*(line_b / 2 + (to_Bs / baseMVA))
      

  #     L["params"]["y"] = Dict()
  #     L["params"]["y"]["re"] = 0.0
  #     # L["params"]["y"]["im"] = imag(Ys)
  #     L["params"]["y_shunt_km"] = +1y_shunt_km # -1 +1 -1 +1
  #     L["params"]["y_shunt_mk"] = +1y_shunt_mk # +1 +1 -1 -1
  #     # L["params"]["y"] = Dict()
  #     # L["params"]["y"]["re"] = 0.0
  #     L["params"]["y"]["im"] = -opfmodeldata[:Y][ff, tt]
  #     # L["params"]["y_shunt_km"] = complex(0.0, opfmodeldata[:YshI][ff])
  #     # L["params"]["y_shunt_mk"] = complex(0.0, opfmodeldata[:YshI][tt])
  #   elseif line_type == "PiModelLineTapratio"
  #     L["params"]["y"] = Dict()
  #     L["params"]["y"]["re"] = 0.0
  #     L["params"]["y"]["im"] = -opfmodeldata[:Y][ff, tt]
  #     L["params"]["y_shunt_km"] = complex(0.0, opfmodeldata[:YshI][ff])
  #     L["params"]["y_shunt_mk"] = complex(0.0, opfmodeldata[:YshI][tt])
  #     tap = (opfmodeldata[:lines][l].ratio == 0) ? (1.0) : (opfmodeldata[:lines][l].ratio)
  #     tap *= exp(opfmodeldata[:lines][l].angle * pi/180 * im)
  #     L["params"]["t_shunt_km"] = tap
  #     L["params"]["t_shunt_mk"] = tap
  #   else
  #     L["params"]["Y"] = Dict()
  #     L["params"]["Y"]["re"] = 0.0
  #     L["params"]["Y"]["im"] = -opfmodeldata[:Y][ff, tt]
  #   end

  #   # add to `nodes` array
  #   push!(out["lines"], L)
  # end