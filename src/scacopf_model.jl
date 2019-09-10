function scacopf_model(opfdata::OPFData, options::Dict=DefaultOptions(),
                       adjustments::Dict=DefaultAdjustments(), contingencies::Dict=Dict())
  ## setup
  opfmodeldata = get_opfmodeldata(opfdata, options, adjustments)
  nbus = length(opfmodeldata[:buses]); nline = length(opfmodeldata[:lines]); ngen = length(opfmodeldata[:generators])
  R     = opfdata.bus_ref
  G     = filter(x -> x ∉ R, findall(.!isempty.(opfmodeldata[:BusGenerators])))
  L     = findall(isempty.(opfmodeldata[:BusGenerators]))
  not_R = deleteat!(collect(1:nbus), R)

  ## base model
  M = acopf_model(opfdata, options, adjustments)
  m = M.m

  ## add contingency security constraints
  if !isempty(contingencies)
    for c_id in keys(contingencies)
      #
      # contingency (line removal)
      #
      c = contingencies[c_id];
      @assert(c.c_type == :line)
      rl = remove_line!(opfdata, c.asset.id);
      c_opfmodeldata = get_opfmodeldata(opfdata, options, adjustments)

      #
      # secondary variables
      #
      # for b in [rl.to; rl.from]
      Pg_ = @variable(m, [i=R], basename="Pg_$c_id")
      Qg_ = @variable(m, [i=[G;R]], basename="Qg_$c_id")
      Vm_ = @variable(m, [i=L], basename="Vm_$c_id")
      Va_ = @variable(m, [i=not_R], basename="Va_$c_id")
      JuMP.registerobject(m, Symbol("Pg_$(c_id)"), Pg_, "Pg_$(c_id)")
      JuMP.registerobject(m, Symbol("Qg_$(c_id)"), Qg_, "Qg_$(c_id)")
      JuMP.registerobject(m, Symbol("Vm_$(c_id)"), Vm_, "Vm_$(c_id)")
      JuMP.registerobject(m, Symbol("Va_$(c_id)"), Va_, "Va_$(c_id)")

      ## composite containers
      Pg    = Array{Variable,1}(undef, ngen)
      Qg    = Array{Variable,1}(undef, ngen)
      Vm    = Array{Variable,1}(undef, nbus)
      Va    = Array{Variable,1}(undef, nbus)
      for gi in zip([R], first.(opfdata.BusGenerators[R]))
        g = gi[1]
        i = gi[2]
        Pg[i] = m[Symbol("Pg_$(c_id)")][g]
      end
      for gi in zip(G, first.(opfdata.BusGenerators[G]))
        g = gi[1]
        i = gi[2]
        Pg[i] = m[Symbol("Pg")][i] # since in acopf Pg is indexed by generator not bus
      end
      for gi in zip([G; R], [first.(opfdata.BusGenerators[G]); first.(opfdata.BusGenerators[R])])
        g = gi[1]
        i = gi[2]
        Qg[i] = m[Symbol("Qg_$(c_id)")][g]
      end
      for i in [G; R]
        Vm[i] = m[:Vm][i]
      end
      for i in L
        Vm[i] = m[Symbol("Vm_$(c_id)")][i]
      end
      for i in [R]
        Va[i] = m[:Va][i]
      end
      for i in [L; G]
        Va[i] = m[Symbol("Va_$(c_id)")][i]
      end
      JuMP.registerobject(m, Symbol("Pg_$(c_id)_container"), Pg, "Pg_$(c_id)_container")
      JuMP.registerobject(m, Symbol("Qg_$(c_id)_container"), Qg, "Qg_$(c_id)_container")
      JuMP.registerobject(m, Symbol("Vm_$(c_id)_container"), Vm, "Vm_$(c_id)_container")
      JuMP.registerobject(m, Symbol("Va_$(c_id)_container"), Va, "Va_$(c_id)_container")

      #
      # power flow balance
      #
      for b in 1:nbus
        ## real part
        add_p_constraint!(m, c_opfmodeldata, b, c_id)
        ## imaginary part
        add_q_constraint!(m, c_opfmodeldata, b, c_id)
      end

      #
      # branch/lines flow limits
      #
      for l in deleteat!(collect(1:nline), c_id)
        if options[:current_rating]
          ## current
          add_line_current_constraint!(m, c_opfmodeldata, l, c_id)
        else
          ## apparent power (to & from)
          add_line_power_constraint!(m, c_opfmodeldata, l, c_id)
        end
      end

      #
      # contingency (line reinstatement)
      #
      reinstate_line!(opfdata, rl.id, rl);
    end
  end
  return OPFModel(m, :InitData, :SC)
end

function get_operating_points(opfmodel::JuMP.Model, contingencies::Dict)
  ops = Dict[]
  for c_id in keys(contingencies)
    D = Dict()
    D[:Pg] = deepcopy(getvalue(opfmodel[Symbol("Pg_$(c_id)_container")]))
    D[:Qg] = deepcopy(getvalue(opfmodel[Symbol("Qg_$(c_id)_container")]))
    D[:Vm] = deepcopy(getvalue(opfmodel[Symbol("Vm_$(c_id)_container")]))
    D[:Va] = deepcopy(getvalue(opfmodel[Symbol("Va_$(c_id)_container")]))
    push!(ops, D)
  end
  return ops
end
function get_operating_points(opfmodel::OPFModel, contingencies::Dict)
  @assert(opfmodel.kind == :SC)
  return get_operating_points(opfmodel.m, contingencies)
end


# function set_n1_limits(opfdata::OPFData, options::Dict=DefaultOptions(), contingencies::Dict,
#                        adjustments::Dict=DefaultAdjustments())
#     ## setup
#     c0_opfdata = deepcopy(opfdata)
#     c0_options = deepcopy(options)
#     c0_options[:feasibility] = true
#     c = get_contingencies(c0_opfdata, c0_options)
#
#     infeas = true
#     while infeas
#         ## base dispatch model
#         c0_dp, c0_m = get_dispatch_point(c0_opfdata, c0_options, adjustments)
#
#         ## contincency models
#         c_opfdata  = OPFData[]
#         cc_opfdata = deepcopy(c0_opfdata);
#         for c_id in keys(c)
#             cc = c[c_id];
#             @assert(cc.c_type == :line)
#
#             ## remove line for contingency
#             rl = remove_line!(cc_opfdata, cc.asset.id);
#
#             ## dispatch contingency
#             push!(c_opfdata, deepcopy(cc_opfdata));
#
#             ## reinstate contingency line
#             reinstate_line!(cc_opfdata, rl.id, rl);
#         end
#
#         ## solve subcases
#         get_dp = idx -> get_dispatch_point(c_opfdata[idx], c0_options, adjustments)
#         res = pmap(get_dp, 1:5)
#
#         ## get contingency dispatch points
#         c_dp = [x[1] for x in res]
#         c_m  = [x[2] for x in res]
#
#         ## update ratings
#         for dp in c_dp
#             c_flowmag2s = get_flowmag2s(dp, c0_opfdata, c0_options)
#             c_ratings   = get_ratings(c_flowmag2s, c0_opfdata.baseMVA)
#             c0_opfdata.lines.rateA .= max.(c_ratings, c0_opfdata.lines.rateA)
#         end
#
#         status = [x.status for x in c_m]
#         infeas = prod(status .!= :Optimal)
#     end
# end
#
#
# ## more like SCACOPF
# function set_n1_limits(opfdata::OPFData, options::Dict=DefaultOptions(), contingencies::Dict,
#                        adjustments::Dict=DefaultAdjustments())
#     ## setup
#     c0_opfdata = deepcopy(opfdata)
#     c0_options = deepcopy(options)
#     c0_options[:feasibility] = true
#     c = get_contingencies(c0_opfdata, c0_options)
#
#     infeas = true
#     while infeas
#         ## base dispatch model
#         # c0_m  = acopf_model(c0_opfdata, c0_options, adjustments)
#         # c0_dp =
#         c0_dp, c0_m = get_dispatch_point(c0_opfdata, c0_options, adjustments)
#         # add_sc_constraint!(c0_m, sc_data)
#
#         ## contincency models
#         c_opfdata  = OPFData[]
#         cc_opfdata = deepcopy(c0_opfdata);
#         for c_id in keys(c)
#             cc = c[c_id];
#             @assert(cc.c_type == :line)
#
#             ## remove line for contingency
#             rl = remove_line!(cc_opfdata, cc.asset.id);
#
#             ## dispatch contingency
#             push!(c_opfdata, deepcopy(cc_opfdata));
#
#             ## reinstate contingency line
#             reinstate_line!(cc_opfdata, rl.id, rl);
#         end
#
#         ## solve subcases
#         get_dp = idx -> get_dispatch_point(c_opfdata[idx], c0_options, adjustments)
#         res = pmap(get_dp, 1:5)
#
#         ## get contingency dispatch points
#         c_dp = [x[1] for x in res]
#         c_m  = [x[2] for x in res]
#
#         ## update ratings
#         for dp in c_dp
#             c_flowmag2s = get_flowmag2s(dp, c0_opfdata, c0_options)
#             c_ratings   = get_ratings(c_flowmag2s, c0_opfdata.baseMVA)
#             c0_opfdata.lines.rateA .= max.(c_ratings, c0_opfdata.lines.rateA)
#         end
#
#         ## check re-dispatch
#         [norm(x[:Pg] - c0_dp[:Pg]) for x in c_dp]
#
#         status = [x.status for x in c_m]
#         infeas = prod(status .== :Optimal)
#     end
# end