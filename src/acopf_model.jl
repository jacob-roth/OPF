function acopf_model(opfdata::OPFData, options::Dict=DefaultOptions(), adjustments::Dict=DefaultAdjustments())
  #
  # model
  #
  opfmodeldata = get_opfmodeldata(opfdata, options, adjustments)
  opfmodel = Model(solver=IpoptSolver(print_level=options[:print_level], tol=options[:tol]))
  nbus = length(opfmodeldata[:buses]); nline = length(opfmodeldata[:lines]); ngen = length(opfmodeldata[:generators])

  ## bound constrained variables
  @variable(opfmodel, opfmodeldata[:generators][i].Pmin <= Pg[i=1:ngen] <= opfmodeldata[:generators][i].Pmax)
  @variable(opfmodel, opfmodeldata[:generators][i].Qmin <= Qg[i=1:ngen] <= opfmodeldata[:generators][i].Qmax)
  @variable(opfmodel, opfmodeldata[:buses][i].Vmin <= Vm[i=1:nbus] <= opfmodeldata[:buses][i].Vmax)

  line_ids = get_line_ids(opfmodeldata)
  for (from_idx, to_idx) in line_ids
    angmin, angmax = get_ang_bounds(opfmodeldata, from_idx, to_idx)
    @variable(opfmodel, max(-pi, angmin * pi / 180) <= Va[from_idx] - Va[to_idx] <= min(pi, angmax * pi / 180))
  end

  ## fix the voltage angle at the reference bus
  setlowerbound(Va[opfdata.bus_ref], opfmodeldata[:buses][opfdata.bus_ref].Va)
  setupperbound(Va[opfdata.bus_ref], opfmodeldata[:buses][opfdata.bus_ref].Va)

  ## objective
  if options[:feasibility] == true
    @NLobjective(opfmodel, Min, 0)
  else
    obj = @NLexpression(opfmodel, sum(opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n-2] * (opfmodeldata[:baseMVA] * Pg[i]) ^ 2
                                    + opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n-1] * (opfmodeldata[:baseMVA] * Pg[i])
                                    + opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n] for i=1:ngen))
    JuMP.registerobject(opfmodel, :obj, obj, "obj")
    @NLobjective(opfmodel, Min, obj)
  end

  #
  # power flow balance
  #
  for b in 1:nbus
    ## real part
    add_p_constraint!(opfmodel, opfmodeldata, b)
    ## imaginary part
    add_q_constraint!(opfmodel, opfmodeldata, b)
  end

  #
  # branch/lines flow limits
  #
  for l in 1:nline
    if options[:current_rating]
      ## current
      add_line_current_constraint!(opfmodel, opfmodeldata, l)
    else
      ## apparent power (to & from)
      add_line_power_constraint!(opfmodel, opfmodeldata, l)
    end
  end

  if options[:print_level] >= 1
    @printf("Buses: %d  Lines: %d  Generators: %d \n", nbus, nline, ngen)
  end
  return OPFModel(opfmodel, :InitData, :D)
end

function get_line_ids(opfmodeldata::Dict)
    num_lines = length(opfmodeldata[:lines])
    from_lines = opfmodeldata[:lines].from
    to_lines = opfmodeldata[:lines].to
    return [(from_lines[i], to_lines[i]) for i in 1:num_lines]
end

function get_ang_bounds(opfmodeldata::Dict, from_idx::Int, to_idx::Int)
    all_lines = opfmodeldata[:lines]
    target_line = lines[(lines.from .== from_idx) .& (lines.to .== to_idx)]
    return (target_line.angmin, target_line.angmax)
end
