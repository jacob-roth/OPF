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
  @variable(opfmodel, -pi <= Va[i=1:nbus] <= pi)

  if options[:pw_angle_limits] == true
    for line in opfmodeldata[:lines]
      angmin, angmax = line.angmin, line.angmax
      if angmin ∉ [-360.0, 0.0, 360.0]
        println("adding pairwise min-angle constraint (angmin = $angmin)")
        @constraint(opfmodel, angmin * pi / 180 <= Va[line.from] - Va[line.to])
      end
      if angmax ∉ [-360.0, 0.0, 360.0]
        println("adding pairwise max-angle constraint (angmax = $angmax)")
        @constraint(opfmodel, Va[line.from] - Va[line.to] <= angmax * pi / 180)
     end
    end
  end

  ## fix the voltage angle at the reference bus
  if options[:slack0] == true
    setlowerbound(Va[opfdata.bus_ref], 0.0)
    setupperbound(Va[opfdata.bus_ref], 0.0)
  else
    setlowerbound(Va[opfdata.bus_ref], opfmodeldata[:buses][opfdata.bus_ref].Va)
    setupperbound(Va[opfdata.bus_ref], opfmodeldata[:buses][opfdata.bus_ref].Va)
  end

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
