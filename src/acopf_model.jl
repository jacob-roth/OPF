function acopf_model(opfdata::OPFData, options::Dict=DefaultOptions(), adjustments::Dict=DefaultAdjustments())
  #
  # model
  #
  opfmodeldata = get_opfmodeldata(opfdata, options, adjustments)
  opfmodel     = Model(solver=IpoptSolver(print_level=options[:print_level], tol=options[:tol]))
  nbus = length(opfmodeldata[:buses]); nline = length(opfmodeldata[:lines]); ngen = length(opfmodeldata[:generators])

  ## bound constrained variables
  @variable(opfmodel, opfmodeldata[:generators][i].Pmin <= Pg[i=1:ngen] <= opfmodeldata[:generators][i].Pmax)
  @variable(opfmodel, opfmodeldata[:generators][i].Qmin <= Qg[i=1:ngen] <= opfmodeldata[:generators][i].Qmax)
  @variable(opfmodel, opfmodeldata[:buses][i].Vmin <= Vm[i=1:nbus] <= opfmodeldata[:buses][i].Vmax)
  @variable(opfmodel, -pi <= Va[1:nbus] <= pi)
  @variable(opfmodel, 0 <= Ps[i=1:nbus] <= opfmodeldata[:buses][i].Pd) # real power shed
  @variable(opfmodel, 0 <= Qs[i=1:nbus] <= opfmodeldata[:buses][i].Qd) # reactive power shed

  ## fix the voltage angle at the reference bus
  setlowerbound(Va[opfdata.bus_ref], opfmodeldata[:buses][opfdata.bus_ref].Va)
  setupperbound(Va[opfdata.bus_ref], opfmodeldata[:buses][opfdata.bus_ref].Va)

  ## objective
  if options[:shed_load]
    @NLobjective(opfmodel, Min, sum((Ps[b] + Qs[b]) for b in 1:nbus))
  else
    setlowerbound.(Ps, 0); setupperbound.(Ps, 0)
    setlowerbound.(Qs, 0); setupperbound.(Qs, 0)
    if options[:feasibility] == true
      @NLobjective(opfmodel, Min, 0)
    else
      @NLobjective(opfmodel, Min, sum( opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n-2]*(opfmodeldata[:baseMVA]*Pg[i])^2
                                      +opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n-1]*(opfmodeldata[:baseMVA]*Pg[i])
                                      +opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n  ] for i=1:ngen))
    end
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
    @printf("Buses: %d  Lines: %d  Generators: %d\n", nbus, nline, ngen)
  end
  return OPFModel(opfmodel, :InitData, :D)
end