function acopf_model(opfdata, options::Dict=DefaultOptions(), adjustments::Dict=DefaultAdjustments())
  # parse options
  lossless       = options[:lossless]
  current_rating = options[:current_rating]
  remove_Bshunt  = options[:remove_Bshunt]
  remove_tap     = options[:remove_tap]
  print_level    = options[:print_level]
  feasibility    = options[:feasibility]
  Pg_hi          = adjustments[:Pg_hi]
  Pg_lo          = adjustments[:Pg_lo]
  Qg_hi          = adjustments[:Qg_hi]
  Qg_lo          = adjustments[:Qg_lo]
  Vm_hi          = adjustments[:Vm_hi]
  Vm_lo          = adjustments[:Vm_lo]
  if lossless && !current_rating
    println("warning: lossless assumption requires `current_rating` instead of `power_rating`\n")
    current_rating = true
  end

  # shortcuts for compactness
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
  busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
  nbus = length(buses); nline = length(lines); ngen = length(generators)

  # branch admitances
  YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(lines, buses, baseMVA;
                                                      lossless=lossless, remove_Bshunt=remove_Bshunt, remove_tap=remove_tap)
  Y = computeAdmittanceMatrix(opfdata, options)

  #
  # model
  #
  if "19" ∈ split(string(Pkg.installed()["JuMP"]), ".")
    opfmodel = Model(with_optimizer(Ipopt.Optimizer))
  else
    opfmodel = Model(solver=IpoptSolver(print_level=print_level))
  end

  ## bound constraints
  Pg_hi_val = zeros(ngen); Pg_lo_val = zeros(ngen); Qg_hi_val = zeros(ngen); Qg_lo_val = zeros(ngen);
  for i in 1:ngen
    Pg_hi_val[i] = generators[i].Pmax - (Pg_hi[:v] * (i ∈ Pg_hi[:i]))
    Pg_lo_val[i] = generators[i].Pmin + (Pg_lo[:v] * (i ∈ Pg_lo[:i]))
    Qg_hi_val[i] = generators[i].Qmax - (Qg_hi[:v] * (i ∈ Qg_hi[:i]))
    Qg_lo_val[i] = generators[i].Qmin + (Qg_lo[:v] * (i ∈ Qg_lo[:i]))
  end
  Vm_hi_val = zeros(nbus); Vm_lo_val = zeros(nbus)
  for i in 1:nbus
    Vm_hi_val[i] = buses[i].Vmax - (Vm_hi[:v] * (i ∈ Vm_hi[:i]))
    Vm_lo_val[i] = buses[i].Vmin + (Vm_lo[:v] * (i ∈ Vm_lo[:i]))
  end
  # println("Pg upper bounds: $Pg_hi_val")
  # println("Pg lower bounds: $Pg_lo_val")
  println("Qg adjusted upper bounds: "); display(Qg_hi_val[Qg_hi[:i]])
  println("Qg adjusted lower bounds: "); display(Qg_lo_val[Qg_lo[:i]])
  # println("Vm upper bounds: $Vm_hi_val")
  # println("Vm lower bounds: $Vm_lo_val")

  ## variables
  @variable(opfmodel, Pg_lo_val[i] <= Pg[i=1:ngen] <= Pg_hi_val[i])
  @variable(opfmodel, Qg_lo_val[i] <= Qg[i=1:ngen] <= Qg_hi_val[i])
  @variable(opfmodel, Vm_lo_val[i] <= Vm[i=1:nbus] <= Vm_hi_val[i])
  @variable(opfmodel,          -pi <= Va[i=1:nbus] <= pi)

  #fix the voltage angle at the reference bus
  if "19" ∈ split(string(Pkg.installed()["JuMP"]), ".")
    fix(Va[opfdata.bus_ref], buses[opfdata.bus_ref].Va; force = true)
  else
    setlowerbound(Va[opfdata.bus_ref], buses[opfdata.bus_ref].Va)
    setupperbound(Va[opfdata.bus_ref], buses[opfdata.bus_ref].Va)
  end

  if feasibility == true
    @NLobjective(opfmodel, Min, 0)
  else
    @NLobjective(opfmodel, Min, sum( generators[i].coeff[generators[i].n-2]*(baseMVA*Pg[i])^2
  			             +generators[i].coeff[generators[i].n-1]*(baseMVA*Pg[i])
  				     +generators[i].coeff[generators[i].n  ] for i=1:ngen))
  end

  #
  # power flow balance
  #
  for b in 1:nbus
    #real part
    @NLconstraint(
      opfmodel,
      ( sum( YffR[l] for l in FromLines[b]) + sum( YttR[l] for l in ToLines[b]) + YshR[b] ) * Vm[b]^2
      + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
      + sum( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
      - ( sum(baseMVA*Pg[g] for g in BusGeners[b]) - buses[b].Pd ) / baseMVA      # Sbus part
      ==0)
    #imaginary part
    @NLconstraint(
      opfmodel,
      ( sum(-YffI[l] for l in FromLines[b]) + sum(-YttI[l] for l in ToLines[b]) - YshI[b] ) * Vm[b]^2
      + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
      + sum( Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
      - ( sum(baseMVA*Qg[g] for g in BusGeners[b]) - buses[b].Qd ) / baseMVA      #Sbus part
      ==0)
  end
  #
  # branch/lines flow limits
  #
  if current_rating
    @constraintref F[1:nline]     ## all lines
  else
    @constraintref F_fr[1:nline]  ## from bus, TODO: won't work in JuMP v0.19
    @constraintref F_to[1:nline]  ## to bus, TODO: won't work in JuMP v0.19
  end
  nlinelim=0
  for l in 1:nline
    if lines[l].rateA!=0 && lines[l].rateA<1.0e10
      nlinelim += 1
      flowmax=(lines[l].rateA/baseMVA)^2

      if current_rating
        #branch current flows
        line = lines[l]
        f = line.from
        t = line.to
        Y_tf = Y[t,f]
        Y_ft = Y[f,t]
        Vm_f = Vm[f]; Va_f = Va[f]
        Vm_t = Vm[t]; Va_t = Va[t]
        Yabs2 = max(abs2(Y_tf), abs2(Y_ft))
        ## NOTE: current from Frank & Rebennack OPF primer: eq 5.11 where turns/tap ratios are accounted for in `Y`
        @NLexpression(opfmodel, current2, (Vm_f^2 + Vm_t^2 - 2 * Vm_f * Vm_t * cos(Va_f - Va_t)) * Yabs2)
        F[l] = @NLconstraint(opfmodel, current2 <= flowmax)
      else
        #branch apparent power limits (from bus)
        Yff_abs2=YffR[l]^2+YffI[l]^2; Yft_abs2=YftR[l]^2+YftI[l]^2
        Yre=YffR[l]*YftR[l]+YffI[l]*YftI[l]; Yim=-YffR[l]*YftI[l]+YffI[l]*YftR[l]
        F_fr[l] = @NLconstraint(opfmodel,
          Vm[busIdx[lines[l].from]]^2 *
                  (
                    Yff_abs2*Vm[busIdx[lines[l].from]]^2 + Yft_abs2*Vm[busIdx[lines[l].to]]^2
                    + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
                  )
                  - flowmax <=0)

        #branch apparent power limits (to bus)
        Ytf_abs2=YtfR[l]^2+YtfI[l]^2; Ytt_abs2=YttR[l]^2+YttI[l]^2
        Yre=YtfR[l]*YttR[l]+YtfI[l]*YttI[l]; Yim=-YtfR[l]*YttI[l]+YtfI[l]*YttR[l]
        F_to[l] = @NLconstraint(opfmodel,
          Vm[busIdx[lines[l].to]]^2 *
                  (
                    Ytf_abs2*Vm[busIdx[lines[l].from]]^2 + Ytt_abs2*Vm[busIdx[lines[l].to]]^2
                    + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
                  )
                  - flowmax <=0)
      end
    end
  end
  if current_rating
    JuMP.registercon(opfmodel, :F, F)
  else
    JuMP.registercon(opfmodel, :F_fr, F_fr)
    JuMP.registercon(opfmodel, :F_to, F_to)
  end

  if print_level >= 1
    @printf("Buses: %d  Lines: %d  Generators: %d\n", nbus, nline, ngen)
    println("Lines with limits  ", nlinelim)
  end
  return OPFModel(opfmodel, :InitData, :D)
end