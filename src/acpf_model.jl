function acpf_model(dispatch_point::Dict, opfdata, options::Dict=DefaultOptions(),
                    tol::Float64=1e-6, unbounded::Bool=true)
  # parse options
  lossless       = options[:lossless]
  current_rating = options[:current_rating]
  remove_Bshunt  = options[:remove_Bshunt]
  remove_tap     = options[:remove_tap]
  print_level    = options[:print_level]
  feasibility    = options[:feasibility]
  if lossless && !current_rating
    println("warning: lossless assumption requires `current_rating` instead of `power_rating`\n")
    current_rating = true
  end

  # shortcuts for compactness
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
  busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
  nbus = length(buses); nline = length(lines); ngen = length(generators)

  # bus IDs
  R     = opfdata.bus_ref
  G     = filter(x -> x ∉ R, findall(.!isempty.(BusGeners)))
  L     = findall(isempty.(BusGeners))
  not_R = deleteat!(collect(1:nbus), R)

  # branch admitances
  YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(lines, buses, baseMVA;
                                                      lossless=lossless, remove_Bshunt=remove_Bshunt, remove_tap=remove_tap)

  #
  # model
  #
  opfmodel = Model(solver=IpoptSolver(print_level=print_level, tol=tol))

  ## Pg constants: Pg @ all generators
  Pg = Array{Union{Float64, JuMP.Variable}}(undef, ngen)
  @variable(opfmodel, Pg_R[i=BusGeners[R]])  # unbounded slack...
  for g in first.(BusGeners[G])
    Pg[g] = dispatch_point[:Pg][g]
  end
  for g in first.(BusGeners[R])
    Pg[g] = Pg_R[g]
    setvalue(Pg[g], dispatch_point[:Pg][g])
  end

  ## Qg variables: Qg @ all non-slack generators
  Qg = Array{Union{Float64, JuMP.Variable}}(undef, ngen)
  @variable(opfmodel, Qg_G[i=first.(BusGeners[G])])  # unbounded generators
  @variable(opfmodel, Qg_R[i=BusGeners[R]])  # unbounded slack...
  for g in first.(BusGeners[G])
    Qg[g] = Qg_G[g]
    setvalue(Qg[g], dispatch_point[:Qg][g])
  end
  for g in first.(BusGeners[R])
    Qg[g] = Qg_R[g]
    setvalue(Qg[g], dispatch_point[:Qg][g])
  end

  ## Vm variables: Vm @ load buses
  Vm = Array{Union{Float64, JuMP.Variable}}(undef, nbus)
  @variable(opfmodel, Vm_L[i=L])  # unbounded loads
  for b in L
    setlowerbound(Vm_L[b], 0.0)
    Vm[b] = Vm_L[b]
    setvalue(Vm[b], dispatch_point[:Vm][b])
  end
  for b in G
    Vm[b] = dispatch_point[:Vm][b]
  end
  for b in R
    Vm[b] = dispatch_point[:Vm][b]
  end

  ## Va @ all non-slack buses
  Va = Array{Union{Float64, JuMP.Variable}}(undef, nbus)
  @variable(opfmodel, -pi <= Va_not_R[i=not_R] <= pi)
  for b in not_R
    Va[b] = Va_not_R[b]
    setvalue(Va[b], dispatch_point[:Va][b])
  end
  for b in R
    Va[b] = dispatch_point[:Va][b]
  end

  #
  # objective
  #
  @NLobjective(opfmodel, Min, 0)

  #
  # power flow balance
  #
  @constraintref Pc[1:nbus]
  @constraintref Qc[1:nbus]
  P = Array{JuMP.NonlinearExpression}(undef, nbus)
  Q = Array{JuMP.NonlinearExpression}(undef, nbus)
  for b in not_R
    #real part
    P[b] = @NLexpression(
      opfmodel,
      ( sum( YffR[l] for l in FromLines[b]) + sum( YttR[l] for l in ToLines[b]) + YshR[b] ) * Vm[b]^2
      + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
      + sum( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
      - ( sum(baseMVA*Pg[g] for g in BusGeners[b]) - buses[b].Pd ) / baseMVA      # Sbus part
      )
      Pc[b] = @NLconstraint(opfmodel, P[b]==0)
    #imaginary part
    Q[b] = @NLexpression(
      opfmodel,
      ( sum(-YffI[l] for l in FromLines[b]) + sum(-YttI[l] for l in ToLines[b]) - YshI[b] ) * Vm[b]^2
      + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
      + sum( Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
      - ( sum(baseMVA*Qg[g] for g in BusGeners[b]) - buses[b].Qd ) / baseMVA      #Sbus part
      )
      Qc[b] = @NLconstraint(opfmodel, Q[b]==0)
  end
  for b in R
    #real part
    P[b] = @NLexpression(
      opfmodel,
      ( sum( YffR[l] for l in FromLines[b]) + sum( YttR[l] for l in ToLines[b]) + YshR[b] ) * Vm[b]^2
      + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
      + sum( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
      - ( sum(baseMVA*Pg[g] for g in BusGeners[b]) - buses[b].Pd ) / baseMVA      # Sbus part
      )
    Pc[b] = @NLconstraint(opfmodel, P[b]==0)
    #imaginary part
    Q[b] = @NLexpression(
      opfmodel,
      ( sum(-YffI[l] for l in FromLines[b]) + sum(-YttI[l] for l in ToLines[b]) - YshI[b] ) * Vm[b]^2
      + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
      + sum( Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
      - ( sum(baseMVA*Qg[g] for g in BusGeners[b]) - buses[b].Qd ) / baseMVA      #Sbus part
      )
      Qc[b] = @NLconstraint(opfmodel, Q[b]==0)
  end
  JuMP.registerobject(opfmodel, :Pg, Pg, "Pg")
  JuMP.registerobject(opfmodel, :Qg, Qg, "Qg")
  JuMP.registerobject(opfmodel, :Vm, Vm, "Vm")
  JuMP.registerobject(opfmodel, :Va, Va, "Va")
  JuMP.registerobject(opfmodel, :P,  P,  "P")
  JuMP.registerobject(opfmodel, :Q,  Q,  "Q")
  JuMP.registercon(opfmodel,    :Pc, Pc)
  JuMP.registercon(opfmodel,    :Qc, Qc)
  if print_level >= 1
    @printf("Buses: %d  Lines: %d  Generators: %d\n", nbus, nline, ngen)
  end
  return OPFModel(opfmodel, :InitData, :D, Dict())
end
#
#
# function compute_p(b, Vm, Va, opfdata, YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI)
#   lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
#   busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
#   nbus = length(buses); nline = length(lines); ngen = length(generators)
#
#   p = (
#         ( reduce(+, YffR[l] for l in FromLines[b]; init=0.0) + reduce(+, YttR[l] for l in ToLines[b]; init=0.0) + YshR[b] ) * Vm[b]^2
#         + reduce(+, Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b]; init=0.0 )
#         + reduce(+, Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]; init=0.0   )
#         - ( reduce(+, baseMVA*Pg[g] for g in BusGeners[b]; init=0.0) - buses[b].Pd ) / baseMVA
#   )
#   return p
# end
#
# v = Vm .* exp.(im .* Va)
# # v = v[C.N.σ_e2i]
# Y = C.YY[C.N.σ_i2e, C.N.σ_i2e]
# S = v .* conj(Y*v)
#
# Pnet = [reduce(+, Pg[g] for g in casedata.opf.BusGenerators[i]; init=0.0) - (buses[i].Pd / 100.0) for i in 1:length(buses)]
# Qnet = [reduce(+, Qg[g] for g in casedata.opf.BusGenerators[i]; init=0.0) - (buses[i].Qd / 100.0) for i in 1:length(buses)]
#
# 1;2;13;22;23;27;3;4;5...
#
# pnet2=0.2779784343168673