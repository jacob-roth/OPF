function add_p_constraint!(opfmodel::JuMP.Model, opfmodeldata::Dict, b::Int64, c::Int64=0)
  Pg0        = opfmodel[:Pg]
  Qg0        = opfmodel[:Qg]
  Vm0        = opfmodel[:Vm]
  Va0        = opfmodel[:Va]
  lines      = opfmodeldata[:lines];
  buses      = opfmodeldata[:buses];
  generators = opfmodeldata[:generators];
  baseMVA    = opfmodeldata[:baseMVA];
  busIdx     = opfmodeldata[:BusIdx];
  FromLines  = opfmodeldata[:FromLines];
  ToLines    = opfmodeldata[:ToLines];
  BusGeners  = opfmodeldata[:BusGenerators];
  YffR       = opfmodeldata[:YffR];
  YffI       = opfmodeldata[:YffI];
  YttR       = opfmodeldata[:YttR];
  YttI       = opfmodeldata[:YttI];
  YftR       = opfmodeldata[:YftR];
  YftI       = opfmodeldata[:YftI];
  YtfR       = opfmodeldata[:YtfR];
  YtfI       = opfmodeldata[:YtfI];
  YshR       = opfmodeldata[:YshR];
  YshI       = opfmodeldata[:YshI];
  Y          = opfmodeldata[:Y];
  nbus  = length(buses); nline = length(lines); ngen = length(generators)
  R     = findall(buses.bustype .== 3)[1]
  G     = filter(x -> x ∉ R, findall(.!isempty.(opfmodeldata[:BusGenerators])))
  L     = findall(isempty.(opfmodeldata[:BusGenerators]))
  not_R = deleteat!(collect(1:nbus), R)

  ## composite
  if c == 0
    Pg = Pg0
    Qg = Qg0
    Vm = Vm0
    Va = Va0
  else
    Pg = opfmodel[Symbol("Pg_$(c)_container")]
    Qg = opfmodel[Symbol("Qg_$(c)_container")]
    Vm = opfmodel[Symbol("Vm_$(c)_container")]
    Va = opfmodel[Symbol("Va_$(c)_container")]
    # Pg    = Array{Variable,1}(undef, ngen)
    # Qg    = Array{Variable,1}(undef, ngen)
    # Vm    = Array{Variable,1}(undef, nbus)
    # Va    = Array{Variable,1}(undef, nbus)
    # Pg = Pg0
    # # for i in first.(BusGeners[G])
    #   # Pg[i] = opfmodel[Symbol("Pg_$(c)")][i]
    # # end
    # # for i in first.(BusGeners[R])
    #   # Pg[i] = Pg0[i]
    # # end
    # for gi in zip([G; R], [first.(BusGeners[G]); first.(BusGeners[R])])
    #   g = gi[1]
    #   i = gi[2]
    #   Qg[i] = opfmodel[Symbol("Qg_$(c)")][g]
    # end
    # for i in [G; R]
    #   Vm[i] = Vm0[i]
    # end
    # for i in L
    #   Vm[i] = opfmodel[Symbol("Vm_$(c)")][i]
    # end
    # for i in [R]
    #   Va[i] = Va0[i]
    # end
    # for i in [L; G]
    #   Va[i] = opfmodel[Symbol("Va_$(c)")][i]
    # end
    # JuMP.registerobject(opfmodel, Symbol("Qg_$(c)_container", Qg, "Qg_$(c)_container")
    # JuMP.registerobject(opfmodel, Symbol("Vm_$(c)_container", Vm, "Vm_$(c)_container")
    # JuMP.registerobject(opfmodel, Symbol("Va_$(c)_container", Qg, "Va_$(c)_container")
  end

  # real part
  P_b = @NLexpression(
    opfmodel,
    ( sum( YffR[l] for l in FromLines[b]) + sum( YttR[l] for l in ToLines[b]) + YshR[b] ) * Vm[b]^2
    + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
    + sum( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
    - ( sum(baseMVA*Pg[g] for g in BusGeners[b]) - buses[b].Pd ) / baseMVA      # Sbus part
    )

  @NLconstraint(opfmodel, P_b==0)
  if c == 0
    Pb = "P$(b)"
  else
    Pb = "c$(c)_P$(b)"
  end
  JuMP.registerobject(opfmodel, Symbol(Pb), P_b, Pb)
end

function add_q_constraint!(opfmodel::JuMP.Model, opfmodeldata::Dict, b::Int64, c::Int64=0)
  Pg0        = opfmodel[:Pg]
  Qg0        = opfmodel[:Qg]
  Vm0        = opfmodel[:Vm]
  Va0        = opfmodel[:Va]
  lines      = opfmodeldata[:lines];
  buses      = opfmodeldata[:buses];
  generators = opfmodeldata[:generators];
  baseMVA    = opfmodeldata[:baseMVA];
  busIdx     = opfmodeldata[:BusIdx];
  FromLines  = opfmodeldata[:FromLines];
  ToLines    = opfmodeldata[:ToLines];
  BusGeners  = opfmodeldata[:BusGenerators];
  YffR       = opfmodeldata[:YffR];
  YffI       = opfmodeldata[:YffI];
  YttR       = opfmodeldata[:YttR];
  YttI       = opfmodeldata[:YttI];
  YftR       = opfmodeldata[:YftR];
  YftI       = opfmodeldata[:YftI];
  YtfR       = opfmodeldata[:YtfR];
  YtfI       = opfmodeldata[:YtfI];
  YshR       = opfmodeldata[:YshR];
  YshI       = opfmodeldata[:YshI];
  Y          = opfmodeldata[:Y];
  nbus  = length(buses); nline = length(lines); ngen = length(generators)
  R     = findall(buses.bustype .== 3)[1]
  G     = filter(x -> x ∉ R, findall(.!isempty.(opfmodeldata[:BusGenerators])))
  L     = findall(isempty.(opfmodeldata[:BusGenerators]))
  not_R = deleteat!(collect(1:nbus), R)

  ## composite
  if c == 0
    Pg = Pg0
    Qg = Qg0
    Vm = Vm0
    Va = Va0
  else
    Pg = opfmodel[Symbol("Pg_$(c)_container")]
    Qg = opfmodel[Symbol("Qg_$(c)_container")]
    Vm = opfmodel[Symbol("Vm_$(c)_container")]
    Va = opfmodel[Symbol("Va_$(c)_container")]
    # Pg    = Array{Variable,1}(undef, ngen)
    # Qg    = Array{Variable,1}(undef, ngen)
    # Vm    = Array{Variable,1}(undef, nbus)
    # Va    = Array{Variable,1}(undef, nbus)
    # Pg = Pg0
    # # for i in first.(BusGeners[G])
    # #   Pg[i] = opfmodel[Symbol("Pg_$(c)")][i]
    # # end
    # # for i in first.(BusGeners[R])
    # #   Pg[i] = Pg0[i]
    # # end
    # for gi in zip([G; R], [first.(BusGeners[G]); first.(BusGeners[R])])
    #   g = gi[1]
    #   i = gi[2]
    #   Qg[i] = opfmodel[Symbol("Qg_$(c)")][g]
    # end
    # for i in [G; R]
    #   Vm[i] = Vm0[i]
    # end
    # for i in L
    #   Vm[i] = opfmodel[Symbol("Vm_$(c)")][i]
    # end
    # for i in [R]
    #   Va[i] = Va0[i]
    # end
    # for i in [L; G]
    #   Va[i] = opfmodel[Symbol("Va_$(c)")][i]
    # end
  end

  # imaginary part
  Q_b = @NLexpression(
    opfmodel,
    ( sum(-YffI[l] for l in FromLines[b]) + sum(-YttI[l] for l in ToLines[b]) - YshI[b] ) * Vm[b]^2
    + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
    + sum( Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
    - ( sum(baseMVA*Qg[g] for g in BusGeners[b]) - buses[b].Qd ) / baseMVA      #Sbus part
    )
  @NLconstraint(opfmodel, Q_b==0)
  if c == 0
    Qb = "Q$(b)"
  else
    Qb = "c$(c)_Q$(b)"
  end
  JuMP.registerobject(opfmodel, Symbol(Qb), Q_b, Qb)
end

function add_line_current_constraint!(opfmodel::JuMP.Model, opfmodeldata::Dict, l::Int64, c::Int64=0)
    Pg0        = opfmodel[:Pg]
    Qg0        = opfmodel[:Qg]
    Vm0        = opfmodel[:Vm]
    Va0        = opfmodel[:Va]
    lines      = opfmodeldata[:lines];
    buses      = opfmodeldata[:buses];
    generators = opfmodeldata[:generators];
    baseMVA    = opfmodeldata[:baseMVA];
    busIdx     = opfmodeldata[:BusIdx];
    FromLines  = opfmodeldata[:FromLines];
    ToLines    = opfmodeldata[:ToLines];
    BusGeners  = opfmodeldata[:BusGenerators];
    YffR       = opfmodeldata[:YffR];
    YffI       = opfmodeldata[:YffI];
    YttR       = opfmodeldata[:YttR];
    YttI       = opfmodeldata[:YttI];
    YftR       = opfmodeldata[:YftR];
    YftI       = opfmodeldata[:YftI];
    YtfR       = opfmodeldata[:YtfR];
    YtfI       = opfmodeldata[:YtfI];
    YshR       = opfmodeldata[:YshR];
    YshI       = opfmodeldata[:YshI];
    Y          = opfmodeldata[:Y];

    ## composite
    if c == 0
        Pg = Pg0
        Qg = Qg0
        Vm = Vm0
        Va = Va0
    else
        Pg = opfmodel[Symbol("Pg_$(c)_container")]
        Qg = opfmodel[Symbol("Qg_$(c)_container")]
        Vm = opfmodel[Symbol("Vm_$(c)_container")]
        Va = opfmodel[Symbol("Va_$(c)_container")]
    end

    line = lines[lines.id .== l][1]
    if line.rateA != 0 && line.rateA < 1.0e10
        flowmax=(line.rateA/baseMVA)^2
        # branch current flows
        f = line.from
        t = line.to
        Y_tf = Y[t,f]
        Y_ft = Y[f,t]
        Vm_f = Vm[f]; Va_f = Va[f]
        Vm_t = Vm[t]; Va_t = Va[t]
        Yabs2 = max(abs2(Y_tf), abs2(Y_ft))
        ## NOTE: current from Frank & Rebennack OPF primer: eq 5.11 where turns/tap ratios are accounted for in `Y`
        F_l = @NLexpression(opfmodel, current2, (Vm_f^2 + Vm_t^2 - 2 * Vm_f * Vm_t * cos(Va_f - Va_t)) * Yabs2 - flowmax)
        @NLconstraint(opfmodel, current2 <= 0)

        if c == 0
            Fl = "F$(l)"
        else
            Fl = "c$(c)_F$(l)"
        end

        JuMP.registerobject(opfmodel, Symbol(Fl), F_l, Fl)
    end
end

function add_line_power_constraint!(opfmodel::JuMP.Model, opfmodeldata::Dict, l::Int64, c::Int64=0)
  Pg         = opfmodel[:Pg]
  Qg         = opfmodel[:Qg]
  Vm         = opfmodel[:Vm]
  Va         = opfmodel[:Va]
  lines      = opfmodeldata[:lines];
  buses      = opfmodeldata[:buses];
  generators = opfmodeldata[:generators];
  baseMVA    = opfmodeldata[:baseMVA];
  busIdx     = opfmodeldata[:BusIdx];
  FromLines  = opfmodeldata[:FromLines];
  ToLines    = opfmodeldata[:ToLines];
  BusGeners  = opfmodeldata[:BusGenerators];
  YffR       = opfmodeldata[:YffR];
  YffI       = opfmodeldata[:YffI];
  YttR       = opfmodeldata[:YttR];
  YttI       = opfmodeldata[:YttI];
  YftR       = opfmodeldata[:YftR];
  YftI       = opfmodeldata[:YftI];
  YtfR       = opfmodeldata[:YtfR];
  YtfI       = opfmodeldata[:YtfI];
  YshR       = opfmodeldata[:YshR];
  YshI       = opfmodeldata[:YshI];
  Y          = opfmodeldata[:Y];
  nbus = length(buses); nline = length(lines); ngen = length(generators)

  if lines[l].rateA!=0 && lines[l].rateA<1.0e10
    flowmax=(lines[l].rateA/baseMVA)^2
    # branch apparent power limits (from bus)
    Yff_abs2=YffR[l]^2+YffI[l]^2; Yft_abs2=YftR[l]^2+YftI[l]^2
    Yre=YffR[l]*YftR[l]+YffI[l]*YftI[l]; Yim=-YffR[l]*YftI[l]+YffI[l]*YftR[l]
    Ff_l = @NLexpression(opfmodel,
      Vm[busIdx[lines[l].from]]^2 *
              (
                Yff_abs2*Vm[busIdx[lines[l].from]]^2 + Yft_abs2*Vm[busIdx[lines[l].to]]^2
                + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
              )
              - flowmax)
    @NLconstraint(opfmodel, Ff_l<=0)

    # branch apparent power limits (to bus)
    Ytf_abs2=YtfR[l]^2+YtfI[l]^2; Ytt_abs2=YttR[l]^2+YttI[l]^2
    Yre=YtfR[l]*YttR[l]+YtfI[l]*YttI[l]; Yim=-YtfR[l]*YttI[l]+YtfI[l]*YttR[l]
    Ft_l = @NLexpression(opfmodel,
      Vm[busIdx[lines[l].to]]^2 *
              (
                Ytf_abs2*Vm[busIdx[lines[l].from]]^2 + Ytt_abs2*Vm[busIdx[lines[l].to]]^2
                + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
              )
              - flowmax)
    @NLconstraint(opfmodel, Ft_l<=0)
    if c == 0
      Ffl = "Ff$(l)"
      Ftl = "Ft$(l)"
    else
      Ffl = "c$(c)_Ff$(l)"
      Ftl = "c$(c)_Ft$(l)"
    end

    JuMP.registerobject(opfmodel, Symbol(Ffl), Ff_l, Ffl)
    JuMP.registerobject(opfmodel, Symbol(Ftl), Ft_l, Ftl)
  end
end

function add_sc_constraint!(opfmodel::JuMP.Model, sc_data::Dict)
  @NLexpression(opfmodel, 1)
end
