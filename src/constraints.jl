function add_S_constraint!(opfmodel::JuMP.Model, data::Dict, options::Dict)
  """
  add apparent power constraints
  """
  ## parse inputs
  #### data
  opfdata = data[:opfdata]
  YffR = data[:YffI]
  YffI = data[:YttR]
  YttR = data[:YttI]
  YttI = data[:YftR]
  YftR = data[:YftI]
  YftI = data[:YtfR]
  YtfR = data[:YtfI]
  YtfI = data[:YshR]
  YshR = data[:YshI]
  YshI = data[:YshI]

  #### options
  lossless       = haskey(options, :lossless)       ? options[:lossless]       : false
  current_rating = haskey(options, :current_rating) ? options[:current_rating] : false
  remove_Bshunt  = haskey(options, :remove_Bshunt)  ? options[:remove_Bshunt]  : false
  remove_tap     = haskey(options, :remove_tap)     ? options[:remove_tap]     : false
  if lossless && !current_rating
    println("warning: lossless assumption requires `current_rating` instead of `power_rating`\n")
    current_rating = true
  end

  ## shortcuts
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
  busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
  nbus = length(buses); nline = length(lines); ngen = length(generators)

  #
  # branch apparent power limits
  #
  @constraintref F_fr[1:nline]  ## from bus, TODO: won't work in JuMP v0.19
  @constraintref F_to[1:nline]  ## to bus, TODO: won't work in JuMP v0.19
  nlinelim=0
  for l in 1:nline
    if lines[l].rateA!=0 && lines[l].rateA<1.0e10
      nlinelim += 1
      flowmax=(lines[l].rateA/baseMVA)^2

      # branch apparent power limits (from bus)
      Yff_abs2=YffR[l]^2+YffI[l]^2; Yft_abs2=YftR[l]^2+YftI[l]^2
      Yre=YffR[l]*YftR[l]+YffI[l]*YftI[l]; Yim=-YffR[l]*YftI[l]+YffI[l]*YftR[l]
      F_fr[l] = @NLconstraint(opfmodel,
	              Vm[BusIdx[lines[l].from]]^2 *
              	( Yff_abs2*Vm[BusIdx[lines[l].from]]^2 + Yft_abs2*Vm[BusIdx[lines[l].to]]^2
              	  + 2*Vm[BusIdx[lines[l].from]]*Vm[BusIdx[lines[l].to]]*(Yre*cos(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]])-Yim*sin(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]]))
              	)
                - flowmax <=0)

      ## skip to-bus current limit: lossless ==> symmetric & current ==> same as from
      ## if current_rating == false, then they are different constraints bc Vm-to vs Vm-from
      # if lossless && current_rating
      #   continue
      # end

      # branch apparent power limits (to bus)
      Ytf_abs2=YtfR[l]^2+YtfI[l]^2; Ytt_abs2=YttR[l]^2+YttI[l]^2
      Yre=YtfR[l]*YttR[l]+YtfI[l]*YttI[l]; Yim=-YtfR[l]*YttI[l]+YtfI[l]*YttR[l]
      F_to[l] = @NLconstraint(opfmodel,
        	      Vm[BusIdx[lines[l].to]]^2 *
                ( Ytf_abs2*Vm[BusIdx[lines[l].from]]^2 + Ytt_abs2*Vm[BusIdx[lines[l].to]]^2
                  + 2*Vm[BusIdx[lines[l].from]]*Vm[BusIdx[lines[l].to]]*(Yre*cos(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]])-Yim*sin(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]]))
                )
                - flowmax <=0)
    end
  end
  JuMP.registercon(opfmodel, :F_fr, F_fr)
  JuMP.registercon(opfmodel, :F_to, F_to)
  nothing
end

function add_i2_constraint!(opfmodel::JuMP.Model, data::Dict, options::Dict)
  """
  add current magnitude constraints
  """
  ## parse inputs
  #### data
  opfdata = data[:opfdata]
  YffR = data[:YffI]
  YffI = data[:YttR]
  YttR = data[:YttI]
  YttI = data[:YftR]
  YftR = data[:YftI]
  YftI = data[:YtfR]
  YtfR = data[:YtfI]
  YtfI = data[:YshR]
  YshR = data[:YshI]
  YshI = data[:YshI]

  #### options
  lossless       = haskey(options, :lossless)       ? options[:lossless]       : false
  current_rating = haskey(options, :current_rating) ? options[:current_rating] : false
  remove_Bshunt  = haskey(options, :remove_Bshunt)  ? options[:remove_Bshunt]  : false
  remove_tap     = haskey(options, :remove_tap)     ? options[:remove_tap]     : false
  if lossless && !current_rating
    println("warning: lossless assumption requires `current_rating` instead of `power_rating`\n")
    current_rating = true
  end

  ## shortcuts
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
  busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
  nbus = length(buses); nline = length(lines); ngen = length(generators)

  #
  # branch flow or branch apparent power limits
  #
  @constraintref F_fr[1:nline]  ## from bus, TODO: won't work in JuMP v0.19
  @constraintref F_to[1:nline]  ## to bus, TODO: won't work in JuMP v0.19
  nlinelim=0
  for l in 1:nline
    if lines[l].rateA!=0 && lines[l].rateA<1.0e10
      nlinelim += 1
      flowmax=(lines[l].rateA/baseMVA)^2

      # branch apparent power limits (from bus)
      Yff_abs2=YffR[l]^2+YffI[l]^2; Yft_abs2=YftR[l]^2+YftI[l]^2
      Yre=YffR[l]*YftR[l]+YffI[l]*YftI[l]; Yim=-YffR[l]*YftI[l]+YffI[l]*YftR[l]
      F_fr[l] = @NLconstraint(opfmodel,
	              Vm[BusIdx[lines[l].from]]^2 *
              	( Yff_abs2*Vm[BusIdx[lines[l].from]]^2 + Yft_abs2*Vm[BusIdx[lines[l].to]]^2
              	  + 2*Vm[BusIdx[lines[l].from]]*Vm[BusIdx[lines[l].to]]*(Yre*cos(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]])-Yim*sin(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]]))
              	)
                - flowmax <=0)

      ## skip to-bus current limit: lossless ==> symmetric & current ==> same as from
      ## if current_rating == false, then they are different constraints bc Vm-to vs Vm-from
      # if lossless && current_rating
      #   continue
      # end

      # branch apparent power limits (to bus)
      Ytf_abs2=YtfR[l]^2+YtfI[l]^2; Ytt_abs2=YttR[l]^2+YttI[l]^2
      Yre=YtfR[l]*YttR[l]+YtfI[l]*YttI[l]; Yim=-YtfR[l]*YttI[l]+YtfI[l]*YttR[l]
      F_to[l] = @NLconstraint(opfmodel,
        	      Vm[BusIdx[lines[l].to]]^2 *
                ( Ytf_abs2*Vm[BusIdx[lines[l].from]]^2 + Ytt_abs2*Vm[BusIdx[lines[l].to]]^2
                  + 2*Vm[BusIdx[lines[l].from]]*Vm[BusIdx[lines[l].to]]*(Yre*cos(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]])-Yim*sin(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]]))
                )
                - flowmax <=0)
    end
  end
  JuMP.registercon(opfmodel, :F_fr, F_fr)
  JuMP.registercon(opfmodel, :F_to, F_to)
  nothing
end