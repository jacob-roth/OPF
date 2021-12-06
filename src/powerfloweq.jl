function powerfloweqs(Vm,Va,opfmodeldata)
  FromLines = opfmodeldata[:FromLines]
  YttI = opfmodeldata[:YttI]
  YshI = opfmodeldata[:YshI]
  YtfR = opfmodeldata[:YtfR]
  YffR = opfmodeldata[:YffR]
  YtfI = opfmodeldata[:YtfI]
  YttR = opfmodeldata[:YttR]
  nbus = opfmodeldata[:nbus]
  YffI = opfmodeldata[:YffI]
  generators = opfmodeldata[:generators]
  BusIdx = opfmodeldata[:BusIdx]
  YftI = opfmodeldata[:YftI]
  buses = opfmodeldata[:buses]
  Y = opfmodeldata[:Y]
  lines = opfmodeldata[:lines]
  YshR = opfmodeldata[:YshR]
  YftR = opfmodeldata[:YftR]
  ToLines = opfmodeldata[:ToLines]
  baseMVA = opfmodeldata[:baseMVA]
  BusGenerators = opfmodeldata[:BusGenerators]
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
  busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
  nbus = length(buses); nline = length(lines); ngen = length(generators)

  nbus = length(buses)
  R     = opfdata.bus_ref
  G     = filter(x -> x ∉ R, findall(.!isempty.(BusGeners)))
  L     = findall(isempty.(BusGeners))
  not_R = deleteat!(collect(1:nbus), R)
  
  P = zeros(nbus)
  Q = zeros(nbus)

  sumreduceempty(s) = isempty(s) ? 0.0 : reduce(+,s)

  for b in not_R
    #real part
    P[b] = (
      ( sumreduceempty( YffR[l] for l in FromLines[b]) + sumreduceempty([YttR[l] for l in ToLines[b]]) + YshR[b]) * Vm[b]^2
      + sumreduceempty( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
      + sumreduceempty( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
      # - ( sumreduceempty(baseMVA*Pg[g] for g in BusGeners[b]) - buses[b].Pd ) / baseMVA      # Sbus part
      )
    #imaginary part
    Q[b] = ( sumreduceempty(-YffI[l] for l in FromLines[b]) + sumreduceempty(-YttI[l] for l in ToLines[b]) - YshI[b] ) * Vm[b]^2
      + sumreduceempty( Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
      + sumreduceempty( Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
      # - ( sumreduceempty(baseMVA*Qg[g] for g in BusGeners[b]) - buses[b].Qd ) / baseMVA      #Sbus part
  end
  for b in R
    #real part
    P[b] =  ( sumreduceempty( YffR[l] for l in FromLines[b]) + sumreduceempty( YttR[l] for l in ToLines[b]) + YshR[b] ) * Vm[b]^2
      + sumreduceempty( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
      + sumreduceempty( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
      # - ( sumreduceempty(baseMVA*Pg[g] for g in BusGeners[b]) - buses[b].Pd ) / baseMVA      # Sbus part
    #imaginary part
    Q[b] = ( sumreduceempty(-YffI[l] for l in FromLines[b]) + sumreduceempty(-YttI[l] for l in ToLines[b]) - YshI[b] ) * Vm[b]^2
      + sumreduceempty( Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
      + sumreduceempty( Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
      # - ( sumreduceempty(baseMVA*Qg[g] for g in BusGeners[b]) - buses[b].Qd ) / baseMVA      #Sbus part
  end
  return P, Q
end