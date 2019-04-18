function ccacopf_model(opf_data, options::Dict=Dict(), data::Dict=Dict())
  # parse options
  lossless = haskey(options, :lossless) ? options[:lossless] : false
  current_rating = haskey(options, :current_rating) ? options[:current_rating] : false
  epsilon = haskey(options, :epsilon) ? options[:epsilon] : 0.05
  if lossless && !current_rating
    println("warning: lossless assumption requires `current_rating` instead of `power_rating`\n")
    current_rating = true
  end

  # shortcuts for compactness
  lines = opf_data.lines; buses = opf_data.buses; generators = opf_data.generators; baseMVA = opf_data.baseMVA
  busIdx = opf_data.BusIdx; FromLines = opf_data.FromLines; ToLines = opf_data.ToLines; BusGeners = opf_data.BusGenerators;
  nbus = length(buses); nline = length(lines); ngen = length(generators); nload = length(findall(buses.bustype .== 1))
  # @assert(nload + ngen == nbus); NOT assert bc some buses can have more than one generator...

  # branch admitances
  YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(lines, buses, baseMVA)
  Y = computeAdmittanceMatrix(opfdata)
  G = real.(Y); g_row, g_col, g_val = findnz(G)
  B = imag.(Y); b_row, b_col, b_val = findnz(B)

  # index sets
  J_idx = jac_idx(opfdata)  ## follows R-G-L ordering for variables, parameters, and PF equations
  J_f_idx = idx_f(J_idx[:eqns], false);  nf = length(J_f_idx)
  J_x_idx = idx_x(J_idx[:pars]       );  nx = length(J_x_idx)
  J_y_idx = idx_y(J_idx[:vars], false);  ny = length(J_y_idx)
  @assert(ny == nf)

  # parse data
  Σ = haskey(data, :Sigma) ? data[:Sigma] : Matrix(Diagonal(ones(nx)))

  #
  # model
  #
  if "19" ∈ split(string(Pkg.installed()["JuMP"]), ".")
    opfmodel = Model(with_optimizer(Ipopt.Optimizer))
  else
    opfmodel = Model(solver=IpoptSolver(print_level=0))
  end

  @variable(opfmodel, generators[i].Pmin <= Pg[i=1:ngen] <= generators[i].Pmax)
  @variable(opfmodel, generators[i].Qmin <= Qg[i=1:ngen] <= generators[i].Qmax)
  @variable(opfmodel, buses[i].Vmin <= Vm[i=1:nbus] <= buses[i].Vmax)
  @variable(opfmodel, -pi <= Va[1:nbus] <= pi)

  ## assumes no buses have generator and load (NOTE: why was this again?)
  @variable(opfmodel, buses[i].Pd/baseMVA <= Pd[i=1:nbus] <= buses[i].Pd/baseMVA)
  @variable(opfmodel, buses[i].Qd/baseMVA <= Qd[i=1:nbus] <= buses[i].Qd/baseMVA)

  ## CC variables
  @variable(opfmodel, Gamma[i=1:nx, j=1:ny])
  @variable(opfmodel, zeta[i=1:nx, j=1:ny])

  #fix the voltage angle at the reference bus
  if "19" ∈ split(string(Pkg.installed()["JuMP"]), ".")
    fix(Va[opf_data.bus_ref], buses[opf_data.bus_ref].Va; force = true)
  else
    setlowerbound(Va[opf_data.bus_ref], buses[opf_data.bus_ref].Va)
    setupperbound(Va[opf_data.bus_ref], buses[opf_data.bus_ref].Va)
  end

  @NLobjective(opfmodel, Min, sum( generators[i].coeff[generators[i].n-2]*(baseMVA*Pg[i])^2
			             +generators[i].coeff[generators[i].n-1]*(baseMVA*Pg[i])
				     +generators[i].coeff[generators[i].n  ] for i=1:ngen))

  #
  # power flow balance
  #
  #real part
  @NLconstraint(opfmodel, P[b=1:nbus],
    ( sum( YffR[l] for l in FromLines[b]) + sum( YttR[l] for l in ToLines[b]) + YshR[b] ) * Vm[b]^2
    + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
    + sum( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
    - ( sum(baseMVA*Pg[g] for g in BusGeners[b]) - sum(baseMVA*Pd[l] for l in busIdx[b]) ) / baseMVA      # Sbus part
    ==0)
  #imaginary part
  @NLconstraint(opfmodel, Q[b=1:nbus],
    ( sum(-YffI[l] for l in FromLines[b]) + sum(-YttI[l] for l in ToLines[b]) - YshI[b] ) * Vm[b]^2
    + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
    + sum( Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
    - ( sum(baseMVA*Qg[g] for g in BusGeners[b]) - sum(baseMVA*Qd[l] for l in busIdx[b]) ) / baseMVA      #Sbus part
    ==0)
  #
  # branch/lines flow limits
  #
  @constraintref F_fr[1:nline]  ## from bus, TODO: won't work in JuMP v0.19
  @constraintref F_to[1:nline]  ## to bus, TODO: won't work in JuMP v0.19
  nlinelim=0
  for l in 1:nline
    if lines[l].rateA!=0 && lines[l].rateA<1.0e10
      nlinelim += 1
      flowmax=(lines[l].rateA/baseMVA)^2

      #branch apparent power limits (from bus)
      Yff_abs2=YffR[l]^2+YffI[l]^2; Yft_abs2=YftR[l]^2+YftI[l]^2
      Yre=YffR[l]*YftR[l]+YffI[l]*YftI[l]; Yim=-YffR[l]*YftI[l]+YffI[l]*YftR[l]
      if current_rating == true
        F_fr[l] = @NLconstraint(opfmodel,
  	              1.0 *
                	( Yff_abs2*Vm[busIdx[lines[l].from]]^2 + Yft_abs2*Vm[busIdx[lines[l].to]]^2
                	  + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
                	)
                  - flowmax <=0)
      else
      F_fr[l] = @NLconstraint(opfmodel,
	              Vm[busIdx[lines[l].from]]^2 *
              	( Yff_abs2*Vm[busIdx[lines[l].from]]^2 + Yft_abs2*Vm[busIdx[lines[l].to]]^2
              	  + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
              	)
                - flowmax <=0)
      end
      #branch apparent power limits (to bus)
      Ytf_abs2=YtfR[l]^2+YtfI[l]^2; Ytt_abs2=YttR[l]^2+YttI[l]^2
      Yre=YtfR[l]*YttR[l]+YtfI[l]*YttI[l]; Yim=-YtfR[l]*YttI[l]+YtfI[l]*YttR[l]
      if current_rating == true
        F_to[l] = @NLconstraint(opfmodel,
          	      1.0 *
                  ( Ytf_abs2*Vm[busIdx[lines[l].from]]^2 + Ytt_abs2*Vm[busIdx[lines[l].to]]^2
                    + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
                  )
                  - flowmax <=0)
      else
        F_to[l] = @NLconstraint(opfmodel,
          	      Vm[busIdx[lines[l].to]]^2 *
                  ( Ytf_abs2*Vm[busIdx[lines[l].from]]^2 + Ytt_abs2*Vm[busIdx[lines[l].to]]^2
                    + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
                  )
                  - flowmax <=0)
      end
    end
  end
  JuMP.registercon(opfmodel, :F_fr, F_fr)
  JuMP.registercon(opfmodel, :F_to, F_to)

  #
  # jacobian: http://schevalier.com/wp-content/uploads/2017/02/Power-Flow-and-Covariance-Matrix.pdf
  #
  #### components
  dP_dVm = Array{JuMP.NonlinearExpression,2}(undef, nbus, nbus)
  dP_dVa = Array{JuMP.NonlinearExpression,2}(undef, nbus, nbus)
  dQ_dVm = Array{JuMP.NonlinearExpression,2}(undef, nbus, nbus)
  dQ_dVa = Array{JuMP.NonlinearExpression,2}(undef, nbus, nbus)
  for q = 1:nbus # P, Q; equations
    for b = 1:nbus # Vm, Va; buses
     h = busIdx[mod1(b, nbus)]
     k = busIdx[mod1(q, nbus)]
     if h == k
       IDX = Y[h,:].nzind
       dP_dVa[b, q] = @NLexpression(opfmodel, -((Vm[h] * sum(Vm[kk] * ( G[h,kk] * sin(Va[h]-Va[kk]) - B[h,kk] * cos(Va[h]-Va[kk]) ) for kk in IDX))      ) - (B[h,h] * Vm[h]^2))
       dP_dVm[b, q] = @NLexpression(opfmodel,  ((Vm[h] * sum(Vm[kk] * ( G[h,kk] * cos(Va[h]-Va[kk]) + B[h,kk] * sin(Va[h]-Va[kk]) ) for kk in IDX))/Vm[h]) + (G[h,h] * Vm[h]))
       dQ_dVa[b, q] = @NLexpression(opfmodel,  ((Vm[h] * sum(Vm[kk] * ( G[h,kk] * cos(Va[h]-Va[kk]) + B[h,kk] * sin(Va[h]-Va[kk]) ) for kk in IDX))      ) - (G[h,h] * Vm[h]^2))
       dQ_dVm[b, q] = @NLexpression(opfmodel,  ((Vm[h] * sum(Vm[kk] * ( G[h,kk] * sin(Va[h]-Va[kk]) - B[h,kk] * cos(Va[h]-Va[kk]) ) for kk in IDX))/Vm[h]) - (B[h,h] * Vm[h]))
     else
       dP_dVa[b, q] = @NLexpression(opfmodel,  Vm[h] * Vm[k] * ( G[h,k] * sin(Va[h]-Va[k]) - B[h,k] * cos(Va[h]-Va[k]) ))
       dP_dVm[b, q] = @NLexpression(opfmodel,  Vm[h]         * ( G[h,k] * cos(Va[h]-Va[k]) + B[h,k] * sin(Va[h]-Va[k]) ))
       dQ_dVa[b, q] = @NLexpression(opfmodel, -Vm[h] * Vm[k] * ( G[h,k] * cos(Va[h]-Va[k]) + B[h,k] * sin(Va[h]-Va[k]) ))
       dQ_dVm[b, q] = @NLexpression(opfmodel,  Vm[h]         * ( G[h,k] * sin(Va[h]-Va[k]) - B[h,k] * cos(Va[h]-Va[k]) ))
     end
    end
  end

  #### aggregate components
  Z_bb = zeros(nbus, nbus)
  I_bgen = zeros(nbus, ngen)
  for (i,j,v) in zip(generators.bus, collect(1:ngen), ones(ngen))
    I_bgen[i,j] = v
  end
  I_bb = Matrix(Diagonal(ones(nbus)))
  Z_bg = zeros(nbus, ngen)
  J = Array{Union{Float64,JuMP.NonlinearExpression}}(undef, 2nbus, 2ngen+2nbus+2nbus)
  J .= [ -I_bgen    Z_bg     dP_dVm   dP_dVa   I_bb   Z_bb;
          Z_bg     -I_bgen   dQ_dVm   dQ_dVa   Z_bb   I_bb    ]

  #### partition aggregated Jacobian
  dF_dx = J[J_f_idx, J_x_idx]
  dF_dy = J[J_f_idx, J_y_idx]

  #
  # Gamma constraint
  #
  @constraintref Gamma_constraint[1:ny, 1:nx]
  for i = 1:nf
    for j = 1:ny
      Gamma_constraint[i, j] = @NLconstraint(opfmodel, sum(dF_dy[i,k] * Gamma[k,j] for k = 1:ny) + dF_dx[i, j] == 0)
    end
  end
  JuMP.registercon(opfmodel, :Gamma_constraint, Gamma_constraint)







  # J = Array{JuMP.NonlinearExpression,2}(undef, 2nbus, 2ngen+2nbus+2nbus)
  ## x = (Pd_RGL, Qd_RGL)
  ## y = (Vm_L, Va_GL)
  ## f = (P_GL, Q_L)
  f_idx = [bus_RGL[:G]; bus_RGL[:L]]
  dFdx = Array{JuMP.NonlinearExpression,2}(undef, nbus, nbus)
  dFdy = Array{JuMP.NonlinearExpression,2}(undef, nbus, nbus)

  Z_bb = spzeros(nbus, nbus)
  Z_bb = spzeros(nbus, nbus)
  I_bgen = spzeros(nbus, ngen)
  for (i,j,v) in zip(opfdata.generators.bus, collect(1:ngen), ones(ngen))
    I_bgen[i,j] = v
  end
  Z_bg = spzeros(nbus, ngen)

  JJ = [ -I_bgen    Z_bg     dPdVm'   dPdVa'   I      Z_bb;
          Z_bg     -I_bgen   dQdVm'   dQdVa'   Z_bb   I    ]
  J_numerical = J_numerical[1:60, :]
  JJ - J_numerical




for q = 1:nbus # P, Q; equations
  for b = 1:nbus # Vm, Va; buses
   h = busIdx[mod1(b, nbus)]
   k = busIdx[mod1(q, nbus)]
   if h == k
     IDX = Y[h,:].nzind
     P = Vm[h] * sum(Vm[kk] * ( G[h,kk] * cos(Va[h]-Va[kk]) + B[h,kk] * sin(Va[h]-Va[kk]) ) for kk in IDX)
     Q = Vm[h] * sum(Vm[kk] * ( G[h,kk] * sin(Va[h]-Va[kk]) - B[h,kk] * cos(Va[h]-Va[kk]) ) for kk in IDX)
     dP_dVa = -Q       - B[h,h] * Vm[h]^2
     dP_dVm =  P/Vm[h] + G[h,h] * Vm[h]
     dQ_dVa =  P       - G[h,h] * Vm[h]^2
     dQ_dVm =  Q/Vm[h] - B[h,h] * Vm[h]
   else
     dP_dVa =  Vm[h] * Vm[k] * ( G[h,k] * sin(Va[h]-Va[k]) - B[h,k] * cos(Va[h]-Va[k]) )
     dP_dVm =  Vm[h]         * ( G[h,k] * cos(Va[h]-Va[k]) + B[h,k] * sin(Va[h]-Va[k]) )
     dQ_dVa = -Vm[h] * Vm[k] * ( G[h,k] * cos(Va[h]-Va[k]) + B[h,k] * sin(Va[h]-Va[k]) )
     dQ_dVm =  Vm[h]         * ( G[h,k] * sin(Va[h]-Va[k]) - B[h,k] * cos(Va[h]-Va[k]) )
   end
   J[q, b]           = dP_dVm
   J[q, nbus+b]      = dP_dVa
   J[nbus+q, b]      = dQ_dVm
   J[nbus+q, nbus+b] = dQ_dVa
  end
end
dPdVm = J[1:nbus, 1:nbus]
dPdVa = J[1:nbus, (nbus+1):(2nbus)]
dQdVm = J[(nbus+1):(2nbus), 1:nbus]
dQdVa = J[(nbus+1):(2nbus), (nbus+1):(2nbus)]

J_numerical__ = J_numerical[:, (2ngen+1):(end)]
J_numerical__ = J_numerical__[:, 1:2nbus]
J_numerical__[1:nbus,1:nbus] - dPdVm'
J_numerical__[1:nbus,(nbus+1):(2nbus)] - dPdVa'
J_numerical__[(nbus+1):(2nbus), 1:nbus] - dQdVm'
J_numerical__[(nbus+1):(2nbus), (nbus+1):(2nbus)] - dQdVa'

          G_1_1 = ( reduce(+, YffR[l] for l in FromLines[h]; init=0) + reduce(+, YttR[l] for l in ToLines[h]; init=0) + YshR[h] )
          G_1_2 = YftR[1]
          reduce(+, Vm[h]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[h]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[h]-Va[busIdx[lines[l].to]]  )) for l in FromLines[h]; init=0 )
          reduce(+, Vm[h]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[h]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[h]-Va[busIdx[lines[l].from]])) for l in ToLines[h]  ; init=0 )
          dP_dVa =  (-1.0)      * ( Vm[h] * sum( Vm[kk] * ( G[h,kk] * sin(Va[h]-Va[k]) - B[h,kk] * cos(Va[h]-Va[kk])) for kk in IDX) ) - B[h,h]*Vm[h]^2

          ( sum( YffR[l] for l in FromLines[b]) + sum( YttR[l] for l in ToLines[b]) + YshR[b] ) * Vm[b]^2
          + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
          + sum( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )

          for l in ToLines[b]
            println(YffR[l])
          end
          for l in FromLines[b]
            println(YffR[l])
          end

          ## JR
          dP_dVm = 0
          for kk in IDX
            global val = Vm[h] * Vm[kk] * ( G[h,kk] * cos(Va[h]-Va[kk]) + B[h,kk] * sin(Va[h]-Va[kk]) )
            global dP_dVm += val
            println(dP_dVm, ", ", val)
          end
          # val = G[h,h]*Vm[h]^2
          # dP_dVm += val
          println(dP_dVm, ", ", val)
          Vm[h] * sum(Vm[kk] * ( G[h,kk] * cos(Va[h]-Va[kk]) + B[h,kk] * sin(Va[h]-Va[kk]) ) for kk in IDX)

          ## opf
          dP_dVm = 0
          for l in FromLines[h]
            global val = Vm[h]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[h]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[h]-Va[busIdx[lines[l].to]]  ))
            global dP_dVm += val
            println(dP_dVm, ", " , val)
          end
          val = sum(YffR[l] for l in FromLines[h])*Vm[h]^2
          dP_dVm += val
          println(dP_dVm, ", " , val)


  @printf("Buses: %d  Lines: %d  Generators: %d\n", nbus, nline, ngen)
  println("Lines with limits  ", nlinelim)

  return OPFModel(opfmodel, :InitData, :S)
end