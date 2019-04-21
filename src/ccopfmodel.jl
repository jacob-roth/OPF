function ccacopf_model(opfdata, options::Dict=Dict(), data::Dict=Dict())
  # parse options
  lossless       = haskey(options, :lossless)       ? options[:lossless]       : false
  current_rating = haskey(options, :current_rating) ? options[:current_rating] : false
  epsilon_Vm     = haskey(options, :epsilon_Vm)     ? options[:epsilon_Vm]     : 0.05
  epsilon_Va     = haskey(options, :epsilon_Va)     ? options[:epsilon_Va]     : 0.05
  epsilon_Qg     = haskey(options, :epsilon_Qg)     ? options[:epsilon_Qg]     : 0.05
  γ              = haskey(options, :gamma)          ? options[:gamma]          : 1.0
  print_level    = haskey(options, :print_level)    ? options[:print_level]    : 0
  if lossless && !current_rating
    println("warning: lossless assumption requires `current_rating` instead of `power_rating`\n")
    current_rating = true
  end

  # shortcuts for compactness
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
  BusIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
  nbus = length(buses); nline = length(lines); ngen = length(generators); nload = length(findall(buses.bustype .== 1))
  # @assert(nload + ngen == nbus); NOT assert bc some buses can have more than one generator...

  # branch admitances
  YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(lines, buses, baseMVA)
  Y = computeAdmittanceMatrix(opfdata)
  G = real.(Y); g_row, g_col, g_val = findnz(G)
  B = imag.(Y); b_row, b_col, b_val = findnz(B)

  # parse data
  Σ_d    = haskey(data, :Sigma_d) ? data[:Sigma_d] : Matrix(Diagonal(ones(2nbus)))
  Va_min = haskey(data, :Va_min)  ? data[:Va_min]  : -pi * ones(nbus)
  Va_max = haskey(data, :Va_max)  ? data[:Va_max]  :  pi * ones(nbus)
  z = Normal(0,1)

  #
  # model
  #
  if "19" ∈ split(string(Pkg.installed()["JuMP"]), ".")
    opfmodel = Model(with_optimizer(Ipopt.Optimizer))
  else
    opfmodel = Model(solver=IpoptSolver(print_level=print_level))
  end

  #
  # variables
  #
  ## all variables
  @variable(opfmodel, generators[i].Pmin  <= Pg[i=1:ngen] <= generators[i].Pmax)
  @variable(opfmodel, generators[i].Qmin  <= Qg[i=1:ngen] <= generators[i].Qmax)
  @variable(opfmodel, buses[i].Vmin       <= Vm[i=1:nbus] <= buses[i].Vmax)
  @variable(opfmodel, Va_min[i]           <= Va[i=1:nbus] <= Va_max[i])
  @variable(opfmodel, buses[i].Pd/baseMVA <= Pd[i=1:nbus] <= buses[i].Pd/baseMVA)
  @variable(opfmodel, buses[i].Qd/baseMVA <= Qd[i=1:nbus] <= buses[i].Qd/baseMVA)
  ## partition variables
  b_RGL_idx, g_RGL_idx = RGL_idx(opfdata)
  #### unknown
  x = [Vm[b_RGL_idx[:L]]; Va[b_RGL_idx[:G]]; Va[b_RGL_idx[:L]]; Qg[g_RGL_idx[:G]]];
  #### control
  u = [Pg[g_RGL_idx[:G]]; Vm[b_RGL_idx[:G]]; Vm[b_RGL_idx[:R]]];
  #### parameter
  p = Va[b_RGL_idx[:R]];
  #### uncertainty
  d = [Pd; Qd];
  #### aggregate "known"
  y = [u; p; d];
  #### dims
  nx = length(x); nu = length(u); np = length(p); nd = length(d); ny = length(y)
  xidx = [xx.col for xx in x]  ## index in model `z`
  uidx = [xx.col for xx in u]  ## index in model `z`
  pidx = [xx.col for xx in p]  ## index in model `z`
  didx = [xx.col for xx in d]  ## index in model `z`
  yidx = [xx.col for xx in y]  ## index in model `z`
  Fidx = [b_RGL_idx[:L]; b_RGL_idx[:G]; nbus .+ b_RGL_idx[:L]; nbus .+ b_RGL_idx[:G]]  ## index in 2nbus equations

  ## CC variables
  @variable(opfmodel, Gamma[i=1:nx, j=1:ny])
  @variable(opfmodel, zeta[i=1:nx, j=1:ny])

  #fix the voltage angle at the reference bus
  if "19" ∈ split(string(Pkg.installed()["JuMP"]), ".")
    fix(Va[opfdata.bus_ref], buses[opfdata.bus_ref].Va; force = true)
  else
    setlowerbound(Va[opfdata.bus_ref], buses[opfdata.bus_ref].Va)
    setupperbound(Va[opfdata.bus_ref], buses[opfdata.bus_ref].Va)
  end

  @NLobjective(opfmodel, Min, sum( generators[i].coeff[generators[i].n-2]*(baseMVA*Pg[i])^2
			                            +generators[i].coeff[generators[i].n-1]*(baseMVA*Pg[i])
				                          +generators[i].coeff[generators[i].n  ] for i=1:ngen)
                            + sum(zeta[i,j]^2 for i=1:nx for j=1:ny))

  #
  # power flow balance
  #
  #real part
  @NLconstraint(opfmodel, P[b=1:nbus],
    ( sum( YffR[l] for l in FromLines[b]) + sum( YttR[l] for l in ToLines[b]) + YshR[b] ) * Vm[b]^2
    + sum( Vm[b]*Vm[BusIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[BusIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[BusIdx[lines[l].to]]  )) for l in FromLines[b] )
    + sum( Vm[b]*Vm[BusIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[BusIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[BusIdx[lines[l].from]])) for l in ToLines[b]   )
    - ( sum(baseMVA*Pg[g] for g in BusGeners[b]) - sum(baseMVA*Pd[l] for l in BusIdx[b]) ) / baseMVA      # Sbus part
    ==0)
  #imaginary part
  @NLconstraint(opfmodel, Q[b=1:nbus],
    ( sum(-YffI[l] for l in FromLines[b]) + sum(-YttI[l] for l in ToLines[b]) - YshI[b] ) * Vm[b]^2
    + sum( Vm[b]*Vm[BusIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[BusIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[BusIdx[lines[l].to]]  )) for l in FromLines[b] )
    + sum( Vm[b]*Vm[BusIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[BusIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[BusIdx[lines[l].from]])) for l in ToLines[b]   )
    - ( sum(baseMVA*Qg[g] for g in BusGeners[b]) - sum(baseMVA*Qd[l] for l in BusIdx[b]) ) / baseMVA      # Sbus part
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
                	( Yff_abs2*Vm[BusIdx[lines[l].from]]^2 + Yft_abs2*Vm[BusIdx[lines[l].to]]^2
                	  + 2*Vm[BusIdx[lines[l].from]]*Vm[BusIdx[lines[l].to]]*(Yre*cos(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]])-Yim*sin(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]]))
                	)
                  - flowmax <=0)
      else
      F_fr[l] = @NLconstraint(opfmodel,
	              Vm[BusIdx[lines[l].from]]^2 *
              	( Yff_abs2*Vm[BusIdx[lines[l].from]]^2 + Yft_abs2*Vm[BusIdx[lines[l].to]]^2
              	  + 2*Vm[BusIdx[lines[l].from]]*Vm[BusIdx[lines[l].to]]*(Yre*cos(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]])-Yim*sin(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]]))
              	)
                - flowmax <=0)
      end
      #branch apparent power limits (to bus)
      Ytf_abs2=YtfR[l]^2+YtfI[l]^2; Ytt_abs2=YttR[l]^2+YttI[l]^2
      Yre=YtfR[l]*YttR[l]+YtfI[l]*YttI[l]; Yim=-YtfR[l]*YttI[l]+YtfI[l]*YttR[l]
      if current_rating == true
        F_to[l] = @NLconstraint(opfmodel,
          	      1.0 *
                  ( Ytf_abs2*Vm[BusIdx[lines[l].from]]^2 + Ytt_abs2*Vm[BusIdx[lines[l].to]]^2
                    + 2*Vm[BusIdx[lines[l].from]]*Vm[BusIdx[lines[l].to]]*(Yre*cos(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]])-Yim*sin(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]]))
                  )
                  - flowmax <=0)
      else
        F_to[l] = @NLconstraint(opfmodel,
          	      Vm[BusIdx[lines[l].to]]^2 *
                  ( Ytf_abs2*Vm[BusIdx[lines[l].from]]^2 + Ytt_abs2*Vm[BusIdx[lines[l].to]]^2
                    + 2*Vm[BusIdx[lines[l].from]]*Vm[BusIdx[lines[l].to]]*(Yre*cos(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]])-Yim*sin(Va[BusIdx[lines[l].from]]-Va[BusIdx[lines[l].to]]))
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
  dP_dPg = zeros(nbus, ngen); dP_dQg = zeros(nbus, ngen)
  dP_dVm = Array{Union{Float64,JuMP.NonlinearExpression},2}(undef, nbus, nbus)
  dP_dVa = Array{Union{Float64,JuMP.NonlinearExpression},2}(undef, nbus, nbus)
  dP_dPd = zeros(nbus, nbus); dP_dQd = zeros(nbus, nbus)

  dQ_dPg = zeros(nbus, ngen); dQ_dQg = zeros(nbus, ngen)
  dQ_dVm = Array{Union{Float64,JuMP.NonlinearExpression},2}(undef, nbus, nbus)
  dQ_dVa = Array{Union{Float64,JuMP.NonlinearExpression},2}(undef, nbus, nbus)
  dQ_dPd = zeros(nbus, nbus); dQ_dQd = zeros(nbus, nbus)

  for i = 1:nbus # P, Q; equations
    for k = 1:nbus # Vm, Va; buses
     i = BusIdx[mod1(i, nbus)]
     k = BusIdx[mod1(k, nbus)]
     if i == k
       IDX = Y[i,:].nzind
       if !isempty(BusGeners[i]); dP_dPg[i, k] = -1.0; end
       dP_dPd[i, k] = 1.0
       dP_dVa[i, k] = Y[i,k] == 0 ? 0.0 : @NLexpression(opfmodel, -((Vm[i] * sum(Vm[kk] * ( G[i,kk] * sin(Va[i]-Va[kk]) - B[i,kk] * cos(Va[i]-Va[kk]) ) for kk in IDX))      ) - (B[i,i] * Vm[i]^2))
       dP_dVm[i, k] = Y[i,k] == 0 ? 0.0 : @NLexpression(opfmodel,  ((Vm[i] * sum(Vm[kk] * ( G[i,kk] * cos(Va[i]-Va[kk]) + B[i,kk] * sin(Va[i]-Va[kk]) ) for kk in IDX))/Vm[i]) + (G[i,i] * Vm[i]))
       dQ_dVa[i, k] = Y[i,k] == 0 ? 0.0 : @NLexpression(opfmodel,  ((Vm[i] * sum(Vm[kk] * ( G[i,kk] * cos(Va[i]-Va[kk]) + B[i,kk] * sin(Va[i]-Va[kk]) ) for kk in IDX))      ) - (G[i,i] * Vm[i]^2))
       dQ_dVm[i, k] = Y[i,k] == 0 ? 0.0 : @NLexpression(opfmodel,  ((Vm[i] * sum(Vm[kk] * ( G[i,kk] * sin(Va[i]-Va[kk]) - B[i,kk] * cos(Va[i]-Va[kk]) ) for kk in IDX))/Vm[i]) - (B[i,i] * Vm[i]))
       if !isempty(BusGeners[i]); dQ_dQg[i, k] = -1.0; end
       dQ_dQd[i, k] = 1.0
     else
       dP_dVa[i, k] = Y[i,k] == 0 ? 0.0 : @NLexpression(opfmodel,  Vm[i] * Vm[k] * ( G[i,k] * sin(Va[i]-Va[k]) - B[i,k] * cos(Va[i]-Va[k]) ))
       dP_dVm[i, k] = Y[i,k] == 0 ? 0.0 : @NLexpression(opfmodel,  Vm[i]         * ( G[i,k] * cos(Va[i]-Va[k]) + B[i,k] * sin(Va[i]-Va[k]) ))
       dQ_dVa[i, k] = Y[i,k] == 0 ? 0.0 : @NLexpression(opfmodel, -Vm[i] * Vm[k] * ( G[i,k] * cos(Va[i]-Va[k]) + B[i,k] * sin(Va[i]-Va[k]) ))
       dQ_dVm[i, k] = Y[i,k] == 0 ? 0.0 : @NLexpression(opfmodel,  Vm[i]         * ( G[i,k] * sin(Va[i]-Va[k]) - B[i,k] * cos(Va[i]-Va[k]) ))
     end
    end
  end
  #### aggregate components
  J = Array{Union{Float64,JuMP.NonlinearExpression}}(undef, 2nbus, 2ngen+2nbus+2nbus)
  J .= [ dP_dPg   dP_dQg   dP_dVm   dP_dVa   dP_dPd   dP_dQd;
         dQ_dPg   dQ_dQg   dQ_dVm   dQ_dVa   dQ_dPd   dQ_dQd ]
  #### partition aggregated Jacobian ([eval(Meta.parse(x)) for x in opfmodel.colNames[xidx]])
  dF_dx = J[Fidx, xidx]
  # dF_du = J[Fidx, uidx]
  dF_dy = J[Fidx, yidx]
  JuMP.registercon(opfmodel, :dF_dx, dF_dx)
  JuMP.registercon(opfmodel, :dF_dy, dF_dy)

  #
  # Gamma constraint
  #
  @constraintref Gamma_constraint[1:nx, 1:ny]
  for i = 1:nx
    for j = 1:ny
      Gamma_constraint[i, j] = @NLconstraint(opfmodel, sum(dF_dx[i, k] * Gamma[k, j] for k = 1:nx) + dF_dy[i, j] + zeta[i,j] == 0)
    end
  end
  JuMP.registercon(opfmodel, :Gamma_constraint, Gamma_constraint)

  #
  # Σ_x (!! NOTE: only computing the diagonal !!)
  #
  Σ_x = Array{JuMP.NonlinearExpression,1}(undef, nx)
  for k = 1:nx
    ## index of d elements in y
    d_mask = [(yy ∈ d) for yy in y]
    d_offset = sum([(yy ∉ d) for yy in y])
    didx = collect(1:ny)[d_mask]
    Σ_x[k] = @NLexpression(opfmodel, sum(Gamma[k,i] * Gamma[k,j] * Σ_d[i-d_offset, j-d_offset] for i in didx for j in didx))
  end

  #
  # single (Bonferroni) chance constraints (min)
  #
  @constraintref cc_Vm_max[1:length(b_RGL_idx[:L])]
  for i in eachindex(b_RGL_idx[:L])
    V = Vm[b_RGL_idx[:L]][i]
    idx = findall(V .== x)[1]
    η = 1.0 - (γ * epsilon_Vm) / nx
    q = quantile(z, η)
    cc_Vm_max[i] = @NLconstraint(opfmodel, Σ_x[idx] <= (buses.Vmax[idx] - V) / q)
  end
  JuMP.registercon(opfmodel, :cc_Vm_max, cc_Vm_max)

  @printf("Buses              : %d\n", nbus)
  @printf("Lines              : %d\n", nline)
  @printf("Generators         : %d\n", ngen)
  # @printf("Chance-constraints : %d\n", length(cc_Vm_max))
  println("Lines with limits  : ", nlinelim)

  return OPFModel(opfmodel, :InitData, :S)
end