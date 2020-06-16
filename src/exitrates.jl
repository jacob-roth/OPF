## -----------------------------------------------------------------------------
## Solve ACOPF with transition rate constraints for line current limits
## -----------------------------------------------------------------------------
function acopf_solve_exitrates(opfmodel::JuMP.Model, casedata, options::Dict=DefaultOptions(), adjustments::Dict=DefaultAdjustments(), warm_point_given=false, other=Dict())
    ## setup
    opfdata, physdata = casedata.opf, casedata.phys
    opfmodeldata = get_opfmodeldata(casedata, options, adjustments)
    opfmodeldata[:Y] = imag.(opfmodeldata[:Y])
    other[:opfmodeldata] = opfmodeldata
    lines        = opfmodeldata[:lines]
    busIdx       = opfmodeldata[:BusIdx]
    nonLoadBuses = opfmodeldata[:nonLoadBuses]

    ## solution vector
    solution = get_initial_point(opfmodel, opfdata, warm_point_given)


    ## Start
    iter = 0
    status = :IpoptInit
    minimize_cost = !options[:shed_load]
    lines_with_added_rate_constraints = Set{Int}()
    while iter < options[:iterlim]
        iter += 1

        ## set initial point
        set_initial_point!(opfmodel, solution)

        ## solve
        status = solve(opfmodel)
        println("\nITER   = $(iter) ")
        println("STATUS = $(string(status))\n")
        println()
        if status != :Optimal
            break
        end

        ## Get optimal values
        solution = get_optimal_values(opfmodel, opfmodeldata)

        ## loop over lines and check exit rates
        pl = deepcopy(options[:print_level])
        options[:print_level] = 0
        if options[:parallel]
            updated = compute_add_exitrates_parallel(solution, opfmodel, opfmodeldata, lines_with_added_rate_constraints, options)
        else
            updated = compute_add_exitrates_serial(solution, opfmodel, opfmodeldata, lines_with_added_rate_constraints, options)
        end
        options[:print_level] = pl

        ## Fix load shedding and change objective to minimizing cost
        if !updated
            minimize_cost && break
            minimize_cost = true
            if isinf(options[:VOLL])
                setlowerbound.(opfmodel[:Ps], solution[:Ps])
                setupperbound.(opfmodel[:Ps], solution[:Ps])
                setlowerbound.(opfmodel[:Qs], solution[:Qs])
                setupperbound.(opfmodel[:Qs], solution[:Qs])
                @NLobjective(opfmodel, Min, sum( opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n-2]*(opfmodeldata[:baseMVA]*opfmodel[:Pg][i])^2
                                            +opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n-1]*(opfmodeldata[:baseMVA]*opfmodel[:Pg][i])
                                            +opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n  ] for i=1:length(opfmodeldata[:generators])))
            else
                @NLobjective(opfmodel, Min, (sum( opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n-2]*(opfmodeldata[:baseMVA]*opfmodel[:Pg][i])^2
                                            +opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n-1]*(opfmodeldata[:baseMVA]*opfmodel[:Pg][i])
                                            +opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n  ] for i=1:length(opfmodeldata[:generators]))/options[:VOLL]) +
                                            sum((opfmodel[:Ps][b]+opfmodel[:Qs][b]) for b in 1:opfmodeldata[:nbus]))
            end
        end
    end

    ## Recompute all exitrates using exact calculation
    pl = deepcopy(options[:print_level])
    options[:print_level] = 0
    if options[:parallel]
        compute_exitrate_exact_all_parallel(solution, opfmodeldata, options, other)
    else
        compute_exitrate_exact_all_serial(solution, opfmodeldata, options, other)
    end
    options[:print_level] = pl

    return (opfmodel, status), other
end
function acopf_solve_exitrates(M, casedata, options::Dict=DefaultOptions(), adjustments::Dict=DefaultAdjustments(), warm_point_given=false)
    opfm = acopf_solve_exitrates(M.m, casedata, options, adjustments, warm_point_given, M.other)
    return OPFModel(opfm[1]..., M.kind, M.other)
end

## -----------------------------------------------------------------------------
## Helper functions for control flow of the main loop
## -----------------------------------------------------------------------------
function compute_add_exitrates_parallel(solution::Dict, opfmodel::JuMP.Model, opfmodeldata::Dict, lines_with_added_rate_constraints, options::Dict)
    updated = false
    lines = opfmodeldata[:lines]
    lines_to_check = symdiff(1:length(lines), lines_with_added_rate_constraints)

    ## compute the exit rate of lines in parallel
    exitrate = SharedVector{Float64}(length(lines))
    #tot = length(lines_to_check)
    function ratecalc(l)
        exit_point = compute_exitrate_kkt(l, solution, opfmodeldata, options)
        if exit_point !== nothing
            exitrate[l] = exit_point[:prefactor] * exit_point[:expterm]
            #println("$l/$tot = ", exitrate[l])
        end
    end
    pmap(ratecalc, lines_to_check)

    # @sync for l in lines_to_check
    #     function ratecalc()
    #         exit_point = compute_exitrate_kkt(l, solution, opfmodeldata, options)
    #         if exit_point !== nothing
    #             exitrate[l] = exit_point[:prefactor] * exit_point[:expterm]
    #         end
    #     end
    #     @async @spawn ratecalc()
    # end

    ## re-compute the exit rates of violated lines (to get starting points)
    for l in lines_to_check
        ## check if exit rate is within threshold
        if exitrate[l] <= options[:ratelimit]
            continue
        end

        ## add exit rate constraint
        updated = true
        exit_point = compute_exitrate_kkt(l, solution, opfmodeldata, options)
        push!(lines_with_added_rate_constraints, l)
        println("adding constraint for line $l")
        add_exitrate_constraint!(l, exit_point, opfmodel, opfmodeldata, options)
    end

    return updated
end

function compute_add_exitrates_serial(solution::Dict, opfmodel::JuMP.Model, opfmodeldata::Dict, lines_with_added_rate_constraints, options::Dict)
    updated = false
    for l in 1:length(opfmodeldata[:lines])
        ## Skip if we have already checked this line
        (l in lines_with_added_rate_constraints) && continue

        ## compute the exit rate of this line
        exit_point = compute_exitrate_kkt(l, solution, opfmodeldata, options)
        (exit_point == nothing) && continue
        exitrate = exit_point[:prefactor] * exit_point[:expterm]

        ## check if exit rate is within threshold
        if exitrate <= options[:ratelimit]
            continue
        end

        ## add exit rate constraint
        updated = true
        push!(lines_with_added_rate_constraints, l)
        println("adding constraint for line $l")
        add_exitrate_constraint!(l, exit_point, opfmodel, opfmodeldata, options)
    end
    return updated
end

function compute_exitrate_exact_all_parallel(solution::Dict, opfmodeldata::Dict, options::Dict, result::Dict)
    lines = opfmodeldata[:lines]
    rates = SharedVector{Float64}(length(lines))
    prefactors = SharedVector{Float64}(length(lines))
    expterms = SharedVector{Float64}(length(lines))
    rates_kkt = SharedVector{Float64}(length(lines))

    function ratecalc(l)
        exit_point = compute_exitrate_kkt(l, solution, opfmodeldata, options)
        if exit_point !== nothing
            rates_kkt[l] = exit_point[:prefactor] * exit_point[:expterm]
        end

        ep2 = compute_exitrate_exact(l, solution, opfmodeldata, options)  ## NOTE: JR - USE THIS ONE
        if ep2 !== nothing
            rates[l] = ep2[:prefactor] * ep2[:expterm]
            prefactors[l] = ep2[:prefactor]
            expterms[l] = ep2[:expterm]
        end
    end
    pmap(ratecalc, 1:length(lines))

    # @sync for l in 1:length(lines)
    #     function ratecalc()
    #         exit_point = compute_exitrate_kkt(l, solution, opfmodeldata, options)
    #         if exit_point !== nothing
    #             rates_kkt[l] = exit_point[:prefactor] * exit_point[:expterm]
    #         end
    #
    #         ep2 = compute_exitrate_exact(l, solution, opfmodeldata, options)  ## NOTE: JR - USE THIS ONE
    #         if ep2 !== nothing
    #             rates[l] = ep2[:prefactor] * ep2[:expterm]
    #             prefactors[l] = ep2[:prefactor]
    #             expterms[l] = ep2[:expterm]
    #         end
    #     end
    #     @async @spawn ratecalc()
    # end

    for l in eachindex(lines)
        @printf("Compare ---> line %4d: exact = %10.2e, log(approx/exact) = %10.2e\n", l, rates[l], abs(log(rates_kkt[l]/rates[l])))
    end

    result[:rates] = rates
    result[:prefactors] = prefactors
    result[:expterms] = expterms
end

function compute_exitrate_exact_all_serial(solution::Dict, opfmodeldata::Dict, options::Dict, result::Dict)
    lines = opfmodeldata[:lines]
    rates = zeros(length(lines))
    prefactors = zeros(length(lines))
    expterms = zeros(length(lines))
    for l in eachindex(lines)
        exit_point = compute_exitrate_kkt(l, solution, opfmodeldata, options)
        (exit_point == nothing) && continue
        exitrate = exit_point[:prefactor] * exit_point[:expterm]

        ep2 = compute_exitrate_exact(l, solution, opfmodeldata, options)  ## NOTE: JR - USE THIS ONE
        if ep2 != nothing
            rates[l] = ep2[:prefactor] * ep2[:expterm]
            prefactors[l] = ep2[:prefactor]
            expterms[l] = ep2[:expterm]
            @printf("Compare ---> line %4d: exact = %10.2e, log(approx/exact) = %10.2e\n", l, rates[l], abs(log(exitrate/rates[l])))
        end
    end

    result[:rates] = rates
    result[:prefactors] = prefactors
    result[:expterms] = expterms
end

## -----------------------------------------------------------------------------
## helpful functions
## -----------------------------------------------------------------------------
function get_initial_point(opfmodel::JuMP.Model, opfdata, warm_point_given=false)
    #
    # initial point - needed especially for pegase cases
    #
    if warm_point_given == false
        Pg0,Qg0,Vm0,Va0 = acopf_initialPt_IPOPT(opfdata)
    else
        Pg0 = warm_point_given[:Pg]
        Qg0 = warm_point_given[:Qg]
        Vm0 = warm_point_given[:Vm]
        Va0 = warm_point_given[:Va]
    end

    solution = Dict()
    solution[:Pg] = Pg0
    solution[:Qg] = Qg0
    solution[:Vm] = Vm0
    solution[:Va] = Va0
    solution[:Ps] = zeros(length(opfdata.buses))
    solution[:Qs] = zeros(length(opfdata.buses))

    return solution
end

function set_initial_point!(opfmodel::JuMP.Model, warm_point::Dict)
    setvalue(getindex(opfmodel, :Pg), warm_point[:Pg])
    setvalue(getindex(opfmodel, :Qg), warm_point[:Qg])
    setvalue(getindex(opfmodel, :Vm), warm_point[:Vm])
    setvalue(getindex(opfmodel, :Va), warm_point[:Va])
    setvalue(getindex(opfmodel, :Ps), warm_point[:Ps])
    setvalue(getindex(opfmodel, :Qs), warm_point[:Qs])
end

function get_optimal_values(opfmodel::JuMP.Model, opfmodeldata::Dict)
    solution = Dict()
    NB = length(opfmodeldata[:buses])
    solution[:Pg_full] = zeros(length(opfmodeldata[:buses]))
    solution[:Qg_full] = zeros(length(opfmodeldata[:buses]))
    Pg = getvalue(getindex(opfmodel,:Pg))
    Qg = getvalue(getindex(opfmodel,:Qg))
    for x in enumerate(opfmodeldata[:generators])
        g,gid = x[2],x[1]
        idx = opfmodeldata[:BusIdx][g.bus]
        solution[:Pg_full][idx] = Pg[gid]
        solution[:Qg_full][idx] = Qg[gid]
    end
    solution[:Pg] = getvalue(getindex(opfmodel,:Pg))
    solution[:Qg] = getvalue(getindex(opfmodel,:Qg))
    solution[:Vm] = getvalue(getindex(opfmodel,:Vm))
    solution[:Va] = getvalue(getindex(opfmodel,:Va))
    solution[:Ps] = getvalue(getindex(opfmodel,:Ps))
    solution[:Qs] = getvalue(getindex(opfmodel,:Qs))

    # Get also power injections
    BusGeners = opfmodeldata[:BusGenerators]
    buses     = opfmodeldata[:buses]
    baseMVA   = opfmodeldata[:baseMVA]
    solution[:Pnet] = [reduce(+, solution[:Pg][g] for g in BusGeners[i]; init=0.0) - ((buses[i].Pd - solution[:Ps][i]) / baseMVA) for i in 1:length(buses)]
    solution[:Qnet] = [reduce(+, solution[:Qg][g] for g in BusGeners[i]; init=0.0) - ((buses[i].Qd - solution[:Qs][i]) / baseMVA) for i in 1:length(buses)]

    # Get also voltages in the reduced space
    solution[:Vmr]  = getvalue(getindex(opfmodel,:Vm))
    solution[:Var]  = getvalue(getindex(opfmodel,:Va))
    for i in opfmodeldata[:nonLoadBuses][end:-1:1]
        splice!(solution[:Vmr], i)
    end
    splice!(solution[:Var], opfmodeldata[:bus_ref])

    # Finally get energy function information
    solution[:H_xbar] =          H([solution[:Vmr]; solution[:Var]]; solution=solution, opfmodeldata=opfmodeldata)
    solution[:grad_H] =  ∇H_direct([solution[:Vm];  solution[:Va]];  solution=solution, opfmodeldata=opfmodeldata)
    solution[:hess_H] = ∇2H_direct([solution[:Vm];  solution[:Va]];  solution=solution, opfmodeldata=opfmodeldata)

    return solution
end

function write_optimal_values(file::String, optimal_values::Dict)
    for k in keys(optimal_values)
        if !isdir("$(file)")
            mkdir("$(file)")
        end
        if isa(optimal_values[k], String)
            write("$(file)$(string(k)).csv", optimal_values[k], '\n')
        else
            writedlm("$(file)$(string(k)).csv", optimal_values[k])
        end
        # open("$(file)$(string(k)).csv", "w") do io
        #     if isa(optimal_values[k], String)
        #         write(io, optimal_values[k], '\n')
        #     else
        #         writedlm(io, optimal_values[k])
        #     end
        # end
    end
end

## -----------------------------------------------------------------------------
## Autodiff calculations
## -----------------------------------------------------------------------------
function H(xr; solution, opfmodeldata)
    Vm = xr[1:length(solution[:Vmr])]
    Va = xr[length(solution[:Vmr])+1:end]
    splice!(Va, opfmodeldata[:bus_ref]:opfmodeldata[:bus_ref]-1, solution[:Va][opfmodeldata[:bus_ref]])
    for i in opfmodeldata[:nonLoadBuses]
        splice!(Vm, i:i-1, solution[:Vm][i])
    end
    VMcosθ = Vm .* cos.(Va)
    VMsinθ = Vm .* sin.(Va)
    return -0.5 * ((VMcosθ' * opfmodeldata[:Y] * VMcosθ) + (VMsinθ' * opfmodeldata[:Y] * VMsinθ)) - (solution[:Pnet]' * Va) - (solution[:Qnet]' * log.(Vm))
end

function ∇H(xr; solution, opfmodeldata)
    g = xr -> H(xr; solution=solution, opfmodeldata=opfmodeldata)
    return ForwardDiff.gradient(g, xr)
end

function ∇2H(xr; solution, opfmodeldata)
    g = xr -> H(xr; solution=solution, opfmodeldata=opfmodeldata)
    return ForwardDiff.hessian(g, xr)
end

function h(xr; solution, opfmodeldata, i, j, flowmax)
    Vm = xr[1:length(solution[:Vmr])]
    Va = xr[length(solution[:Vmr])+1:end]
    splice!(Va, opfmodeldata[:bus_ref]:opfmodeldata[:bus_ref]-1, solution[:Va][opfmodeldata[:bus_ref]])
    for b in opfmodeldata[:nonLoadBuses]
        splice!(Vm, b:b-1, solution[:Vm][b])
    end
    return 0.5*(Vm[i]^2 + Vm[j]^2 - (2*Vm[i]*Vm[j]*cos(Va[i]-Va[j])) - flowmax)
end

function ∇h(xr; solution, opfmodeldata, i, j, flowmax)
    g = xr -> h(xr; solution=solution, opfmodeldata=opfmodeldata, i=i, j=j, flowmax=flowmax)
    return ForwardDiff.gradient(g, xr)
end

function ∇2h(xr; solution, opfmodeldata, i, j, flowmax)
    g = xr -> h(xr; solution=solution, opfmodeldata=opfmodeldata, i=i, j=j, flowmax=flowmax)
    return ForwardDiff.hessian(g, xr)
end

## -----------------------------------------------------------------------------
## Explicit/direct gradient and Hessian calculations
## -----------------------------------------------------------------------------
function ∇H_direct(xfullspace; solution, opfmodeldata)
    g = zeros(opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1)
    VM = xfullspace[1:opfmodeldata[:nbus]]
    VA = xfullspace[opfmodeldata[:nbus]+1:2opfmodeldata[:nbus]]
    for n in 1:opfmodeldata[:nbus]
        nv = opfmodeldata[:linindex_VMr][n]
        nθ = opfmodeldata[:linindex_VAr][n]
        if nv > 0
            g[nv] = -(+(VM[n]*opfmodeldata[:Y][n,n])
                      +reduce(+, opfmodeldata[:Y][n,k]*VM[k]*cos(VA[n] - VA[k]) for k in 1:opfmodeldata[:nbus] if k != n && opfmodeldata[:Y][n,k] != 0)
                      +(solution[:Qnet][n]/VM[n]))
        end
        if nθ > 0
            g[nθ] = -solution[:Pnet][n]+reduce(+, opfmodeldata[:Y][n,k]*VM[n]*VM[k]*sin(VA[n] - VA[k]) for k in 1:opfmodeldata[:nbus] if k != n && opfmodeldata[:Y][n,k] != 0)
        end
    end
    return g
end

function ∇2H_direct_old(xfullspace; solution, opfmodeldata)
    #A = Array{typeof(zero(xfullspace[1])), 2}(undef, opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1, opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1);
    A = zeros(opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1, opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1)
    VM = xfullspace[1:opfmodeldata[:nbus]]
    VA = xfullspace[opfmodeldata[:nbus]+1:2opfmodeldata[:nbus]]
    i, j = 0, 0
    for n in 1:opfmodeldata[:nbus]
        (n in opfmodeldata[:nonLoadBuses]) && continue
        j = i; i += 1
        for b in 1:opfmodeldata[:nbus]
            (b in opfmodeldata[:nonLoadBuses] || b < n) && continue
            j += 1
            A[i, j] = (b == n) ? (-opfmodeldata[:Y][n,n] + (solution[:Qnet][n] / VM[n]^2)) :
                                 (-opfmodeldata[:Y][n,b] * cos(VA[n] - VA[b]))
        end
        for b in 1:opfmodeldata[:nbus]
            (b == opfmodeldata[:bus_ref]) && continue
            j += 1
            A[i, j] = (b == n) ? reduce(+, opfmodeldata[:Y][n,k]*VM[k]*sin(VA[n] - VA[k]) for k in 1:opfmodeldata[:nbus]) :
                                 (-opfmodeldata[:Y][n,b] * VM[b] * sin(VA[n] - VA[b]))
        end
    end
    for n in 1:opfmodeldata[:nbus]
        (n == opfmodeldata[:bus_ref]) && continue
        j = i; i += 1
        for b in 1:opfmodeldata[:nbus]
            (b == opfmodeldata[:bus_ref] || b < n) && continue
            j += 1
            A[i, j] = (b == n) ? reduce(+, VM[k]*VM[n]*opfmodeldata[:Y][n,k]*cos(VA[n] - VA[k]) for k in 1:opfmodeldata[:nbus] if k != n) :
                                 (-opfmodeldata[:Y][n,b] * VM[n] * VM[b] * cos(VA[n] - VA[b]))
        end
    end
    A = Symmetric(A)
    return A
end

function ∇2H_direct(xfullspace; solution, opfmodeldata)
    #A = Array{typeof(zero(xfullspace[1])), 2}(undef, opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1, opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1);
    A = zeros(opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1, opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1)
    VM = xfullspace[1:opfmodeldata[:nbus]]
    VA = xfullspace[opfmodeldata[:nbus]+1:2opfmodeldata[:nbus]]
    for n in 1:opfmodeldata[:nbus]
        nv = opfmodeldata[:linindex_VMr][n]
        nθ = opfmodeldata[:linindex_VAr][n]
        if nv > 0
            A[nv,nv] = -opfmodeldata[:Y][n,n] + (solution[:Qnet][n] / VM[n]^2)
            A[nv,nθ] = reduce(+, opfmodeldata[:Y][n,k]*VM[k]*sin(VA[n] - VA[k]) for k in 1:opfmodeldata[:nbus] if opfmodeldata[:Y][n,k] != 0)
        end
        if nθ > 0
            A[nθ,nθ] = reduce(+, VM[k]*VM[n]*opfmodeldata[:Y][n,k]*cos(VA[n] - VA[k]) for k in 1:opfmodeldata[:nbus] if k != n && opfmodeldata[:Y][n,k] != 0)
            for b in 1:opfmodeldata[:nbus]
                bv = opfmodeldata[:linindex_VMr][b]
                bθ = opfmodeldata[:linindex_VAr][b]
                if opfmodeldata[:Y][n,b] != 0
                    if nv > 0 && bv > 0 && b > n
                        A[nv,bv] = -opfmodeldata[:Y][n,b] * cos(VA[n] - VA[b])
                    end
                    if nv > 0 && bθ > 0 && b != n
                        A[nv,bθ] = -opfmodeldata[:Y][n,b] * VM[b] * sin(VA[n] - VA[b])
                    end
                    if bθ > 0 && b > n
                        A[nθ,bθ] = -opfmodeldata[:Y][n,b] * VM[n] * VM[b] * cos(VA[n] - VA[b])
                    end
                end
            end
        end
    end
    A = Symmetric(A)
    return A
end

function ∇h_direct(xfullspace; solution, opfmodeldata, i, j, flowmax)
    g = zeros(opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1)
    a = opfmodeldata[:linindex_VMr][i]
    b = opfmodeldata[:linindex_VMr][j]
    c = opfmodeldata[:linindex_VAr][i]
    d = opfmodeldata[:linindex_VAr][j]
    VM = xfullspace[1:opfmodeldata[:nbus]]
    VA = xfullspace[opfmodeldata[:nbus]+1:2opfmodeldata[:nbus]]
    (a > 0) && (g[a] = +VM[i] - (VM[j]*cos(VA[i] - VA[j])))
    (b > 0) && (g[b] = +VM[j] - (VM[i]*cos(VA[i] - VA[j])))
    (c > 0) && (g[c] = +VM[i]*VM[j]*sin(VA[i] - VA[j]))
    (d > 0) && (g[d] = -VM[i]*VM[j]*sin(VA[i] - VA[j]))
    return g
end

function ∇2h_direct(xfullspace; solution, opfmodeldata, i, j, flowmax)
    A = zeros(opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1, opfmodeldata[:nloads] + opfmodeldata[:nbus] - 1);
    a = opfmodeldata[:linindex_VMr][i]
    b = opfmodeldata[:linindex_VMr][j]
    c = opfmodeldata[:linindex_VAr][i]
    d = opfmodeldata[:linindex_VAr][j]
    VM = xfullspace[1:opfmodeldata[:nbus]]
    VA = xfullspace[opfmodeldata[:nbus]+1:2opfmodeldata[:nbus]]
    (a > 0) && (A[a,a] = 1)
    (b > 0) && (A[b,b] = 1)
    (c > 0) && (A[c,c] = VM[i]*VM[j]*cos(VA[i] - VA[j]))
    (d > 0) && (A[d,d] = VM[i]*VM[j]*cos(VA[i] - VA[j]))
    if a > 0 && b > 0
        A[a,b] = -cos(VA[i] - VA[j])
        A[b,a] = A[a,b]
    end
    if a > 0 && c > 0
        A[a,c] = VM[j]*sin(VA[i] - VA[j])
        A[c,a] = A[a,c]
    end
    if a > 0 && d > 0
        A[a,d] = -VM[j]*sin(VA[i] - VA[j])
        A[d,a] = A[a,d]
    end
    if b > 0 && c > 0
        A[b,c] = VM[i]*sin(VA[i] - VA[j])
        A[c,b] = A[b,c]
    end
    if b > 0 && d > 0
        A[b,d] = -VM[i]*sin(VA[i] - VA[j])
        A[d,b] = A[b,d]
    end
    if c > 0 && d > 0
        A[c,d] = -VM[i]*VM[j]*cos(VA[i] - VA[j])
        A[d,c] = A[c,d]
    end
    return A
end

## -----------------------------------------------------------------------------
## Main constraints
## -----------------------------------------------------------------------------
function add_kkt_stationarity_constraint!(opfmodel::JuMP.Model, opfmodeldata::Dict, i::Int, j::Int, Vmbar, Vabar, Qsbar, Vmstar, Vastar, Kstar)
    nbus = opfmodeldata[:nbus]
    buses = opfmodeldata[:buses]
    nonLoadBuses = opfmodeldata[:nonLoadBuses]
    bus_ref = opfmodeldata[:bus_ref]
    Y = opfmodeldata[:Y]
    S = opfmodeldata[:S]

    Vmi = (i in nonLoadBuses) ? Vmbar[i] : Vmstar[i]
    Vmj = (j in nonLoadBuses) ? Vmbar[j] : Vmstar[j]
    Vai = (i == bus_ref     ) ? Vabar[i] : Vastar[i]
    Vaj = (j == bus_ref     ) ? Vabar[j] : Vastar[j]

    #
    # Vm component
    #
    for n in 1:nbus
        n in nonLoadBuses && continue
        if n == i
            rhs = @NLexpression(opfmodel, Kstar*(Vmi - Vmj*cos(Vai-Vaj)))
        elseif n == j
            rhs = @NLexpression(opfmodel, Kstar*(Vmj - Vmi*cos(Vai-Vaj)))
        else
            rhs = 0.0
        end
        @NLconstraint(opfmodel,
                 -(Vmstar[n] - Vmbar[n])*(Y[n,n]+ ((buses[n].Qd - Qsbar[n]) / (Vmbar[n]^2 * opfmodeldata[:baseMVA])))
            - sum((Vmstar[b] - Vmbar[b])* Y[n,b]*         cos(Vabar[n] - Vabar[b]) for b in 1:nbus if Y[n,b] != 0 && b ∉ nonLoadBuses && b != n)
            + sum((Vastar[n] - Vabar[n])* Y[n,k]*Vmbar[k]*sin(Vabar[n] - Vabar[k]) for k in 1:nbus if Y[n,k] != 0)
            - sum((Vastar[b] - Vabar[b])* Y[n,b]*Vmbar[b]*sin(Vabar[n] - Vabar[b]) for b in 1:nbus if Y[n,b] != 0 && b != bus_ref && b != n)
            - rhs == 0
        )
    end
    #
    # Va component
    #
    for n in 1:nbus
        n == bus_ref && continue
        if n == i
            rhs = @NLexpression(opfmodel, Kstar*Vmi*Vmj*sin(Vai-Vaj))
        elseif n == j
            rhs = @NLexpression(opfmodel, Kstar*Vmi*Vmj*sin(Vaj-Vai))
        else
            rhs = 0.0
        end
        @NLconstraint(opfmodel,
              sum((Vmstar[n] - Vmbar[n])*Y[n,k]*Vmbar[k]*         sin(Vabar[n] - Vabar[k]) for k in 1:nbus if Y[n,k] != 0 && n ∉ nonLoadBuses)
            + sum((Vmstar[b] - Vmbar[b])*Y[n,b]*Vmbar[n]*         sin(Vabar[n] - Vabar[b]) for b in 1:nbus if Y[b,n] != 0 && b ∉ nonLoadBuses && b != n)
            + sum((Vastar[n] - Vabar[n])*Y[n,k]*Vmbar[k]*Vmbar[n]*cos(Vabar[n] - Vabar[k]) for k in 1:nbus if Y[n,k] != 0 && k != n)
            - sum((Vastar[b] - Vabar[b])*Y[n,b]*Vmbar[n]*Vmbar[b]*cos(Vabar[n] - Vabar[b]) for b in 1:nbus if Y[n,b] != 0 && b != bus_ref && b != n)
            - rhs == 0)
    end
end

function add_exitrate_constraint!(l::Int, exit_point::Dict, opfmodel::JuMP.Model, opfmodeldata::Dict, options::Dict)
    nbus         = opfmodeldata[:nbus]
    buses        = opfmodeldata[:buses]
    nonLoadBuses = opfmodeldata[:nonLoadBuses]
    bus_ref      = opfmodeldata[:bus_ref]
    baseMVA      = opfmodeldata[:baseMVA]
    Y            = opfmodeldata[:Y]
    S            = opfmodeldata[:S]
    line         = opfmodeldata[:lines][l]
    nrow         = 2nbus - length(nonLoadBuses) - 1
    # flowmax      = (options[:constr_limit_scale]*line.rateA/(abs(1.0/(line.x*im))*baseMVA))^2
    i            = opfmodeldata[:BusIdx][line.from]
    j            = opfmodeldata[:BusIdx][line.to]
    flowmax      = options[:constr_limit_scale]*(line.rateA^2 / (opfmodeldata[:baseMVA]^2 * abs2(1.0 / line.x)))
    Vm           = opfmodel[:Vm]
    Va           = opfmodel[:Va]
    Qs           = opfmodel[:Qs]

    (i in nonLoadBuses && j in nonLoadBuses) && return nothing


    #
    # Define variables and Lagrange multiplier
    #
    VMstar = @variable(opfmodel, [b=1:nbus; b ∉ nonLoadBuses])
    VAstar = @variable(opfmodel, [b=1:nbus; b != bus_ref])
    Kstar  = @variable(opfmodel)
    for b in 1:nbus
        (b ∉ nonLoadBuses) && setvalue(VMstar[b], exit_point[:Vm][b])
        (b != bus_ref)     && setvalue(VAstar[b], exit_point[:Va][b])
        (b != bus_ref)     && setlowerbound(VAstar[b], -pi)
        (b != bus_ref)     && setupperbound(VAstar[b], +pi)
    end
    setlowerbound(Kstar, 0)
    setvalue(Kstar, exit_point[:K])

    Vmi = (i in nonLoadBuses) ? Vm[i] : VMstar[i]
    Vmj = (j in nonLoadBuses) ? Vm[j] : VMstar[j]
    Vai = (i == bus_ref)      ? Va[i] : VAstar[i]
    Vaj = (j == bus_ref)      ? Va[j] : VAstar[j]

    #
    # Determine type of line and the LDL factorization of its constraint hessian
    #
    a = opfmodeldata[:linindex_VMr][i]
    b = opfmodeldata[:linindex_VMr][j]
    c = opfmodeldata[:linindex_VAr][i]
    d = opfmodeldata[:linindex_VAr][j]
    type = "nonload-nonload"
    if a > 0 && b > 0
        type = "load-load"
        ncol = 3
        L    = Array{JuMP.NonlinearExpression,2}(undef, nrow, ncol)
        D    = Array{JuMP.NonlinearExpression,2}(undef, ncol, ncol)
        L   .= @NLexpression(opfmodel, 0)
        get_LDL_factors_load_load!(L, D, Vmi, Vai, Vmj, Vaj, a, c, b, d; model = opfmodel)
    else
        ncol = 2
        L    = Array{JuMP.NonlinearExpression,2}(undef, nrow, ncol)
        D    = Array{JuMP.NonlinearExpression,2}(undef, ncol, ncol)
        L   .= @NLexpression(opfmodel, 0)
        if a > 0 && d > 0
            type = "load-gen"
            get_LDL_factors_load_gen!(L, D, Vmi, Vai, Vmj, Vaj, a, c, b, d; model = opfmodel)
        elseif a > 0
            type = "load-slack"
            get_LDL_factors_load_slk!(L, D, Vmi, Vai, Vmj, Vaj, a, c, b, d; model = opfmodel)
        elseif b > 0 && c > 0
            type = "gen-load"
            get_LDL_factors_load_gen!(L, D, Vmj, Vaj, Vmi, Vai, b, d, a, c; model = opfmodel)
        elseif b > 0
            type = "slack-load"
            get_LDL_factors_load_slk!(L, D, Vmj, Vaj, Vmi, Vai, b, d, a, c; model = opfmodel)
        else
            @assert false
        end
    end
    if options[:print_level] >= 1
        println("Adding exit rate constraint for line ", l, " = (", i, ",", j, ") of type = ", type)
    end

    #
    # KKT conditions - primal feasibility
    #
    # @NLconstraint(opfmodel, baseMVA*(Vmi^2 + Vmj^2 - (2*Vmi*Vmj*cos(Vai-Vaj)) - flowmax) == 0)
    @NLconstraint(opfmodel, baseMVA*(flowmax - (Vmi^2 + Vmj^2 - (2*Vmi*Vmj*cos(Vai-Vaj)))) == 0)

    #
    # KKT conditions - stationarity
    #
    add_kkt_stationarity_constraint!(opfmodel, opfmodeldata, i, j, opfmodel[:Vm], opfmodel[:Va], opfmodel[:Qs], VMstar, VAstar, Kstar)


    #
    # Define variables Z such that hess_H * Z = L
    #   where L is from the LDL factor of constraint hessian hess_h
    #
    ZVM = @variable(opfmodel, [b=1:nbus, 1:ncol; b ∉ nonLoadBuses])
    ZVA = @variable(opfmodel, [b=1:nbus, 1:ncol; b != bus_ref])
    for col in 1:ncol
        # Vm component
        for n in 1:nbus
            n in nonLoadBuses && continue
            setvalue(ZVM[n,col], exit_point[:Z][n,col])
            @NLconstraint(opfmodel,
                     -ZVM[n,col]*(Y[n,n]+ ((buses[n].Qd - Qs[n]) / (Vm[n]^2 * baseMVA)))
                - sum(ZVM[b,col]* Y[n,b]*      cos(Va[n] - Va[b]) for b in 1:nbus if Y[n,b] != 0 && b ∉ nonLoadBuses && b != n)
                + sum(ZVA[n,col]* Y[n,k]*Vm[k]*sin(Va[n] - Va[k]) for k in 1:nbus if Y[n,k] != 0)
                - sum(ZVA[b,col]* Y[n,b]*Vm[b]*sin(Va[n] - Va[b]) for b in 1:nbus if Y[n,b] != 0 && b != bus_ref && b != n)
                - L[opfmodeldata[:linindex_VMr][n],col] == 0)
        end
        # Va component
        for n in 1:nbus
            n == bus_ref && continue
            setvalue(ZVA[n,col], exit_point[:Z][nbus+n,col])
            @NLconstraint(opfmodel,
                  sum(ZVM[n,col]*Y[n,k]*Vm[k]*      sin(Va[n] - Va[k]) for k in 1:nbus if Y[n,k] != 0 && n ∉ nonLoadBuses)
                + sum(ZVM[b,col]*Y[n,b]*Vm[n]*      sin(Va[n] - Va[b]) for b in 1:nbus if Y[b,n] != 0 && b ∉ nonLoadBuses && b != n)
                + sum(ZVA[n,col]*Y[n,k]*Vm[k]*Vm[n]*cos(Va[n] - Va[k]) for k in 1:nbus if Y[n,k] != 0 && k != n)
                - sum(ZVA[b,col]*Y[n,b]*Vm[n]*Vm[b]*cos(Va[n] - Va[b]) for b in 1:nbus if Y[n,b] != 0 && b != bus_ref && b != n)
                - L[opfmodeldata[:linindex_VAr][n],col] == 0)
        end
    end

    #
    # Compute ||∇h||^2 and (x⋆ - x̄)'*∇h
    #
    norm_squared_grad_h =
        @NLexpression(opfmodel,
            ((a > 0 ? 1.0 : 0.0)*( Vmi-(Vmj*cos(Vai - Vaj)))^2) +
            ((b > 0 ? 1.0 : 0.0)*( Vmj-(Vmi*cos(Vai - Vaj)))^2) +
            ((c > 0 ? 1.0 : 0.0)*( Vmi* Vmj*sin(Vai - Vaj))^2) +
            ((d > 0 ? 1.0 : 0.0)*(-Vmi* Vmj*sin(Vai - Vaj))^2)
            ## !NOTE! can access grad_H_xstarᵀ S grad_H_xstar elements of S with:
            ##        S[a]
            ##        S[b]
            ##        S[c]
            ##        S[d]
        )
    norm_squared_grad_h_weighted_S =
        @NLexpression(opfmodel,
            ((a > 0 ? S[a] : 0.0)*( Vmi-(Vmj*cos(Vai - Vaj)))^2) +
            ((b > 0 ? S[b] : 0.0)*( Vmj-(Vmi*cos(Vai - Vaj)))^2) +
            ((c > 0 ? S[c] : 0.0)*( Vmi* Vmj*sin(Vai - Vaj))^2) +
            ((d > 0 ? S[d] : 0.0)*(-Vmi* Vmj*sin(Vai - Vaj))^2)
        )
    xstar_minus_xbar_times_grad_h =
        @NLexpression(opfmodel,
            ((Vmi - Vm[i])*( Vmi-(Vmj*cos(Vai - Vaj)))) +
            ((Vmj - Vm[j])*( Vmj-(Vmi*cos(Vai - Vaj)))) +
            ((Vai - Va[i])*( Vmi* Vmj*sin(Vai - Vaj))) +
            ((Vaj - Va[j])*(-Vmi* Vmj*sin(Vai - Vaj)))
        )

    #
    # 1. Compute det(W) where W = I - KDT and T = L'Z
    #       note that T is symmetric
    #       T[m,n] = L[a,m]*ZVM[i,n] + L[b,m]*ZVM[j,n] + L[c,m]*ZVA[i,n] + L[d,m]*ZVA[j,n]
    #
    #
    # 2. Compute quad = det(-D)*J'C(X)J
    #       where
    #           J    = L'(x⋆ - x̄)
    #           J[n] = L[a,n]*(Vmi - Vm[i]) + L[b,n]*(Vmj - Vm[j]) + L[c,n]*(Vai - Va[i]) + L[d,n]*(Vaj - Va[j])
    #           C(X) = cofactor matrix of X = -inv(D) + KT (T defined as above)
    #
    T = Array{JuMP.NonlinearExpression, 2}(undef, ncol, ncol)
    J = Array{JuMP.NonlinearExpression, 1}(undef, ncol)
    for n=1:ncol, m=n:ncol
        if a > 0 && b > 0
            T[m,n] = @NLexpression(opfmodel, (L[a,m]*ZVM[i,n])    +(L[b,m]*ZVM[j,n])    +(L[c,m]*ZVA[i,n])      +(L[d,m]*ZVA[j,n])     ); (m != n) && continue
            J[m]   = @NLexpression(opfmodel, (L[a,m]*(Vmi-Vm[i])) +(L[b,m]*(Vmj-Vm[j])) +(L[c,m]*(Vai - Va[i])) +(L[d,m]*(Vaj - Va[j])))
        elseif a > 0 && d > 0
            T[m,n] = @NLexpression(opfmodel, (L[a,m]*ZVM[i,n])                          +(L[c,m]*ZVA[i,n])      +(L[d,m]*ZVA[j,n])     ); (m != n) && continue
            J[m]   = @NLexpression(opfmodel, (L[a,m]*(Vmi-Vm[i]))                       +(L[c,m]*(Vai - Va[i])) +(L[d,m]*(Vaj - Va[j])))
        elseif a > 0
            T[m,n] = @NLexpression(opfmodel, (L[a,m]*ZVM[i,n])                          +(L[c,m]*ZVA[i,n])                             ); (m != n) && continue
            J[m]   = @NLexpression(opfmodel, (L[a,m]*(Vmi-Vm[i]))                       +(L[c,m]*(Vai - Va[i]))                        )
        elseif b > 0 && c > 0
            T[m,n] = @NLexpression(opfmodel,                      +(L[b,m]*ZVM[j,n])    +(L[c,m]*ZVA[i,n])      +(L[d,m]*ZVA[j,n])     ); (m != n) && continue
            J[m]   = @NLexpression(opfmodel,                      +(L[b,m]*(Vmj-Vm[j])) +(L[c,m]*(Vai - Va[i])) +(L[d,m]*(Vaj - Va[j])))
        elseif b > 0
            T[m,n] = @NLexpression(opfmodel,                      +(L[b,m]*ZVM[j,n])                            +(L[d,m]*ZVA[j,n])     ); (m != n) && continue
            J[m]   = @NLexpression(opfmodel,                      +(L[b,m]*(Vmj-Vm[j]))                         +(L[d,m]*(Vaj - Va[j])))
        end
    end
    if type == "load-load"
        detW = @NLexpression(opfmodel,
                    -prod(-1 + (Kstar*D[n,n]*T[n,n]) for n in 1:3) -
                    (2*Kstar^3*D[1,1]*D[2,2]*D[3,3]*T[2,1]*T[3,1]*T[3,2]) -
                    (Kstar^2*
                        sum(T[m,n]^2*D[m,m]*D[n,n]*(1 - (Kstar*D[6-m-n,6-m-n]*T[6-m-n,6-m-n]))
                            for n=1:3 for m=n+1:3)
                    )
                )
        #
        # C is the (symmetric) cofactor matrix times det(-D)
        #
        C = Array{JuMP.NonlinearExpression, 2}(undef, 3, 3)
        for n=1:3
            C[n,n] = @NLexpression(opfmodel, D[n,n]*(
                            sum(D[e,e]*D[f,f]*Kstar^2*T[e,f]^2 for f in 1:3 for e in f+1:3 if e != n && f != n) -
                            prod(1 - (Kstar*D[m,m]*T[m,m]) for m in 1:3 if m != n)
                        )
                      )
        end
        for n=1:3, m=n+1:3
            C[m,n] = @NLexpression(opfmodel, (-1)^(m+n)*Kstar*D[m,m]*D[n,n]*(T[m,n] -
                            (
                                Kstar*D[6-m-n,6-m-n]*((T[m,n]*T[6-m-n,6-m-n]) -
                                    prod(T[e,f] for f in 1:3 for e in f+1:3 if e != m || f != n))
                            )
                        )
                      )
        end
        quad = @NLexpression(opfmodel,
                       sum(C[n,n]*(J[n])^2 for n in 1:2) +
                    (2*C[2,1]*prod(J[n]    for n in 1:2)) +
                    (2*C[3,2]*prod(J[n]    for n in 2:3)) +
                    (2*C[3,1]*prod(J[n]    for n in [1,3]))
                )
    else
        #
        # By exploiting symmetry of T, we have
        # det(I - KDT) = (1 - KD[1,1]T[1,1])(1 - KD[2,2]T[2,2]) - K^2 T[2,1]^2 D[1,1]D[2,2]
        #
        detW = @NLexpression(opfmodel,
                    ((1 - (Kstar*D[1,1]*T[1,1]))*(1 - (Kstar*D[2,2]*T[2,2]))) -
                    (Kstar^2*D[1,1]*D[2,2]*T[2,1]^2)
                )
        #
        # Again by explpoiting symmetry of T, we have
        # quad =    (P[1]^2*(-D[1,1] + Kstar*D[1,1]*D[2,2]*T[2,2])) +
        #         + (P[2]^2*(-D[2,2] + Kstar*D[1,1]*D[2,2]*T[1,1]))
        #         - 2*Kstar*T[2,1]*D[1,1]*D[2,2]*P[1]*P[2]
        #
        quad = @NLexpression(opfmodel,
                    sum(
                        (-D[n,n]+(Kstar*D[1,1]*D[2,2]*T[3-n,3-n])) * (J[n])^2
                        for n in 1:2)
                    -(
                        2*Kstar*D[1,1]*D[2,2]*T[2,1]*prod(J[n] for n in 1:2)
                     )
                )
    end


    #
    # the denominator in the expression for the failure rate
    #
    denom = @NLexpression(opfmodel,
                    (xstar_minus_xbar_times_grad_h*norm_squared_grad_h*detW)
                    -quad
                )

    #
    # The main constraint written in log form
    #
    #@NLconstraint(opfmodel, denom >= 0)
    if options[:psd_constraint]
        #--------------------------------------------------------------------------------
        # Enforcing minimization of lower-level problem is equivalent to enforcing
        # that the Hessian of its Lagrangian at the exit point is positive definite.
        #
        # ASSUMING that Hessian of the energy function at x_bar is positive definite,
        # this can be shown to be equivalent to enforcing the spectral radius of the matrix
        # X = KDT, where T is as defined above,
        # to be less than 1.
        #
        # However, even though X is a small kxk matrix, it is nonsymmetric
        # and there is no easy way to enforce this spectral constraint in general.
        #
        # When k=2, this is equivalent to enforcing:
        # tr(X) - det(X) <= 1
        #
        # When k=3, there are 2 strategies:
        #       (a) Implement a relaxed version of the constraint: tr(X) = sum of eigenvalues < k
        #       (b) Try and decompose hess h = LDL^T, where D > 0. However, this might not always be possible.
        #           In this case, we can replace the condition on X with that on (K*sqrt(D)*T*sqrt(D))
        #           which is symmetric and hence we can use Sylvester’s criterion of principal minors > 0.
        #
        # We adopt strategy (a) for simplicity
        #--------------------------------------------------------------------------------
        if ncol == 2
            @NLconstraint(opfmodel, sum(Kstar*D[m,m]*T[m,m] for m=1:ncol) -
                (Kstar^2*D[1,1]*D[2,2]*((T[1,1]*T[2,2]) - T[2,1]^2))
                <= 1
            )
        else
            @NLconstraint(opfmodel, sum(Kstar*D[m,m]*T[m,m] for m=1:ncol) <= ncol)
        end
    end
    if options[:high_temp_adj]
        @NLconstraint(opfmodel,
                        log(((Kstar^1.5)*(norm_squared_grad_h_weighted_S)/(denom^2 + 1.0e-16)^0.25) +
                            (2options[:temperature]/(Kstar*xstar_minus_xbar_times_grad_h))) -
                        (Kstar*xstar_minus_xbar_times_grad_h/2options[:temperature]) -
                        log(options[:ratelimit]*sqrt(2*pi*options[:temperature])/options[:damping])
                        <= 0
                    )
    else
        @NLconstraint(opfmodel,
                        1.5log(Kstar) + log(norm_squared_grad_h_weighted_S) -
                        (Kstar*xstar_minus_xbar_times_grad_h/2options[:temperature]) -
                        0.5log(sqrt(denom^2 + 1.0e-16)) -
                        log(options[:ratelimit]*sqrt(2*pi*options[:temperature])/options[:damping])
                        <= 0
                    )
    end
end

## -----------------------------------------------------------------------------
## transition rate calculations
## -----------------------------------------------------------------------------
function compute_exitrate_exact(l::Int, xbar::Dict, opfmodeldata::Dict, options::Dict)
    @assert options[:lossless]
    @assert options[:current_rating]

    buses        = opfmodeldata[:buses]
    line         = opfmodeldata[:lines][l]
    nonLoadBuses = opfmodeldata[:nonLoadBuses]
    bus_ref      = opfmodeldata[:bus_ref]
    Y            = opfmodeldata[:Y]
    S            = opfmodeldata[:S]
    # flowmax      = options[:constr_limit_scale]*(line.rateA/(abs(1.0/(line.x*im))*opfmodeldata[:baseMVA]))^2
    i            = opfmodeldata[:BusIdx][line.from]
    j            = opfmodeldata[:BusIdx][line.to]
    flowmax      = options[:constr_limit_scale]*(line.rateA^2 / (opfmodeldata[:baseMVA]^2 * abs2(1.0 / line.x)))
    nbus         = length(buses)
    nrow         = 2nbus - length(nonLoadBuses) - 1
    VMbar        = xbar[:Vm]
    VAbar        = xbar[:Va]
    baseMVA      = opfmodeldata[:baseMVA]



    ## check if exit rate formula is applicable
    if ((i in nonLoadBuses && j in nonLoadBuses) || line.rateA==0 || line.rateA>=1.0e10)
        return nothing
    end



    # --------------------
    # solve constrained energy minimization problem
    # --------------------
    em = Model(solver = IpoptSolver(print_level=options[:print_level]))

    #
    # Variables
    #
    @variable(em, Vm[1:nbus])
    @variable(em, -pi <= Va[1:nbus] <= pi)
    setlowerbound(Va[bus_ref], VAbar[bus_ref])
    setupperbound(Va[bus_ref], VAbar[bus_ref])
    setlowerbound(Vm[bus_ref], VMbar[bus_ref])
    setupperbound(Vm[bus_ref], VMbar[bus_ref])
    for b in 1:nbus
        if b in nonLoadBuses
            setlowerbound(Vm[b], VMbar[b])
            setupperbound(Vm[b], VMbar[b])
        end
    end

    #
    # Set objective to minimize energy function
    #
    @NLobjective(em, Min,
        -0.5*sum(Y[m,n]*Vm[m]*Vm[n]*cos(Va[m] - Va[n]) for m in 1:nbus for n in 1:nbus if Y[m,n] != 0)
        - sum(xbar[:Pnet][m]*Va[m] + xbar[:Qnet][m]*log(Vm[m]) for m in 1:nbus))

    #
    # Line failure constraint
    #
    # @NLconstraint(em, baseMVA*(Vm[i]^2 + Vm[j]^2 - (2*Vm[i]*Vm[j]*cos(Va[i]-Va[j])) - flowmax) == 0)
    @NLconstraint(em, baseMVA*( flowmax - (Vm[i]^2 + Vm[j]^2 - (2*Vm[i]*Vm[j]*cos(Va[i]-Va[j]))) ) == 0)


    #
    # initial value
    #
    setvalue(getindex(em, :Vm), VMbar)
    setvalue(getindex(em, :Va), VAbar)


    #
    # solve model
    #
    status = solve(em)
    if status != :Optimal
        if options[:print_level] >= 1
            println("warning: unable to compute EXACT exit rate for line ", l, " = (", i, ",", j, "). ipopt returned ", status)
        end
        return nothing
    end


    #
    # Construct solution in reduced space
    #
    VMstar = getvalue(getindex(em,:Vm))
    VAstar = getvalue(getindex(em,:Va))
    VMstar_r = copy(VMstar)
    VAstar_r = copy(VAstar)
    for i in opfmodeldata[:nonLoadBuses][end:-1:1]
        splice!(VMstar_r, i)
    end
    splice!(VAstar_r, opfmodeldata[:bus_ref])

    ##
    ##
    exit_point = Dict()
    exit_point[:Vm]  = VMstar
    exit_point[:Va]  = VAstar
    exit_point[:Vmr] = VMstar_r
    exit_point[:Var] = VAstar_r
    exit_point[:Pnet]= xbar[:Pnet]
    exit_point[:Qnet]= xbar[:Qnet]
    ##

    H_xstar      =          H([VMstar_r; VAstar_r]; solution=exit_point, opfmodeldata=opfmodeldata)
    grad_H_xstar =  ∇H_direct([VMstar  ; VAstar  ]; solution=exit_point, opfmodeldata=opfmodeldata)
    hess_H_xstar = ∇2H_direct([VMstar  ; VAstar  ]; solution=exit_point, opfmodeldata=opfmodeldata)
    grad_h_xstar =  ∇h_direct([VMstar  ; VAstar  ]; solution=exit_point, opfmodeldata=opfmodeldata, i=i, j=j, flowmax=flowmax)
    hess_h_xstar = ∇2h_direct([VMstar  ; VAstar  ]; solution=exit_point, opfmodeldata=opfmodeldata, i=i, j=j, flowmax=flowmax)


    # calculate rate
    Kstar      = norm(grad_H_xstar) / norm(grad_h_xstar)
    tempM      = hess_H_xstar - (Kstar * hess_h_xstar)
    z1         = exp(logabsdet(xbar[:hess_H])[1] - logabsdet(tempM)[1])
    z2         = grad_H_xstar' * (tempM\grad_H_xstar)
    z3         = sqrt(abs(z1/z2))
    kappa      = Kstar^3/z3
    prefactor  = z3 * sum(S.*grad_H_xstar.^2) * (options[:damping]/sqrt(2*pi*options[:temperature]) )
    energydiff = H_xstar - xbar[:H_xbar]
    if options[:high_temp_adj]
        prefactor += options[:temperature]/energydiff
    end
    expterm    = max(exp(-energydiff/options[:temperature]), eps(0.0))
    caputil    = 100.0 * options[:constr_limit_scale] * sqrt(abs(VMbar[i]^2 + VMbar[j]^2 - 2*VMbar[i]*VMbar[j]*cos(VAbar[i]-VAbar[j])))/sqrt(flowmax)
    if isnan(kappa)
        (options[:print_level] >= 1) && println("warning: unable to compute EXACT exit rate for line ", l, " = (", i, ",", j, "). kappa = ", kappa, " (NaN)")
        return nothing
    end

    ##
    ##
    exit_point[:K]          = Kstar
    exit_point[:prefactor]  = prefactor
    exit_point[:energydiff] = energydiff
    exit_point[:expterm]    = expterm
    exit_point[:kappa]      = kappa
    exit_point[:caputil]    = caputil
    ##
    ##

    return exit_point
end

function compute_exitrate_kkt(l::Int, xbar::Dict, opfmodeldata::Dict, options::Dict)
    @assert options[:lossless]
    @assert options[:current_rating]

    buses        = opfmodeldata[:buses]
    line         = opfmodeldata[:lines][l]
    nonLoadBuses = opfmodeldata[:nonLoadBuses]
    bus_ref      = opfmodeldata[:bus_ref]
    Y            = opfmodeldata[:Y]
    baseMVA      = opfmodeldata[:baseMVA]
    S            = opfmodeldata[:S]
    # flowmax      = (options[:constr_limit_scale]*line.rateA/(abs(1.0/(line.x*im))*baseMVA))^2
    i            = opfmodeldata[:BusIdx][line.from]
    j            = opfmodeldata[:BusIdx][line.to]
    flowmax      = options[:constr_limit_scale]*(line.rateA^2 / (opfmodeldata[:baseMVA]^2 * abs2(1.0 / line.x)))
    nbus         = length(buses)
    nrow         = 2nbus - length(nonLoadBuses) - 1
    VMbar        = xbar[:Vm]
    VAbar        = xbar[:Va]



    ## check if exit rate formula is applicable
    if ((i in nonLoadBuses && j in nonLoadBuses) || line.rateA==0 || line.rateA>=1.0e10)
        return nothing
    end


    # --------------------
    # solve KKT system
    # --------------------
    em = Model(solver = IpoptSolver(print_level=options[:print_level],max_iter=1000))

    #
    # Variables
    #
    @variable(em, Vm[1:nbus])
    @variable(em, Va[1:nbus])
    @variable(em, K >= 0)
    setlowerbound(Va[bus_ref], VAbar[bus_ref])
    setupperbound(Va[bus_ref], VAbar[bus_ref])
    setlowerbound(Vm[bus_ref], VMbar[bus_ref])
    setupperbound(Vm[bus_ref], VMbar[bus_ref])
    for b in 1:nbus
        if b in nonLoadBuses
            setlowerbound(Vm[b], VMbar[b])
            setupperbound(Vm[b], VMbar[b])
        end
    end


    #
    # KKT - Primal feasibility
    #
    # @NLconstraint(em, baseMVA*(Vm[i]^2 + Vm[j]^2 - (2*Vm[i]*Vm[j]*cos(Va[i]-Va[j])) - flowmax) == 0)
    @NLconstraint(em, baseMVA*(flowmax - (Vm[i]^2 + Vm[j]^2 - (2*Vm[i]*Vm[j]*cos(Va[i]-Va[j])))) == 0)


    #
    # Set objective function to be 2nd-order Taylor expansion
    #
    lin_vm = opfmodeldata[:linindex_VMr]
    lin_va = opfmodeldata[:linindex_VAr]
    idx_vm = [n for n in 1:nbus if n ∉ nonLoadBuses]
    idx_va = [n for n in 1:nbus if n != bus_ref]
    @objective(em, Min,
                    sum( xbar[:hess_H][lin_vm[a],lin_vm[b]]*(Vm[a] - VMbar[a])*(Vm[b] - VMbar[b])
                            for a in idx_vm for b in idx_vm
                                if xbar[:hess_H][lin_vm[a],lin_vm[b]] != 0) +
                    sum( xbar[:hess_H][lin_va[a],lin_va[b]]*(Va[a] - VAbar[a])*(Va[b] - VAbar[b])
                            for a in idx_va for b in idx_va
                                if xbar[:hess_H][lin_va[a],lin_va[b]] != 0) +
                    sum(2xbar[:hess_H][lin_vm[a],lin_va[b]]*(Vm[a] - VMbar[a])*(Va[b] - VAbar[b])
                            for a in idx_vm for b in idx_va
                                if xbar[:hess_H][lin_vm[a],lin_va[b]] != 0)
    )

    #
    # initial value
    #
    setvalue(getindex(em, :Vm), VMbar)
    setvalue(getindex(em, :Va), VAbar)


    #
    # solve model
    #
    status = solve(em)
    if status != :Optimal
        if options[:print_level] >= 1
            println("warning: unable to compute KKT exit rate for line ", l, " = (", i, ",", j, "). ipopt returned ", status)
        end
        return nothing
    end


    # Set initial value
    Vmopt = getvalue(getindex(em,:Vm))
    Vaopt = getvalue(getindex(em,:Va))
    setvalue(getindex(em, :Vm), Vmopt)
    setvalue(getindex(em, :Va), Vaopt)



    # reset objective
    @objective(em, Min, 0)


    #
    # KKT - stationarity
    #
    add_kkt_stationarity_constraint!(em, opfmodeldata, i, j, VMbar, VAbar, xbar[:Qs], Vm, Va, K)


    # solve model
    status = solve(em)
    if status != :Optimal
        if options[:print_level] >= 1
            println("warning: unable to compute KKT exit rate for line ", l, " = (", i, ",", j, "). ipopt returned ", status)
        end
        return nothing
    end

    #
    # calculate rate now
    #
    VMstar = getvalue(getindex(em,:Vm))
    VAstar = getvalue(getindex(em,:Va))
    Kstar  = getvalue(getindex(em,:K))

    ## Construct solution in reduced space
    VMstar_r = copy(VMstar)
    VAstar_r = copy(VAstar)
    xstar_minus_xbar = [VMstar; VAstar] - [VMbar; VAbar]
    splice!(xstar_minus_xbar, opfmodeldata[:nbus] + opfmodeldata[:bus_ref])
    for i in opfmodeldata[:nonLoadBuses][end:-1:1]
        splice!(VMstar_r, i)
        splice!(xstar_minus_xbar, i)
    end
    splice!(VAstar_r, opfmodeldata[:bus_ref])


    ##
    ##
    exit_point = Dict()
    exit_point[:Vm]  = VMstar
    exit_point[:Va]  = VAstar
    exit_point[:Vmr] = VMstar_r
    exit_point[:Var] = VAstar_r
    exit_point[:K]   = Kstar
    ##


    #
    # Calculate constraint gradient
    #
    grad_h = ∇h_direct([VMstar; VAstar]; solution=exit_point, opfmodeldata=opfmodeldata, i=i, j=j, flowmax=flowmax)
    xstar_minus_xbar_times_grad_h = xstar_minus_xbar'*grad_h

    #
    # Get LDL' factorization of constraint Hessian
    #
    a = opfmodeldata[:linindex_VMr][i]
    b = opfmodeldata[:linindex_VMr][j]
    c = opfmodeldata[:linindex_VAr][i]
    d = opfmodeldata[:linindex_VAr][j]
    if a > 0 && b > 0         ## load-load
        L = zeros(nrow, 3)
        D = Diagonal{Float64}(I, 3)
        get_LDL_factors_load_load!(L, D, VMstar[i], VAstar[i], VMstar[j], VAstar[j], a, c, b, d)
    else
        L = zeros(nrow, 2)
        D = Diagonal{Float64}(I, 2)
        if a > 0 && d > 0     ## load-gen
            get_LDL_factors_load_gen!(L, D, VMstar[i], VAstar[i], VMstar[j], VAstar[j], a, c, b, d)
        elseif a > 0          ## load-slack
            get_LDL_factors_load_slk!(L, D, VMstar[i], VAstar[i], VMstar[j], VAstar[j], a, c, b, d)
        elseif b > 0 && c > 0 ## gen-load
            get_LDL_factors_load_gen!(L, D, VMstar[j], VAstar[j], VMstar[i], VAstar[i], b, d, a, c)
        elseif b > 0          ## slack-load
            get_LDL_factors_load_slk!(L, D, VMstar[j], VAstar[j], VMstar[i], VAstar[i], b, d, a, c)
        else                  ## nonload-nonload ==> should not happen
            @assert false
        end
    end



    ##
    ## Just for sanity check -- comment out after debugging
    ##
    #=
    xstar_minus_xbar_times_grad_h_check =
        ((VMstar[i] - VMbar[i])*( VMstar[i]-(VMstar[j]*cos(VAstar[i] - VAstar[j])))) +
        ((VMstar[j] - VMbar[j])*( VMstar[j]-(VMstar[i]*cos(VAstar[i] - VAstar[j])))) +
        ((VAstar[i] - VAbar[i])*( VMstar[i]* VMstar[j]*sin(VAstar[i] - VAstar[j]))) +
        ((VAstar[j] - VAbar[j])*(-VMstar[i]* VMstar[j]*sin(VAstar[i] - VAstar[j])))
    norm_squared_grad_h_check =
        ((a > 0 ? 1.0 : 0.0)*( VMstar[i]-(VMstar[j]*cos(VAstar[i] - VAstar[j])))^2) +
        ((b > 0 ? 1.0 : 0.0)*( VMstar[j]-(VMstar[i]*cos(VAstar[i] - VAstar[j])))^2) +
        ((c > 0 ? 1.0 : 0.0)*( VMstar[i]* VMstar[j]*sin(VAstar[i] - VAstar[j]))^2) +
        ((d > 0 ? 1.0 : 0.0)*(-VMstar[i]* VMstar[j]*sin(VAstar[i] - VAstar[j]))^2)
    norm_squared_grad_h_weighted_S_check =
        ((a > 0 ? S[a] : 0.0)*( VMstar[i]-(VMstar[j]*cos(VAstar[i] - VAstar[j])))^2) +
        ((b > 0 ? S[b] : 0.0)*( VMstar[j]-(VMstar[i]*cos(VAstar[i] - VAstar[j])))^2) +
        ((c > 0 ? S[c] : 0.0)*( VMstar[i]* VMstar[j]*sin(VAstar[i] - VAstar[j]))^2) +
        ((d > 0 ? S[d] : 0.0)*(-VMstar[i]* VMstar[j]*sin(VAstar[i] - VAstar[j]))^2)
    hess_h = ∇2h_direct([VMstar; VAstar]; solution=exit_point, opfmodeldata=opfmodeldata, i=i, j=j, flowmax=flowmax)
    @assert maximum(abs.((L*D*L') - hess_h)) <= 1e-8
    @assert abs(norm_squared_grad_h_check - norm(grad_h)^2) <= 1e-8
    @assert abs(norm_squared_grad_h_weighted_S_check - sum(S.*grad_h.^2)) <= 1e-8
    @assert abs(xstar_minus_xbar_times_grad_h_check - xstar_minus_xbar_times_grad_h) <= 1e-8
    =#
    ##
    ##


    inv_hess_H_times_L       = xbar[:hess_H] \ L
    Proj                     = compute_projection_matrix(nbus, nonLoadBuses, bus_ref)
    Z_initial                = Proj' * inv_hess_H_times_L ## used in main optimization as an initial point
    X                        = -inv(D) + (Kstar*L'*inv_hess_H_times_L)
    L_times_xstar_minus_xbar = L'*xstar_minus_xbar
    kappa                    = det(-D) * det(X) * (## we expect this term to be close to zero O(10^-6 - 10^-2)
                                xstar_minus_xbar_times_grad_h -
                                    (L_times_xstar_minus_xbar' * (X \ L_times_xstar_minus_xbar))
                               )
    if kappa < 0
        (options[:print_level] >= 1) && println("warning: forcing kappa > 0 for line ", l, " = (", i, ",", j, "). kappa = ", kappa)
        kappa = abs(kappa)
    end
    if isnan(kappa)
        (options[:print_level] >= 1) && println("warning: unable to compute KKT exit rate for line ", l, " = (", i, ",", j, "). kappa = ", kappa, " (negative)")
        return nothing
    end
    prefactor  = Kstar^(1.5) * sum(S.*grad_h.^2) * (options[:damping]/sqrt(2*pi*options[:temperature]) ) / sqrt(kappa)
    energydiff = Kstar * xstar_minus_xbar_times_grad_h/2.0
    if options[:high_temp_adj]
        prefactor += options[:temperature]/energydiff
    end
    expterm    = max(exp(-energydiff/options[:temperature]), eps(0.0))
    ## line capacity utilization
    caputil    = 100.0 * options[:constr_limit_scale] * sqrt(abs(VMbar[i]^2 + VMbar[j]^2 - 2*VMbar[i]*VMbar[j]*cos(VAbar[i]-VAbar[j])))/sqrt(flowmax)

    if prefactor < 0 || isnan(prefactor)
        (options[:print_level] >= 1) && println("warning: unable to compute KKT exit rate for line ", l, " = (", i, ",", j, "). prefactor = ", prefactor, " (negative)")
        return nothing
    end


    ##
    ##
    exit_point[:prefactor]  = prefactor
    exit_point[:energydiff] = energydiff
    exit_point[:expterm]    = expterm
    exit_point[:Z]          = Z_initial
    exit_point[:kappa]      = kappa
    exit_point[:caputil]    = caputil
    ##
    ##

    return exit_point
end

## returns a matrix P of size (nreduced x nbus) such that
## for any A of size (nreduced x nreduced)
## P'AP is an expanded matrix of size (nbus x nbus)
## with zeros for entries that are not in nreduced
function compute_projection_matrix(nbus, nonLoadBuses, bus_ref)
    Proj = zeros(2nbus - length(nonLoadBuses) - 1, 2nbus)
    b = 1
    for i in 1:nbus
        if i ∉ nonLoadBuses
            Proj[b,i] = 1
            b += 1
        end
    end
    for i in 1:nbus
        if i != bus_ref
            Proj[b,i+nbus] = 1
            b += 1
        end
    end
    return Proj
end

## NOTE: Assume the line is (i,j) where i is a load and j can be any other bus type
##          V1,    θ1,    V2,    θ2 = Vmstar[i], Vastar[i], Vmstar[j], Vastar[j]
##       idxV1, idxθ1, idxV2, idxθ2 = ^their linear indices in the reduced vector
function get_LDL_factors_load_load!(L, D, V1, θ1, V2, θ2, idxV1, idxθ1, idxV2, idxθ2; model = nothing)
    if model == nothing
        D[1,1] = 1
        D[2,2] = 1
        D[3,3] = -1
        L[idxV1,1] = 1
        L[idxV2,1] = -cos(θ1 - θ2)
        L[idxV2,2] = sin(θ1 - θ2)
        L[idxθ1,1] = V2*sin(θ1 - θ2)
        L[idxθ1,2] = V1+ (V2*cos(θ1 - θ2))
        L[idxθ1,3] = sqrt(V1^2 + V2^2 + (V1*V2*cos(θ1 - θ2)))
        L[idxθ2,:] = -L[idxθ1,:]
    else
        D[1,1] = @NLexpression(model, 1)
        D[2,2] = @NLexpression(model, 1)
        D[3,3] = @NLexpression(model, -1)
        L[idxV1,1] = @NLexpression(model, 1)
        L[idxV2,1] = @NLexpression(model, -cos(θ1 - θ2))
        L[idxV2,2] = @NLexpression(model, sin(θ1 - θ2))
        for idx in [idxθ1, idxθ2]
            coef = (idx == idxθ1) ? 1 : -1
            L[idx,1] = @NLexpression(model, coef*V2*sin(θ1 - θ2))
            L[idx,2] = @NLexpression(model, coef*(V1+ (V2*cos(θ1 - θ2))))
            L[idx,3] = @NLexpression(model, coef*sqrt(V1^2 + V2^2 + (V1*V2*cos(θ1 - θ2))))
        end
    end
end
function get_LDL_factors_load_slk!(L, D, V1, θ1, V2, θ2, idxV1, idxθ1, idxV2, idxθ2; model = nothing)
    if model == nothing
        D[1,1] = 1
        D[2,2] = (V1*V2*cos(θ1 - θ2)) - (V2*sin(θ1 - θ2))^2
        L[idxV1,1] = 1
        L[idxθ1,1] = V2*sin(θ1 - θ2)
        L[idxθ1,2] = 1
    else
        D[1,1] = @NLexpression(model, 1)
        D[2,2] = @NLexpression(model, (V1*V2*cos(θ1 - θ2)) - (V2*sin(θ1 - θ2))^2)
        L[idxV1,1] = @NLexpression(model, 1)
        L[idxθ1,1] = @NLexpression(model, V2*sin(θ1 - θ2))
        L[idxθ1,2] = @NLexpression(model, 1)
    end
end
function get_LDL_factors_load_gen!(L, D, V1, θ1, V2, θ2, idxV1, idxθ1, idxV2, idxθ2; model = nothing)
    get_LDL_factors_load_slk!(L, D, V1, θ1, V2, θ2, idxV1, idxθ1, idxV2, idxθ2; model = model)
    if model == nothing
        L[idxθ2,1] = -L[idxθ1,1]
        L[idxθ2,2] = -L[idxθ1,2]
    else
        L[idxθ2,1] = @NLexpression(model, -L[idxθ1,1])
        L[idxθ2,2] = @NLexpression(model, -L[idxθ1,2])
    end
end
