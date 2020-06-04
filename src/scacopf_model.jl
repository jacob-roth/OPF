function scacopf_model(opfdata::OPFData, options::Dict=DefaultOptions(), adjustments::Dict=DefaultAdjustments(), contingencies::Dict=Dict())

    ## setup
    opfmodeldata = get_opfmodeldata(opfdata, options, adjustments)
    nbus = length(opfmodeldata[:buses])
    nline = length(opfmodeldata[:lines])
    ngen = length(opfmodeldata[:generators])
    R = opfdata.bus_ref
    G = filter(x -> x ∉ R, findall(.!isempty.(opfmodeldata[:BusGenerators])))
    L = findall(isempty.(opfmodeldata[:BusGenerators]))
    not_R = deleteat!(collect(1:nbus), R)

    ## base model
    M = acopf_model(opfdata, options, adjustments)
    m = M.m

    ## objective expression of each contingency
    obj_cs = Array{Union{Float64,JuMP.NonlinearExpression},1}(undef, length(keys(contingencies)))

    ## add contingency security constraints
    if !isempty(contingencies)
        for c_id in keys(contingencies)
            #
            # contingency (line removal)
            #
            c = contingencies[c_id];
            @assert(c.c_type == :line)
            rl = remove_line!(opfdata, c_id);
            c_opfmodeldata = get_opfmodeldata(opfdata, options, adjustments)

            #
            # secondary variables
            #
            Pg_ = @variable(m, [i=1:ngen], basename="Pg_$c_id")
            Qg_ = @variable(m, [i=1:ngen], basename="Qg_$c_id")
            Pg  = m[:Pg]
            Vm  = m[:Vm]
            for i in 1:ngen
                ## bound
                @constraint(m, opfmodeldata[:generators].Pmin[i] <= Pg_[i] <= opfmodeldata[:generators].Pmax[i])
                ## ramp
                @constraint(m, Pg[i] - Pg_[i] <= options[:ramp_pct] * opfmodeldata[:generators].Pmax[i])
                @constraint(m, Pg_[i] - Pg[i] <= options[:ramp_pct] * opfmodeldata[:generators].Pmax[i])
            end
            for i in 1:ngen
                ## bound
                @constraint(m, opfmodeldata[:generators].Qmin[i] <= Qg_[i] <= opfmodeldata[:generators].Qmax[i])
            end
            Vm_ = @variable(m, [i=1:nbus], basename="Vm_$c_id")
            for i in 1:nbus
                @constraint(m, (1 - options[:ctg_Vm_adj])*opfmodeldata[:buses].Vmin[i] <= Vm_[i] <= (1 + options[:ctg_Vm_adj])*opfmodeldata[:buses].Vmax[i])
            end
            Va_ = @variable(m, [i=1:nbus], basename="Va_$c_id")
            for i in 1:nbus
                @constraint(m, -pi <= Va_[i] <= pi)
            end
            JuMP.registerobject(m, Symbol("Pg_$(c_id)"), Pg_, "Pg_$(c_id)")
            JuMP.registerobject(m, Symbol("Qg_$(c_id)"), Qg_, "Qg_$(c_id)")
            JuMP.registerobject(m, Symbol("Vm_$(c_id)"), Vm_, "Vm_$(c_id)")
            JuMP.registerobject(m, Symbol("Va_$(c_id)"), Va_, "Va_$(c_id)")

            # ## composite containers
            # Pg = Array{Variable,1}(undef, ngen)
            # Qg = Array{Variable,1}(undef, ngen)
            # Vm = Array{Variable,1}(undef, nbus)
            # Va = Array{Variable,1}(undef, nbus)
            # for gi in zip([R], first.(opfdata.BusGenerators[R]))
            #     g = gi[1]
            #     i = gi[2]
            #     Pg[i] = m[Symbol("Pg_$(c_id)")][g]
            # end
            # for gi in zip(G, first.(opfdata.BusGenerators[G]))
            #     g = gi[1]
            #     i = gi[2]
            #     Pg[i] = m[Symbol("Pg")][i] # since in acopf Pg is indexed by generator not bus
            # end
            # for gi in zip([G; R], [first.(opfdata.BusGenerators[G]); first.(opfdata.BusGenerators[R])])
            #     g = gi[1]
            #     i = gi[2]
            #     Qg[i] = m[Symbol("Qg_$(c_id)")][g]
            # end
            # for i in [G; R]
            #     Vm[i] = m[:Vm][i]
            # end
            # for i in L
            #     Vm[i] = m[Symbol("Vm_$(c_id)")][i]
            # end
            # for i in [R]
            #     Va[i] = m[:Va][i]
            # end
            # for i in [L; G]
            #     Va[i] = m[Symbol("Va_$(c_id)")][i]
            # end
            # JuMP.registerobject(m, Symbol("Pg_$(c_id)_container"), Pg, "Pg_$(c_id)_container")
            # JuMP.registerobject(m, Symbol("Qg_$(c_id)_container"), Qg, "Qg_$(c_id)_container")
            # JuMP.registerobject(m, Symbol("Vm_$(c_id)_container"), Vm, "Vm_$(c_id)_container")
            # JuMP.registerobject(m, Symbol("Va_$(c_id)_container"), Va, "Va_$(c_id)_container")

            #
            # secondary objective
            #
            if options[:ctg_feas] == true
                obj_c = 0.0
            else
                obj_c = @NLexpression(m, sum(opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n-2] * (opfmodeldata[:baseMVA] * Pg_[i]) ^ 2 + opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n-1] * (opfmodeldata[:baseMVA] * Pg_[i]) +
                        opfmodeldata[:generators][i].coeff[opfmodeldata[:generators][i].n] for i=1:ngen))
            end
            c_idx = first(findall(c_id .== collect(keys(contingencies))))
            obj_cs[c_idx] = obj_c

            #
            # power flow balance
            #
            islanding_buses = get_islanding_buses(c_opfmodeldata, options)
            for b in 1:nbus
                if b ∉ islanding_buses
                    ## real part
                    add_p_constraint!(m, c_opfmodeldata, b, c_id)
                    ## imaginary part
                    add_q_constraint!(m, c_opfmodeldata, b, c_id)
                end
            end

            #
            # branch/lines flow limits
            #
            for l in filter(x -> x != c_id, keys(contingencies))
                if options[:current_rating]
                    ## current
                    add_line_current_constraint!(m, c_opfmodeldata, options, l, c_id)
                else
                    ## apparent power (to & from)
                    add_line_power_constraint!(m, c_opfmodeldata, l, c_id)
                end
            end

            #
            # contingency (line reinstatement)
            #
            reinstate_line!(opfdata, rl.id, rl);
        end
    end

    #
    # full objective function
    #
    objective_full = m[:obj]
    if options[:sc_total_obj] == true
        for c_id in keys(contingencies)
            c_idx = first(findall(c_id .== collect(keys(contingencies))))
            objective_full = @NLexpression(m, objective_full + obj_cs[c_idx])
        end
    end
    @NLobjective(m, Min, objective_full)

    return OPFModel(m, :InitData, :SC, Dict())
end