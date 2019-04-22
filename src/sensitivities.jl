"""
## `get_Gamma`: get sensitivities wrt fixed/uncertain dispatch `d` space
### arguments:
    - `z::AbstractArray`: full data `z := (x, y) = (Pg; Qg; Vm; Va; Pd; Qd)` from an `OPFModel`
### returns:
    - `Γ::AbstractArray`: dx/dd
"""
function get_Gamma(z::AbstractArray, data::Dict, Gamma_type::Symbol=:ew, epsilon::Float64=1e-5)
    @assert(haskey(data, :m_idx))
    m_idx         = data[:m_idx]
    Y             = data[:Y]
    Pnet          = data[:Pnet]
    Qnet          = data[:Qnet]
    BusIdx        = data[:BusIdx]
    BusGenerators = data[:BusGenerators]

    if Gamma_type ∈ [:vec, :num, :ew]
        if Gamma_type == :num
            J = jac_z(z, data, Gamma_type)
        else
            J,_,_ = jac_z(z, data, Gamma_type)
        end
        dF_dx = J[m_idx[:F], m_idx[:x]]
        dF_dy = J[m_idx[:F], m_idx[:y]]
        dF_dd = J[m_idx[:F], m_idx[:d]]
        Γ = dF_dx \ -dF_dd
        return Γ
    elseif Gamma_type == :fd
        return get_Gamma_fd(z, Y, Pnet, Qnet, m_idx, BusIdx, BusGenerators, epsilon)
    end
end

"""
## `get_Gamma_fd`: get finite difference approximation of `Gamma`
"""
function get_Gamma_fd(z::AbstractArray, Y::AbstractArray,
                      Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict,
                      BusIdx::Dict, BusGenerators,
                      epsilon::Float64, verb::Int64=1)
    ## setup
    x = deepcopy(z[m_idx[:x]])
    nz = length(z)
    nx = length(x)
    ny = length(m_idx[:y])
    nd = length(m_idx[:d])
    nbus = length(Pnet); @assert(nbus == length(Qnet)); @assert(nbus == size(Y,1))
    #### output
    Γ_fd = fill(NaN, nx, nd)
    xs = fill(NaN, nx, nd)

    ## compute
    for i in 1:nd
        ## overwrite
        x_ = deepcopy(x)
        z_ = deepcopy(z)
        Pnet_ = deepcopy(Pnet)
        Qnet_ = deepcopy(Qnet)

        ## perturb
        if i <= nbus
            Pnet_[i] -= epsilon
        elseif i > nbus
            Qnet_[mod1(i, nbus)] -= epsilon
        end
        if verb >= 1
            println("---------- bus index $(mod1(i, nbus)) ----------")
            println("Pnet diff: ", Pnet-Pnet_)
            println("Qnet diff: ", Qnet-Qnet_)
        end

        ## solve
        PFE_wrap!(F, x) = PFE!(F, x, z_, Y, Pnet_, Qnet_, m_idx)
        PFE_J_wrap!(J, x) = PFE_J!(J, x, z_, Y, Pnet_, Qnet_, m_idx, BusIdx, BusGenerators)
        sol = nlsolve(PFE_wrap!, PFE_J_wrap!, x_, iterations=5_000, ftol=epsilon/100)
        @info sol.f_converged == true
        xs[:,i] = sol.zero
        if verb >= 1
            println("     NaNs:  ", sum(isnan.((xs[:,i] - x) ./ epsilon)))
            println()
        end
        Γ_fd[:,i] = (xs[:,i] - x) / epsilon
    end
    return Γ_fd
end
