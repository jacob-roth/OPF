## -----------------------------------------------------------------------------
## complex, vectorized
## -----------------------------------------------------------------------------
function PF(opfmodel_z::AbstractArray, Y::AbstractArray, z_idx::Dict)
    ## setup
    Vm = opfmodel_z[z_idx[:Vm]]
    Va = opfmodel_z[z_idx[:Va]]
    ngen = length(z_idx[:Pg]); @assert(ngen == length(z_idx[:Qg]))
    nbus = length(z_idx[:Vm]); @assert(nbus == length(z_idx[:Va]))
    @assert(length(opfmodel_z) == 2ngen + 2nbus + 2nbus)

    # branch admitances
    Vtilde = Vm .* exp.(im .* Va)
    Itilde = Y*Vtilde
    Stilde = Vtilde .* conj.(Itilde)
    return Vtilde, Itilde, Stilde
end
function PF(y_RGL::AbstractArray, data::Dict)
    ## setup
    Y = data[:Y]
    y_RGL_idx = data[:y_RGL_idx]
    z_idx = data[:z_idx]
    opfmodel_z = deepcopy(data[:opfmodel_z])
    ngen = length(z_idx[:Pg]); @assert(ngen == length(z_idx[:Qg]))
    nbus = length(z_idx[:Vm]); @assert(nbus == length(z_idx[:Va]))
    @assert(length(opfmodel_z) == 2ngen + 2nbus + 2nbus)

    opfmodel_z[y_RGL_idx] = y_RGL
    return PF(opfmodel_z, Y, z_idx)
end
function PFE_RGL!(F::AbstractArray, y_RGL::AbstractArray, data::Dict)
    ## setup
    @assert(length(F) == length(y_RGL))
    f_RGL_idx = data[:f_RGL_idx]

    Vtilde, Itilde, Stilde = PF(deepcopy(y_RGL), data)
    F_P = real.(Stilde) - data[:Pnet]
    F_Q = imag.(Stilde) - data[:Qnet]
    FF = [F_P; F_Q]
    F .= FF[f_RGL_idx]
    nothing
end

## -----------------------------------------------------------------------------
## real, entrywise
## -----------------------------------------------------------------------------
function PF_real(opfmodel_z::AbstractArray, Y::AbstractArray, z_idx::Dict)
    ## setup
    Vm = opfmodel_z[z_idx[:Vm]]
    Va = opfmodel_z[z_idx[:Va]]
    ngen = length(z_idx[:Pg]); @assert(ngen == length(z_idx[:Qg]))
    nbus = length(z_idx[:Vm]); @assert(nbus == length(z_idx[:Va]))
    @assert(length(opfmodel_z) == 2ngen + 2nbus + 2nbus)
    G = real.(Y)
    B = imag.(Y)
    F = fill(NaN, 2nbus)

    # branch admitances
    for i = 1:nbus
        h = mod1(i, nbus)
        IDX = Y[h,:].nzind
        P = Vm[h] * sum(Vm[kk] * ( G[h,kk] * cos(Va[h]-Va[kk]) + B[h,kk] * sin(Va[h]-Va[kk]) ) for kk in IDX)
        Q = Vm[h] * sum(Vm[kk] * ( G[h,kk] * sin(Va[h]-Va[kk]) - B[h,kk] * cos(Va[h]-Va[kk]) ) for kk in IDX)
        F[i] = P
        F[i+nbus] = Q
    end
    return F
end
function PF_real(y_RGL::AbstractArray, data::Dict)
    ## setup
    Y = data[:Y]
    y_RGL_idx = data[:y_RGL_idx]
    z_idx = data[:z_idx]
    opfmodel_z = deepcopy(data[:opfmodel_z])
    ngen = length(z_idx[:Pg]); @assert(ngen == length(z_idx[:Qg]))
    nbus = length(z_idx[:Vm]); @assert(nbus == length(z_idx[:Va]))
    @assert(length(opfmodel_z) == 2ngen + 2nbus + 2nbus)

    opfmodel_z[y_RGL_idx] = y_RGL
    return PF_real(opfmodel_z, Y, z_idx)
end
function PFE_RGL_real!(F::AbstractArray, y_RGL::AbstractArray, data::Dict)
    ## setup
    @assert(length(F) == length(y_RGL))
    f_RGL_idx = data[:f_RGL_idx]
    z_idx = data[:z_idx]
    nbus = length(z_idx[:Vm]); @assert(nbus == length(z_idx[:Va]))

    FF = PF_real(deepcopy(y_RGL), data)
    F_P = FF[1:nbus] - data[:Pnet]
    F_Q = FF[(nbus+1):(2nbus)] - data[:Qnet]
    FF = [F_P; F_Q]
    F .= FF[f_RGL_idx]
    nothing
end