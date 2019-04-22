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
function P_i(Vm::AbstractArray, Va::AbstractArray, i::Int64, Y::SparseMatrixCSC{Complex{Float64},Int64})
    IDX = Y[i,:].nzind
    P = Vm[i] * sum(Vm[kk] * ( real(Y[i,kk]) * cos(Va[i]-Va[kk]) + imag(Y[i,kk]) * sin(Va[i]-Va[kk]) ) for kk in IDX)
    return P
end
function Q_i(Vm::AbstractArray, Va::AbstractArray, i::Int64, Y::SparseMatrixCSC{Complex{Float64},Int64})
    IDX = Y[i,:].nzind
    Q = Vm[i] * sum(Vm[kk] * ( real(Y[i,kk]) * sin(Va[i]-Va[kk]) - imag(Y[i,kk]) * cos(Va[i]-Va[kk]) ) for kk in IDX)
    return Q
end
function P_i(Vm::AbstractArray, Va::AbstractArray, i::Int64, Y::Array{Complex{Float64},2})
    P = Vm[i] * sum(Vm[kk] * ( real(Y[i,kk]) * cos(Va[i]-Va[kk]) + imag(Y[i,kk]) * sin(Va[i]-Va[kk]) ) for kk in 1:size(Y,1))
    return P
end
function Q_i(Vm::AbstractArray, Va::AbstractArray, i::Int64, Y::Array{Complex{Float64},2})
    Q = Vm[i] * sum(Vm[kk] * ( real(Y[i,kk]) * sin(Va[i]-Va[kk]) - imag(Y[i,kk]) * cos(Va[i]-Va[kk]) ) for kk in 1:size(Y,1))
    return Q
end
function PF(Vm::AbstractArray, Va::AbstractArray, Y::AbstractArray)
    @assert(length(Vm) == length(Va)); nbus = length(Vm)
    FP = fill(NaN, nbus); FQ = fill(NaN, nbus)
    for i = 1:nbus
        FP[i] = P_i(Vm, Va, i, Y)
        FQ[i] = Q_i(Vm, Va, i, Y)
    end
    return FP, FQ
end
function PFE(Vm::AbstractArray, Va::AbstractArray, Y::AbstractArray,
             Pnet::AbstractArray, Qnet::AbstractArray)
    FP, FQ = PF(Vm, Va, Y)
    rP = FP - Pnet
    rQ = FQ - Qnet
    return rP, rQ
end
function PFE(x::AbstractArray, z::AbstractArray, Y::AbstractArray,
             Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict)
    z[m_idx[:x]] .= x
    Vm = z[m_idx[:Vm]]
    Va = z[m_idx[:Va]]
    return PFE(Vm, Va, Y, Pnet, Qnet)
end
function PFE!(F::AbstractArray, x::AbstractArray,
              z::AbstractArray, Y::AbstractArray,
              Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict)
    @assert(length(F) == length(x))
    rP, rQ = PFE(x, z, Y, Pnet, Qnet, m_idx)
    F .= [rP; rQ][m_idx[:F]]
end
function PFE_J!(J::AbstractArray, x::AbstractArray,
                z::AbstractArray, Y::AbstractArray,
                Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict,
                BusIdx::Dict, BusGenerators)
    _,_, dF = jac_z_alg(z, Y, BusIdx, BusGenerators, m_idx, false)
    J .= dF[:dF_dx]
end


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