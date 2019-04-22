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
function PFE(z::AbstractArray, Y::AbstractArray,
             Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict)
    Vm = z[m_idx[:Vm]]
    Va = z[m_idx[:Va]]
    return PFE(Vm, Va, Y, Pnet, Qnet)
end
function PFE(x::AbstractArray, z::AbstractArray, Y::AbstractArray,
             Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict)
    z[m_idx[:x]] .= x
    Vm = z[m_idx[:Vm]]
    Va = z[m_idx[:Va]]
    return PFE(Vm, Va, Y, Pnet, Qnet)
end
function PFE!(F::AbstractArray, z::AbstractArray, Y::AbstractArray,
              Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict)
    @assert(length(F) == length(m_idx[:x]))
    rP, rQ = PFE(z, Y, Pnet, Qnet, m_idx)
    F .= [rP; rQ][m_idx[:F]]
end
function PFE!(F::AbstractArray, x::AbstractArray,
              z::AbstractArray, Y::AbstractArray,
              Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict)
    @assert(length(F) == length(x))
    z[m_idx[:x]] .= x
    rP, rQ = PFE(x, z, Y, Pnet, Qnet, m_idx)
    F .= [rP; rQ][m_idx[:F]]
end
function PFE_J!(J::AbstractArray, x::AbstractArray,
                z::AbstractArray, Y::AbstractArray,
                Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict,
                BusIdx::Dict, BusGenerators)
    _,_, dF = jac_z_alg_ew(z, Y, BusIdx, BusGenerators, m_idx)
    J .= dF[:dF_dx]
end

## -----------------------------------------------------------------------------
## complex, vectorized
## -----------------------------------------------------------------------------
function PF_vec(Vm::AbstractArray, Va::AbstractArray, Y::AbstractArray)
    Vtilde = Vm .* exp.(im .* Va)
    Itilde = Y*Vtilde
    Stilde = Vtilde .* conj.(Itilde)
    return Vtilde, Itilde, Stilde
end
function PFE_vec(Vm::AbstractArray, Va::AbstractArray, Y::AbstractArray,
                 Pnet::AbstractArray, Qnet::AbstractArray)
    ## setup
    Vtilde, Itilde, Stilde = PF_vec(Vm, Va, Y)
    rP = real.(Stilde) - Pnet
    rQ = imag.(Stilde) - Qnet
    return rP, rQ
end
function PFE_vec(z::AbstractArray, Y::AbstractArray,
                 Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict)
    ## setup
    Vm = z[m_idx[:Vm]]
    Va = z[m_idx[:Va]]
    return PFE_vec(Vm, Va, Y, Pnet, Qnet)
end
function PFE_vec(x::AbstractArray, z::AbstractArray, Y::AbstractArray,
                 Pnet::AbstractArray, Qnet::AbstractArray, m_idx::Dict)
    ## setup
    z[m_idx[:x]] .= x
    Vm = z[m_idx[:Vm]]
    Va = z[m_idx[:Va]]
    return PFE_vec(z, Y, Pnet, Qnet, m_idx)
end