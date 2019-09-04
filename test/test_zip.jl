## -----------------------------------------------------------------------------
## (α, β, γ) = (1, 0, 0)
## μ_P = 0, μ_Q = 0
## -----------------------------------------------------------------------------
m      = OPF.acopf_model(opfdata)
m      = OPF.acopf_solve(m, opfdata)
m_eval = setup(m.m);               ## deterministic model evaluator
m_zbar = deepcopy(m_eval.last_x);  ## deterministic model equilibrium z̄
OPF.acopf_outputAll(m, opfdata)

z      = acopf_zip_model(opfdata)
z      = acopf_solve(z, opfdata)
z_eval = setup(z.m);               ## deterministic model evaluator
z_zbar = deepcopy(z_eval.last_x);  ## deterministic model equilibrium z̄
acopf_outputAll(z, opfdata)

@testset "zip base equality" begin
    @test(norm(m_zbar - z_zbar) < tol)
end # testset

## -----------------------------------------------------------------------------
## range over α, β, γ
## μ_P = 0, μ_Q = 0
## -----------------------------------------------------------------------------
ϵ = [1e-5; 1e-4; 1e-3; 1e-2; 1e-1; 5e-1; 1.0-1e-1; 1.0-1e-2; 1.0-1e-3; 1.0-1e-4; 1.0-1e-5]
α = 1.0 .- ϵ
β = (1.0 .- α) ./ 2.0
γ = (1.0 .- α) ./ 2.0
@assert( α .+ β .+ γ == ones(length(ϵ)))
obj = zeros(length(ϵ))
opt = Array{Symbol}(undef, length(ϵ))
for i in eachindex(α)
    D = DefaultOptions()
    D[:zip][:alpha] = α[i]
    D[:zip][:beta]  = β[i]
    D[:zip][:gamma] = γ[i]
    z      = acopf_zip_model(opfdata, D)
    z      = acopf_solve(z, opfdata)
    z_eval = setup(z.m);               ## deterministic model evaluator
    z_zbar = deepcopy(z_eval.last_x);  ## deterministic model equilibrium z̄
    acopf_outputAll(z, opfdata, D)
    obj[i] = getobjectivevalue(z.m)
    opt[i] = z.status
end
# semilogx(ϵ, obj)

@testset "zip optimality" begin
for i in eachindex(ϵ)
    @test(opt[i] == :Optimal)
end
end # testset