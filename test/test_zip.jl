## -----------------------------------------------------------------------------
## (α, β, γ) = (1, 0, 0)
## μ_P = 0, μ_Q = 0
## -----------------------------------------------------------------------------
m      = OPF.acopf_model(opfdata)
m      = OPF.acopf_solve(m, opfdata)
m_eval = setup(m.m);               ## deterministic model evaluator
m_zbar = deepcopy(m_eval.last_x);  ## deterministic model equilibrium z̄
OPF.acopf_outputAll(m, opfdata)

z      = OPF.acopf_zip_model(opfdata)
z      = OPF.acopf_solve(z, opfdata)
z_eval = setup(z.m);               ## deterministic model evaluator
z_zbar = deepcopy(z_eval.last_x);  ## deterministic model equilibrium z̄
OPF.acopf_outputAll(z, opfdata)

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
obj     = zeros(length(ϵ))
opt     = Array{Symbol}(undef, length(ϵ))
ang_mu  = zeros(length(ϵ))
mag_mu  = zeros(length(ϵ))
Pg_mu   = zeros(length(ϵ))
ang_sig = zeros(length(ϵ))
mag_sig = zeros(length(ϵ))
Pg_sig  = zeros(length(ϵ))

for i in eachindex(α)
    D = OPF.DefaultOptions()
    D[:zip][:alpha] = α[i]
    D[:zip][:beta]  = β[i]
    D[:zip][:gamma] = γ[i]
    z      = OPF.acopf_zip_model(opfdata, D)
    z      = OPF.acopf_solve(z, opfdata)
    z_eval = setup(z.m);               ## deterministic model evaluator
    z_zbar = deepcopy(z_eval.last_x);  ## deterministic model equilibrium z̄
    OPF.acopf_outputAll(z, opfdata, D)
    obj[i] = getobjectivevalue(z.m)
    opt[i] = z.status
    ang_mu[i]  = mean(abs.(getvalue(z.m[:Va])))
    mag_mu[i]  = mean(getvalue(z.m[:Vm]))
    Pg_mu[i]   = mean(getvalue(z.m[:Pg]))
    ang_sig[i] = std(abs.(getvalue(z.m[:Va])))
    mag_sig[i] = std(getvalue(z.m[:Vm]))
    Pg_sig[i]  = std(getvalue(z.m[:Pg]))
end

if plotting
    using PyPlot
    rcParams = PyPlot.matplotlib[:rcParams]
    prop_cycle = rcParams["axes.prop_cycle"]
    colors = prop_cycle[:by_key]()["color"]

    fig = subplots()
    ax  = subplot(111)
    ax[:semilogx](ϵ, obj, label = "objective value")
    ax[:set_xlabel]("ϵ")
    ax[:set_ylabel]("objective value")
    ax[:legend]()

    fig = subplots()
    ax  = subplot(111)
    ax[:loglog](ϵ, ang_mu, label="mean(θ)", color=colors[1], alpha=0.5)
    ax[:loglog](ϵ, ang_sig, label="std(θ)", color=colors[1], alpha=0.5, linestyle=":")
    ax[:loglog](ϵ, mag_mu, label="mean(V)", color=colors[2], alpha=0.5)
    ax[:loglog](ϵ, mag_sig, label="std(V)", color=colors[2], alpha=0.5, linestyle=":")
    ax[:loglog](ϵ, Pg_mu, label="mean(Pg)", color=colors[3], alpha=0.5)
    ax[:loglog](ϵ, Pg_sig, label="std(Pg)", color=colors[3], alpha=0.5, linestyle=":")
    ax[:set_xlabel]("ϵ")
    ax[:set_ylabel]("objective value")
    ax[:legend]()
end

@testset "zip optimality" begin
for i in eachindex(ϵ)
    @test(opt[i] == :Optimal)
end
end # testset