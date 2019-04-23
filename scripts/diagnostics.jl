using PyCall, PyPlot
@pyimport matplotlib.gridspec as gspec
@pyimport matplotlib.colors as matcolors
@pyimport matplotlib.cm as cmx
using LaTeXStrings, Printf

## data
data = Dict()
data[:Sigma_d]           = Matrix(Diagonal(ones(2nbus)))
data[:Va_min]            = -pi * ones(nbus)
data[:Va_max]            =  pi * ones(nbus)

## options
options = Dict()
options[:lossless]       = false
options[:current_rating] = false
options[:epsilon_Vm]     = 0.01
options[:gamma]          = 1.0
options[:relax_Gamma]    = false
options[:Gamma_type]     = :d
options[:xtilde]         = false
options[:print_level]    = 5

function test_input(input_type::Symbol, opfdata::OPFData, data::Dict, options::Dict, name::Symbol, vals::AbstractArray)
    if input_type == :data
        data_ = deepcopy(data)
    elseif input_type == :options
        opts_ = deepcopy(options)
    end
    nexp  = length(vals)
    cost  = Array{Any}(undef, nexp)
    xs    = Array{Any}(undef, nexp)
    Pgs   = Array{Any}(undef, nexp)
    Qgs   = Array{Any}(undef, nexp)
    Vms   = Array{Any}(undef, nexp)
    Vas   = Array{Any}(undef, nexp)
    Pds   = Array{Any}(undef, nexp)
    Qds   = Array{Any}(undef, nexp)
    m_idx = OPF.model_idx(opfdata, options[:xtilde])
    for i in eachindex(vals)
        if input_type == :data
            data_[name] = vals[i]
            println("~~~~ $(data_[name]) = $(vals[i]) ~~~~")
            cm = OPF.cc_acopf_model(opfdata, options, data_)
        elseif input_type == :options
            opts_[name] = vals[i]
            println("~~~~ $(opts_[name]) = $(vals[i]) ~~~~")
            cm = OPF.cc_acopf_model(opfdata, opts_, data)
        end
        cm = OPF.cc_acopf_solve(cm, opfdata, options)
        cost[i] = cm.m.objVal
        v       = OPF.get_values(cm)
        xs[i]   = v[:z][m_idx[:x]]
        Pgs[i]  = v[:Pg]
        Qgs[i]  = v[:Qg]
        Vms[i]  = v[:Vm]
        Vas[i]  = v[:Va]
        Pds[i]  = v[:Pd]
        Qds[i]  = v[:Qd]
    end
    out = Dict()
    out[:cost] = cost
    out[:xs]   = xs
    out[:Pgs]  = Pgs
    out[:Qgs]  = Qgs
    out[:Vms]  = Vms
    out[:Vas]  = Vas
    out[:Pds]  = Pds
    out[:Qds]  = Qds
    return out
end

vals = [0.1; 0.05; 0.01; 0.005; 0.001; 0.0005; 0.0001]
epsilon_Vm = test_input(:options, opfdata, data, options, :epsilon_Vm, vals)
out = epsilon_Vm

σ = [0.01; 0.1; 0.2; 0.5; 1; 1.5; 2; 2.5]
vals = [Matrix(Diagonal(s * ones(2nbus))) for s in σ]
Sigma_d = test_input(:data, opfdata, data, options, :Sigma_d, vals)
out = Sigma_d


fig = figure(figsize=(15,5))
# z_components = [:Pgs; :Qgs; :Vms; :Vas; :Pds; :Qds]
# z_components_names = [raw"$Pg$", raw"$Q_g$", raw"$V_m$", raw"$V_a$", raw"$P_d$", raw"$Q_d$"]
z_components = [:Qgs; :Vms; :Vas]
z_components_names = [raw"$Pg$", raw"$V_m$", raw"$V_a$"]
nzc = length(z_components)
gs = gspec.GridSpec(1, nzc, width_ratios=[1 for i in 1:nzc], height_ratios=[1])
slc = (lo,hi) -> pycall(pybuiltin("slice"), PyObject, lo, hi)
for i = 1:nzc
    println("~~~ component = $(i) ~~~")
    ax = subplot(get(gs, (0,i-1)))
    # ax[:semilogx](vals, out[z_components[i]], label=z_components_names[i])
    ax[:semilogx](σ, out[z_components[i]], label=z_components_names[i])
    ax[:legend]()
    ax[:set_ylabel](z_components_names[i])
    # ax[:set_xlabel](raw"$\epsilon_{Vm}$")
    ax[:set_xlabel](raw"$\sigma$")
end
fig[:tight_layout]()

fig = figure(figsize=(10,10))
gs = gspec.GridSpec(1, 1, width_ratios=[1], height_ratios=[1])
slc = (lo,hi) -> pycall(pybuiltin("slice"), PyObject, lo, hi)
ax = subplot(get(gs, (0,0)))
# ax[:semilogx](vals, out[:cost], label="cost")
ax[:semilogx](σ, out[:cost], label="cost")
ax[:legend]()
ax[:set_ylabel]("objective function")
ax[:set_xlabel](raw"$\sigma$")
fig[:tight_layout]()

