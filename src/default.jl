## -----------------------------------------------------------------------------
## default options
## -----------------------------------------------------------------------------
function DefaultOptions()
  D                  = Dict()
  D[:lossless]       = false   ## lossless transmission lines?
  D[:loss_scale]     = 1.0     ## scaling of real component of admittance elements
  D[:current_rating] = false   ## use |i| for line limits?
  D[:remove_Bshunt]  = false   ## remove susceptance from shunt?
  D[:remove_tap]     = false   ## remove tap adjustments?
  D[:print_level]    = 0       ## 0: no printing
  D[:feasibility]    = false   ## feasibility problem only (no OPF objective)?
  D[:sc_total_obj]   = false   ## minimize cost under contingency cases in security-constrained acopf
  D[:rampup]         = 2.0     ## Pg coupling up
  D[:rampdn]         = 0.5     ## Pg coupling down
  D[:ctg_Vm_adj]     = 0.05    ## allowable voltage increase/decrease from Vmax/Vmin in contingency case
  zip0 = Dict(:alpha => 1.0,   ## α + β + γ = 1
              :beta  => 0.0,
              :gamma => 0.0,
              :V0    => 1.0)   ## nominal voltage
  sol0 = Dict(:mu_P  => 0.0,   ## P solar penetration (mu < 0)
              :I_P   => false, ## P irradience -> active power forecast
              :mu_Q  => 0.0,   ## Q solar penetration (mu < 0)
              :I_Q   => false) ## Q irradience -> reactive power forecast
  D[:zip]            = zip0
  D[:sol]            = sol0
  D[:tol]            = 1e-9
  #
  # Options for enforcing transition rate limits
  #
  D[:damping]        = 1.0     ## damping coefficient
  D[:temperature]    = 0.0001  ## temperature
  D[:constr_limit_scale] = 1.2 ## line limits will be multiplied by this factor
  D[:shed_load]      = true    ## shed load?
  D[:ratelimit]      = 10.0    ## max upper bound on the transition rates
  D[:iterlim]        = 1e9     ## max number of rounds of adding constraints
  D[:VOLL]           = Inf     ## Value of lost load (can be Inf)
  D[:psd_constraint] = false   ## Enforce convexity of KKT subproblem
  D[:high_temp_adj]  = false   ## High temperature adjustment to prefactor
  D[:slack0]         = true    ## fix slack reference angle to be zero
  D[:pw_angle_limits]= false   ## pairwise angle difference constraints
  D[:ramp_pct]       = 0.1     ## |Pg - Pg_contingency| <= ramp_pct * Pg_max
  D[:ctg_feas]       = true    ## contingency objective: true = 0 or false = sum(Pg cost)
  D[:parallel] = (nprocs() > 1)## evaluate exit rates in parallel
  return D
end

function DefaultAdjustments()
    D = Dict()
    D[:Pg_hi] = Dict(:v => 0.0, :i => Int64[])  ## dec scale bound-constraints at bus `i` as (1 - v)
    D[:Pg_lo] = Dict(:v => 0.0, :i => Int64[])  ## inc scale bound-constraints at bus `i` as (1 + v)
    D[:Qg_hi] = Dict(:v => 0.0, :i => Int64[])  ## dec scale bound-constraints at bus `i` as (1 - v)
    D[:Qg_lo] = Dict(:v => 0.0, :i => Int64[])  ## inc scale bound-constraints at bus `i` as (1 + v)
    D[:Vm_hi] = Dict(:v => 0.0, :i => Int64[])  ## dec scale bound-constraints at bus `i` as (1 - v)
    D[:Vm_lo] = Dict(:v => 0.0, :i => Int64[])  ## inc scale bound-constraints at bus `i` as (1 + v)
    return D
end

function DefaultLoading();    return 0.50; end  # 0.85
function DefaultAdjPF();      return 0.0;  end
function DefaultFeasTol();    return 1e-6; end
function DefaultBuffer();     return 0.25; end
function DefaultMaxIter();    return 10;   end
function DefaultMaxLimit();   return 1e7;  end
function DefaultVPct();       return 1.0;  end
function DefaultTightenAmt(); return 1e-3; end
