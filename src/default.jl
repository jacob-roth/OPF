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
  D[:pw_angle_limits]= false   ## pairwise angle bound limits
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
