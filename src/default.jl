## -----------------------------------------------------------------------------
## default options
## -----------------------------------------------------------------------------
function DefaultOptions()
  D                  = Dict()
  D[:lossless]       = false   ## lossless transmission lines?
  D[:current_rating] = false   ## use |i| for line limits?
  D[:remove_Bshunt]  = false   ## remove susceptance from shunt?
  D[:remove_tap]     = false   ## remove tap adjustments?
  D[:print_level]    = 0       ## 0: no printing
  D[:feasibility]    = false   ## feasibility problem only (no OPF objective)?
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
  return D
end