options = DefaultOptions()
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true

opfdata, dp, scm = get_n1_limits("case118", path, options, 0.45)
# opfdata, dp, scm = get_n1_limits("case118", path, options, 0.57)  ## this works too...
opfdata, dp, scm = get_n1_limits("case30", path, options, 0.58)
acopf_outputAll(scm, opfdata, options)

# casedata = load_case("30-files/"* "n1-lowdamp", "/Users/jakeroth/Desktop/planning-large-deviation/data/cases", other=true)
# dm = acopf_model(casedata.opf)
# dm = acopf_solve(dm, casedata.opf)
# V̄a = deepcopy(getvalue(getindex(dm.m, :Va)))
# V̄a .-= V̄a[casedata.opf.bus_ref]
# V̄m = deepcopy(getvalue(getindex(dm.m, :Vm)))
# P̄g = deepcopy(getvalue(getindex(dm.m, :Pg)))
# Q̄g = deepcopy(getvalue(getindex(dm.m, :Qg)))