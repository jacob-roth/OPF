options = DefaultOptions()
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true

opfdata = load_case("case30", path; other=false)
opfmodel = acopf_model(opfdata, options)
opfmodel = acopf_solve(opfmodel, opfdata)
acopf_outputAll(opfmodel, opfdata, options)
# V̄a = deepcopy(getvalue(getindex(dm.m, :Va)))
# V̄a .-= V̄a[casedata.opf.bus_ref]
# V̄m = deepcopy(getvalue(getindex(dm.m, :Vm)))
# P̄g = deepcopy(getvalue(getindex(dm.m, :Pg)))
# Q̄g = deepcopy(getvalue(getindex(dm.m, :Qg)))