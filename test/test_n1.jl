options = DefaultOptions()
options[:current_rating] = true
options[:lossless]       = true
options[:remove_Bshunt]  = true
options[:remove_tap]     = true

opfdata, dp, scm = get_n1_limits("case118", path, options, 0.45)
acopf_outputAll(scm, opfdata, options)
# d  = setup(scm.m);
# x  = MathProgBase.getsolution(scm.m.internalModel);
# F1 = d.subexpressions_as_julia_expressions[236+1]
# F2 = d.subexpressions_as_julia_expressions[236+2]
# F3 = d.subexpressions_as_julia_expressions[236+3]
# flowmax = (opfdata.lines.rateA ./ 100).^2
# flowmax1 = (new_ratings[1]/100)^2
# flowmax2 = (new_ratings[2]/100)^2
# flowmax3 = (new_ratings[3]/100)^2
#
# e = [getvalue(scm.m[Symbol("F$(i)")]) for i in 1:186]
# (e + flowmax) ./ flowmax
#
# (e - opfdata.lines.rateA) ./ opfdata.lines.rateA