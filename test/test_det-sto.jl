## deterministic acopf
dm = OPF.acopf_model(opfdata)
dm = OPF.acopf_solve(dm, opfdata)
dm_eval = setup(dm.m);               ## deterministic model evaluator
dm_zbar = deepcopy(dm_eval.last_x);  ## deterministic model equilibrium z̄
OPF.acopf_outputAll(dm, opfdata)

## stochastic acopf
sm = OPF.s_acopf_model(opfdata)
sm = OPF.acopf_solve(sm, opfdata)
sm_eval = setup(sm.m);               ## stochastic model evaluator
sm_zbar = deepcopy(sm_eval.last_x);  ## stochastc model equilibrium z̄
OPF.acopf_outputAll(sm, opfdata)

@testset "indexing" begin
@test norm([ getvalue(dm.m[:Pg]);
             getvalue(dm.m[:Qg]);
             getvalue(dm.m[:Vm]);
             getvalue(dm.m[:Va]) ] - dm_zbar) <= tol
@test norm([ getvalue(sm.m[:Pg]);
             getvalue(sm.m[:Qg]);
             getvalue(sm.m[:Vm]);
             getvalue(sm.m[:Va]);
             getvalue(sm.m[:Pd]);
             getvalue(sm.m[:Qd]) ] - sm_zbar) <= tol
@test norm([ getvalue(getindex(dm.m, :Pg));
             getvalue(getindex(dm.m, :Qg));
             getvalue(getindex(dm.m, :Vm));
             getvalue(getindex(dm.m, :Va)) ] - dm_zbar) <= tol
@test norm([ getvalue(getindex(sm.m, :Pg));
             getvalue(getindex(sm.m, :Qg));
             getvalue(getindex(sm.m, :Vm));
             getvalue(getindex(sm.m, :Va));
             getvalue(getindex(sm.m, :Pd));
             getvalue(getindex(sm.m, :Qd)) ] - sm_zbar) <= tol
end # testset
