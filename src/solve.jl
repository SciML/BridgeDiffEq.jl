function solve(prob::Union{DiffEqBase.AbstractODEProblem{uType, tType, isinplace},
                           DiffEqBase.AbstractSDEProblem{uType, tType, isinplace}},
               alg::AlgType,
               timeseries = [], ts = [], ks = [];
               verbose = true,
               dt = nothing,
               timeseries_errors = true,
               callback = nothing,
               kwargs...) where {uType, tType, isinplace, AlgType <: BridgeAlgorithm}
    if dt == nothing
        error("dt required for fixed timestep methods.")
    end

    if verbose
        warned = !isempty(kwargs) && check_keywords(alg, kwargs, warnlist)
        warned && warn_compat()
    end

    if callback != nothing || :callback in keys(prob.kwargs)
        error("Bridge is not compatible with callbacks.")
    end

    u0 = prob.u0
    p = prob.p

    if isinplace
        f = (t, u) -> (du = similar(u); prob.f(du, u, p, t); Diagonal(du))
        if typeof(prob) <: DiffEqBase.AbstractSDEProblem
            g = (t, u) -> (du = similar(u); prob.g(du, u, p, t, u); Diagonal(du))
        end
    else
        f = (t, u) -> prob.f(u, p, t)
        if typeof(prob) <: DiffEqBase.AbstractSDEProblem
            if typeof(u0) <: Number
                g = (t, u) -> prob.g(u, p, t)
            else
                g = (t, u) -> Diagonal(prob.g(u, p, t))
            end
        end
    end

    t = prob.tspan[1]:dt:prob.tspan[2]
    if typeof(prob) <: DiffEqBase.AbstractSDEProblem
        W = Bridge.samplepath(t, zero(u0))
        samp = Bridge.sample(t, Bridge.Wiener{typeof(u0)}())
        if typeof(alg) <: BridgeEuler
            u = Bridge.solve!(Bridge.Euler(), W, u0, samp, (f, g))
        elseif typeof(alg) <: BridgeHeun
            u = Bridge.solve!(Bridge.StochasticHeun(), W, u0, samp, (f, g))
        elseif typeof(alg) <: BridgeSRK && typeof(u0) <: Number
            u = Bridge.solve!(Bridge.StochasticRungeKutta(), W, u0, samp, (f, g))
        else
            error("BridgeSRK is not compatible with non Number types")
        end
    else # ODE
        samp = Bridge.SamplePath(t, Vector{typeof(u0)}(undef, length(t)))
        W = nothing
        if typeof(alg) <: BridgeR3
            u = Bridge.solve!(Bridge.R3(), samp, u0, f)
        elseif typeof(alg) <: BridgeBS3
            u, _ = Bridge.solve!(Bridge.BS3(), samp, u0, f)
        end
    end

    DiffEqBase.build_solution(prob, alg, u.tt, u.yy,
                              W = W,
                              timeseries_errors = timeseries_errors,
                              retcode = :Success)
end
