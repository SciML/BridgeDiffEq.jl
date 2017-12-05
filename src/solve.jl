function solve{uType,tType,isinplace,AlgType<:BridgeAlgorithm}(
    prob::Union{AbstractODEProblem{uType,tType,isinplace},
                AbstractSDEProblem{uType,tType,isinplace}},
    alg::AlgType,
    timeseries=[],ts=[],ks=[];
    verbose=true,
    save_start=true, dt = nothing,
    timeseries_errors=true,
    callback=nothing,kwargs...)

    if dt == nothing
        error("dt required for fixed timestep methods.")
    end

    if verbose
        warned = !isempty(kwargs) && check_keywords(alg, kwargs, warnlist)
        warned && warn_compat()
    end

    if prob.callback != nothing || callback != nothing
        error("Bridge is not compatible with callbacks.")
    end

    u0 = prob.u0

    if isinplace
        f = (t,u) -> (du = similar(u); prob.f(t,u,du); Diagonal(du))
        if typeof(prob) <: AbstractSDEProblem
            g = (t,u) -> (du = similar(u); prob.g(t,u,du); Diagonal(du))
        end
    else
        f = prob.f
        if typeof(prob) <: AbstractSDEProblem
            if typeof(u0) <: Number
                g = (t,u) -> prob.g(t,u)
            else
                g = (t,u) -> Diagonal(prob.g(t,u))
            end
        end
    end

    t = prob.tspan[1]:dt:prob.tspan[2]
    if typeof(prob) <: AbstractSDEProblem
        W = Bridge.samplepath(t, zero(u0))
        samp = Bridge.sample(t, Bridge.Wiener{typeof(u0)}())
        if typeof(alg) <: BridgeEuler
            u = Bridge.solve!(Bridge.Euler(), W, u0, samp, (f,g))
        elseif typeof(alg) <: BridgeHeun
            u = Bridge.solve!(Bridge.StochasticHeun(), W, u0, samp, (f,g))
        elseif typeof(alg) <: BridgeSRK && typeof(u0) <: Number
            u = Bridge.solve!(Bridge.StochasticRungeKutta(), W, u0, samp, (f,g))
        else
            error("BridgeSRK is not compatible with non Number types")
        end
    else # ODE
        samp = Bridge.SamplePath(t,Vector{typeof(u0)}(length(t)))
        W = nothing
        if typeof(alg) <: BridgeR3
            u = Bridge.solve!(Bridge.R3(),samp,u0,f)
        elseif typeof(alg) <: BridgeBS3
            u,_ = Bridge.solve!(Bridge.BS3(),samp,u0,f)
        end
    end
    #=
    if save_start
        start_idx = 1
        ts = prob.tspan[1]:dt:prob.tspan[end]
    else
        start_idx = 2
        ts = (prob.tspan[1]+dt):dt:prob.tspan[end]
    end

    if typeof(u0) <: Union{AbstractArray}
        _timeseries = Vector{uType}(0)
        for i=start_idx:size(sol.q, 2)
            push!(_timeseries, reshape(view(sol.q, :, i-1)', sizeu))
        end
    elseif typeof(u0) <: Union{Tuple}
        _timeseries = Vector{typeof(u0[1])}(0)
        for i=start_idx:size(sol.q, 2)
            push!(_timeseries, reshape(view(sol.q, :, i-1)', sizeu))
        end
    else
        _timeseries = vec(sol.q)
    end
    =#

    build_solution(prob,  alg, u.tt, u.yy,
                   W = W,
                   timeseries_errors = timeseries_errors,
                   retcode = :Success)
end
