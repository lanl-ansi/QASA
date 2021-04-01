using JuMP, Ipopt

function mean_quantum(h, u)
    if h^2 + u^2 == 0
        return h
    else
        return h*tanh(sqrt(h^2 + u^2))/sqrt(h^2 + u^2)
    end
end
function mean_noise_quantum(h, β, βb, βγ, βη)
    mean_upper = mean_quantum(β*h + βb + βη, βγ*h)
    mean_lower = mean_quantum(β*h + βb - βη, βγ*h)
    return (mean_upper + mean_lower) / 2
end
function mle_noise_quantum(p_down, h, spin)
    n = length(p_down)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)
    @variable(model, β >= 0)
    @variable(model, βb)
    @variable(model, βγ >= 0)
    @variable(model, βη >= 0)

    register(model, :mean_noise_quantum, 5, mean_noise_quantum, autodiff = true)
    @NLexpression(model, log_likelihood,
        sum(p_down[i] * log(1+mean_noise_quantum(h[i], β, βb, βγ, βη)) + ((1 - p_down[i]) * log(1 - mean_noise_quantum(h[i], β, βb, βγ, βη))) for i in 1:n))

    @NLobjective(model, Max, log_likelihood)
    optimize!(model)
    if termination_status(model) != MOI.LOCALLY_SOLVED
        throw(ErrorException("Solver status: $(termination_status(model))"))
    end
    println("Finished $(spin)")
    return value(β), value(βb), value(βγ), value(βη)
end
function parallel_mle(p, h, spin)
    try
        return mle_noise_quantum(p, h, spin)
    catch e
        println(e.msg)
        println("Skipping $(spin)")
        return [NaN, NaN, NaN, NaN]
    end
end
