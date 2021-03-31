using JuMP, Ipopt, DataFrames, CSV, ArgParse, Distributed

SPIN_TABLE_FILENAME="spin_table.csv"
PARAM_FILENAME="model_parameters.csv"
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
function mle_noise_quantum(p_down, h)
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
        throw(ErrorException("Solver status: $(termination_status(model))."))
    end
    return value(β), value(βb), value(βγ), value(βη)
end
function main(args)
    println("Loading spin table")
    csv_df = DataFrame(CSV.File("$(args["directory"])/$(SPIN_TABLE_FILENAME)"))

    col_names = names(csv_df)
    col_names = col_names[3:ncol(csv_df)]
    spin_names=[]
    β_list = []
    b_list = []
    γ_list = []
    η_list = []
    h_in = csv_df.h
    println("Starting parameter reconstruction via MLE")
    for i in 1:length(col_names)
        p = (csv_df[:, Symbol(col_names[i])])/csv_df."samples"[i]
        try
            β, βb, βγ, βη = mle_noise_quantum(p, h_in)
            push!(β_list, β)
            push!(b_list, βb/β)
            push!(γ_list, βγ/β)
            push!(η_list, βη/β)
            push!(spin_names, col_names[i])
            println("Finished $(col_names[i])")
        catch e
            println(e.msg)
            println("Skipping spin: $(col_names[i])")
        end
    end
    spin_names = [split(c,"_")[end] for c in string.(spin_names)]
    result_df = DataFrame(spin = spin_names, β = β_list, b = b_list, γ = γ_list, η = η_list)
    println("Writing parameters to output csv file")
    CSV.write("$(args["directory"])/$(PARAM_FILENAME)", result_df)
end
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--directory", "-d"
            help = "a directory to read and write csv files"
            required = true
    end

    return parse_args(s)
end


if isinteractive() == false
    main(parse_commandline())
end
