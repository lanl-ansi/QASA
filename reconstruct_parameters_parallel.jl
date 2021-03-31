using JuMP, Ipopt, DataFrames, CSV, ArgParse, Distributed
SPIN_TABLE_FILENAME = "spin_table.csv"
PARAM_FILENAME="model_parameters.csv"


@everywhere begin
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
end
function main(args)
    directory = args["directory"]
    println("Loading spin table")

    csv_df = DataFrame(CSV.File("$(directory)/$(SPIN_TABLE_FILENAME)"))

    col_names = names(csv_df)
    col_names = col_names[3:ncol(csv_df)]
    h_list = []
    h_in = csv_df.h
    p_list = []
    for i in 1:length(col_names)
        p = (csv_df[:, Symbol(col_names[i])])/csv_df."samples"[i]
        push!(p_list, p)
        push!(h_list, h_in)
    end
    println("Starting parameter reconstruction via MLE")

    parameters = pmap(parallel_mle, p_list, h_list, col_names)
    β_list = [p[1] for p in parameters if !isnan(p[1])]
    b_list = [p[2]/p[1] for p in parameters if !isnan(p[2])]
    γ_list = [p[3]/p[1] for p in parameters if !isnan(p[3])]
    η_list = [p[4]/p[1] for p in parameters if !isnan(p[4])]
    col_names = [col_names[i] for i in 1:length(col_names) if !isnan(parameters[i][1])]

    col_names = [split(c,"_")[end] for c in string.(col_names)]
    result_df = DataFrame(spin = col_names, β = β_list, b = b_list, γ = γ_list, η = η_list)
    println("Writing parameters to output csv file")
    CSV.write("$(directory)/$(PARAM_FILENAME)", result_df)
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
