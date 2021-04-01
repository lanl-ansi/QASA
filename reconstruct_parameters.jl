using DataFrames, CSV, ArgParse, Distributed
SPIN_TABLE_FILENAME = "spin_table.csv"
PARAM_FILENAME="model_parameters.csv"

@everywhere include("mle.jl")

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
