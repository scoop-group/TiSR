
# import Pkg
# Pkg.develop(path=".")

using TiSR
using DataFrames
using CSV

# # ==================================================================================================
# # preparation
# # ==================================================================================================
# # load data # ------------------------------------------------------------------------------------
# file_path = "path_to_your_data_file.txt"
#
# df = CSV.read(file_path, DataFrame)
#
# # set variables for algorithm
# data_matr = Matrix(df)

# or with synthetic data # -------------------------------------------------------------------------
data_matr = rand(100, 3)
data_matr[:, end] .= 3.0 .* (data_matr[:, 1] .* 5.0 .+ data_matr[:, 2]) .^ 7.0 + exp.(data_matr[:, 1] .* 5.0 .+ data_matr[:, 2])
# -> 3 * (v1 * 5 + v2)^7 + exp(v1 * 5 + v2)

# ==================================================================================================
# options -> specify some custom settings
# ==================================================================================================
# # custom functions
# pow(x, y) = abs(x)^y
# pow2(x) = x^2

ops, data                       = Options(
    data_matr,
    binops                      = (+, -, *, /, ^),
    unaops                      = (exp, log, sqrt),
    general                     = general_params(
        t_lim                   = 60 * 60.0,
        multithreading          = true,
    ),
    grammar                     = grammar_params(
        max_compl               = 30,
    ),
    fitting                     = fitting_params(
        max_iter                = 10,
    ),
);

# ==================================================================================================
# main generational loop
# ==================================================================================================
hall_of_fame, population, prog_dict, stop_msg = generational_loop(data, ops);

# ==================================================================================================
# with starting population
# ==================================================================================================
# -> variables must be v1, v2, ... and functions must be available in the function set

# start_pop = [
#     "1.0 * v1^1.0 * v2^1.0"
# ]

# hall_of_fame, population, prog_dict, stop_msg = generational_loop(data, ops, start_pop);

# ==================================================================================================
# hot start with previous population
# ==================================================================================================
# start_pop = vcat(hall_of_fame, population)
# hall_of_fame, population, prog_dict, stop_msg = generational_loop(data, ops, start_pop);

# ==================================================================================================
# Inspect the results
# ==================================================================================================
df_hall_of_fame = TiSR.convert_to_dataframe(hall_of_fame, ops, sort_by=:mare)
show(
    df_hall_of_fame[:, [:eq_rounded, :mare, :compl]],
    truncate = maximum(length, df_hall_of_fame.eq_rounded)
)

# show the Pareto front # --------------------------------------------------------------------------
using UnicodePlots

scatterplot(
    df_hall_of_fame[!, :compl],
    df_hall_of_fame[!, :mare],
    xlabel="complexity",
    ylabel="mean rel. dev.",
    yscale=:log10,
    marker=:o,
    unicode_exponent = false,
)

# ==================================================================================================
# write the results to a fwf, csv, or a excel
# ==================================================================================================
save_to_fwf(hall_of_fame, ops)

# save_to_csv(hall_of_fame, ops)
# save_to_excel(hall_of_fame, population, prog_dict, ops)
