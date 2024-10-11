
# import Pkg
# Pkg.develop(path=".")
# Pkg.resolve()
# Pkg.instantiate()

using TiSR

# # ==================================================================================================
# # preparation
# # ==================================================================================================
# using DataFrames
# using CSV
#
# # load data # ------------------------------------------------------------------------------------
# file_path = "measurements_methanol_reduced_T150_P100_rho1000.txt"
#
# df = CSV.read(file_path, DataFrame)
#
# # set variables for algorithm
# data_matr = Matrix(df[!, 1:end-1])

# # or with synthetic data # -------------------------------------------------------------------------
# # -> 3 * (v1 * 5 + v2)^7 + exp(v1 * 5 + v2) 
# data_matr = rand(100, 3)
# data_matr[:, end] .= 3.0 .* (data_matr[:, 1] .* 5.0 .+ data_matr[:, 2]) .^ 7.0 + exp.(data_matr[:, 1] .* 5.0 .+ data_matr[:, 2])

# Netwons gravity 
data_matr = rand(100, 9)
data_matr[:, 1:2] .*= 1000
data_matr[:, 3:8] .-= 0.5 
data_matr[:, 3:8] .*= 100
data_matr[:, end] .= @. (
                        1e-5 * data_matr[:, 1] * data_matr[:, 2] / ( 
                          (data_matr[:, 3] - data_matr[:, 4])^2 
                        + (data_matr[:, 5] - data_matr[:, 6])^2 
                        + (data_matr[:, 7] - data_matr[:, 8])^2
                       )^0.5
                     )


split_inds = [collect(1:10), collect(10:100)]

# ==================================================================================================
# options -> specify some custom settings, where the default setting is unsatisfactory
# ==================================================================================================
pow(x, y) = abs(x)^y
pow2(x) = x^2

ops, data                          = Options(
    data_matr,
    data_split                     = data_split_params(
        # parts                      = [0.5, 0.5],
        split_inds                 = split_inds,
    ),
    p_binops                       = (1.0, 1.0, 1.0, 1.0, 1.0),#0.0,
    binops                         = (+,   -,   *,   /,   pow),#^,
    p_unaops                       = (1.0, 1.0, 0.0, 0.0, 1.0,  1.0),
    unaops                         = (exp, log, sin, cos, pow2, sqrt),
    general                        = general_params(
        n_gens                     = typemax(Int64),
        t_lim                      = 60 * 10.0,
        multithreading             = true,
        # remove_doubles_across_islands = true,
        migration_interval         = 100,
        island_extinction_interval = 1000,
    ),
    selection                      = selection_params(
        hall_of_fame_objectives    = [:ms_processed_e, :compl],
        selection_objectives       = [:ms_processed_e, :minus_abs_spearman, :compl, :age],
    ),
    grammar                        = grammar_params(
        max_compl                  = 30,
    ),
    fitting                        = fitting_params(
        max_iter                   = 10,
        # early_stop_iter          = 5,
        # lasso_factor             = 1e-3,
    ),
    # mutation                     = mutation_params(;
    #     p_crossover              = 5.0,
    #     p_point                  = 0.5,
    #     p_insert                 = 0.2,
    #     p_hoist                  = 0.5,
    #     p_subtree                = 0.2,
    #     p_drastic_simplify       = 0.1,
    #     p_insert_times_param     = 0.1,
    #     p_add_term               = 0.1,
    #     p_simplify               = 0.0,
    # )
);

# ==================================================================================================
# create a staring population
# ==================================================================================================
# -> variables must be v1, v2, ... and functions must be available in the function set

# start_pop = [
#     "5.0 * v1 + v2^3"
#     "5.0 * log(v1) + v2^3"
# ]

# ==================================================================================================
# main generational loop
# ==================================================================================================
hall_of_fame, population, prog_dict, stop_msg = generational_loop(data, ops);

# hall_of_fame, population, prog_dict, stop_msg = generational_loop(data, ops, start_pop);

# hot start with previous population # -------------------------------------------------------------
start_pop = vcat(hall_of_fame, population)
hall_of_fame, population, prog_dict, stop_msg = generational_loop(data, ops, start_pop);

# Inspect the results # ---------------------------------------------------------------------------
df_hall_of_fame = TiSR.convert_to_dataframe(hall_of_fame, ops, sort_by="max_are")
show(
    df_hall_of_fame[:, [:eqs_orig_rounded, :mare, :max_are, :weighted_compl, :n_params]],
    truncate = maximum(length, df_hall_of_fame.eqs_orig_rounded)
)

# show the Pareto front # --------------------------------------------------------------------------
using UnicodePlots

scatterplot(
    getfield.(hall_of_fame, :compl),
    getfield.(hall_of_fame, :mare),
    xlabel="complexity",
    ylabel="mean rel. dev.",
    yscale=:log10
)

# ==================================================================================================
# write the results to a fwf, csv, or a excel
# ==================================================================================================
save_to_fwf(hall_of_fame, ops)

# save_to_csv(hall_of_fame, ops)
# save_to_excel(hall_of_fame, population, prog_dict, ops)


node = TiSR.string_to_node("log(1.0)", ops)



length(ops.data_descript.split_inds[2])


