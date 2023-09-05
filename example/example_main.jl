
# ==================================================================================================
# load TiSR module
# ==================================================================================================

using TiSR

# ==================================================================================================
# preparation
# ==================================================================================================
# using DataFrames
# using CSV

# load data # ------------------------------------------------------------------------------------
# file_path = ""
#
# df = CSV.read(file_path, DataFrame)
#
# # reduce
# maxes = maximum(Matrix(df), dims=1)
# df ./= maxes
#
# # set variables for algorithm
# data = Matrix(df)

# or with synthetic data # -------------------------------------------------------------------------
data_matr = rand(1000, 3)
data_matr[:, end] .= 3.0 .* (data_matr[:, 1] .* 5.0 .+ data_matr[:, 2]) .^ 7.0 + exp.(data_matr[:, 1] .* 5.0 .+ data_matr[:, 2])
# -> 3 * (v1 * 5 + v2)^7 + exp(v1 * 5 + v2)

# prepare remainder for settings # -----------------------------------------------------------------
fit_weights = 1 ./ data_matr[:, end] # weights to minimize relative deviation
arbitrary_name = ""
parts = [0.8, 0.2]

# ==================================================================================================
# options -> specify some custom settings, where the default setting is unsatisfactory
# ==================================================================================================

ops, data = Options(
    data_descript=data_descript(
        data_matr;
        arbitrary_name = arbitrary_name,
        parts          = parts,
        fit_weights    = fit_weights
    ),
    general=general_params(
        n_gens          = typemax(Int64),
        pop_size        = 200,
        max_compl       = 30,
        pow_abs_param   = true,
        prevent_doubles = 1e-2,
        t_lim           = 60 * 2.0,
        multithreadding = true,
    ),
    fitting=fitting_params(
        early_stop_iter = 5,
        max_iter        = 30,
        lasso_factor    = 1e-7,
    ),
);

# ==================================================================================================
# create a staring population
# ==================================================================================================
# -> variables must be v1, v2, ... and functions must be available in the function set

# using SymbolicUtils
#
# @syms v1 v2 v3 v4 # ....
#
# start_pop = [
#     5.0 * v1 + v2^3,
#     5.0 * log(v1) + v2^3,
# ]
#
# start_pop = Node[string_to_node(eq, ops) for eq in start_pop]

# ==================================================================================================
# main generational loop
# ==================================================================================================
hall_of_fame, population, prog_dict = generational_loop(data, ops);

# hall_of_fame, population, prog_dict = generational_loop(data, ops, start_pop=start_pop);

# hot start with previous population (age information is lost) # -----------------------------------
# start_pop = vcat(hall_of_fame["node"], population["node"])
# start_pop = vcat(hall_of_fame["eqs_trees"], population["eqs_trees"], start_pop)
# hall_of_fame, population, prog_dict = generational_loop(data, ops, start_pop = start_pop);

# Inspect the results # ---------------------------------------------------------------------------

col = "mare"
perm = sortperm(hall_of_fame[col])
hall_of_fame[col][perm]
hall_of_fame["compl"][perm]
hall_of_fame["node"][perm]#[1:5]

# show the Pareto front # --------------------------------------------------------------------------

using UnicodePlots

scatterplot(
    hall_of_fame["compl"],
    hall_of_fame["mare"],
    xlabel="complexity",
    ylabel="mean rel. dev."
)

# ==================================================================================================
# write pareto optimal ones to excel
# ==================================================================================================
write_to_excel(hall_of_fame, population, prog_dict, ops)



