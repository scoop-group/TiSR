
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

# set variables for algorithm
# data_matr = Matrix(df)

# or with synthetic data # -------------------------------------------------------------------------
# -> 3 * (v1 * 5 + v2)^7 + exp(v1 * 5 + v2) 
data_matr = rand(100, 3)
data_matr[:, end] .= 3.0 .* (data_matr[:, 1] .* 5.0 .+ data_matr[:, 2]) .^ 7.0 + exp.(data_matr[:, 1] .* 5.0 .+ data_matr[:, 2])

# # Netwons gravity 
# data_matr = rand(1000, 9)
# data_matr[:, 1:2] .*= 1000
# data_matr[:, 3:8] .-= 0.5 
# data_matr[:, 3:8] .*= 100
# data_matr[:, end] .= @. (
#                         1e-5 * data_matr[:, 1] * data_matr[:, 2] / ( 
#                           (data_matr[:, 3] - data_matr[:, 4])^2 
#                         + (data_matr[:, 5] - data_matr[:, 6])^2 
#                         + (data_matr[:, 7] - data_matr[:, 8])^2
#                        )^0.5
#                      )

# prepare remainder for settings # -----------------------------------------------------------------
fit_weights = 1 ./ data_matr[:, end] # weights to minimize relative deviation
parts = [0.8, 0.2]

# ==================================================================================================
# options -> specify some custom settings, where the default setting is unsatisfactory
# ==================================================================================================

# TODO: improve asserts and testing of illegal_dict

pow_abs(x, y) = abs(x)^y
pow2(x) = x^2

ops, data                              = Options(
    data_matr,
    fit_weights                        = fit_weights,
    p_binops                           = (1.0, 1.0, 1.0, 1.0, 0.0, 1.0),
    binops                             = (+,   -,   *,   /,   ^, pow_abs),
    p_unaops                           = (1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0),
    unaops                             = (exp, log, sin, cos, abs, pow2, sqrt),
    parts                              = parts,
    general                            = general_params(
        n_gens                         = typemax(Int64),
        t_lim                          = 60 * 10.0,
        pop_size                       = 500,
        remove_doubles_sigdigits       = 3,
        remove_doubles_across_islands  = true,
        multithreadding                = true,
        always_drastic_simplify        = 1e-7,
        adaptive_compl_increment       = 10,
        # callback = (hall_of_fame, population, ops) -> any(i.compl < 30 && i.mare < 1e-5 for i in hall_of_fame)
    ),
    selection                          = selection_params(
        population_niching_sigdigits   = 3,
        selection_objectives           = [:ms_processed_e, :minus_abs_spearman, :compl, :age],
    ),
    grammar                            = grammar_params(
        max_compl                      = 30,
        illegal_dict = Dict(
            "pow_abs" => (lef = (),             rig = ("+", "-", "*", "/", "^", "pow_abs", "pow2", "sqrt", "exp", "log", "sin", "cos", "VAR")),
            "^"       => (lef = (),             rig = ("+", "-", "*", "/", "^", "pow_abs", "pow2", "sqrt", "exp", "log", "sin", "cos", "VAR")),
            "sin"     => (lef = ("sin", "cos"), rig = ()),
            "cos"     => (lef = ("sin", "cos"), rig = ()),
            "exp"     => (lef = ("exp", "log"), rig = ()),
            "log"     => (lef = ("exp", "log"), rig = ()),
        )
    ),
    fitting                            = fitting_params(
        early_stop_iter                = 0,
        max_iter                       = 15,
        lasso_factor                   = 1e-7,
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
hall_of_fame, population, prog_dict, stop_msg = generational_loop(data, ops);

# hall_of_fame, population, prog_dict = generational_loop(data, ops, start_pop=start_pop);

# hot start with previous population (age information is lost) # -----------------------------------
# start_pop = vcat(hall_of_fame["node"], population["node"])
# start_pop = vcat(hall_of_fame["eqs_trees"], population["eqs_trees"], start_pop)
# hall_of_fame, population, prog_dict = generational_loop(data, ops, start_pop = start_pop);

# Inspect the results # ---------------------------------------------------------------------------

col = "mare"
perm = sortperm(hall_of_fame[col])
hall_of_fame[col][perm]
hall_of_fame["mare"][perm]
hall_of_fame["node"][perm]#[1:5]


for n in sort(string.(population["node"]), by=length)
    println(n)
end


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

