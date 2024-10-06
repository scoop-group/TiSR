
# ==================================================================================================
# load TiSR module
# ==================================================================================================

# import Pkg
# Pkg.develop(path=".")
# Pkg.resolve()
# Pkg.instantiate()

using TiSR

# ==================================================================================================
# preparation
# ==================================================================================================
using DataFrames
using CSV

# load data # ------------------------------------------------------------------------------------
file_path = "measurements_methanol_reduced_T150_P100_rho1000.txt"

df = CSV.read(file_path, DataFrame)

# set variables for algorithm
data_matr = Matrix(df[!, 1:end-1])

# # or with synthetic data # -------------------------------------------------------------------------
# # -> 3 * (v1 * 5 + v2)^7 + exp(v1 * 5 + v2) 
# data_matr = rand(100, 3)
# data_matr[:, end] .= 3.0 .* (data_matr[:, 1] .* 5.0 .+ data_matr[:, 2]) .^ 7.0 + exp.(data_matr[:, 1] .* 5.0 .+ data_matr[:, 2])

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
parts = [1.0]

# ==================================================================================================
# options -> specify some custom settings, where the default setting is unsatisfactory
# ==================================================================================================
pow_abs(x, y) = abs(x)^y
pow2(x) = x^2

ops, data                             =  Options(
    data_matr,
    fit_weights                       =  fit_weights,
    p_binops                          =  (1.0, 1.0, 1.0, 1.0, 1.0, 0.0),
    binops                            =  (+,   -,   *,   /,   ^, pow_abs),
    p_unaops                          =  (1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0),
    unaops                            =  (exp, log, sin, cos, abs, pow2, sqrt),
    parts                             =  parts,
    general                           =  general_params(
        n_gens                        =  typemax(Int64),
        t_lim                         =  60 * 10.0,
        island_extinction_interval    =  1000,            # TODO: experiment
        remove_doubles_sigdigits      =  3,               # TODO: experiment
        remove_doubles_across_islands =  false,            # TODO: experiment
        multithreading                =  true,
        always_drastic_simplify       =  1e-7,
        adaptive_compl_increment      =  20, # TODO: broken # TODO: add island extiction interval parameter and experiment
        # callback                    =  (hall_of_fame, population, ops) -> any(i.compl < 30 && i.mare < 1e-5 for i in hall_of_fame)
    ),
    selection                         =  selection_params(
        population_niching_sigdigits  =  3,
        selection_objectives          =  [:ms_processed_e, :minus_abs_spearman, :compl], # , :age # TODO: experimetn age
        hall_of_fame_objectives       =  [:ms_processed_e, :compl],
    ),
    grammar                           =  grammar_params(
        max_compl                     =  30,
        max_nodes_per_term            =  Inf,
        illegal_dict                  =  Dict(
            # # "pow_abs"             => (lef = (),              rig = ("+", "-", "*", "/", "^", "pow_abs", "pow2", "sqrt", "exp", "log", "sin", "cos", "VAR")),
            # # "^"                   => (lef = (),              rig = ("+", "-", "*", "/", "^", "pow_abs", "pow2", "sqrt", "exp", "log", "sin", "cos", "VAR")),
            "pow_abs"                 => (lef = (),              rig = ("+", "-", "*", "/", "^", "pow_abs", "pow2", "sqrt", "exp", "log", "sin", "cos")),
            "^"                       => (lef = (),              rig = ("+", "-", "*", "/", "^", "pow_abs", "pow2", "sqrt", "exp", "log", "sin", "cos")),
            "/"                       => (lef = (),              rig = ("+", "-")),
            "sin"                     => (lef = ("sin",  "cos"), rig = ()),
            "cos"                     => (lef = ("sin",  "cos"), rig = ()),
            "exp"                     => (lef = ("exp",  "log"), rig = ()),
            "log"                     => (lef = ("exp",  "log"), rig = ()),
            "sqrt"                    => (lef = ("sqrt", "log"), rig = ()),
        )
    ),
    fitting                           =  fitting_params(
        early_stop_iter               =  0,
        max_iter                      =  15,
        lasso_factor                  =  1e-5,
    ),
    mutation                          =  mutation_params(;
        p_crossover                   =  5.0,
        p_point                       =  0.5,
        p_insert_times_param          =  0.5,
        p_drastic_simplify            =  0.2,
        p_insert                      =  0.2,
        p_hoist                       =  0.2,
        p_subtree                     =  0.2,
        p_add_term                    =  0.1,
        p_simplify                    =  0.1,
    )
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
df_hall_of_fame = TiSR.convert_to_dataframe(hall_of_fame, ops, sort_by="mare")
show(
    df_hall_of_fame[:, [:eqs_orig_rounded, :mare, :max_are, :compl]],
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
# write the results to a fwf or a excel
# ==================================================================================================
save_to_fwf(hall_of_fame, ops)

# save_to_excel(hall_of_fame, population, prog_dict, ops)







