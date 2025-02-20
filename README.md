# Thermodynamics-informed Symbolic Regression - A Tool for the Thermodynamic Equation of State Development

 [![DOI](https://zenodo.org/badge/685443077.svg)](https://zenodo.org/badge/latestdoi/685443077)

This code base is based on the work of the project "Machine Learning and Optimal Experimental Design for Thermodynamic Property Modeling" in the DFG priority program [SPP2331](https://chemengml.org/).
It implements an adapted NSGA-II genetic algorithm for symbolic regression and is aimed at the development of thermodynamic equations of state.
The project and the code are a work-in-progress and there are many more features planned.
For more details on past and future uses, a paper accompanying this code is available at [arxiv](https://arxiv.org/abs/2309.02805).

# Getting started

This package is currently not registered at the Julia package registry.

To get started, install Julia 1.11 and run the following in a Julia REPL:

```julia
import Pkg
Pkg.add(url="https://github.com/scoop-group/TiSR#main")
```

A more detailed example on how to use it is provided in [example_main.jl](example/example_main.jl).
However, you may also start with the following minimal example:

```julia
using TiSR

# create synthetic data
data_matr = rand(100, 3)
data_matr[:, end] .= 3.0 .* (data_matr[:, 1] .* 5.0 .+ data_matr[:, 2]) .^ 7.0 + exp.(data_matr[:, 1] .* 5.0 .+ data_matr[:, 2])
# -> 3 * (v1 * 5 + v2)^7 + exp(v1 * 5 + v2)

ops, data =  Options(
    data_matr,
    general = general_params(
        t_lim     =  60 * 10.0,
    ),
    grammar = grammar_params(
        max_compl =  30,
    ),
);

# start the equation search
hall_of_fame, population, prog_dict, stop_msg = generational_loop(data, ops);

# Inspect the results
df_hall_of_fame = TiSR.convert_to_dataframe(hall_of_fame, ops, sort_by=:max_are)
show(
    df_hall_of_fame[:, [:eq_rounded, :mare, :compl]],
    truncate = maximum(length, df_hall_of_fame.eq_rounded)
)
```

# Available settings and their default values

In the following code snippet, all available, exposed settings are shown with their default values.
A brief explanation is provided for most settings as a comment.

```julia
ops, data = Options(
    data_matr;                                     # -> nxm matrix containing the n data points, m-1 variables and the output
    fit_weights    = abs.(1 ./ data_matr[:, end]), # -> weights for the data fitting -> residual .* weight
    binops         = (+,   -,   *,   /,   ^  ),    # -> binary function set to choose from
    unaops         = (exp, log, sin, cos, abs),    # -> unary function set to choose from
    data_split     = data_split_params(;
        parts      = [1.0],                     # -> how the data should be split. The first part is used for fitting, the second is used for early stopping, while the first AND the second part is used for the selection measures. The third part is used for the test-selection measures. The third part is used for the test-selection measures.
        split_inds = nothing,                   # -> rather then splitting randomly and automatically, the indices for each split can be specified in a Vector{Vector{Int64}}. If specified, parts cannot be specified.
    ),
    general = general_params(
        n_gens::Int64                          = typemax(Int64),                                           # -> number of generations to conduct
        t_lim::Float64                         = 60. * 5.,                                                 # -> time limit for the algorithm
        pop_size::Int64                        = 600,                                                      # -> number of individuals selected for next generation / population size
        parent_selection::Bool                 = true,                                                     # -> wheather to conduct parent selection or to use all individuals in the population as parents
        num_islands::Int64                     = 12,                                                       # -> number of parallel islands
        children_ratio::Float64                = 1.0,                                                      # -> the ratio of children that should be generated in each generation 0 ... 2
        migration_interval::Int64              = 200,                                                      # -> generation interval, in which an individual is moved to other islands. (ring topology)
        island_extinction_interval::Int64      = 5000,                                                     # -> interval in which all individuals from one islands are distributed across all other islands and the extiction islands starts from scratch. -> typemax(Int64) is off; 1000 ... 10000
        migrate_after_extinction_prob::Float64 = 1.0,                                                      # -> probability that an individual migrates to another islands, if its island goes extinct. -> 0 ... 1
        fitting_island_function::Function      = isle -> floor(isle / 2) % 2 == 0,                         # -> function to determine on which islands fitting is conducted. Must take an integer and return a bool
        hall_of_fame_migration_interval::Int64 = 1000,                                                     # -> interval in which a random individual from the hall of fame is returned to a random island
        always_drastic_simplify::Float64       = 1e-8,                                                     # -> for individuals with parameters smaller than `always_drastic_simplify` a copy is created, those parameters removed with some probability, and simplified accordingly. -> 0 is off; 1e-10 ... 1e-6
        remove_doubles_sigdigits::Int64        = 3,                                                        # -> removes individuals in an island if their MAE and MSE rouned to `remove_doubles_sigdigits` digits are the same. The one with the lowest complexity is retained. -> 0 is off; 2 ... 5
        remove_doubles_across_islands::Bool    = false,                                                    # -> same as remove_doubles_sigdigits, but across islands
        max_age::Int64                         = typemax(Int64),                                           # -> maximal age after which individuals are removed from the popoulation
        n_refitting::Int64                     = 1,                                                        # -> how many individuals from the hall_of_fame are copied and fitted again
        adaptive_compl_increment::Int64        = 100,                                                      # -> highest complexity in the hall of fame + `adaptive_compl_increment` is the highest allowed complexity for a new individual; -> Inf is off; 5 ... 10
        callback::Function                     = (hall_of_fame, population, gen, prog_dict, ops) -> false, # -> a function, which is executed in each iteration and allows more flexible termination. If the function returns true, the execution is terminated. For example, the following stops the equation search, if one individual in the hall of fame has a complexity lower than 30 and a mean absolute relative deviation of lower then 1e-5: `(hall_of_fame, population, gen, prog_dict, ops) -> any(i.measures[:compl] < 30 && i.measures[:mare] < 1e-5 for i in hall_of_fame)`
        multithreading::Bool                   = false,                                                    # -> whether to use multithreading for the fitting (most expensive). Not always faster -> depends on how expensive fitting is for the problem at hand. Also, for this to apply, Julia needs to be started with more threads, like `julia -t 4`.
        print_progress::Bool                   = true,                                                     # -> whether to print the elapsed time and some other KPIs.
        plot_hall_of_fame::Bool                = true,                                                     # -> whether to plot the hall of fame
        print_hall_of_fame::Bool               = true,                                                     # -> whether to print some of the individuals in the hall of fame. For this, `plot_hall_of_fame` must also be true
    ),
    measures = measure_params(
        additional_measures::Dict{Symbol, Function} = Dict(                 # -> specify a Dict{Symbol, Function} containing name and function pairs that calculate custom measures. TiSR offers some additional ones, which all start with `TiSR.get_measure_...`
            :one_minus_abs_spearman => get_measure_one_minus_abs_spearman,
            :mare                   => get_measure_mare,
            :max_are                => get_measure_max_are,
        )
    ),
    selection = selection_params(
        hall_of_fame_objectives::Vector{Symbol}    = [:ms_processed_e, :compl],                          # -> objectives for the hall_of_fame
        selection_objectives::Vector{Symbol}       = [:ms_processed_e, :one_minus_abs_spearman, :compl], # -> objectives for the Pareto-optimal selection part of selection
        hall_of_fame_niching_sigdigits::Int64      = 2,                                                  # -> number of significant digits to round hall_of_fame_objectives for hall_of_fame selection after their normalization. -> 2 ... 5
        population_niching_sigdigits::Int64        = 3,                                                  # -> number of significant digits to round selection_objectives for population selection after their normalization. -> 2 ... 5
        ratio_pareto_tournament_selection::Float64 = 0.5,                                                # -> ratio to which the selection is conducted using the Pareto-optimal selection vs. tournament selection
        tournament_size::Int64                     = 5,                                                  # -> tournament size
    ),
    fitting = fitting_params(
        max_iter::Int64                   = 10,                 # -> maximum iterations for parameter fitting. -> 10 ... 50 ==> biggest time consumer <==
        early_stop_iter::Int64            = 0,                  # -> how many iterations to account for early stopping regularization. to use, the data needs to be partitioned into at least 2 parts. The early stopping evaluation is performed on the second partition. -> 0 is off; 4 ... 10
        t_lim::Float64                    = Inf,                # -> time limit for parameter fitting of individual. -> Inf is off; 0.1 ... 0.5
        rel_f_tol_5_iter::Float64         = 1e-2 * 0.01,        # -> relative tolerance for parameter fitting. considered converged if relative improvement over 5 iterations is smaller. -> 0 is off; 1e-2 * 1.0 ... 1e-2 * 0.01
        lasso_factor::Float64             = 0.0,                # -> factor for the lasso regularization. pushing parameter values to 0. -> 0 is off; 1e-8 ... 1e-4
        pre_residual_processing::Function = (x, ind, ops) -> x, # -> processing of the equation output before the residual is calculated. The inds refer for the indices of the current residuals, which may be used to slice some data in the function like "(x, inds, ops) -> x ./= data[end][inds]"
        residual_processing::Function     = (x, ind, ops) -> x, # -> processing of the residuals. NOT an inplace function. The inds refer for the indices of the current residuals, which may be used to slice some data in the function like "(x, inds, ops) -> x ./ data[end][inds]"
    ),
    grammar = grammar_params(
        max_compl::Int64                           = 30,                        # -> max allowed complexity
        min_compl::Int64                           = 3,                         # -> min allowed complexity. -> 2 ... 3
        init_tree_depth::Int64                     = 4,                         # -> maximal initial tree depth. -> 3 ... 6
        max_nodes_per_term::Int64                  = typemax(Int64),            # -> maximal number of nodes per top-level term. All terms above are trimmed until they satisfy this threshold
        weighted_compl_dict::Dict{String, Float64} = Dict{String, Float64}(),   # -> weights for weighted_compl calculation. For any that are not included, 1.0 is assumed. The weights for variables and parameters "VAR" and "PARAM" may be used. An example is shown below at (2).
        bank_of_terms::Vector{String}              = String[],                  # -> specify terms that can be added via the add term mutation, whose relativ probability is set via the p_add_from_bank_of_terms parameter
        illegal_dict::Dict                         = Dict(),                    # -> check for illegal nestings in existing nodes. For it, ops.illegal_dict needs to be specified like below. Variables can be specified using "VAR" and parameters with "PARAM". An example is shown below at (1).
        custom_check_legal::Function               = (node, data, ops) -> true, # -> specify a custom function, which checks the legality of nodes. Must return true or false
    ),
    mutation = mutation_params(                    #|-> probabilites for the various mutations (don't need to add up to 1)
        p_crossover::Float64              = 10.0,  #|
        p_point::Float64                  = 1.0,   #|
        p_insert::Float64                 = 1.0,   #|
        p_hoist::Float64                  = 1.0,   #|
        p_subtree::Float64                = 0.5,   #|
        p_drastic_simplify::Float64       = 0.1,   #|-> remove parameter nodes with small values and simplify accordingly
        p_insert_times_param::Float64     = 0.1,   #|
        p_add_term::Float64               = 0.1,   #|
        p_simplify::Float64               = 0.1,   #|-> simplify with SymbolicUtils
        p_add_from_bank_of_terms::Float64 = 0.0,   #|-> probability to add a term from the provided bank_of_terms
    ),
    meta_data = Dict(), # -> can be used to provide data for use in, for example, user-defined functions like the callback or pre_residual_processing
)

# (1) example illegal_dict
# illegal_dict = Dict(
#     "^" =>   (lef = (),             rig = ("+", "-", "*", "/", "VAR", )),
#     "/" =>   (lef = (),             rig = ("-", "+")),
#     "log" => (lef = ("log", "exp"), rig = ()),
#     "exp" => (lef = ("exp", "log"), rig = ()),
#     "cos" => (lef = ("sin", "cos"), rig = ())
# )

# (2) example weighted_compl_dict
# weighted_compl_dict =  Dict(
#     "PARAM" => 1.2, "VAR"  => 1.0,
#     "+"     => 1.2, "-"    => 1.4,
#     "*"     => 1.0, "/"    => 1.6,
#     "^"     => 3.0, "pow"  => 3.5,
#     "neg"   => 1.4,
#     "sqrt"  => 2.0, "pow2" => 2.0,
#     "sin"   => 3.0, "cos"  => 3.0, "tanh" => 3.0,
#     "exp"   => 3.0, "log"  => 3.0,
# ),
```

# Funding

This work was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) within the Priority Programme SPP 2331: "Machine Learning in Chemical Engineering" - project no. 466528284 - HE 6077/14-1 and RI 2482/10-1

