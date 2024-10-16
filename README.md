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

ops, data         =  Options(
    data_matr,
    general       =  general_params(
        t_lim     =  60 * 10.0,
    ),
    grammar       =  grammar_params(
        max_compl =  30,
    ),
);

# start the equation search
hall_of_fame, population, prog_dict, stop_msg = generational_loop(data, ops);

# Inspect the results
df_hall_of_fame = TiSR.convert_to_dataframe(hall_of_fame, ops, sort_by="max_are")
show(
    df_hall_of_fame[:, [:eqs_simpl_rounded, :mare, :max_are, :compl]],
    truncate = maximum(length, df_hall_of_fame.eqs_orig_rounded)
)

```

# Available settings and their default values

In the following code snippet, all available, exposed settings are shown with their default values.
A brief explanation is provided for most settings as a comment.

```julia
ops, data = Options(
    data_matr;                                  # -> nxm matrix containing the n data points, m-1 variables and the output
    fit_weights = abs.(1 ./ data_matr[:, end]), # -> weights for the data fitting -> residual .* weight
    binops      = (+,   -,   *,   /,   ^  ),    # -> binary function set to choose from
    unaops      = (exp, log, sin, cos, abs),    # -> unary function set to choose from
    data_split  = data_split_params(;
        parts      = [1.0],                     # -> how the data should be splitted. The first part is used for fitting, the second is used for early stopping, while (currently) all are used to calculate the fit quality measures for selection. 
        split_inds = nothing,                   # -> rather then splitting randomly automatically, the indices for each split can be specified in a Vector{Vector{Int64}}. If specified, parts cannot be specified.
    ),
    general     = general_params(
        n_gens                          = typemax(Int64),                           # -> number of generations to conduct
        t_lim                           = 60. * 5.,                                 # -> time limit for the algorithm
        pop_size                        = 600,                                      # -> number of individuals selected for next generation
        num_islands                     = 12,                                       # -> number of parallel islands
        migration_interval              = 200,                                      # -> generation interval, in which an individual is moved to other islands. (ring topology)
        island_extinction_interval      = 5000,                                     # -> interval in which all individuals from one islands are distributed across all other islands and the extiction islands starts from scratch. -> typemax(Int64) is off; 1000 ... 10000
        migrate_after_extinction_prob   = 1.0,                                      # -> probability that an individual migrates to another islands, if its island goes extinct. -> 0 ... 1
        fitting_island_function         = isle -> floor(isle / 2) % 2 == 0,         # -> function to determine on which islands fitting is conducted
        hall_of_fame_migration_interval = 1000,                                     # -> interval in which a random individual from the hall of fame is returned to a random island
        always_drastic_simplify         = 1e-8,                                     # -> for individuals with parameters smaller than `always_drastic_simplify` a copy is created, those parameters removed, and simplified accordingly. -> 0 is off; 1e-10 ... 1e-6
        remove_doubles_sigdigits        = 3,                                        # -> removes individuals in an island if their  MAE and MSE rouned to `remove_doubles_sigdigits` digits are the same. The one with the lowest complexity is retained. -> 0 is off; 2 ... 5
        max_age                         = pop_size / num_islands,                   # -> maximal age after which individuals are removed from the popoulation
        adaptive_compl_increment        = Inf,                                      # -> highest complexity in the hall of fame + `adaptive_compl_increment` is the highest allowed complexity for a new individual; -> Inf is off; 5 ... 10
        callback                        = (hall_of_fame, population, ops) -> false, # -> a function, which is executed in each iteration and allows more flexible termination. For example, the following stops the equation search, if one individual in the hall of fame has a complexity lower than 30 and a mean absolute relative deviation of lower then 1e-5: `(hall_of_fame, population, ops) -> any(i.compl < 30 && i.mare < 1e-5 for i in hall_of_fame)`
        multithreading                  = false,                                    # -> whether to use multithreading for the fitting (most expensive). Not always faster -> depends on how expensive fitting is for the problem at hand. Also, for this to apply, Julia needs to be started with more threads, like `julia -t 4`.
        print_progress                  = true,                                     # -> whether to print the elapsed time and some other KPIs.
        plot_hall_of_fame               = true,                                     # -> whether to plot the hall of fame
        print_hall_of_fame              = true,                                     # -> whether to print some of the individuals in the hall of fame. For this, `plot_hall_of_fame` must also be true
    ),
    selection   = selection_params(
        hall_of_fame_objectives           = [:ms_processed_e, :compl],                      # -> objectives for the hall_of_fame
        selection_objectives              = [:ms_processed_e, :minus_abs_spearman, :compl], # -> objectives for the Pareto-optimal selection part of selection
        hall_of_fame_niching_sigdigits    = 2,                                              # -> number of significant digits to round hall_of_fame_objectives for hall_of_fame selection. -> 2 ... 5
        population_niching_sigdigits      = 3,                                              # -> number of significant digits to round selection_objectives for population selection. -> 2 ... 5
        ratio_pareto_tournament_selection = 0.5,                                            # -> ratio to which the selection is conducted using the Pareto-optimal selection vs. tournament selection
        tournament_size                   = 5,                                              # -> tournament size
    ),
    fitting     = fitting_params(
        max_iter                 = 10,                 # -> maximum iterations for parameter fitting. -> 10 ... 100 ==> biggest time consumer <==
        early_stop_iter          = 0,                  # -> how many iterations to account for early stopping regularization. to use, the data needs to be partitioned into at least 2 parts. The early stopping evaluation is performed on the second partition. -> 0 is off; 4 ... 10
        t_lim                    = Inf,                # -> time limit for parameter fitting of individual. -> Inf is off; 0.1 ... 0.5
        rel_f_tol_5_iter         = 1e-2 * 0.01,        # -> relative tolerance for parameter fitting. considered converged if relative improvement over 5 iterations is smaller. -> 0 is off; 1e-2 * 1.0 ... 1e-2 * 0.01
        lasso_factor             = 0.0,                # -> factor for the lasso regularization. pushing parameter values to 0. -> 0 is off; 1e-8 ... 1e-4
        pre_residual_processing! = (x, ind, ops) -> x, # -> processing of the equation output before the residual is calculated. Must be an inplace function. The inds refer for the indices of the current residuals, which may be used to slice some data in the function like "(x, inds) -> x ./= data[end][inds]"
        residual_processing      = (x, ind, ops) -> x, # -> processing of the residuals. NOT an inplace function. The inds refer for the indices of the current residuals, which may be used to slice some data in the function like "(x, inds) -> x ./ data[end][inds]"
    ),
    grammar     = grammar_params(
        max_compl           = 30,     # -> max allowed complexity.
        min_compl           = 3,      # -> min allowed complexity. -> 2 ... 3
        init_tree_depth     = 4,      # -> maximal initial tree depth. -> 3 ... 6
        max_nodes_per_term  = Inf,    # -> maximal number of nodes per top-level term. All terms above are trimmed until they satisfy this threshold.
        illegal_dict        = Dict(), # -> check for illegal nestings in existing nodes. For it, ops.illegal_dict needs to be specified like below. Variables can be specified using "VAR" and parameters with "PARAM". An example is shown below at (1).
        weighted_compl_dict = Dict(), # -> weights for weighted_compl calculation. For any that are not included, 1.0 is assumed. The weights for variables and parameters "VAR" and "PARAM" may be used. An example is shown below at (2).
    ),
    mutation    = mutation_params(    # -> probabilites for the various mutations (don't need to add up to 1)
        p_crossover          = 10.0,  #|
        p_point              = 1.0,   #|
        p_insert             = 1.0,   #|
        p_hoist              = 1.0,   #|
        p_subtree            = 0.5,   #|
        p_drastic_simplify   = 0.1,   #|-> remove parameter nodes with values smaller than 1e-1
        p_insert_times_param = 0.1,   #|
        p_add_term           = 0.1,   #|
        p_simplify           = 0.1,   #|-> simplify with SymbolicUtils
    ),
    meta_data = Dict() # -> can be used to provide data for use in, for example, user-defined functions like the callback or pre_residual_processing!
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
#     "PARAM" => 1.5, "VAR" => 1.0,
#     "+"     => 1.0, "-"   => 1.5,
#     "*"     => 2.0, "/"   => 2.5, "^"   => 3.0,
#     "exp"   => 2.0, "log" => 2.0, "sin" => 2.0, "cos" => 2.0,
# )
```

# Funding

This work was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) within the Priority Programme SPP 2331: "Machine Learning in Chemical Engineering" - project no. 466528284 - HE 6077/14-1 and RI 2482/10-1

