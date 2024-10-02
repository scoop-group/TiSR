# Thermodynamics-informed Symbolic Regression - A Tool for the Thermodynamic Equation of State Development

 [![DOI](https://zenodo.org/badge/685443077.svg)](https://zenodo.org/badge/latestdoi/685443077)

This code base is based on the work of the project "Machine Learning and Optimal Experimental Design for Thermodynamic Property Modeling" in the DFG priority program [SPP2331](https://chemengml.org/).
It implements a NSGA-II genetic algorithm for symbolic regression and is aimed at the development of thermodynamic equations of state.
The project and the code are a work-in-progress and there are many more features planned.
For more details on past and future uses, a paper accompanying this code is available at [arxiv](https://arxiv.org/abs/2309.02805).

# Getting started

This package is currently not registered at the Julia package registry. 

To get started, install Julia 1.9 and run the following in a Julia REPL:

```julia
import Pkg
Pkg.add(url="https://github.com/scoop-group/TiSR#main")
```

A more detailed example on how to use it is provided in [example_main.jl](example/example_main.jl).
However, you may also start with the following minimal example:

```julia
using TiSR

# create synthetic data
data_matr = rand(1000, 3)
data_matr[:, end] .= 3.0 .* (data_matr[:, 1] .* 5.0 .+ data_matr[:, 2]) .^ 7.0 + exp.(data_matr[:, 1] .* 5.0 .+ data_matr[:, 2])
# -> 3 * (v1 * 5 + v2)^7 + exp(v1 * 5 + v2)

# make some custom settings
fit_weights = 1 ./ data_matr[:, end] # weights to minimize relative deviation
parts = [0.8, 0.2]

ops, data = Options(
    data_descript=data_descript(
        data_matr;
        parts          = parts,
        fit_weights    = fit_weights
    ),
    general=general_params(
        n_gens          = typemax(Int64),
        pop_size        = 200,
        max_compl       = 30,
        pow_abs_param   = true,
        remove_doubles_sigdigits = 1e-2,
        t_lim           = 60 * 2.0,                  # will run for 2 minutes
        multithreading = true,
    ),
    fitting=fitting_params(
        early_stop_iter = 5,
        max_iter        = 15,
    ),
);

# start the equation search
hall_of_fame, population, prog_dict = generational_loop(data, ops);

# inspect the results
col = "mare" # mean relative error
perm = sortperm(hall_of_fame[col])

hall_of_fame[col][perm]
hall_of_fame["compl"][perm]
hall_of_fame["node"][perm]

```

# Available settings and their default values

In the following code snippet, all available, exposed settings are shown with their default values.
A brief explanation is provided for most settings as a comment.
For a demo on how to use the code, please find the example_main.jl file.

```julia
ops, data = Options(
    data_descript = data_descript(
        data_matr;                                  # -> nxm matrix containing the n data points, m-1 variables and the output
        parts          = [1.0],                     # -> how to split the data. e.g. [1.0] -> (no split) or [0.8, 0.2]
        fit_weights    = ones(size(data_matr, 1)),  # -> weights for the data fitting -> residual .* weight
    ),

    general=general_params(
        n_gens                         = typemax(Int64), # -> number of generations to conduct
        t_lim                          = 60. * 5.,       # -> time limit for the algorithm
        pop_size                       = 100,            # -> number of individuals selected for next generation
        num_islands                    = 4,              # -> numer of parallel islands
        migration_interval             = 50,             # -> generation interval, in which an indivudual is moved to other islands. (ring topology)
        init_tree_depth                = 4,              # -> initial tree depth (with full method)
        max_compl                      = 50,             # -> max allowed complexity. Individuals exceeding it, are trimmed with repeated hoist mutations rather than removed.
        pow_abs_param                  = false,          # -> allows only parameter terminals as power -> (x - 5)^3 is allowed, but 3^(x - 5) is not.
        remove_doubles_sigdigits       = 2,              # -> filters indivuduals in an island if their rounded MAE and MSE are the same -> 0 is off; 2 ... 5
        remove_doubles_across_islands  = false,          # -> depends on the remove_doubles_sigdigits setting above. Applies it across islands  -> true / false
        multithreading                 = false           # -> whether to use multithreading for the most expensive computations. Not always faster  -> depends on how expensive fitting is for the problem at hand.
    ),

    selection=selection_params(
        hall_of_fame_objectives           = [:ms_processed_e, :compl, :mare],          # -> objectives for the hall_of_fame
        selection_objectives              = [:ms_processed_e, :compl, :age],           # -> objectives for the Pareto-optimal selection part of selection
        tournament_selection_fitness      = [(1.0, :ms_processed_e), (1e-5, :compl)],  # -> how to calculate the fitness for the tournament selection. Must be a vector as shown in the example. It is then calculated as a sum weighted be the specified scalars.
        ratio_pareto_tournament_selection = 0.7,                                       # -> ratio to which the selection should be using the Pareto-optimal selection vs. tournament selection
        tournament_size                   = 5,                                         # -> tournament size
    ),

    fitting=fitting_params(
        max_iter                 = 20,             # -> maximum iterations for parameter fitting (10 ... 100) ==> biggest time consumer <==
        early_stop_iter          = 0,              # -> how many iterations to account for early stopping regularization (0 -> turned off) (to use, the data needs to be partitioned into at least 2 parts. The early stopping evaluation is performed on the second partition.)
        t_lim                    = Inf,            # -> time limit for parameter fitting of individual (0.1 ... 5)
        rel_f_tol_5_iter         = 1e-2 * 0.01,    # -> relative tolerance for parameter fitting of individual -> considered converged if improvement over 5 iterations is smaller than this
        lasso_factor             = 1e-7,           # -> factor for the lasso regularization -> pushing parameter values to 0 (0 -> turned off) (1e-7 ... 1e-4)
        pre_residual_processing! = (x, ind) -> x,  # -> processing of the equation output before the residual is calculated. Must be an inplace function. The inds refer for the indices of the current residuals, which may be used to slice some data in the function like "(x, inds) -> x ./= data[end][inds]"
        residual_processing      = (x, ind) -> x,  # -> processing of the residuals. NOT an inplace function. The inds refer for the indices of the current residuals, which may be used to slice some data in the function like "(x, inds) -> x ./ data[end][inds]"
    ),

    mutation=mutation_params(         # -> probabilites for the various mutations (don't need to add up to 1)
        p_crossover        = 4.0,     #|
        p_point            = 0.5,     #|
        p_innergrow        = 0.0,     #|-> a new directed-acyclic-graph connection
        p_insert           = 0.2,     #|
        p_hoist            = 0.2,     #|
        p_subtree          = 0.2,     #|
        p_add_term         = 0.1,     #|
        p_simplify         = 0.5,     #|-> simplify with SymbolicUtils
        p_drastic_simplify = 0.5,     #|-> remove parameter nodes with values smaller than 1e-4
    ),

    binops         = (  +,   -,   *,   /,   ^),  # -> binary function set to choose from
    p_binops       = (1.0, 1.0, 1.0, 1.0, 1.0),  # -> probabilites for selection of each binary functions (same length as provided binops) (dont need to add up to 1, adjusted accordingly)
    unaops         = (exp, log, sin, cos, abs),  # -> unary function set to choose from
    p_unaops       = (1.0, 1.0, 1.0, 1.0, 1.0),  # -> probabilites for selection of each unary functions (same length as provided unaops) (dont need to add up to 1, adjusted accordingly)

    illegal_dict = Dict(), # -> May be set to a dict to specify nestings of unary functions that are removed, if they appear.
                           #    For example: illegal_dict = Dict(:sin => (sin, cos),
);                         #                                     :cos => (sin, cos))
```

# Funding

This work was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) within the Priority Programme SPP 2331: "Machine Learning in Chemical Engineering" - project no. 466528284 - HE 6077/14-1 and RI 2482/10-1

