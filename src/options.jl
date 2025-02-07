
""" The Options struct contains all settings. For most of its fields, there are functions with
    default arguments to make sure all fields are filled and potentially adapted to make sense.
    There are also some safeguards, preventing certain settings combinations or warn the user,
    which may not cover all cases.
"""
struct Options{A, B, C, D, F, G, H, I, J, K}
    general::A
    unaops::B
    binops::C
    measures::D

    selection::F
    fitting::G
    mutation::H
    grammar::I
    data_descript::J
    meta_data::K

    function Options(
        data::Matrix;                                             # -> nxm matrix containing the n data points, m-1 variables and the output
        fit_weights::Vector{Float64} = abs.(1 ./ data[:, end]),   # -> weights for the data fitting -> residual .* weight
        binops                       = (+,   -,   *,   /,   ^  ), # -> binary function set to choose from
        unaops                       = (exp, log, sin, cos, abs), # -> unary function set to choose from
        data_split                   = data_split_params(),
        general                      = general_params(),
        measures                     = measure_params(),
        selection                    = selection_params(),
        fitting                      = fitting_params(),
        mutation                     = mutation_params(),
        grammar                      = grammar_params(),
        meta_data                    = Dict(),                  # -> can be used to provide data for use in, for example, user-defined functions like the callback, pre_residual_processing, or measures
    )

        # make sure all required measres are calculated # ------------------------------------------
        @assert all(s in keys(measures) for s in selection.selection_objectives)    "some selection_objectives are not specified in the meausres"
        @assert all(s in keys(measures) for s in selection.hall_of_fame_objectives) "some hall_of_fame_objectives are not specified in the meausres"

        # prepare data and its params --------------------------------------------------------------
        data_vect, data_descript = prepare_data(data, fit_weights)

        split_inds = _data_split_params(data, data_split...)
        data_descript = (data_descript..., split_inds = split_inds)

        # ------------------------------------------------------------------------------------------
        @assert length(unaops) <= 128 "TiSR cannot handle more than 128 functions in the unaops out of the box"
        @assert length(binops) <= 128 "TiSR cannot handle more than 128 functions in the binops out of the box"

        # function set asserts ---------------------------------------------------------------------
        @assert ("+" in string.(binops)) "+ is required in the function set for some mutations to work"
        @assert ("*" in string.(binops)) "* is required in the function set for some mutations to work"
        :^ in Symbol.(binops) && @warn   "^ is only valid for positive bases. Otherwise provide a pow(abs(x), y) function with its own drawbacks."

        #-------------------------------------------------------------------------------------------
        if fitting.early_stop_iter != 0
            @assert !isempty(data_descript.split_inds[2]) "second data split must contain data for early stopping"
        end

        if !isempty(grammar.weighted_compl_dict)
            :weighted_compl in selection.selection_objectives    || @warn "weighted_compl_dict specified, but :weighted_compl not in selection_objectives"
            :weighted_compl in selection.hall_of_fame_objectives || @warn "weighted_compl_dict specified, but :weighted_compl not in hall_of_fame_objectives"
        end

        # ------------------------------------------------------------------------------------------
        @assert !xor(isempty(grammar.bank_of_terms), mutation.mut_probs[9] == 0) "for p_add_from_bank_of_terms != 0, bank_of_terms cannot be empty and vv"

        # resulting parameters
        n_pareto_select_per_isle = ceil(Int64, general.pop_per_isle * selection.ratio_pareto_tournament_selection)
        selection = (;selection..., n_pareto_select_per_isle = n_pareto_select_per_isle)

        global __operators = (unaops = unaops, binops = binops) # required for Base.show() Nodes

        return new{typeof(general),
                   typeof(unaops),
                   typeof(binops),
                   typeof(measures),
                   typeof(selection),
                   typeof(fitting),
                   typeof(mutation),
                   typeof(grammar),
                   typeof(data_descript),
                   typeof(meta_data)}(
                       general,
                       unaops,
                       binops,
                       measures,
                       selection,
                       fitting,
                       mutation,
                       grammar,
                       data_descript,
                       meta_data
        ), data_vect
    end
end

# ==================================================================================================
# User facing function with default parameters, precalculatons, and checks
# ==================================================================================================
""" Function to specify either parts for data splitting or split_inds directy.
"""
data_split_params(;
    parts      = nothing, # -> how the data should be split. The first part is used for fitting, the second is used for early stopping, while the first AND the second part is used for the selection measures. The third part is used for the test-selection measures. The third part is used for the test-selection measures.
    split_inds = nothing, # -> rather then splitting randomly and automatically, the indices for each split can be specified in a Vector{Vector{Int64}}. If specified, parts cannot be specified.
) = (parts, split_inds)
_data_split_params(data, parts::Nothing, split_inds::Nothing) = _data_split_params(data, nothing, [collect(1:size(data, 1))])
_data_split_params(data, parts, split_inds) = throw("parts and split_inds cannot both be specified")

function _data_split_params(data, parts::Vector{Float64}, split_inds::Nothing)

    # parts assertions
    @assert length(parts) > 0  "parts must have at least one float entry"
    @assert sum(parts) == 1.0  "sum of the parts should be 1.0"
    @assert all(parts .>= 0)   "data split parts have to be >= 0"
    length(parts) < 4 || @warn "data split seems off. 1. part is used for fitting, 2. part is used for early stopping, while 1&2 are used for selection. All further parts are not used, but can be useful, like for a user defined callback."

    # create the split inds according to specified
    eachind = collect(1:size(data, 1))
    shuffle!(eachind)

    start_inds = cumsum([ceil(Int64, size(data, 1) * p) for p in parts])
    pushfirst!(start_inds, 0)

    split_inds = [eachind[start_inds[i]+1:min(start_inds[i+1], size(data, 1))] for i in 1:length(parts)]
    sort!.(split_inds)

    return _data_split_params(data, nothing, split_inds)
end

function _data_split_params(data, parts::Nothing, split_inds::Vector{Vector{Int64}})
    eachind = reduce(vcat, split_inds)

    allunique(eachind)                       || @warn "in split_inds, some data are part of different splits"
    length(unique(eachind)) == size(data, 1) || @warn "some data are not part of any split"

    while length(split_inds) < 3
        push!(split_inds, Int64[])
    end

    push!(split_inds, sort(vcat(split_inds[1], split_inds[2])))

    return split_inds
end

function prepare_data(data, fit_weights)
    data_vect = [data[:, i] for i in 1:size(data, 2)]

    # prepare data # ---------------------------------------------------------------------------
    @assert eltype(data) <: AbstractFloat "data must be float type "
    @assert size(data, 2) > 1             "data has only one column"
    @assert size(data, 2) <= 128 + 1      "TiSR cannot handle more than 128 variables out of the box"
    !any(length(unique(d)) == 1 for d in data_vect) || @warn "data containts columns, which are constant"

    @assert length(fit_weights) == size(data, 1) "fit_weights seems to have too many or to few entries"
    @assert all(fit_weights .>= 0.0)             "fit_weights must be larger than 0"
    @assert all(isfinite, fit_weights)           "some non finite values in fit_weights"

    return data_vect, (
        data_type   = eltype(data_vect[1]),
        n_vars      = length(data_vect) - 1,
        n_points    = length(data_vect[1]),
        fit_weights = fit_weights,
    )
end

""" Returns a NamedTuple of the general parameters to be appeded to the Options.
    The default parameters can be adapted according the requirements. Some resulting
    parameters are also calculated and included.
"""
function general_params(;
    n_gens::Int64                          = typemax(Int64),                                           # -> number of generations to conduct
    t_lim::Float64                         = 60. * 5.,                                                 # -> time limit for the algorithm
    pop_size::Int64                        = 600,                                                      # -> number of individuals selected for next generation / population size
    parent_selection::Bool                 = true,                                                     # -> wheather to conduct parent selection or to use all individuals in the population as parents
    num_islands::Int64                     = 12,                                                       # -> number of parallel islands
    children_ratio::Float64                = 0.5,                                                      # -> the ratio of children that should be generated in each generation 0 ... 2
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
)
    @assert num_islands > 0                             "num_islands should be at least 1                          "
    @assert migration_interval > 0                      "migration_interval should be at least 1                   "
    @assert always_drastic_simplify >= 0                "always_drastic_simplify must be >= 0                      "
    @assert adaptive_compl_increment > 0                "adaptive_compl_increment must be larger 0                 "
    @assert island_extinction_interval > 0              "island_extinction_interval must be > 0                    "
    @assert max_age > 1                                 "max_age must be > 1                                       "
    @assert hall_of_fame_migration_interval > 0         "hall_of_fame_migration_interval must be larger 0          "
    @assert 0.0 <= migrate_after_extinction_prob <= 1.0 "migrate_after_extinction_prob must be between 0 and 1     "
    @assert 0.0 <= children_ratio <= 2.0                "children_ratio should be inbetween 0.0 and 2.0            "
    @assert 0 <= n_refitting <= pop_size / num_islands  "n_refitting should be between 0 and pop_size / num_islands"

    if print_hall_of_fame
        @assert plot_hall_of_fame "for print_hall_of_fame, plot_hall_of_fame must be true"
    end

    remove_doubles_sigdigits > 1     || @warn "a low remove_doubles_sigdigits may filter non-equal individuals "
    remove_doubles_sigdigits < 6     || @warn "a low remove_doubles_sigdigits may not detect equal individuals "
    always_drastic_simplify < 1e-3   || @warn "always_drastic_simplify seems high                              "
    adaptive_compl_increment > 4     || @warn "adaptive_compl_increment should be >= 5                         "
    island_extinction_interval > 500 || @warn "island_extinction_interval seems small                          "
    max_age > 10                     || @warn "max_age seems small                                             "

    island_extinction_interval == migration_interval == typemax(Int64) && @warn "island_extinction_interval & migration_interval should not both be off"

    if multithreading
        Threads.nthreads() > 1 || @warn "To acually utilize multithreading, Julia must be started with the desired number of threads, i.e., `julia -t 4`"
    end

    # resulting parameters
    pop_per_isle = ceil(Int64, pop_size / num_islands)
    n_children = max(2, round(Int64, children_ratio * pop_per_isle))

    return (
        n_gens                          = n_gens,
        pop_size                        = pop_size,
        parent_selection                = parent_selection,
        n_children                      = n_children,
        pop_per_isle                    = pop_per_isle,
        num_islands                     = num_islands,
        migration_interval              = migration_interval,
        hall_of_fame_migration_interval = hall_of_fame_migration_interval,
        island_extinction_interval      = island_extinction_interval,
        migrate_after_extinction_prob   = migrate_after_extinction_prob,
        fitting_island_function         = fitting_island_function,
        always_drastic_simplify         = always_drastic_simplify,
        remove_doubles_sigdigits        = remove_doubles_sigdigits,
        remove_doubles_across_islands   = remove_doubles_across_islands,
        max_age                         = max_age,
        n_refitting                     = n_refitting,
        t_lim                           = t_lim,
        multithreading                  = multithreading,
        adaptive_compl_increment        = adaptive_compl_increment,
        callback                        = callback,
        print_progress                  = print_progress,
        plot_hall_of_fame               = plot_hall_of_fame,
        print_hall_of_fame              = print_hall_of_fame,
    )
end

""" Specify the fit quality or complexity measures that should be calculated for
    all individuals. Those will be tracked, printed, logged, and can be used as
    selection_objetives or hall_of_fame_objectives. The measures "ms_processed_e",
    "compl", "mse", and "mae" are always included. The following functions are
    provided out of the box: ... TODO User-specified should should take 6
    positional arguments: residual, residual_relative, prediction, data, node, ops.
"""
function measure_params(;
    additional_measures::Dict{Symbol, Function} = Dict(                 # -> specify a Dict{Symbol, Function} containing name and function pairs that calculate custom measures. TiSR offers some additional ones, which all start with `TiSR.get_measure_...`
        :one_minus_abs_spearman => get_measure_one_minus_abs_spearman,
        :mare                   => get_measure_mare,
        :max_are                => get_measure_max_are,
    )
)
    :ms_processed_e in keys(additional_measures) && @warn "ms_processed_e is overwritten in measures"
    :compl          in keys(additional_measures) && @warn "compl is overwritten in measures         "
    :mse            in keys(additional_measures) && @warn "mse is overwritten in measures           "
    :mae            in keys(additional_measures) && @warn "mae is overwritten in measures           "

    additional_measures[:ms_processed_e] = get_measure_ms_processed_e
    additional_measures[:compl]          = get_measure_compl
    additional_measures[:mse]            = get_measure_mse
    additional_measures[:mae]            = get_measure_mae

    return additional_measures
end

function selection_params(;
    hall_of_fame_objectives::Vector{Symbol}    = [:ms_processed_e, :compl],                          # -> objectives for the hall_of_fame
    selection_objectives::Vector{Symbol}       = [:ms_processed_e, :one_minus_abs_spearman, :compl], # -> objectives for the Pareto-optimal selection part of selection
    normalize_objectives::Bool                 = false,                                              # -> whether to normalize the objectives for niching and crowding distance
    hall_of_fame_niching_sigdigits::Int64      = 2,                                                  # -> number of significant digits to round hall_of_fame_objectives for hall_of_fame selection after their normalization. -> 2 ... 5
    population_niching_sigdigits::Int64        = 3,                                                  # -> number of significant digits to round selection_objectives for population selection after their normalization. -> 2 ... 5
    ratio_pareto_tournament_selection::Float64 = 0.5,                                                # -> ratio to which the selection is conducted using the Pareto-optimal selection vs. tournament selection
    tournament_size::Int64                     = 5,                                                  # -> tournament size
)
    @assert tournament_size > 1                             "tournament size must be greater than 1"
    @assert 0.0 <= ratio_pareto_tournament_selection <= 1.0 "ratio_pareto_tournament_selection must be between 0.0 and 1.0"

    @assert hall_of_fame_niching_sigdigits > 0 "hall_of_fame_niching_sigdigits must be larger than 0"
    @assert population_niching_sigdigits   > 0 "population_niching_sigdigits must be larger than 0"

    0 < hall_of_fame_niching_sigdigits < 5 || @warn "hall_of_fame_niching_sigdigits should be between 0 and 5"
    0 < population_niching_sigdigits   < 5 || @warn "population_niching_sigdigits should be between 0 and 5"

    return (
        hall_of_fame_objectives           = hall_of_fame_objectives,
        selection_objectives              = selection_objectives,
        normalize_objectives              = normalize_objectives,
        hall_of_fame_niching_sigdigits    = hall_of_fame_niching_sigdigits,
        population_niching_sigdigits      = population_niching_sigdigits,
        tournament_size                   = tournament_size,
        ratio_pareto_tournament_selection = ratio_pareto_tournament_selection
    )
end

function fitting_params(;
    max_iter::Int64                   = 10,                 # -> maximum iterations for parameter fitting. -> 10 ... 50 ==> biggest time consumer <==
    early_stop_iter::Int64            = 0,                  # -> how many iterations to account for early stopping regularization. to use, the data needs to be partitioned into at least 2 parts. The early stopping evaluation is performed on the second partition. -> 0 is off; 4 ... 10
    t_lim::Float64                    = Inf,                # -> time limit for parameter fitting of individual. -> Inf is off; 0.1 ... 0.5
    rel_f_tol_5_iter::Float64         = 1e-2 * 0.01,        # -> relative tolerance for parameter fitting. considered converged if relative improvement over 5 iterations is smaller. -> 0 is off; 1e-2 * 1.0 ... 1e-2 * 0.01
    lasso_factor::Float64             = 0.0,                # -> factor for the lasso regularization. pushing parameter values to 0. -> 0 is off; 1e-8 ... 1e-4
    pre_residual_processing::Function = (x, ind, ops) -> x, # -> processing of the equation output before the residual is calculated. The inds refer for the indices of the current residuals, which may be used to slice some data in the function like "(x, inds, ops) -> x ./= data[end][inds]"
    residual_processing::Function     = (x, ind, ops) -> x, # -> processing of the residuals. NOT an inplace function. The inds refer for the indices of the current residuals, which may be used to slice some data in the function like "(x, inds, ops) -> x ./ data[end][inds]"
)
    t_lim > 1e-1             || @warn "fitting t_lim may be too low"
    max_iter >= 5            || @warn "max_iter may be too low"
    lasso_factor < 1.0       || @warn "lasso_factor seems to large"
    0 < early_stop_iter < 5  && @warn "early stopping may be too strict -> higher values may produce better results"

    @assert max_iter >= early_stop_iter "early_stop_iter should be smaller than max_iter"
    @assert 0 <= rel_f_tol_5_iter < 1.0 "rel_f_tol_5_iter must smaller than 1.0 and larger or equal to 0"
    @assert lasso_factor >= 0           "lasso factor must be >= 0"

    return (
        max_iter                = max_iter,
        early_stop_iter         = early_stop_iter,
        rel_f_tol_5_iter        = rel_f_tol_5_iter,
        t_lim                   = t_lim,
        lasso_factor            = lasso_factor,
        pre_residual_processing = pre_residual_processing,
        residual_processing     = residual_processing,
    )
end

""" Returns the equation grammar related parameters.
"""
function grammar_params(;
    max_compl::Int64                           = 30,                        # -> max allowed complexity
    min_compl::Int64                           = 3,                         # -> min allowed complexity. -> 2 ... 3
    init_tree_depth::Int64                     = 4,                         # -> maximal initial tree depth. -> 3 ... 6
    max_nodes_per_term::Int64                  = typemax(Int64),            # -> maximal number of nodes per top-level term. All terms above are trimmed until they satisfy this threshold
    weighted_compl_dict::Dict{String, Float64} = Dict{String, Float64}(),   # -> weights for weighted_compl calculation. For any that are not included, 1.0 is assumed. The weights for variables and parameters "VAR" and "PARAM" may be used. An example is shown below at (2).
    bank_of_terms::Vector{String}              = String[],                  # -> specify terms that can be added via the add term mutation, whose relativ probability is set via the p_add_from_bank_of_terms parameter
    illegal_dict::Dict                         = Dict(),                    # -> check for illegal nestings in existing nodes. For it, ops.illegal_dict needs to be specified like below. Variables can be specified using "VAR" and parameters with "PARAM". An example is shown below at (1).
    custom_check_legal::Function               = (node, data, ops) -> true, # -> specify a custom function, which checks the legality of nodes. Must return true or false
)
    @assert max_compl > 3             "max_compl must be larger than 3          "
    @assert 0 < min_compl < max_compl "min_compl must be between 1 and max_compl"
    @assert init_tree_depth > 2       "init_tree_depth should be 3 or higher    "
    @assert max_nodes_per_term > 1    "max_nodes_per_term must be larger than 1 "

    if !isempty(illegal_dict)
        @assert illegal_dict isa Dict                  "illegal_dict is not formatted correctly"
        for (k, v) in illegal_dict
            @assert k isa String                       "illegal_dict is not formatted correctly"
            @assert haskey(v, :lef) && haskey(v, :rig) "illegal_dict is not formatted correctly"
            for v_ in v.lef
                @assert v_ isa String                  "illegal_dict is not formatted correctly"
            end
            for v_ in v.rig
                @assert v_ isa String                  "illegal_dict is not formatted correctly"
            end
        end
    end

    @assert weighted_compl_dict isa Dict{String, Float64} "weighted_compl_dict is not formatted correctly"

    max_compl > 100     && @warn "a high max_compl may lead to high calculation times"
    init_tree_depth > 6 && @warn "a high init_tree_depth may lead to high calculation times"

    return (
        illegal_dict        = illegal_dict,
        weighted_compl_dict = weighted_compl_dict,
        init_tree_depth     = init_tree_depth,
        max_compl           = max_compl,
        min_compl           = min_compl,
        max_nodes_per_term  = max_nodes_per_term,
        bank_of_terms       = bank_of_terms,
        custom_check_legal  = custom_check_legal,
    )
end

""" Returns a Tuple of the normalized cumulative sum of probabilities of the various
    mutations to be appended to the Options. The input values don't need to add
    up to 1, since they are normalized here.
"""
function mutation_params(;                     #|-> probabilites for the various mutations (don't need to add up to 1)
    p_crossover::Float64              = 10.0,  #|
    p_point::Float64                  = 1.0,   #|
    p_insert::Float64                 = 1.0,   #|
    p_hoist::Float64                  = 1.0,   #|
    p_subtree::Float64                = 1.0,   #|
    p_drastic_simplify::Float64       = 0.1,   #|-> remove parameter nodes with small values and simplify accordingly
    p_insert_times_param::Float64     = 0.1,   #|
    p_add_term::Float64               = 0.1,   #|
    p_simplify::Float64               = 0.1,   #|-> simplify with SymbolicUtils
    p_add_from_bank_of_terms::Float64 = 0.0,   #|-> probability to add a term from the provided bank_of_terms
    p_multiple_mutations::Float64     = 0.5,   # -> probability for more than one mutation
)
    @assert all(p >= 0 for p in (
        p_crossover, p_point, p_insert, p_hoist, p_subtree, p_drastic_simplify,
        p_insert_times_param, p_add_term, p_simplify, p_add_from_bank_of_terms
    )) "all mutation probabilities must be >= 0"

    p_crossover *= 0.5 # because it consumes a second one, if hit

    @assert p_multiple_mutations < 1         "p_multiple_mutations must be < 1"
    0 <= p_multiple_mutations < 0.9 || @warn "p_multiple_mutations should be between 0 and 0.9"

    # TODO: reduce remainder by p_multiple_mutations

    mut_probs = (
        p_insert,
        p_point,
        p_add_term,
        p_insert_times_param,
        p_hoist,
        p_subtree,
        p_drastic_simplify,
        p_simplify,
        p_add_from_bank_of_terms,
        p_crossover,
    )

    mut_probs = mut_probs ./ sum(mut_probs)
    return (
        mut_probs = mut_probs,
        p_multiple_mutations = p_multiple_mutations,
    )
end
