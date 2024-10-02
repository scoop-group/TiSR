
# ==================================================================================================
# Options struct
# ==================================================================================================
""" The Options struct contains all settings. For most of its fields, there are functions with
    default arguments to make sure all fields are filled and potentially adapted to make sense.
    There are also some safeguards, preventing certain settings combinations or warn the user,
    which may not cover all cases.
"""
struct Options{A, B, C, D, E, F, G, H, I, J}
    general::A
    unaops::B
    binops::C
    p_unaops::D
    p_binops::E
    selection::F
    fitting::G
    mutation::H
    grammar::I
    data_descript::J

    function Options(
        data::Matrix;
        fit_weights::Vector   = ones(size(data, 1)),
        parts::Vector         = [1.0],
        general               = general_params(),
        p_unaops              = (1.0, 1.0, 1.0, 1.0, 1.0),
        unaops                = (exp, log, sin, cos, abs),
        p_binops              = (1.0, 1.0, 1.0, 1.0, 1.0),
        binops                = (+,   -,   *,   /,   ^  ),
        selection             = selection_params(),
        fitting               = fitting_params(),
        mutation              = mutation_params(),
        grammar               = grammar_params(),
    )
        # convert to vectors of vectors
        data_vect = [data[:, i] for i in 1:size(data, 2)]

        # prepare data # ---------------------------------------------------------------------------
        @assert eltype(data) <: AbstractFloat "data must be float type "
        @assert size(data, 2) > 1             "data has only one column"
        !any(length(unique(d)) == 1 for d in data_vect) || @warn "data containts rows, which are constant"

        @assert length(fit_weights) == size(data, 1) "fit_weights seems to have too many or to few entries"
        @assert all(fit_weights .>= 0.0)             "fit_weights must be larger than 0"

        # data splitting # -------------------------------------------------------------------------
        @assert length(parts) > 0  "parts must have at least one float entry"
        @assert sum(parts) == 1.0  "sum of the parts should be 1.0"
        @assert all(parts .>= 0)   "data split parts have to be >= 0"
        length(parts) < 3 || @warn "only two of the data split parts are used, while all data is used to calculate the final measures" # TODO: use 3rd place to calculate test measures

        # split inds # -----------------------------------------------------------------------------
        eachind = collect(1:size(data, 1))
        shuffle!(eachind) # TODO: don't shuffle if not needed?

        start_inds = cumsum([ceil(Int64, size(data_vect, 1) * p) for p in parts])
        pushfirst!(start_inds, 1)
        split_inds = [eachind[start_inds[i]:start_inds[i+1]] for i in 1:length(parts)]

        data_descript = (
            data_type             = eltype(data_vect[1]),
            n_vars                = length(data_vect) - 1,
            n_points              = length(data_vect[1]),
            split_inds            = split_inds,
            fit_weights           = fit_weights,
        )

        @assert data_descript.n_vars <= 128 "TiSR cannot handle more than 128 variables out of the box"
        @assert length(unaops) <= 128       "TiSR cannot handle more than 128 functions in the unaops out of the box"
        @assert length(binops) <= 128       "TiSR cannot handle more than 128 functions in the binops out of the box"

        # function set asserts ---------------------------------------------------------------------
        @assert ("+" in string.(binops)) "+ is required in the function set for some mutations to work"
        @assert ("*" in string.(binops)) "* is required in the function set for some mutations to work"

        #-------------------------------------------------------------------------------------------
        @assert length(unaops) == length(p_unaops) "unaops tuple and its wights must be the same length"
        @assert length(binops) == length(p_binops) "binops tuple and its wights must be the same length"

        :^ in Symbol.(binops) && @warn "^ is only valid for positive bases. Otherwise provide a pow(abs(x), y) function with its own drawbacks."

        if fitting.early_stop_iter != 0
            @assert length(data_descript.split_inds) > 1 "Data split required for early stopping"
        end

        # relative reference offset
        any(abs(d) < 0.1 for d in data[end]) && @warn "some target data < 0.1 -> 0.1 is used as lower bound for the relative measures like mare"

        # resulting parameters
        n_pareto_select_per_isle = ceil(Int64, general.pop_per_isle * selection.ratio_pareto_tournament_selection)
        selection = (;selection..., n_pareto_select_per_isle = n_pareto_select_per_isle)

        global __operators = (unaops = unaops, binops = binops) # required for Base.show() Nodes

        return new{typeof(general),
                   typeof(unaops),
                   typeof(binops),
                   typeof(p_unaops),
                   typeof(p_binops),
                   typeof(selection),
                   typeof(fitting),
                   typeof(mutation),
                   typeof(grammar),
                   typeof(data_descript)}(
                       general,
                       unaops,
                       binops,
                       p_unaops,
                       p_binops,
                       selection,
                       fitting,
                       mutation,
                       grammar,
                       data_descript
        ), data_vect
    end
end

# ==================================================================================================
# User facing function with default parameters, precalculatons, and checks
# ==================================================================================================
""" Returns a NamedTuple of the general parameters to be appeded to the Options.
    The default parameters can be adapted according the requirements. Some resulting
    parameters are also calculated and included.
"""
function general_params(;
    n_gens                        = typemax(Int64),
    t_lim                         = 60. * 5.,
    pop_size                      = 500,
    num_islands                   = 10,
    migration_interval            = 30,
    always_drastic_simplify       = 1e-7,
    remove_doubles_sigdigits      = 3,
    remove_doubles_across_islands = false,
    multithreading               = false,
    adaptive_compl_increment      = Inf,
    callback                      = (hall_of_fame, population, ops) -> false,
    print_progress                = true,
    plot_hall_of_fame             = true,
)
    @assert num_islands > 0              "num_islands should be at least 1                "
    @assert migration_interval > 0       "migration_interval should be at least 1         "
    @assert always_drastic_simplify >= 0 "always_drastic_simplify must be >= 0            "
    @assert adaptive_compl_increment > 0 "adaptive_compl_increment must be larger 0       "
    @assert callback isa Function        "callback must be a function                     "

    remove_doubles_sigdigits > 1   || @warn "a low remove_doubles_sigdigits may filter non-equal individuals "
    remove_doubles_sigdigits < 6   || @warn "a low remove_doubles_sigdigits may not detect equal individuals "
    always_drastic_simplify < 1e-3 || @warn "always_drastic_simplify seems high                              "
    adaptive_compl_increment > 4   || @warn "adaptive_compl_increment should be >= 5                  "

    if multithreading
        Threads.nthreads() > 1 || @warn "To acually utilize multithreading, Julia must be started with the desired number of threads, i.e., `julia -t 4`" 
    end

    # resulting parameters
    pop_per_isle = ceil(Int64, pop_size / num_islands)

    return (
        n_gens                        = n_gens,
        pop_size                      = pop_size,
        pop_per_isle                  = pop_per_isle,
        num_islands                   = num_islands,
        migration_interval            = migration_interval,
        always_drastic_simplify       = always_drastic_simplify,
        remove_doubles_sigdigits      = remove_doubles_sigdigits,
        remove_doubles_across_islands = remove_doubles_across_islands,
        t_lim                         = t_lim,
        multithreading                = multithreading,
        adaptive_compl_increment      = adaptive_compl_increment,
        callback                      = callback,
        print_progress                = print_progress,
        plot_hall_of_fame             = plot_hall_of_fame,
    )
end

function selection_params(;
    hall_of_fame_objectives           = [:ms_processed_e, :compl, :mare],
    selection_objectives              = [:ms_processed_e, :minus_abs_spearman, :compl, :age],
    hall_of_fame_niching_sigdigits    = 2,
    population_niching_sigdigits      = 3,
    tournament_selection_fitness      = [(1.0, :ms_processed_e), (1e-5, :compl)],
    ratio_pareto_tournament_selection = 0.7,
    tournament_size                   = 5,
)
    @assert tournament_size > 1                             "tournament size must be greater than 1"
    @assert 0.0 <= ratio_pareto_tournament_selection <= 1.0 "ratio_pareto_tournament_selection must be between 0.0 and 1.0"
    @assert typeof(tournament_selection_fitness) == Vector{Tuple{Float64, Symbol}}

    @assert hall_of_fame_niching_sigdigits > 0 "hall_of_fame_niching_sigdigits must be larger than 0"
    @assert population_niching_sigdigits   > 0 "population_niching_sigdigits must be larger than 0"

    0 < hall_of_fame_niching_sigdigits < 5 || @warn "hall_of_fame_niching_sigdigits should be between 0 and 5"
    0 < population_niching_sigdigits    < 5 || @warn "population_niching_sigdigits should be between 0 and 5"

    return (
        hall_of_fame_objectives           = hall_of_fame_objectives,
        selection_objectives              = selection_objectives,
        hall_of_fame_niching_sigdigits    = hall_of_fame_niching_sigdigits,
        population_niching_sigdigits      = population_niching_sigdigits,
        tournament_selection_fitness      = tournament_selection_fitness,
        tournament_size                   = tournament_size,
        ratio_pareto_tournament_selection = ratio_pareto_tournament_selection
    )
end

function fitting_params(;
    max_iter                 = 20,
    early_stop_iter          = 0,
    t_lim                    = Inf,
    rel_f_tol_5_iter         = 1e-2 * 0.01,
    lasso_factor             = 1e-7,
    pre_residual_processing! = (x, ind) -> x,
    residual_processing      = (x, ind) -> x,
)
    t_lim > 1e-1             || @warn "fitting t_lim may be too low"
    max_iter >= 5            || @warn "max_iter may be too low"
    lasso_factor < 1.0       || @warn "lasso_factor seems to large"
    0 < early_stop_iter < 5  && @warn "early stopping may be too strict -> higher values may produce better results"

    @assert max_iter >= early_stop_iter "early_stop_iter should be smaller than max_iter"
    @assert 0 <= rel_f_tol_5_iter < 1.0 "rel_f_tol_5_iter must smaller than 1.0 and larger or equal to 0"
    @assert lasso_factor >= 0           "lasso factor must be >= 0"

    return (
        max_iter                 = max_iter,
        early_stop_iter          = early_stop_iter,
        rel_f_tol_5_iter         = rel_f_tol_5_iter,
        t_lim                    = t_lim,
        lasso_factor             = lasso_factor,
        pre_residual_processing! = pre_residual_processing!,
        residual_processing      = residual_processing,
    )
end

""" Returns the equation grammar related parameters.
"""
function grammar_params(;
    illegal_dict       = Dict(),
    max_compl          = 30,
    min_compl          = 2,
    max_nodes_per_term = Inf,
    init_tree_depth    = 4,
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

    max_compl > 100         && @warn "a high max_compl may lead to high calculation times"
    init_tree_depth > 6     && @warn "a high init_tree_depth may lead to high calculation times"

    return (
        illegal_dict       = illegal_dict,
        init_tree_depth    = init_tree_depth,
        max_compl          = max_compl,
        min_compl          = min_compl,
        max_nodes_per_term = max_nodes_per_term,
    )
end

""" Returns a Tuple of the normalized cumulative sum of probabilities of the various
    mutations to be appended to the Options. The input values don't need to add
    up to 1, since they are normalized here.
"""
function mutation_params(;
    p_crossover          = 5.0,
    p_point              = 0.5,
    p_insert             = 0.2,
    p_hoist              = 0.2,
    p_subtree            = 0.2,
    p_drastic_simplify   = 0.2,
    p_insert_times_param = 0.1,
    p_add_term           = 0.1,
    p_simplify           = 0.1,
)
    p_crossover *= 0.5 # because it consumes a second one, if hit

    mut_params = (
        p_insert,
        p_point,
        p_add_term,
        p_insert_times_param,
        p_hoist,
        p_subtree,
        p_drastic_simplify,
        p_simplify,
        p_crossover,
    )

    return cumsum(mut_params ./ sum(mut_params)) # TODO: change to do that where its needed
end

