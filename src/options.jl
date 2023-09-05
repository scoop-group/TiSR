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
    illegal_dict::I
    data_descript::J

    function Options(;
        general       = general_params(),
        p_unaops      = (1.0, 1.0, 1.0, 0.0, 0.0),
        unaops        = (exp, log, sin, cos, abs),
        p_binops      = (1.0, 1.0, 1.0, 1.0, 1.0),
        binops        = (+,   -,   *,   /,   ^  ),
        selection     = selection_params(),
        fitting       = fitting_params(),
        mutation      = mutation_params(),
        illegal_dict  = Dict(),
        data_descript = data_descript(),
    )

        data_descript, data = data_descript

        @assert length(unaops) == length(p_unaops) "unaops tuple and its wights must be the same length"
        @assert length(binops) == length(p_binops) "binops tuple and its wights must be the same length"

        :^ in Symbol.(binops) && @warn "^ is only valid for positive bases. Otherwise provide a pow(abs(x), y) function with its own drawbacks."

        if fitting.early_stop_iter != 0
            @assert length(data_descript.split_inds) > 1 "Data split required for early stopping"
        end

        # resulting parameters
        n_pareto_select_per_isle = ceil(Int64, general.pop_per_isle * selection.ratio_pareto_tournament_selection)
        selection = (;selection..., n_pareto_select_per_isle = n_pareto_select_per_isle)

        global __operators = (unaops = unaops, binops = binops) # required to Base.show() Nodes

        return new{typeof(general),
                   typeof(unaops),
                   typeof(binops),
                   typeof(p_unaops),
                   typeof(p_binops),
                   typeof(selection),
                   typeof(fitting),
                   typeof(mutation),
                   typeof(illegal_dict),
                   typeof(data_descript)}(
                       general,
                       unaops,
                       binops,
                       p_unaops,
                       p_binops,
                       selection,
                       fitting,
                       mutation,
                       illegal_dict,
                       data_descript
                   ), data
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
    n_gens                         = typemax(Int64),
    t_lim                          = 60. * 5.,
    pop_size                       = 100,
    num_islands                    = 4,
    migration_interval             = 50,
    n_migrations                   = 5,
    init_tree_depth                = 4,
    max_compl                      = 50,
    pow_abs_param                  = false,
    prevent_doubles                = 1e-5,
    prevent_doubles_across_islands = false,
    multithreadding                = false
)
    @assert num_islands >= 1        "num_islands should be at least 1"
    @assert migration_interval >= 1 "migration_interval should be at least 1"
    @assert n_migrations >= 1       "num_islands should be at least 1"
    @assert init_tree_depth > 2     "init_tree_depth should be 3 or higher"

    init_tree_depth > 6     && @warn "a high init_tree_depth may lead to high calculation times"
    prevent_doubles > 1e-3  && @warn "a high prevent_doubles may filter non-equal individuals"
    prevent_doubles < 1e-14 && @warn "a low prevent_doubles may not detect equal individuals"
    max_compl > 100         && @warn "a high max_compl may lead to high calculation times"

    # resulting parameters
    pop_per_isle = ceil(Int64, pop_size / num_islands)

    nam_tup = (
        pop_size                       = pop_size,
        pop_per_isle                   = pop_per_isle,
        n_gens                         = n_gens,
        num_islands                    = num_islands,
        migration_interval             = migration_interval,
        n_migrations                   = n_migrations,
        init_tree_depth                = init_tree_depth,
        max_compl                      = max_compl,
        pow_abs_param                  = pow_abs_param,
        prevent_doubles                = prevent_doubles,
        prevent_doubles_across_islands = prevent_doubles_across_islands,
        t_lim                          = t_lim,
        multihreadding                 = multithreadding,
    )
    return nam_tup
end

function selection_params(;
    hall_of_fame_objectives           = [:ms_processed_e, :compl, :mare],
    selection_objectives              = [:ms_processed_e, :compl, :age],
    tournament_selection_fitness      = [(1.0, :ms_processed_e), (1e-5, :compl)],
    ratio_pareto_tournament_selection = 0.7,
    tournament_size                   = 5,
)

    @assert 2   <= tournament_size                          "tournament size must be greater than 1"
    @assert 0.0 <= ratio_pareto_tournament_selection <= 1.0 "ratio_pareto_tournament_selection must be between 0.0 and 1.0"
    @assert typeof(tournament_selection_fitness) == Vector{Tuple{Float64, Symbol}}

    sel_params = (
        hall_of_fame_objectives           = hall_of_fame_objectives,
        selection_objectives              = selection_objectives,
        tournament_selection_fitness      = tournament_selection_fitness,
        tournament_size                   = tournament_size,
        ratio_pareto_tournament_selection = ratio_pareto_tournament_selection
    )

    return sel_params
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

    t_lim < 1e-1             && @warn "t_lim may be too low"
    max_iter < 5             && @warn "max_iter may be too low"
    early_stop_iter == 0     && @warn "no early stopping may lead to overfitting"
    0 < early_stop_iter <= 2 && @warn "early stopping may be too strict -> higher values may produce better results"
    lasso_factor < 1.0       && @warn "lasso_factor seems to large"

    @assert max_iter >= early_stop_iter "early_stop_iter should be smaller than max_iter"
    @assert 0 <= rel_f_tol_5_iter < 1.0 "rel_f_tol_5_iter must smaller than 1.0 and larger or equal to 0"
    @assert lasso_factor >= 0           "lasso factor must me >= 0"

    nam_tup = (
        max_iter                 = max_iter,
        early_stop_iter          = early_stop_iter,
        rel_f_tol_5_iter         = rel_f_tol_5_iter,
        t_lim                    = t_lim,
        lasso_factor             = lasso_factor,
        pre_residual_processing! = pre_residual_processing!,
        residual_processing      = residual_processing,
    )
    return nam_tup
end

""" Returns a Tuple of the normalized cumulative sum of probabilities of the various
    mutations to be appended to the Options. The input values don't need to add
    up to 1, since they are normalized here.
"""
function mutation_params(;
    p_crossover        = 4.0,
    p_point            = 0.5,
    p_innergrow        = 0.0,
    p_insert           = 0.2,
    p_hoist            = 0.2,
    p_subtree          = 0.2,
    p_add_term         = 0.1,
    p_simplify         = 0.5,
    p_drastic_simplify = 0.5,
)
    p_crossover *= 0.5 # because it searches for a second one, if hit

    mut_params = (
        p_insert,
        p_point,
        p_add_term,
        p_hoist,
        p_innergrow,
        p_subtree,
        p_drastic_simplify,
        p_simplify,
        p_crossover,
    )

    return cumsum(mut_params ./ sum(mut_params))
end

""" Performs a random data split according to parts variable. Then, the data_descript function
    below is called.
"""
function data_descript(
    data::Matrix;
    arbitrary_name = "",
    parts::Vector  = [1.0],
    fit_weights    = ones(size(data, 1)),
)
    eachind = collect(1:size(data, 1))

    @assert sum(parts) == 1.0    "sum of the parts should be 1.0"
    length(parts) > 1 && parts[1] < parts[2] && @warn "more test data than training data"
    @assert length(fit_weights) == size(data, 1) "size of fit_weights don't match data"

    # make a random data split
    shuffle!(eachind)
    start_inds = cumsum([round(Int64, size(data, 1) * p) for p in parts])
    pushfirst!(start_inds, 1)
    split_inds = [eachind[start_inds[i]:start_inds[i+1]] for i in 1:length(parts)]

    return data_descript(
        data,
        split_inds,
        arbitrary_name = arbitrary_name,
        fit_weights    = fit_weights
    )
end

""" Takes the data prepares it for the algorithm. The split_inds variable is a vector of vectors
    of Int64, which specifies which rows to use for fitting (split_inds[1]) and which to use for
    Early Stopping (split_inds[2]) if turned on.
"""
function data_descript(
    data::Matrix,
    split_inds::Vector{Vector{Int64}};
    arbitrary_name = "",
    fit_weights    = ones(size(data, 1)),
)
    @assert eltype(data) <: AbstractFloat "Data is not of float type"

    # convert to vectors of vectors
    data = [data[:, i] for i in 1:size(data, 2)]

    @assert !any(length(unique(d)) == 1 for d in data)                  "data containts rows, which are constant"
    @assert length(unique(reduce(vcat, split_inds))) == length(data[1]) "split_inds seems to have too many or to few entries"
    @assert length(fit_weights) == length(data[1])                      "fit_weights seems to have too many or to few entries"
    @assert typeof(arbitrary_name) <: String                            "arbitrary_name must be a string"

    @assert length(data) > 1                                            "data has only one column"

    length(split_inds) > 2 && @warn   "Only the first and the second data split are used. The remainder are ignored during fitting but the final measures are calculated for all data."

    data_type = eltype(data[1])
    n_vars    = length(data) - 1
    n_points  = length(data[1])

    nam_tup = (
        arbitrary_name        = arbitrary_name,
        data_type             = data_type,
        split_inds            = split_inds,
        n_vars                = n_vars,
        n_points              = n_points,
        fit_weights           = fit_weights,
    )
    return nam_tup, data
end

