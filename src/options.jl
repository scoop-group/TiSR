
""" The Options struct contains all settings. For most of its fields, there are functions with
    default arguments to make sure all fields are filled and potentially adapted to make sense.
    There are also some safeguards, preventing certain settings combinations or warn the user,
    which may not cover all cases.
"""
struct Options{A, B, C, D, F, G, H, I, J, K}#, L}
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
    # dynam_expr::L

    function Options(
        data::Matrix;                                               # -> nxm matrix containing the n data points, m-1 variables and the output
        fit_weights::Vector{Float64} = inv.(abs.(data[:, end])),    # -> weights for the data fitting -> residual vector .* weights vector
        binops                       = (+,   -,   *,   /,   ^  ),   # -> binary function set
        unaops                       = (exp, log, sin, cos, abs),   # -> unary function set
        data_split                   = data_split_params(),
        general                      = general_params(),
        measures                     = measure_params(),
        selection                    = selection_params(),
        fitting                      = fitting_params(),
        mutation                     = mutation_params(),
        grammar                      = grammar_params(),
        meta_data                    = Dict(),                      # -> can be used to provide data for use in, for example, user-defined functions like the callback, pre_residual_processing, or measures
    )

        # make sure all required measres are calculated # ------------------------------------------
        @assert all(s in keys(measures) for s in selection.selection_objectives)    "some selection_objectives are not specified in the meausres"
        @assert all(s in keys(measures) for s in selection.hall_of_fame_objectives) "some hall_of_fame_objectives are not specified in the meausres"

        # prepare data and its params --------------------------------------------------------------
        data_vect, data_descript = prepare_data(data, fit_weights, general.replace_inf)

        split_inds = _data_split_params(data, data_split...)
        data_descript = (data_descript..., split_inds = split_inds)

        # ------------------------------------------------------------------------------------------
        @assert length(unaops) <= 128 "TiSR cannot handle more than 128 functions in the unaops out of the box"
        @assert length(binops) <= 128 "TiSR cannot handle more than 128 functions in the binops out of the box"

        # function set asserts ---------------------------------------------------------------------
        @assert ("+" in string.(binops)) "+ is required in the function set for some mutations to work"
        @assert ("*" in string.(binops)) "* is required in the function set for some mutations to work"
        :^ in Symbol.(binops) && @warn   "^ is only valid for positive bases. Otherwise provide a pow(abs(x), y) function with its own drawbacks."

        # some inter-subsection parameter assertions # ---------------------------------------------
        if fitting.early_stop_iter != 0
            @assert !isempty(data_descript.split_inds[2]) "second data split must contain data for early stopping"
        end

        if :weighted_compl in selection.selection_objectives || :weighted_compl in selection.hall_of_fame_objectives
            "VAR"   in keys(grammar.weighted_compl_dict) || @warn "'VAR' not specified in weighted_compl_dict -> 3.0 is assumed"
            "PARAM" in keys(grammar.weighted_compl_dict) || @warn "'PARAM' not specified in weighted_compl_dict -> 3.0 is assumed"
            for op in binops
                string(op) in keys(grammar.weighted_compl_dict) || @warn "$op not specified in weighted_compl_dict -> 3.0 is assumed"
            end
            for op in unaops
                string(op) in keys(grammar.weighted_compl_dict) || @warn "$op not specified in weighted_compl_dict -> 3.0 is assumed"
            end
        end

        @assert !xor(isempty(grammar.bank_of_terms), mutation.mut_probs[10] == 0) "for p_add_from_bank_of_terms != 0, bank_of_terms cannot be empty and vv"

        # test measures requrie data split # -------------------------------------------------------
        test_funcs = [
            get_measure_max_ae_test, get_measure_mae_test, get_measure_mse_test, get_measure_one_minus_r2_test,
            get_measure_one_minus_abs_spearman_test, get_measure_mare_test, get_measure_q75_are_test,
            get_measure_max_are_test, get_measure_ms_processed_e_test
        ]

        if any(t in values(measures) for t in test_funcs)
            @assert !isempty(data_descript.split_inds[3]) "for the test measures, the data split in position 3 must not be empty"
        end

        # resulting parameters # -------------------------------------------------------------------
        n_pareto_select_per_isle = ceil(Int64, general.pop_per_isle * selection.ratio_pareto_tournament_selection)
        selection = (;selection..., n_pareto_select_per_isle = n_pareto_select_per_isle)

        global __operators = (unaops = unaops, binops = binops) # required for Base.show() Nodes

        # # dynamic expressions # --------------------------------------------------------------------
        # operators      = DynamicExpressions.OperatorEnum(; binary_operators=collect(binops), unary_operators=collect(unaops))
        # variable_names = ["v$i" for i in 1:data_descript.n_vars]
        # dynam_expr     = (operators = operators, variable_names = variable_names)

        return new{typeof(general),
                   typeof(unaops),
                   typeof(binops),
                   typeof(measures),
                   typeof(selection),
                   typeof(fitting),
                   typeof(mutation),
                   typeof(grammar),
                   typeof(data_descript),
                   typeof(meta_data),
                   # typeof(dynam_expr),
                   }(
                       general,
                       unaops,
                       binops,
                       measures,
                       selection,
                       fitting,
                       mutation,
                       grammar,
                       data_descript,
                       meta_data,
                       # dynam_expr,
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

function prepare_data(data, fit_weights, replace_inf)
    data_vect = [data[:, i] for i in 1:size(data, 2)]

    # prepare data # ---------------------------------------------------------------------------
    @assert eltype(data) <: AbstractFloat "data must be float type "
    @assert size(data, 2) > 1             "data has only one column"
    @assert size(data, 2) <= 128 + 1      "TiSR cannot handle more than 128 variables out of the box"
    !any(length(unique(d)) == 1 for d in data_vect) || @warn "data containts columns, which are constant"


    @assert length(fit_weights) == size(data, 1) "fit_weights must have the same length as the data points"

    if any(!isfinite, fit_weights)
        println("some fit_weights are nonfinite, replacing those with ops.general.replace_inf")
        replace!(fit_weights, Inf => replace_inf)
        @assert all(isfinite, fit_weights) "fit_weights still nonfinite after replacing Infs, there might be -Inf or NaN"
    end

    @assert all(fit_weights .> 0.0)              "fit_weights must be larger than 0"

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
    n_gens::Int64                          = typemax(Int64),                                                    # -> number of generations to conduct
    t_lim::Float64                         = 60. * 5.,                                                          # -> time limit for the algorithm
    pop_size::Int64                        = 1000,                                                              # -> number of individuals selected for next generation / population size
    parent_selection::Bool                 = true,                                                              # -> wheather to conduct parent selection or to use all individuals in the population as parents
    num_islands::Int64                     = 20,                                                                # -> number of parallel islands
    children_ratio::Float64                = 1.0,                                                               # -> the ratio of children that should be generated in each generation 0 ... 2
    migration_interval::Int64              = 200,                                                               # -> generation interval, in which an individual is moved to other islands. (ring topology)
    island_extinction_interval::Int64      = 1000,                                                              # -> interval in which all individuals from one islands are distributed across all other islands and the extiction islands starts from scratch. -> typemax(Int64) is off; 1000 ... 10000
    island_extinction_rotation::Bool       = false,                                                             # -> whether the island extinctions should rotate as opposed to random islands
    migrate_after_extinction_prob::Float64 = 1.0,                                                               # -> probability that an individual migrates to another islands, if its island goes extinct. -> 0 ... 1
    migrate_after_extinction_dist::Int64   = 4,                                                                 # -> maximum relative distance an individual from an extinct island can propagate to a new island in case it survives. -> 0.2 ... 0.5
    fitting_island_function::Function      = isle -> floor(isle / 2) % 2 == 0,                                  # -> function to determine on which islands fitting is conducted. Must take an integer and return a bool
    hall_of_fame_migration_interval::Int64 = 500,                                                               # -> interval in which a random individual from the hall of fame is returned to a random island
    max_age::Int64                         = 2 * round(Int64, pop_size / num_islands),                          # -> maximal age after which individuals are removed from the popoulation
    n_refitting::Int64                     = 1,                                                                 # -> how many individuals from the hall_of_fame are copied and fitted again
    adaptive_compl_increment::Int64        = 100,                                                               # -> highest complexity in the hall of fame + `adaptive_compl_increment` is the highest allowed complexity for a new individual; -> Inf is off; 5 ... 10
    seen_reject_prob::Float64              = 0.9,                                                               # -> probability to reject an expression, if it is present in the forgetting expression log, i.e., has been seen recently or frequently
    seen_forget_interval::Int64            = 100,                                                               # -> generation interval in which the expression_log forgets
    seen_merge_interval::Int64             = 100,                                                               # -> generation interval in which the expression_logs of each island are merged
    seen_gc_length::Int64                  = 100_000,                                                           # -> number of entries in the expression_log above which an additional forgetting is applied
    replace_inf::Float64                   = 1e100,                                                             # -> the value infs measures should be replaced with. The inverse is also used as the guarding offset g -> 1/(0 + g)
    callback::Function                     = (hall_of_fame, population, gen, t_since, prog_dict, ops) -> false, # -> a function, which is executed in each iteration and allows more flexible termination. If the function returns true, the execution is terminated. For example, the following stops the equation search, if one individual in the hall of fame has a complexity lower than 30 and a mean absolute relative deviation of lower then 1e-5: `(hall_of_fame, population, gen, prog_dict, ops) -> any(i.measures[:compl] < 30 && i.measures[:mare] < 1e-5 for i in hall_of_fame)`
    multithreading::Bool                   = false,                                                             # -> whether to use multithreading for the fitting (most expensive). Not always faster -> depends on how expensive fitting is for the problem at hand. Also, for this to apply, Julia needs to be started with more threads, like `julia -t 4`.
    print_progress::Bool                   = true,                                                              # -> whether to print the elapsed time and some other KPIs.
    plot_hall_of_fame::Bool                = true,                                                              # -> whether to plot the hall of fame
    print_hall_of_fame::Bool               = true,                                                              # -> whether to print some of the individuals in the hall of fame. For this, `plot_hall_of_fame` must also be true
)
    @assert num_islands > 0                                       "num_islands should be at least 1                                    "
    @assert migration_interval > 0                                "migration_interval should be at least 1                             "
    @assert adaptive_compl_increment > 0                          "adaptive_compl_increment must be larger 0                           "
    @assert island_extinction_interval > 0                        "island_extinction_interval must be > 0                              "
    @assert max_age > 1                                           "max_age must be > 1                                                 "
    @assert hall_of_fame_migration_interval > 0                   "hall_of_fame_migration_interval must be larger 0                    "
    @assert 0.0 <= migrate_after_extinction_prob <= 1.0           "migrate_after_extinction_prob must be between 0 and 1               "
    @assert 0.0 <= children_ratio <= 2.0                          "children_ratio should be inbetween 0.0 and 2.0                      "
    @assert 0 <= n_refitting <= pop_size / num_islands            "n_refitting should be between 0 and pop_size / num_islands          "
    @assert 0 <= migrate_after_extinction_dist <= num_islands / 2 "migrate_after_extinction_dist must be between 0 and num_islands/2   "
    @assert !isnan(replace_inf)                                   "replace_inf must not be NaN"
    @assert 0.0 <= seen_reject_prob <= 1.0                        "seen_rejection_probability must be between 0 and 1"
    @assert seen_forget_interval > 0                              "seen_forget_interval must be larger 0"
    @assert seen_merge_interval > 0                               "seen_merge_interval must be larger 0"
    @assert seen_gc_length > 0                                    "seen_gc_length must be larger 0"

    if print_hall_of_fame
        @assert plot_hall_of_fame "for print_hall_of_fame, plot_hall_of_fame must be true"
    end

    adaptive_compl_increment > 4     || @warn "adaptive_compl_increment should be >= 5                         "
    island_extinction_interval > 500 || @warn "island_extinction_interval seems small                          "
    max_age > 10                     || @warn "max_age seems small                                             "

    10 <= seen_forget_interval <= 500 || @warn "seen_forget_interval should be between 10 and 500"
    10 <= seen_merge_interval <= 500  || @warn "seen_merge_interval should be between 10 and 500"
    10_000 <= seen_gc_length <= 500_000 || @warn "seen_gc_length should be between 10_000 and 500_000"

    island_extinction_interval == migration_interval == typemax(Int64) && @warn "island_extinction_interval & migration_interval should not both be off"

    if multithreading
        Threads.nthreads() > 1 || @warn "To acually utilize multithreading, Julia must be started with the desired number of threads, i.e., `julia -t 4`"
    end

    # resulting parameters
    pop_per_isle = ceil(Int64, pop_size / num_islands)
    n_children = max(2, round(Int64, children_ratio * pop_per_isle))
    migrate_after_extinction_num = round(Int64, pop_per_isle * migrate_after_extinction_prob)
    if migrate_after_extinction_num > 0
        @assert migrate_after_extinction_dist > 0 "if migrate_after_extinction_prob is > 0, migrate_after_extinction_dist must be > 0"
    end

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
        island_extinction_rotation      = island_extinction_rotation,
        migrate_after_extinction_num    = migrate_after_extinction_num,
        migrate_after_extinction_dist   = migrate_after_extinction_dist,
        fitting_island_function         = fitting_island_function,
        max_age                         = max_age,
        n_refitting                     = n_refitting,
        t_lim                           = t_lim,
        multithreading                  = multithreading,
        adaptive_compl_increment        = adaptive_compl_increment,
        seen_reject_prob                = seen_reject_prob,
        seen_forget_interval            = seen_forget_interval,
        seen_merge_interval             = seen_merge_interval,
        seen_gc_length                  = seen_gc_length,
        callback                        = callback,
        replace_inf                     = replace_inf,
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
    additional_measures::Dict{Symbol, Function} = Dict(                # -> specify a Dict{Symbol, Function} containing name and function pairs that calculate measures. The function must take 4 positional arguments `prediction::Vector{T}, target::Vector{T}, node, ops` and return a Float. All currently pre-implemented measures are listed below at (2).
        :one_minus_abs_spearman => get_measure_one_minus_abs_spearman,
        :mare                   => get_measure_mare,
        :max_are                => get_measure_max_are,
        :weighted_compl         => get_measure_weighted_compl,
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
    hall_of_fame_objectives::Vector{Symbol}    = [:ms_processed_e, :weighted_compl],                          # -> objectives for the hall_of_fame selection
    selection_objectives::Vector{Symbol}       = [:ms_processed_e, :one_minus_abs_spearman, :weighted_compl], # -> objectives for the population selection
    hall_of_fame_niching_sigdigits::Int64      = 5,                                                           # -> number of significant digits to round hall_of_fame_objectives for hall_of_fame selection after their normalization. -> 2 ... 5; 0 is off
    population_niching_sigdigits::Int64        = 5,                                                           # -> number of significant digits to round selection_objectives for population selection after their normalization. -> 2 ... 5; 0 is off
    ratio_pareto_tournament_selection::Float64 = 0.5,                                                         # -> ratio to which the selection is conducted using the Pareto-optimal selection vs. tournament selection
    tournament_size::Int64                     = 5,                                                           # -> tournament size
)
    @assert tournament_size > 1                             "tournament size must be greater than 1"
    @assert 0.0 <= ratio_pareto_tournament_selection <= 1.0 "ratio_pareto_tournament_selection must be between 0.0 and 1.0"

    0 <= hall_of_fame_niching_sigdigits < 10 || @warn "hall_of_fame_niching_sigdigits should be between 0 and 5"
    0 <= population_niching_sigdigits   < 10 || @warn "population_niching_sigdigits should be between 0 and 5"

    return (
        hall_of_fame_objectives           = hall_of_fame_objectives,
        selection_objectives              = selection_objectives,
        hall_of_fame_niching_sigdigits    = hall_of_fame_niching_sigdigits,
        population_niching_sigdigits      = population_niching_sigdigits,
        tournament_size                   = tournament_size,
        ratio_pareto_tournament_selection = ratio_pareto_tournament_selection
    )
end

function fitting_params(;
    max_iter::Int64                       = 10,                 # -> maximum iterations for parameter fitting. -> 10 ... 50 ==> biggest time consumer <==
    NM_iter::Int64                        = 50,                 # -> maximum iterations for parameter fitting with Nelder-Mead. -> 20 ... 100 ==> biggest time consumer <==
    NM_prob::Float64                      = 0.1,                # -> probabitlity that fitting is conducted with Nelder-Mead as opposed to Levenberg-Marquard -> 0.01 ... 0.2
    early_stop_iter::Int64                = 0,                  # -> how many iterations to account for early stopping regularization. to use, the data needs to be partitioned into at least 2 parts. The early stopping evaluation is performed on the second partition. -> 0 is off; 4 ... 10
    t_lim::Float64                        = Inf,                # -> time limit for parameter fitting of individual. -> Inf is off; 0.1 ... 0.5
    rel_f_tol_5_iter::Float64             = 1e-2 * 0.01,        # -> relative tolerance for parameter fitting. considered converged if relative improvement over 5 iterations is smaller. -> 0 is off; 1e-2 * 1.0 ... 1e-2 * 0.01
    lasso_factor::Float64                 = 0.0,                # -> factor for the lasso regularization. pushing parameter values to 0. -> 0 is off; 1e-8 ... 1e-4
    pre_residual_processing::Function     = (x, ind, ops) -> x, # -> processing of the equation output before the residual is calculated. The inds refer for the indices of the current residuals, which may be used to slice some data in the function like "(x, inds, ops) -> x ./= data[end][inds]"
    residual_processing::Function         = (x, ind, ops) -> x, # -> processing of the residuals. NOT an inplace function. The inds refer for the indices of the current residuals, which may be used to slice some data in the function like "(x, inds, ops) -> x ./ data[end][inds]"
    all_constr_f_select::Vector{Function} = Function[],         # TODO
    eq_constr::Vector{Function}           = Function[],         # TODO
    ineq_constr::Vector{Function}         = Function[],         # TODO
    max_mare_for_constr_fit::Float64      = 0.1,                # -> max mean relative error at which shape constraints should be minimized -> don't cast pearls before swine
    additional_constr_fit_iter::Int64     = 5,                  # -> additional fitting iterations that are conducted with the constraint violation minimization
    constr_tol::Float64                   = 1e-5,               # TODO
)
    @assert NM_iter >= 0                        "NM_iter must be >= 0                                           "
    @assert max_iter >= early_stop_iter         "early_stop_iter should be smaller than max_iter                "
    @assert 0 <= rel_f_tol_5_iter < 1.0         "rel_f_tol_5_iter must smaller than 1.0 and larger or equal to 0"
    @assert lasso_factor >= 0                   "lasso factor must be >= 0                                      "
    @assert 0 <= NM_prob <= 1.0                 "NM_prob must be between 0 and 1                                "
    @assert constr_tol >= 0                     "constr_tol must be >= 0"

    t_lim > 1e-1             || @warn "fitting t_lim may be too low"
    max_iter >= 5            || @warn "max_iter may be too low"
    0 <= NM_iter <= 100      || @warn "0 <= NM_iter <= 100"
    lasso_factor < 1.0       || @warn "lasso_factor seems to large"
    0 < early_stop_iter < 5  && @warn "early stopping may be too strict -> higher values may produce better results"

    # for constraints
    if !isempty(eq_constr) || !isempty(ineq_constr)
        @assert early_stop_iter == 0           "early stopping and constrained fitting not implemented together. Please set early_stop_iter = 0 and reconsider data split"
        @assert 0 <= max_mare_for_constr_fit   "max_mare_for_constr_fit must be larger than 0"
        @assert additional_constr_fit_iter > 0 "additional_constr_fit_iter must be > 0. If the constraints should not be minimized during parameter estimation, leave the eq_constr and the ineq_constr empty and provide only all_constr_f_select"

        length(eq_constr) + length(ineq_constr) <= length(all_constr_f_select) || @warn "it seems that not all eq_constr and ineq_constr also have a respective version in all_constr_f_select. This leads to their minimization during the parameter estimation but their violations are disregarded during the selection."
        rel_f_tol_5_iter == 0                                                  || @warn "rel_f_tol_5_iter does not work with constrained least squares and will be ignored"
        1 <= additional_constr_fit_iter <= 20                                  || @warn "additional_constr_fit_iter should be between 1 and 20"
    end

    return (
        max_iter                   = max_iter,
        NM_iter                    = NM_iter,
        NM_prob                    = NM_prob,
        early_stop_iter            = early_stop_iter,
        rel_f_tol_5_iter           = rel_f_tol_5_iter,
        t_lim                      = t_lim,
        lasso_factor               = lasso_factor,
        pre_residual_processing    = pre_residual_processing,
        residual_processing        = residual_processing,
        all_constr_f_select        = all_constr_f_select,
        eq_constr                  = eq_constr,
        ineq_constr                = ineq_constr,
        max_mare_for_constr_fit    = max_mare_for_constr_fit,
        additional_constr_fit_iter = additional_constr_fit_iter,
        constr_tol                 = constr_tol,
    )
end

""" Returns the equation grammar related parameters.
"""
function grammar_params(;
    max_compl::Int64                           = 30,                        # -> max allowed complexity
    min_compl::Int64                           = 3,                         # -> min allowed complexity. -> 2 ... 3
    init_tree_depth::Int64                     = 5,                         # -> maximal initial tree depth. -> 3 ... 6
    max_nodes_per_term::Int64                  = typemax(Int64),            # -> maximal number of nodes per top-level term. All terms above are trimmed until they satisfy this threshold
    bank_of_terms::Vector{String}              = String[],                  # -> specify terms that can be added via the add term mutation, whose relativ probability is set via the p_add_from_bank_of_terms parameter
    illegal_dict::Dict                         = Dict(),                    # -> check for illegal nestings in existing nodes. For it, ops.illegal_dict needs to be specified like below. Variables can be specified using "VAR" and parameters with "PARAM". An example is shown below at (1).
    custom_check_legal_before_fit::Function    = (node, data, ops) -> true, # -> specify a custom function, which checks the legality of nodes. Must return true or false
    weighted_compl_dict::Dict{String, Float64} = Dict(                      # -> weights for weighted_compl calculation. For any that are not included, 3.0 is assumed. The weights for variables and parameters "VAR" and "PARAM" may be used. An example is shown below at (2).
        "PARAM" => 1.2, "VAR"  => 1.0,
        "+"     => 1.2, "-"    => 1.4,
        "*"     => 1.0, "/"    => 1.6,
        "^"     => 3.0, "pow"  => 3.5,
        "neg"   => 1.4, "abs"  => 2.0,
        "sqrt"  => 2.0, "pow2" => 2.0, "pow3" => 2.1,
        "sin"   => 3.0, "cos"  => 3.0, "tanh" => 3.0,
        "exp"   => 3.0, "log"  => 3.0,
    ),
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
        illegal_dict                  = illegal_dict,
        weighted_compl_dict           = weighted_compl_dict,
        init_tree_depth               = init_tree_depth,
        max_compl                     = max_compl,
        min_compl                     = min_compl,
        max_nodes_per_term            = max_nodes_per_term,
        bank_of_terms                 = bank_of_terms,
        custom_check_legal_before_fit = custom_check_legal_before_fit,
    )
end

""" Returns a Tuple of the normalized cumulative sum of probabilities of the various
    mutations to be appended to the Options. The input values don't need to add
    up to 1, since they are normalized here.
"""
function mutation_params(;                     #|-> probabilites for the various mutations (don't need to add up to 1)
    p_crossover::Float64              = 0.5,   #|
    p_point::Float64                  = 1.0,   #|
    p_point2::Float64                 = 0.3,   #|
    p_insert::Float64                 = 2.0,   #|
    p_hoist::Float64                  = 2.0,   #|
    p_subtree::Float64                = 0.5,   #|
    p_drastic_simplify::Float64       = 1.0,   #|-> remove parameter nodes with small values and simplify accordingly
    p_insert_times_param::Float64     = 0.1,   #|
    p_add_term::Float64               = 0.1,   #|
    p_add_from_bank_of_terms::Float64 = 0.0,   #|-> probability to add a term from the provided bank_of_terms
    p_multiple_mutations::Float64     = 0.5,   # -> probability for more than one mutation
    max_muts_ratio::Float64           = 0.75,  # -> max ratio of multiple mutations wrt to the complexity of the expression
)
    p_simplify = 0.0 #|-> simplify with SymbolicUtils # -> does not work in multithreading, and is slow
    @assert all(p >= 0 for p in (
        p_crossover, p_point, p_point2, p_insert, p_hoist, p_subtree, p_drastic_simplify,
        p_insert_times_param, p_add_term, p_simplify, p_add_from_bank_of_terms
    )) "all mutation probabilities must be >= 0"

    p_crossover *= 0.5 # because it consumes a second one, if hit

    @assert p_multiple_mutations < 1         "p_multiple_mutations must be < 1"
    0 <= p_multiple_mutations < 0.9 || @warn "p_multiple_mutations should be between 0 and 0.9"

    @assert 0.0 <= max_muts_ratio            "max_muts_ratio must be >= 0                 "
    0.0 <= max_muts_ratio <= 1.0    || @warn "max_muts_ratio should be between 0.0 and 2.0"

    mut_probs = (
        p_insert,
        p_point,
        p_point2,
        p_add_term,
        p_insert_times_param,
        p_hoist,
        p_subtree,
        p_drastic_simplify,
        p_simplify,
        p_add_from_bank_of_terms,
        p_crossover,
    )

    multiple_mut_probs = (p_insert, p_point, p_point2, p_hoist)

    mut_probs          = mut_probs          ./ sum(mut_probs)
    multiple_mut_probs = multiple_mut_probs ./ sum(multiple_mut_probs)
    return (
        mut_probs            = mut_probs,
        multiple_mut_probs   = multiple_mut_probs,
        p_multiple_mutations = p_multiple_mutations,
        max_muts_ratio       = max_muts_ratio,
    )
end
