
mutable struct Individual
    node::Node
    valid::Bool
    rank::Int64
    crowding::Float64
    age::Float64

    measures::Dict{Symbol, Float64}

    Individual() = new()
    Individual(node::Node) = new(node)
end

function fit_individual!(indiv, data, ops, cur_max_compl, fit_iter)
    indiv.age = 0.0

    # prepare node -> simplify, trim, reorder # ----------------------------------------------------
    apply_simple_simplifications!(indiv.node, ops)
    trim_to_max_nodes_per_term!(indiv.node, ops)
    apply_simple_simplifications!(indiv.node, ops)
    div_to_mul_param!(indiv.node, ops)
    reorder_add_n_mul!(indiv.node, ops)

    # remove invalids # ----------------------------------------------------------------------------
    if !(ops.grammar.min_compl <= count_nodes(indiv.node) <= min(ops.grammar.max_compl, cur_max_compl + ops.general.adaptive_compl_increment))
        indiv.valid = false
        return
    end

    if !isempty(ops.grammar.illegal_dict) && !is_legal_nesting(indiv.node, ops)
        indiv.valid = false
        return
    end

    if !ops.grammar.custom_check_legal(indiv.node, data, ops)
        indiv.valid = false
        return
    end

    # fitting # ------------------------------------------------------------------------------------
    prediction, valid = fit_n_eval!(indiv.node, data, ops, fit_iter)

    indiv.valid = valid
    indiv.valid || return

    # calculate measures # -------------------------------------------------------------------------
    indiv.measures = Dict(
        m => f(prediction, data[end], indiv.node, ops)
        for (m, f) in ops.measures
    )

    if any(!isfinite(v) for (k, v) in indiv.measures)
        indiv.valid = false
        return
    end
end

# ==================================================================================================
# some helpers
# ==================================================================================================
""" prints or displays an Individual
"""
Base.show(io::IO, indiv::Individual) = println(io, "Indiv($(node_to_string(indiv.node, __operators; sigdigits=3)))")

""" copies an individual.
"""
function Base.copy(indiv::Individual)
    new = Individual()
    for field_ in fieldnames(Individual)
        if isdefined(indiv, field_)
            setfield!(new, field_, copy(getfield(indiv, field_)))
        end
    end
    return new
end

Base.deepcopy(indiv::Individual)::Individual = Base.copy(indiv)

""" make individual type iterable.
"""
Base.iterate(s::Individual, state=1) = state > length(fieldnames(Individual)) ? nothing : ((fieldnames(Individual)[state], getfield(s, state)), state + 1)
Base.length(s::Individual) = length(fieldnames(Individual))

""" The NSGA-II definition of isless.
"""
function Base.isless(i1::Individual, i2::Individual)
    if i1.rank < i2.rank
        return true
    elseif i1.rank > i2.rank
        return false
    else
        if i1.crowding < i2.crowding
            return false
        else
            return true
        end
    end
end

# ==================================================================================================
# different remove doubles -> maybe combine somehow with views?
# ==================================================================================================
""" Removes similar individual from a population isle. The individuals are considered similar,
    if their MSE and MAE rounded to ops.general.remove_doubles_sigdigits significant digits are
    the same. Of those, the one with the least complexity remains in the population.
"""
function remove_doubles!(individs::Vector{Individual}, ops)

    indiv_obj_vals = [
        Float64[
            round(indiv.measures[obj], sigdigits=ops.general.remove_doubles_sigdigits)
            for obj in [:mse, :mae, :compl]
        ]
        for indiv in individs
    ]

    unique_inds = collect(eachindex(individs))

    _remove_doubles_helper!(indiv_obj_vals, unique_inds)

    sort!(unique_inds)
    keepat!(individs, unique_inds)
end

""" Removes similar individual across all islands. The individuals are considered similar,
    if their MSE and MAE rounded to ops.general.remove_doubles_sigdigits significant digits are
    the same. Of those, the one with the least complexity remains in the population.
"""
function remove_doubles_across_islands!(individs::Vector{Vector{Individual}}, ops)

    indiv_obj_vals = [
        Float64[
            round(indiv.measures[obj], sigdigits=ops.general.remove_doubles_sigdigits)
            for obj in [:mse, :mae, :compl]
        ]
        for isle in 1:ops.general.num_islands for indiv in individs[isle]
    ]

    unique_inds = [
        (isle, ind)
        for isle in 1:ops.general.num_islands for ind in eachindex(individs[isle])
    ]

    _remove_doubles_helper!(indiv_obj_vals, unique_inds)

    # apply the unique_inds # ---------------------------------------------------------------------
    unique_inds = [
        [
            unique_inds[i][2]
            for i in 1:length(unique_inds)
            if unique_inds[i][1] == isle
        ] for isle in 1:ops.general.num_islands
    ]

    foreach(isle -> sort!(unique_inds[isle]), 1:ops.general.num_islands)
    foreach(isle -> keepat!(individs[isle], unique_inds[isle]), 1:ops.general.num_islands)
end

""" Helper function while for remove_doubles!, which performs the comparisons and modifies
    the unique_inds.
"""
function _remove_doubles_helper!(indiv_obj_vals, unique_inds)

    perm = sortperm(indiv_obj_vals)
    indiv_obj_vals .= indiv_obj_vals[perm]
    unique_inds    .= unique_inds[perm]

    i = 1
    while i < length(indiv_obj_vals)
        while i < length(indiv_obj_vals) && indiv_obj_vals[i][1:2] == indiv_obj_vals[i+1][1:2]
            deleteat!(unique_inds, i+1)
            deleteat!(indiv_obj_vals, i+1)
        end
        i += 1
    end
end

