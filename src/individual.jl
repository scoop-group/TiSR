
mutable struct Individual   #{T}
    node::Node
    age::Float64
    valid::Bool

    # complexity measures
    compl::Float64
    weighted_compl::Float64
    recursive_compl::Float64
    n_params::Float64

    # fit measures -> how to split in train and test
    ms_processed_e::Float64           # T
    mae::Float64                      # T
    mse::Float64                      # T
    max_ae::Float64                   # T
    minus_r2::Float64                 # T
    minus_abs_spearman::Float64       # T
    mare::Float64                     # T
    q75_are::Float64                  # T
    max_are::Float64                  # T

    Individual() = new()

    function Individual(node, data, ops, cur_max_compl, fit_iter)

        apply_simple_simplifications!(node, ops)

        trim_to_max_nodes_per_term!(node, ops)

        if count_nodes(node) > ops.grammar.max_compl
            target_compl = rand(5:min(
                    cur_max_compl + ops.general.adaptive_compl_increment,
                    ops.grammar.max_compl
            ))
            trim_to_max_compl!(
                node,
                target_compl,
                ops
            )
        end

        apply_simple_simplifications!(node, ops)

        div_to_mul_param!(node, ops)
        reorder_add_n_mul!(node, ops)

        indiv = new(node)

        if count_nodes(node) < ops.grammar.min_compl
            indiv.valid = false
            return indiv
        end

        if !isempty(ops.grammar.illegal_dict) && !is_legal_nesting(node, ops)
            indiv.valid = false
            return indiv
        end

        fit_n_eval!(indiv, data, ops, fit_iter)

        indiv.valid || return indiv

        indiv.age             = 0.0
        indiv.compl           = count_nodes(node)

        if isempty(ops.grammar.weighted_compl_dict)
            indiv.weighted_compl = indiv.compl
        else
            indiv.weighted_compl = get_weighted_compl(node, ops)
        end

        indiv.recursive_compl = recursive_compl(node, ops)
        indiv.n_params        = length(list_of_param_nodes(node))

        return indiv
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
            round(getfield(indiv, obj), sigdigits=ops.general.remove_doubles_sigdigits)
            for obj in [:mse, :mae, :compl]
        ]
        for indiv in individs
    ]

    unique_inds = collect(eachindex(individs))

    _remove_doubles_helper!(indiv_obj_vals, unique_inds)

    sort!(unique_inds)
    keepat!(individs, unique_inds)
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

