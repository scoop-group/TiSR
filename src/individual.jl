
mutable struct Individual
    node::Node
    valid::Bool
    rank::Int64
    crowding::Float64
    age::Float64
    measures::NamedTuple

    Individual() = new()
    Individual(node::Node) = new(node)
end

function fit_individual!(indiv, data, ops, cur_max_compl, fit_iter)
    indiv.age = 0.0

    # prepare node -> simplify, trim, reorder # ----------------------------------------------------
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
    indiv.measures = NamedTuple(
        m => f(prediction, data[end], indiv.node, ops)
        for (m, f) in ops.measures
    )

    if any(!isfinite(v) for v in values(indiv.measures))
        indiv.valid = false
        return
    end
end

""" Prints or displays an Individual.
"""
Base.show(io::IO, indiv::Individual) = println(io, "Indiv($(node_to_string(indiv.node, __operators; sigdigits=3)))")

""" Creates a new individual with same node as the provided one.
"""
fastcopy(indiv::Individual) = Individual(deepcopy(indiv.node))

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
