
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

function fit_individual!(indiv, data, ops, cur_max_compl, fit_iter, expression_log)
    indiv.age = 0.0

    # prepare node -> simplify, trim, reorder # ----------------------------------------------------
    #@timeit to "prepare node" begin
        trim_to_max_nodes_per_term!(indiv.node, ops)
        apply_simple_simplifications!(indiv.node, ops)
        div_to_mul_param!(indiv.node, ops)
        reorder_add_n_mul!(indiv.node, ops)
    #end # @timeit

    # remove invalids # ----------------------------------------------------------------------------
    #@timeit to "remove invalids" begin
        if !(ops.grammar.min_compl <= count_nodes(indiv.node) <= min(ops.grammar.max_compl, cur_max_compl + ops.general.adaptive_compl_increment))
            indiv.valid = false
            return 0
        end

        if !isempty(ops.grammar.illegal_dict) && !is_legal_nesting(indiv.node, ops)
            indiv.valid = false
            return 0
        end

        if !ops.grammar.custom_check_legal(indiv.node, data, ops)
            indiv.valid = false
            return 0
        end
    #end # @timeit

    # check if seen # ------------------------------------------------------------------------------
    #@timeit to "check expression_log" begin
        if ops.general.seen_reject_prob > 0
            #@timeit to "convert for expression_log" begin
                hnode = hash(indiv.node)
            #end # @timeit
            #@timeit to "check and expand expression_log" begin
                expression_log[hnode] = visits = min(get!(expression_log, hnode, 0) + 1, 127)
            #end # @timeit
            if visits > 1 && rand() < ops.general.seen_reject_prob # if rand() > 0.5^visits
                indiv.valid = false
                return 0
            end
        end
    #end # @timeit

    # fitting # ------------------------------------------------------------------------------------
    #@timeit to "fitting" begin
        prediction, valid, eval_counter = fit_n_eval!(indiv.node, data, ops, fit_iter)
    #end # @timeit

    indiv.valid = valid
    indiv.valid || return eval_counter

    # calculate measures # -------------------------------------------------------------------------
    #@timeit to "calc measures" begin
        indiv.measures = NamedTuple(
            m => clamp(
                f(prediction, data[end], indiv.node, ops)::Float64,
                -ops.general.replace_inf,
                ops.general.replace_inf
            )
            for (m, f) in ops.measures
        )
    #end # @timeit

    if any(!isfinite(v) for v in values(indiv.measures))
        indiv.valid = false
    end

    return eval_counter
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
