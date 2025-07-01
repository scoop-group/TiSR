
mutable struct Individual
    node::Node
    valid::Int8
    rank::Int16
    crowding::Float64
    age::Float64
    measures::NamedTuple

    Individual() = new()
    Individual(node::Node) = new(node)
end

function fit_individual!(indiv, data, ops, cur_max_compl, expression_log, isle)
    indiv.age   = 0.0
    indiv.valid = 0

    # prepare node -> simplify, trim, reorder # ----------------------------------------------------
    #@timeit to "prepare node" begin
        apply_simple_simplifications!(indiv.node, ops)
        div_to_mul_param!(indiv.node, ops)
        reorder_add_n_mul!(indiv.node, ops)
    #end # @timeit

    # remove invalids # ----------------------------------------------------------------------------
    #@timeit to "remove invalids" begin
        if count_nodes(indiv.node) > min(ops.grammar.max_compl, cur_max_compl + ops.general.adaptive_compl_increment)
            indiv.valid = 1
            return
        end

        if !isempty(ops.grammar.illegal_dict) && !is_legal_nesting(indiv.node, ops)
            indiv.valid = 2
            return
        end

        if !ops.grammar.custom_check_legal_before_fit(indiv.node, data, ops)
            indiv.valid = 3
            return
        end
    #end # @timeit

    # check if seen # ------------------------------------------------------------------------------
    #@timeit to "check expression_log" begin
        if ops.general.seen_reject_prob > 0
            #@timeit to "convert for expression_log" begin
                node_hash = hash(indiv.node)
                #reject_rate[2] += 1 # DEBUG expression_log
            #end # @timeit
            #@timeit to "check and expand expression_log" begin
                expression_log[node_hash] = visits = min(get!(expression_log, node_hash, 0) + 1, 127)
            #end # @timeit
            if visits > 1 && rand() < ops.general.seen_reject_prob # if rand() > 0.5^visits
                #reject_rate[1] += 1 # DEBUG expression_log
                indiv.valid = 4
                return
            end
        end
    #end # @timeit

    # fitting # ------------------------------------------------------------------------------------
    #@timeit to "fitting" begin
        prediction, valid  = fit_n_eval!(
            indiv.node,
            data,
            ops,
            do_fit    = ops.general.fitting_island_function(isle),
            do_constr = ops.general.constraint_island_function(isle),
        )
    #end # @timeit

    if !valid
        indiv.valid = 5
        return
    end

    # after fitting legal checks # -----------------------------------------------------------------
    if !ops.grammar.custom_check_legal_after_fit(indiv.node, prediction, data, ops)
        indiv.valid = 6
        return
    end

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
        indiv.valid = 7
    end

    return
end

""" Prints or displays an Individual.
"""
Base.show(io::IO, indiv::Individual) = println(io, "Indiv($(node_to_string(indiv.node, __operators; sigdigits=3)))")

""" The NSGA-II definition of isless.
"""
function Base.isless(i1::Individual, i2::Individual)
    if i1.rank != i2.rank
        return i1.rank < i2.rank
    else
        return i1.crowding > i2.crowding
    end
end
