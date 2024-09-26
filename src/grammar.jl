
""" Check for illegal nestings in existing nodes. For it, ops.illegal_dict needs 
    to be specified like thus:

    ops.grammar.illegal_dict = Dict(
        :^ => (lef = (), rig = (+, -, *, /, sin, cos)),
        :/ => (lef = (), rig = (*, /)),
    )
"""
function check_legal_function_nesting(node, ops)
    isempty(ops.grammar.illegal_dict) && return true

    unaweights = Float64[ops.p_unaops...]
    binweights = Float64[ops.p_binops...]

    return check_legal_function_nesting_(node, ops, unaweights, binweights)
end

function check_legal_function_nesting_(node, ops, unaweights, binweights)

    if node.ari == 2
        iszero(binweights[node.ind]) && return false

        cur_fun = Symbol(ops.binops[node.ind])
        ill_nexts = get(ops.grammar.illegal_dict, cur_fun, (lef = (), rig = ()))

        # left side
        ill_nexts_lef = ill_nexts.lef

        unaweights_lef = copy(unaweights)
        binweights_lef = copy(binweights)

        adjust_ops_weights!(unaweights_lef, binweights_lef, ill_nexts_lef, ops)
        check_legal_function_nesting_(node.lef, ops, unaweights_lef, binweights_lef) || return false

        # right side
        ill_nexts_rig = ill_nexts.rig

        unaweights_rig = copy(unaweights)
        binweights_rig = copy(binweights)

        adjust_ops_weights!(unaweights_rig, binweights_rig, ill_nexts_rig, ops)
        return check_legal_function_nesting_(node.rig, ops, unaweights_rig, binweights_rig)

    elseif node.ari == 1
        iszero(unaweights[node.ind]) && return false

        cur_fun = Symbol(ops.unaops[node.ind])
        ill_nexts = get(ops.grammar.illegal_dict, cur_fun, (lef = (),))

        unaweights_lef = copy(unaweights)
        binweights_lef = copy(binweights)

        adjust_ops_weights!(unaweights_lef, binweights_lef, ill_nexts.lef, ops)

        return check_legal_function_nesting_(node.lef, ops, unaweights_lef, binweights_lef)
    end
    return true
end

""" Adjust the weight-vectors for the operator selection probability of the grow_equation and
    check_legal_function_nesting function to prevent illegal nestings according to illegal_dict.
"""
function adjust_ops_weights!(unaweights, binweights, ill_nexts, ops)
    for fun in ill_nexts
        ind = findfirst(isequal(fun), ops.unaops)
        if !isnothing(ind)
            unaweights[ind] = 0.0
        end

        ind = findfirst(isequal(fun), ops.binops)
        if !isnothing(ind)
            binweights[ind] = 0.0
        end
    end
end

""" Return the maximal number of nodes across all top-level terms.
"""
function get_max_nodes_per_term(node, ops)
    if node.ari == 2 && ops.binops[node.ind] in (+, -)
        return max(
            get_max_nodes_per_term(node.lef, ops),
            get_max_nodes_per_term(node.rig, ops)
        )
    else
        return count_nodes(node)
    end
end

""" Trim the node according to max terms per node.
"""
function trim_to_max_nodes_per_term!(node, ops)
    if node.ari == 2 && ops.binops[node.ind] in (+, -)
        trim_to_max_nodes_per_term!(node.lef, ops)
        trim_to_max_nodes_per_term!(node.rig, ops)
    else
        trim_to_max_compl!(node, ops.grammar.max_nodes_per_term, ops)
    end
end
