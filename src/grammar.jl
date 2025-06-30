
""" Check for illegal nestings in existing nodes. For it, ops.illegal_dict needs
    to be specified like thus:

    illegal_dict = Dict(
        "^" =>   (lef = (),             rig = ("+", "-", "*", "/", "VAR", )),
        "/" =>   (lef = (),             rig = ("-", "+")),
        "log" => (lef = ("log", "exp"), rig = ()),
        "exp" => (lef = ("exp", "log"), rig = ()),
        "cos" => (lef = ("cos",),       rig = ())
    )
"""
is_legal_nesting(node, ops)::Bool = _is_legal_nesting(node, ops, String[])::Bool

function _is_legal_nesting(node, ops, illegals)::Bool

    if node.ari == 2
        cur_fun = string(ops.binops[node.ind]) # TODO: maybe avoid conversion by adding str_binops to Options
        cur_fun in illegals && return false

        ill_nexts = get(ops.grammar.illegal_dict, cur_fun, (lef = (), rig = ()))
        illegals_rig = copy(illegals)
        append!(illegals, ill_nexts.lef)

        if !_is_legal_nesting(node.lef, ops, illegals)
            return false
        end

        append!(illegals_rig, ill_nexts.rig)

        return _is_legal_nesting(node.rig, ops, illegals_rig)

    elseif node.ari == 1
        cur_fun = string(ops.unaops[node.ind]) # TODO: maybe avoid conversion by adding str_unaops to Options
        cur_fun in illegals && return false

        ill_nexts = get(ops.grammar.illegal_dict, cur_fun, (lef = (), rig = ()))
        append!(illegals, ill_nexts.lef)

        return _is_legal_nesting(node.lef, ops, illegals)

    elseif node.ari == 0
        return !("VAR" in illegals)

    elseif node.ari == -1
        return !("PARAM" in illegals)

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

