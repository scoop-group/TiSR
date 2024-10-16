# ==================================================================================================
# manual simplifications
# ==================================================================================================
""" Replaces unary function of parameters with parameters.
    Like: abs(param) -> param, cos(param) -> param
"""
function simplify_unary_of_param!(node)
    if node.ari >= 1
        simplify_unary_of_param!(node.lef)
        if node.ari == 2
            simplify_unary_of_param!(node.rig)
        end
    end

    if node.ari == 1 && node.lef.ari == -1
        node.ari = node.lef.ari
        node.val = node.lef.val
    end
end

""" Replaces a binary function of two parameters with a parameter.
    Like: param + param -> param, param * param -> param
"""
function simplify_binary_of_param!(node)
    if node.ari >= 1
        simplify_binary_of_param!(node.lef)
        if node.ari == 2
            simplify_binary_of_param!(node.rig)
        end
    end

    if node.ari == 2 && node.lef.ari == -1 && node.rig.ari == -1
        node.ari = -1
        node.val = node.lef.val
    end
end

""" Reorders the children of + and * in decreasing arity.
"""
function reorder_add_n_mul!(node, ops)
    if node.ari >= 1
        reorder_add_n_mul!(node.lef, ops)
        if node.ari == 2
            reorder_add_n_mul!(node.rig, ops)
        end
    end

    if node.ari == 2 && ops.binops[node.ind] in (+, *)
        if node.rig.ari < node.lef.ari
            node.lef, node.rig = node.rig, node.lef
        elseif node.rig.ari == node.lef.ari && node.lef.ari != -1
            if node.rig.ind < node.lef.ind
                node.lef, node.rig = node.rig, node.lef
            end
        end
    end
end

""" Removes x - x and x / x and replaces by x, where x can be any subtree. Prevents domian errors
    in SymbolicUtils.
"""
function replace_same_subst_n_div!(node, ops)
    if node.ari == 1
        replace_same_subst_n_div!(node.lef, ops)
    elseif node.ari == 2
        if ops.binops[node.ind] in (-, /) && isapprox(node.lef, node.rig)
            copy_node_wo_copy!(node, node.lef)
            replace_same_subst_n_div!(node, ops)
        else
            replace_same_subst_n_div!(node.lef, ops)
            replace_same_subst_n_div!(node.rig, ops)
        end
    end
end

""" Simplifies redunand operations across one tree level.
    Like: (x + param) + param -> x + param, (x * param) / param -> x * param
    Idea from SymbolicRegression.jl
"""
global const simplify_binary_across_1_level_dict = Dict(
    :+ => (+, -),
    :- => (+, -),
    :* => (*, /),
    :/ => (*, /),
)

function simplify_binary_across_1_level!(node, ops)
    if node.ari >= 1
        simplify_binary_across_1_level!(node.lef, ops)
        if node.ari == 2
            simplify_binary_across_1_level!(node.rig, ops)
        end
    end

    if node.ari == 2 && xor(node.lef.ari == -1, node.rig.ari == -1) && xor(node.lef.ari == 2, node.rig.ari == 2)
        node_1 = getfield(node, node.lef.ari == 2 ? :lef : :rig)
        node_1_param = getfield(node, node.lef.ari == 2 ? :rig : :lef)

        par_op_str = Symbol(ops.binops[node.ind])
        chi_op_str = ops.binops[node_1.ind]

        if par_op_str in keys(simplify_binary_across_1_level_dict) && chi_op_str in simplify_binary_across_1_level_dict[par_op_str]
            if xor(node_1.lef.ari == -1, node_1.rig.ari == -1)
                node_2       = getfield(node_1, node_1.lef.ari == -1 ? :rig : :lef)
                node_2_param = getfield(node_1, node_1.lef.ari == -1 ? :lef : :rig)

                setfield!(node, :lef, node_1_param.ari < node_2_param.ari ? node_1_param : node_2_param)
                setfield!(node, :rig, node_2)
            end
        end
    end
end

""" Convert x / param to x * param -> better for drastic simplify.
"""
function div_to_mul_param!(node, ops) # TODO: test
    if node.ari >= 1
        div_to_mul_param!(node.lef, ops)
        if node.ari == 2
            div_to_mul_param!(node.rig, ops)
        end
    end

    if node.ari == 2 && isequal(ops.binops[node.ind], /) && node.rig.ari == -1
        node.ind = findfirst(isequal(*), ops.binops)
        node.rig.val = 1 / node.rig.val
    end
end

function apply_simple_simplifications!(node, ops)
    str1 = node_to_string(node, ops) # TODO: avoid going through strings -> return Bool from simplifes
    while true
        simplify_unary_of_param!(node)
        simplify_binary_of_param!(node)
        simplify_binary_across_1_level!(node, ops)
        replace_same_subst_n_div!(node, ops)
        str2 = node_to_string(node, ops)
        str1 == str2 && break
        str1 = str2
    end
end

# ==================================================================================================
# simplify drasticly -> not neccessarily the same & used as a genetic operation
# ==================================================================================================
""" This simplification changes the node. It is used as a genetic operation and not a
    simplification. It removes parameters that are smaller than threshold.
    Like x + 1e-5 -> x, x * 1e-5 -> 1e-5
    In the latter example, the 1e-5 will again be detected and removed on the level above.
"""
function drastic_simplify!(node, ops; threshold=1e-1, potential=false, prob=0.5)

    pot = false

    if node.ari >= 1
        pot_l = drastic_simplify!(node.lef, ops, threshold=threshold, potential=potential, prob=prob)
        pot = pot || pot_l
        if node.ari == 2
            pot_r = drastic_simplify!(node.rig, ops, threshold=threshold, potential=potential, prob=prob)
            pot = pot || pot_r
        end
    end

    node.ari == 2 || return pot

    op = ops.binops[node.ind]

    if (node.lef.ari == -1 || node.rig.ari == -1)
        node_1 = getfield(node, node.lef.ari != -1 ? :lef : :rig)
        node_1_param = getfield(node, node.lef.ari != -1 ? :rig : :lef)

        if abs(node_1_param.val) < threshold # if parameter 0.0
            if op in (+, -) # x + 0.0 -> x
                potential && return true
                if rand() < prob
                    copy_node_wo_copy!(node, node_1)
                end
            elseif (op == *) # x * 0.0 -> 0.0
                potential && return true
                if rand() < prob
                    copy_node_wo_copy!(node, node_1_param)
                end
            elseif (op == ^) && node.rig.ari == -1 # x^0.0 -> 1.0
                potential && return true
                if rand() < prob
                    node.rig.val = 1.0
                    copy_node_wo_copy!(node, node.rig)
                end
            elseif (op == ^) && node.lef.ari == -1 # 0.0^x -> 0 # TODO: new, test
                potential && return true
                if rand() < prob
                    node.lef.val = 0.0
                    copy_node_wo_copy!(node, node.lef)
                end
            elseif (op == /) && node.lef.ari == -1 # 0.0 / x -> 0.0
                potential && return true
                if rand() < prob
                    copy_node_wo_copy!(node, node.lef)
                end
            end
        elseif abs(1.0 - node_1_param.val) < threshold # if parameter 1.0
            if (op == *) # 1.0 * x -> x
                potential && return true
                if rand() < prob
                    copy_node_wo_copy!(node, node_1)
                end
            elseif (op == /) && node.rig.ari == -1 # x / 1.0 -> x
                if rand() < prob
                    copy_node_wo_copy!(node, node_1)
                end
                potential && return true
            elseif (op == ^) && node.rig.ari == -1 # x^1.0 -> x
                potential && return true
                if rand() < prob
                    copy_node_wo_copy!(node, node.lef)
                end
            elseif (op == ^) && node.lef.ari == -1 # 1.0^x -> 1.0 # TODO: test
                potential && return true
                if rand() < prob
                    copy_node_wo_copy!(node, node.lef)
                end
            end
        elseif (op == /) && node.rig.ari == -1 && abs(1 / node_1_param.val) < threshold
            potential && return true
            node_1_param.val = 0.0
            if rand() < prob
                copy_node_wo_copy!(node, node_1_param)
            end
        end
    end
    return pot || false
end

# ==================================================================================================
# SymbolicUtils simplify -> used as a genetic operation
# ==================================================================================================
""" Converts a node to a SymbolicUtils expression. Partially inspired by SymbolicRegression.jl
"""
function node_to_symbolic(node::Node, ops::Options)
    if node.ari == -1
        return node.val
    elseif node.ari == 0
        return SymbolicUtils.Sym{Real}(Symbol("v$(node.ind)"))
    elseif node.ari == 1
        lef = node_to_symbolic(node.lef, ops)
        return ops.unaops[node.ind](lef)
    elseif node.ari == 2
        lef = node_to_symbolic(node.lef, ops)
        rig = node_to_symbolic(node.rig, ops)
        return ops.binops[node.ind](lef, rig)
    end
end

""" Simplify using SymbolicUtils. If through_polyform=true, the expressions converted into
    polynomial form and further simplified.
"""
function simplify_w_symbolic_utils!(node::Node, ops::Options; use_simplify=false, though_polyform=false)
    sym_eq = node_to_symbolic(node, ops)

    if use_simplify
        sym_eq = SymbolicUtils.simplify(sym_eq)
        if though_polyform
            sym_eq = SymbolicUtils.PolyForm(sym_eq)
            sym_eq = SymbolicUtils.simplify(sym_eq)
        end
    end

    expr = SymbolicUtils.Code.toexpr(sym_eq)
    simp_node = string_to_node(expr, ops)

    reorder_add_n_mul!(simp_node, ops)
    copy_node_wo_copy!(node, simp_node)
end

""" 
"""
function simplify_to_string(node::Node, ops::Options; sigdigits=15)
    sym_eq = node_to_symbolic(node, ops)
    # sym_eq = SymbolicUtils.simplify(sym_eq)
    eq_str = string(sym_eq)
    if sigdigits >= 15
        return eq_str
    else
        return round_equation_string(eq_str, sigdigits=sigdigits)
    end
end

