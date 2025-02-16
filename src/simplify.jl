# ==================================================================================================
# manual simplifications
# ==================================================================================================
""" Replaces unary function of parameters with parameters.
    Like: abs(param) -> param, cos(param) -> param
"""
function simplify_unary_of_param!(node)
    c = c_rig = c_lef = false
    if node.ari >= 1
        c_lef = simplify_unary_of_param!(node.lef)
        if node.ari == 2
            c_rig = simplify_unary_of_param!(node.rig)
        end
    end
    if node.ari == 1 && node.lef.ari == -1
        node.ari = node.lef.ari
        node.val = node.lef.val
        c = true
    end
    return c || c_lef || c_rig # TODO: test
end

""" Replaces a binary function of two parameters with a parameter.
    Like: param + param -> param, param * param -> param
"""
function simplify_binary_of_param!(node)
    c = c_rig = c_lef = false
    if node.ari >= 1
        c_lef = simplify_binary_of_param!(node.lef)
        if node.ari == 2
            c_rig = simplify_binary_of_param!(node.rig)
        end
    end
    if node.ari == 2 && node.lef.ari == -1 && node.rig.ari == -1
        node.ari = -1
        node.val = node.lef.val
        c = true
    end
    return c || c_lef || c_rig # TODO: test
end

""" Reorders the children of + and * in decreasing arity.
"""
function reorder_add_n_mul!(node, ops)
    c = c_rig = c_lef = false
    if node.ari >= 1
        c_lef = reorder_add_n_mul!(node.lef, ops)
        if node.ari == 2
            c_rig = reorder_add_n_mul!(node.rig, ops)
        end
    end
    if node.ari == 2 && ops.binops[node.ind] in (+, *)
        if node.rig.ari < node.lef.ari
            node.lef, node.rig = node.rig, node.lef
            c = true
        elseif node.rig.ari == node.lef.ari && node.lef.ari != -1
            if node.rig.ind < node.lef.ind
                node.lef, node.rig = node.rig, node.lef
                c = true
            end
        end
    end
    return c || c_lef || c_rig
end

""" Removes x - x and x / x and replaces by x, where x can be any subtree. Prevents domian errors
    in SymbolicUtils.
"""
function replace_same_subst_n_div!(node, ops)
    c = c_rig = c_lef = false
    if node.ari == 1
        c_lef = replace_same_subst_n_div!(node.lef, ops)
    elseif node.ari == 2
        if ops.binops[node.ind] in (-, /) && isapprox(node.lef, node.rig)
            copy_node_wo_copy!(node, node.lef)
            replace_same_subst_n_div!(node, ops)
            c = true
        else
            c_lef = replace_same_subst_n_div!(node.lef, ops)
            c_rig = replace_same_subst_n_div!(node.rig, ops)
        end
    end
    return c || c_lef || c_rig
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
    c = c_rig = c_lef = false
    if node.ari >= 1
        c_lef = simplify_binary_across_1_level!(node.lef, ops)
        if node.ari == 2
            c_rig = simplify_binary_across_1_level!(node.rig, ops)
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
                c = true
            end
        end
    end
    return c || c_lef || c_rig
end

""" Convert x / param to x * param -> better for drastic simplify.
"""
function div_to_mul_param!(node, ops) # TODO: test
    c = c_rig = c_lef = false
    if node.ari >= 1
        c_lef = div_to_mul_param!(node.lef, ops)
        if node.ari == 2
            c_rig = div_to_mul_param!(node.rig, ops)
        end
    end
    if node.ari == 2 && isequal(ops.binops[node.ind], /) && node.rig.ari == -1
        node.ind = findfirst(isequal(*), ops.binops)
        node.rig.val = 1 / node.rig.val
        c = true
    end
    return c || c_lef || c_rig
end

function apply_simple_simplifications!(node, ops)
    c1 = simplify_unary_of_param!(node)
    c2 = simplify_binary_of_param!(node)
    c3 = simplify_binary_across_1_level!(node, ops)
    c4 = replace_same_subst_n_div!(node, ops)
    (c1 || c2 || c3 || c4) && apply_simple_simplifications!(node, ops)
end

# ==================================================================================================
# simplify drasticly -> not neccessarily the same & used as a genetic operation
# ==================================================================================================
""" This simplification changes the node. It is used as a genetic operation and not a
    simplification. It removes parameters that are smaller than threshold.
    Like x + 1e-5 -> x, x * 1e-5 -> 1e-5
    In the latter example, the 1e-5 will again be detected and removed on the level above.
"""
function is_drastic_simplifyable(node, ops; threshold=1e-1)
    bool = false
    if node.ari >= 1
        bool = is_drastic_simplifyable(node.lef, ops, threshold=threshold)
        if node.ari == 2 && !bool
            bool = is_drastic_simplifyable(node.rig, ops, threshold=threshold)
        end
    end
    if !bool && node.ari == 2 && (node.lef.ari == -1 || node.rig.ari == -1) # if param either is a parameter
        op = ops.binops[node.ind]
        node_1       = getfield(node, node.lef.ari != -1 ? :lef : :rig)
        node_1_param = getfield(node, node.lef.ari != -1 ? :rig : :lef)
        if abs(node_1_param.val) < threshold           # if parameter 0.0
            if op in (+, -, *, ^, /)
                return true
            end
        elseif abs(1.0 - node_1_param.val) < threshold # if parameter 1.0
            if op in (*, ^) || ((op == /) && node.rig.ari == -1)
                return true
            end
        end
    end
    return bool
end

function get_drastic_simplify_nodes(node, ops; threshold=1e-1)
    nodes = Node[]
    if node.ari >= 1
        append!(nodes, get_drastic_simplify_nodes(node.lef, ops, threshold=threshold))
        if node.ari == 2
            append!(nodes, get_drastic_simplify_nodes(node.rig, ops, threshold=threshold))
        end
    end
    if node.ari == 2 && (node.lef.ari == -1 || node.rig.ari == -1) # if param is 0.0
        op = ops.binops[node.ind]
        node_1       = getfield(node, node.lef.ari != -1 ? :lef : :rig)
        node_1_param = getfield(node, node.lef.ari != -1 ? :rig : :lef)
        if abs(node_1_param.val) < threshold           # if parameter 0.0
            if op in (+, -, *, ^, /)
                return push!(nodes, node)
            end
        elseif abs(1.0 - node_1_param.val) < threshold # if parameter 1.0
            if op in (*, ^) || ((op == /) && node.rig.ari == -1)
                return push!(nodes, node)
            end
        end
    end
    return nodes
end

function drastic_simplify_!(node, ops; threshold=1e-1)
    op = ops.binops[node.ind]
    if (node.lef.ari == -1 || node.rig.ari == -1) # if param is 0.0
        node_1       = getfield(node, node.lef.ari != -1 ? :lef : :rig)
        node_1_param = getfield(node, node.lef.ari != -1 ? :rig : :lef)

        if abs(node_1_param.val) < threshold           # if parameter 0.0
            if op in (+, -)                        # x + 0.0 -> x
                copy_node_wo_copy!(node, node_1)
            elseif (op == *)                       # x * 0.0 -> 0.0
                copy_node_wo_copy!(node, node_1_param)
            elseif (op == ^) && node.rig.ari == -1 # x^0.0 -> 1.0
                node.rig.val = 1.0
                copy_node_wo_copy!(node, node.rig)
            elseif (op == ^) && node.lef.ari == -1 # 0.0^x -> 0
                node.lef.val = 0.0
                copy_node_wo_copy!(node, node.lef)
            elseif (op == /) && node.lef.ari == -1 # 0.0 / x -> 0.0
                copy_node_wo_copy!(node, node.lef)
            end
        elseif abs(1.0 - node_1_param.val) < threshold # if parameter 1.0
            if (op == *)                           # 1.0 * x -> x
                copy_node_wo_copy!(node, node_1)
            elseif (op == /) && node.rig.ari == -1 # x / 1.0 -> x
                copy_node_wo_copy!(node, node_1)
            elseif (op == ^) && node.rig.ari == -1 # x^1.0 -> x
                copy_node_wo_copy!(node, node.lef)
            elseif (op == ^) && node.lef.ari == -1 # 1.0^x -> 1.0
                copy_node_wo_copy!(node, node.lef)
            end
        end
    end
end

function drastic_simplify!(node, ops; threshold=1e-1, full=false) # TODO: can go in infinite loop
    while true
        drastic_nodes = get_drastic_simplify_nodes(node, ops; threshold=threshold)
        shuffle!(drastic_nodes)
        for n in drastic_nodes
            drastic_simplify_!(n, ops; threshold=threshold)
            full || break
        end
        full || break
        is_drastic_simplifyable(node, ops; threshold=threshold) || break
    end
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

""" Simplify using SymbolicUtils.
"""
function simplify_w_symbolic_utils!(node::Node, ops::Options)
    sym_eq = node_to_symbolic(node, ops)

    expr = SymbolicUtils.Code.toexpr(sym_eq)
    simp_node = string_to_node(expr, ops)

    reorder_add_n_mul!(simp_node, ops)
    copy_node_wo_copy!(node, simp_node)
end

