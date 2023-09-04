# ==================================================================================================
# equation struct
# ==================================================================================================
mutable struct Node{T <: Number}
    ari::Int64
    val::T
    ind::Int64
    lef::Node{T}
    rig::Node{T}

    Node(v::T_; val_type::T_=1.0)                         where {T_ <: Number} = new{T_}(-1, v)
    Node(a::Int64, v::Number; val_type::T_=1.0)           where {T_ <: Number} = new{T_}(a, v)
    Node(a::Int64, v::Number, i::Int64; val_type::T_=1.0) where {T_ <: Number} = new{T_}(a, v, i)
    Node(i::Int64; val_type::T_=1.0)                      where {T_ <: Number} = new{T_}(0, 0, i)
    Node(a::Int64, op::Int64; val_type::T_=1.0)           where {T_ <: Number} = new{T_}(a, 0, op)
    Node(op::Int64, l::Node; val_type::T_=1.0)            where {T_ <: Number} = new{T_}(1, 0, op, l)
    Node(op::Int64, l::Node, r::Node; val_type::T_=1.0)   where {T_ <: Number} = new{T_}(2, 0, op, l, r)
end

""" Creates a new node with the specified type.
"""
function convert_node(node::Node, val_T::T) where {T <: Number}
    if node.ari == -1
        return Node(convert(T, node.val), val_type=val_T)
    elseif node.ari == 0
        return Node(node.ind, val_type=val_T)
    elseif node.ari == 1
        lef = convert_node(node.lef, val_T)
        return Node(node.ind, lef, val_type=val_T)
    elseif node.ari == 2
        lef = convert_node(node.lef, val_T)
        rig = convert_node(node.rig, val_T)
        return Node(node.ind, lef, rig, val_type=val_T)
    end
end

function Base.convert(::Type{T1}, node::Node{T2}) where {T1 <: Real, T2 <: Real}
    T1 == T2 && return node
    return convert_node(node, zero(T1))
end

function Base.convert(::Type{ForwardDiff.Dual{T, V, N}}, node::Node{T2}) where {T, V, N, T2<:Real}
    T1 = ForwardDiff.Dual{T, V, N}
    T1 == T2 && return node
    return convert_node(node, zero(T1))
end

# ==================================================================================================
# eval equation
# ==================================================================================================
@inline bad_array(arr) = !isfinite(sum(arr))

""" Eval a tree with the "cautious" approach -> prevent evaluations, which cause error, rather
    than redefining function like pow(abs(x), y). The faulty and each following one is prevented
    and a false is returned as the second return value.
    For now, the following evalulations are prevented:
    log(x <= 0), pow(x < 0, y), x / 0
"""
function eval_equation(node::Node{T}, data::AbstractArray, ops::Options)::Tuple{AbstractArray{T},Bool} where {T <: Real}
    if node.ari == -1
        return fill(node.val, size(data[1])), true
    elseif node.ari == 0
        return convert.(T, data[node.ind]), true
    elseif node.ari == 1
        arr_l, finite = eval_equation(node.lef, data, ops)
        (!finite || bad_array(arr_l)) && return arr_l, false

        if (ops.unaops[node.ind] == log) && any(l <= 0.0 for l in arr_l)
            return arr_l, false
        end

        arr_l .= ops.unaops[node.ind].(arr_l)
        return arr_l, true
    else
        arr_l, finite = eval_equation(node.lef, data, ops)
        (!finite || bad_array(arr_l)) && return arr_l, false

        arr_r, finite = eval_equation(node.rig, data, ops)
        (!finite || bad_array(arr_r)) && return arr_r, false

        if (ops.binops[node.ind] == /) && any(iszero(r) for r in arr_r)
            return arr_l, false
        end

        if (ops.binops[node.ind] == ^) && any(l < 0.0 || (iszero(l) && r < 0.0)
            for (l, r) in zip(arr_l, arr_r)
        )
            return arr_l, false
        end

        arr_l .= ops.binops[node.ind].(arr_l, arr_r)
        return arr_l, true
    end
end

# ==================================================================================================
# equation utilities
# ==================================================================================================
Base.show(io::IO, node::Node) = println(io, node_to_string(node, __operators; sigdigits=3, unique_nodes=true))

""" Convert a node to a string. This is also used to display nodes.
    - unique_nodes = true makes sense for directed acyclic graphs. Then, reused nodes are surrounded
      by curly braces {}
"""
node_to_string(node::Node, ops; sigdigits=15, unique_nodes=false, node_list=Node[]) = node_to_string_(node, ops, sigdigits, unique_nodes, node_list)[1]

function node_to_string_(node::Node, ops, sigdigits, unique_nodes, node_list)

    p_lef, p_rig = "(", ")"
    p_seen_lef, p_seen_rig = "", ""

    if unique_nodes
        if (node in node_list)
            p_seen_lef, p_seen_rig = "{", "}"
        else
            push!(node_list, node)
        end
    end


    if node.ari == 2
        op = string(ops.binops[node.ind])
        lef, node_list = node_to_string_(node.lef, ops, sigdigits, unique_nodes, node_list)
        rig, node_list = node_to_string_(node.rig, ops, sigdigits, unique_nodes, node_list)

        if op in ("+", "-", "*", "/", "^")
            str = p_seen_lef * p_lef * lef * op * rig * p_rig * p_seen_rig
        else
            str = p_seen_lef * op * p_lef * lef * ", " * rig * p_rig * p_seen_rig
        end

    elseif node.ari == 1
        op = string(ops.unaops[node.ind])
        lef, node_list = node_to_string_(node.lef, ops, sigdigits, unique_nodes, node_list)

        str = p_seen_lef * op * p_lef * lef * p_rig * p_seen_rig

    elseif node.ari == 0
        str = "v" * string(node.ind)

    elseif node.ari == -1
        str = string(round(extract_from_dual(node.val), sigdigits=sigdigits))
    end

    return str, node_list
end


""" Copy a node without using deepcopy -> 9/10 taken from SymbolicRegression.jl
"""
function copy_node(node::Node)::Node
    if node.ari == -1
        return Node(copy(node.val))
    elseif node.ari == 0
        return Node(copy(node.ind))
    elseif node.ari == 1
        return Node(copy(node.ind), copy_node(node.lef))
    elseif node.ari == 2
        return Node(copy(node.ind), copy_node(node.lef), copy_node(node.rig))
    end
end

""" Sets all fields of node1 to the values of node2. Is useful to exchange a node without losing
    the references pointing to it.
"""
function copy_node_wo_copy!(node1, node2)
    for field in fieldnames(typeof(node1))
        if isdefined(node2, field)
            setfield!(node1, field, deepcopy(getfield(node2, field)))
        end
    end
end

""" Count all sub-nodes of a node.
"""
function count_nodes(node::Node)::Float64
    if node.ari == 1
        return count_nodes(node.lef) + 1.0
    elseif node.ari == 2
        return count_nodes(node.lef) + count_nodes(node.rig) + 1.0
    else
        return 1.0
    end
end


""" Count all sub-nodes of a node, but do not count the same instance twice
    -> in case of directed acyclic graphs
"""
count_nodes_unique(node::Node; node_list=Node[]) = count_nodes_unique_(node, node_list)[1]

function count_nodes_unique_(node::Node, node_list)::Tuple{Float64,Vector{Node}}
    compl = 1.0
    if !(node in node_list)
        push!(node_list, node)
        if node.ari >= 1
            compl_l, node_list = count_nodes_unique_(node.lef, node_list)
            compl += compl_l
        end

        if node.ari == 2
            compl_r, node_list = count_nodes_unique_(node.rig, node_list)
            compl += compl_r
        end
    end
    return compl, node_list
end

""" Return the maximal depth of the tree. Interrupts in case the minim tree depth is reached.
"""
function maxim_tree_depth(node::Node; depth=0, minim=typemax(Int64))
    depth += 1
    node.ari <= 0 && return depth

    lef = maxim_tree_depth(node.lef, depth=depth, minim=minim)

    if lef < minim && node.ari == 2
        rig = maxim_tree_depth(node.rig, depth=depth, minim=minim)
        depth = max(lef, rig)
    else
        depth = lef
    end
    return depth
end


""" Compares two nodes and returns true if they are approx. same.
"""
function Base.isapprox(node1::Node, node2::Node; rtol=0.0)

    node1.ari == node2.ari || return false

    if node1.ari == -1
        Base.isapprox(node1.val, node2.val, rtol=rtol) || return false

    elseif node1.ari == 0
        node1.ind == node2.ind || return false

    elseif node1.ari >= 1
        node1.ind == node2.ind || return false
        Base.isapprox(node1.lef, node2.lef, rtol=rtol) || return false

        if node1.ari == 2
            node1.ind == node2.ind || return false
            Base.isapprox(node1.rig, node2.rig, rtol=rtol) || return false
        end
    end
    return true
end

""" Create an array containing the parameter nodes of a tree.
"""
function list_of_param_nodes(node; list=Node[])
    if node.ari >= 1
        list = list_of_param_nodes(node.lef; list=list)
        if node.ari == 2
            list = list_of_param_nodes(node.rig; list=list)
        end
    elseif node.ari == -1 && !(node in list)
        push!(list, node)
    end
    return list
end


""" Recursive functions using multiple dispatch to extract a value from arbitrarily
    nested dual numbers.
"""
extract_from_dual(v::AbstractFloat) = v
extract_from_dual(v::Number) = extract_from_dual(v.value)


