
# ==================================================================================================
# equation struct
# ==================================================================================================
mutable struct Node{T <: Number}
    ari::Int8
    val::T
    ind::Int8
    lef::Node{T}
    rig::Node{T}

    Node(v::T_;                         val_type::T_=1.0) where {T_ <: Number} = new{T_}(Int8(-1), T_(v))
    Node(i::Integer;                    val_type::T_=1.0) where {T_ <: Number} = new{T_}(Int8(0),  T_(0), Int8(i))
    Node(a::Integer,  op::Integer;      val_type::T_=1.0) where {T_ <: Number} = new{T_}(Int8(a),  T_(0), Int8(op))
    Node(op::Integer, l::Node;          val_type::T_=1.0) where {T_ <: Number} = new{T_}(Int8(1),  T_(0), Int8(op), l)
    Node(op::Integer, l::Node, r::Node; val_type::T_=1.0) where {T_ <: Number} = new{T_}(Int8(2),  T_(0), Int8(op), l, r)
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

function Base.hash(node::Node, h::UInt)
    if node.ari == -1
        val = sign(node.val) * round(log1p(abs(node.val)), sigdigits=2)
        return hash(node.ari, hash(val, h))
    elseif node.ari == 0
        return hash(node.ari, hash(node.ind, h))
    elseif node.ari == 1
        return hash(node.ari, hash(node.ind, hash(node.lef, h)))
    else
        return hash(node.ari, hash(node.ind, hash(node.lef, hash(node.rig, h))))
    end
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
function eval_equation(node::Node{T}, data::AbstractArray, ops::Options)::Tuple{AbstractArray{T}, Bool} where {T <: Real}
    if node.ari == -1
        return fill(node.val, length(data[1])), true
    elseif node.ari == 0
        return data[node.ind], true
    elseif node.ari == 1
        arr_l, finite = eval_equation(node.lef, data, ops)
        (!finite || bad_array(arr_l)) && return arr_l, false

        if (ops.unaops[node.ind] in (log, log10)) && any(l <= 0.0 for l in arr_l)
            return arr_l, false
        end

        if (ops.unaops[node.ind] == sqrt) && any(l < 0.0 for l in arr_l)
            return arr_l, false
        end

        arr_l = ops.unaops[node.ind].(arr_l)
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

        arr_l = ops.binops[node.ind].(arr_l, arr_r)
        return arr_l, true
    end
end

# ==================================================================================================
# equation utilities
# ==================================================================================================

""" How to display nodes.
"""
Base.show(io::IO, node::Node) = print(io, node_to_string(node, __operators; sigdigits=3))

""" Convert a node to a string. This is also used to display nodes.
"""
node_to_string(node::Node, ops; sigdigits=15) = node_to_string_(node, ops, sigdigits)

function node_to_string_(node::Node, ops, sigdigits)
    if node.ari == 2
        op = string(ops.binops[node.ind])
        lef = node_to_string_(node.lef, ops, sigdigits)
        rig = node_to_string_(node.rig, ops, sigdigits)

        if op in ("+", "-", "*", "/", "^")
            return "($lef$op$rig)"
        else
            return "$op($lef,$rig)"
        end

    elseif node.ari == 1
        op = string(ops.unaops[node.ind])
        lef = node_to_string_(node.lef, ops, sigdigits)
        return "$op($lef)"

    elseif node.ari == 0
        return "v$(string(node.ind))"

    elseif node.ari == -1
        return string(round(node.val, sigdigits=sigdigits))
    end
end

function encode_single_char(node::Node, ops)
    if node.ari == 2
        op_num = node.ind + length(ops.unaops)
        lef = encode_single_char(node.lef, ops)
        rig = encode_single_char(node.rig, ops)
        return "$(Char(96 + op_num))$lef$rig"
    elseif node.ari == 1
        op_num = node.ind
        lef = encode_single_char(node.lef, ops)
        return "$(Char(96 + op_num))$lef"
    elseif node.ari == 0
        return "V$(string(node.ind))"
    elseif node.ari == -1
        return string(round(sign(node.val)) * round(log1p(min(ops.general.replace_inf, abs(node.val))), sigdigits=2))
    end
end

""" Copy a node without using deepcopy -> 9/10 taken from SymbolicRegression.jl
"""
function Base.copy(node::Node)::Node
    if node.ari == -1
        return Node(copy(node.val))
    elseif node.ari == 0
        return Node(copy(node.ind))
    elseif node.ari == 1
        return Node(copy(node.ind), Base.copy(node.lef))
    elseif node.ari == 2
        return Node(copy(node.ind), Base.copy(node.lef), Base.copy(node.rig))
    end
end

Base.deepcopy(node::Node)::Node = Base.copy(node)

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

""" Compares two nodes and returns true if they are same.
"""
Base.:(==)(node1::Node, node2::Node) = Base.isapprox(node1, node2, rtol=0.0)

""" Return an array containing the parameter nodes of a tree.
"""
function list_of_param_nodes(node; list=Node[])
    if node.ari >= 1
        list = list_of_param_nodes(node.lef; list=list)
        if node.ari == 2
            list = list_of_param_nodes(node.rig; list=list)
        end
    elseif node.ari == -1
        push!(list, node)
    end
    return list
end

""" Search and overwrite nodes, which point to nodes that are no longer required.
"""
function clean_trash_nodes!(pop::Vector, null_node)
    for elem in pop
        clean_trash_nodes!(elem, null_node)
    end
end

clean_trash_nodes!(indiv, null_node) = clean_trash_nodes!(indiv.node, null_node)

function clean_trash_nodes!(node::Node{T}, null_node) where {T <: Number}
    if node.ari == 1
        clean_trash_nodes!(node.lef, null_node)
        if isdefined(node, :rig)
            node.rig = null_node
        end
    elseif node.ari == 2
        clean_trash_nodes!(node.lef, null_node)
        clean_trash_nodes!(node.rig, null_node)
    else
        if isdefined(node, :lef)
            node.lef = null_node
        end
        if isdefined(node, :rig)
            node.rig = null_node
        end
    end
end
