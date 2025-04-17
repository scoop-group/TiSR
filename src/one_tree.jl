mutable struct State
    times_visited::Float64
    children::Dict{Symbol, State}
    lock::ReentrantLock

    function State(
        times_visited = 0.,
        children      = Dict{Symbol, State}(),
    )
        return new(times_visited, children, ReentrantLock())
    end
end

expr_to_prefix(expr::Expr) = expr_to_prefix(expr.args) # TODO: maybe just have a prefix option in TiSR?
expr_to_prefix(expr::T) where {T <: Number} = [:p]
expr_to_prefix(expr) = [expr]
function expr_to_prefix(expr::Vector)
    arr = Symbol[]
    for (i, ex) in enumerate(expr)
        if i > 1 && length(expr) - i > 1  # required to convert n-ary to binary
            append!(arr, expr_to_prefix(expr[1]))
        end
        ret = expr_to_prefix(ex)
        append!(arr, ret)
    end
    return arr
end

function feed_the_tree!(state, prefix_eq)
    lock(state.lock) do
        state.times_visited += 1
    end

    if !(prefix_eq[1] in keys(state.children))
        lock(state.lock) do
            state.children[prefix_eq[1]] = State()
        end
    end

    if length(prefix_eq) > 1
        feed_the_tree!(state.children[prefix_eq[1]], prefix_eq[2:end])
    elseif length(prefix_eq) == 1
        lock(state.lock) do
            state.children[prefix_eq[1]].times_visited += 1
        end
    end
end

import Base: in

function Base.in(prefix_eq::Vector{Symbol}, state::State)
    prefix_eq[1] in keys(state.children) || return false
    length(prefix_eq) == 1               && return true
    return Base.in(state.children[prefix_eq[1]], prefix_eq[2:end])
end

function check_and_add!(state, prefix_eq)
    unseen = false
    lock(state.lock) do
        state.times_visited += 1
    end

    if !(prefix_eq[1] in keys(state.children))
        lock(state.lock) do
            state.children[prefix_eq[1]] = State()
        end
        unseen = true
    end

    if length(prefix_eq) > 1
        return unseen | check_and_add!(state.children[prefix_eq[1]], prefix_eq[2:end])
    elseif length(prefix_eq) == 1
        lock(state.lock) do
            state.children[prefix_eq[1]].times_visited += 1
        end
    end
    return unseen
end

function garbage_collect!(state::State; forget_rate = 1)
    for child in keys(state.children)
        state.children[child].times_visited -= forget_rate
        if state.children[child].times_visited <= 0
            delete!(state.children, child)
        else
            garbage_collect!(state.children[child], forget_rate = forget_rate)
        end
    end
end

count_terminals(state::State) = sum(
    isempty(s.children) ? 1 : count_terminals(s) for s in values(state.children)
)

function node_to_string_prefix(node::Node, ops)
    if node.ari == 2
        op = Symbol(ops.binops[node.ind])
        lef = node_to_string_prefix(node.lef, ops)
        rig = node_to_string_prefix(node.rig, ops)
        return [op, lef..., rig...]

    elseif node.ari == 1
        op = Symbol(ops.unaops[node.ind])
        lef = node_to_string_prefix(node.lef, ops)
        return [op, lef...]

    elseif node.ari == 0
        return [Symbol("v$(string(node.ind))")]

    elseif node.ari == -1
        return [:p]
    end
end


