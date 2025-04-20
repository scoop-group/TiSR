mutable struct Trie
    children::Dict{Symbol, Trie}
    times_visited::Int8

    function Trie(
        children      = Dict{Symbol, Trie}(),
        times_visited = Int8(0),
    )
        return new(children, times_visited)
    end
end

function check_and_add!(trie, eq_str)
    if !(eq_str[1] in keys(trie.children))
        trie.children[eq_str[1]] = Trie()
    end
    trie.children[eq_str[1]].times_visited = min((trie.children[eq_str[1]].times_visited + 1), 127)
    if length(eq_str) > 1
        return check_and_add!(trie.children[eq_str[1]], eq_str[2:end])
    end
    return trie.times_visited - 1
end

function garbage_collect!(trie::Trie)
    for child in keys(trie.children)
        # trie.children[child].times_visited -= Int8(1)
        trie.children[child].times_visited ÷= Int8(2)
        if trie.children[child].times_visited <= 0
            delete!(trie.children, child)
        else
            garbage_collect!(trie.children[child])
        end
    end
end

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
        return [Symbol(round(sign(node.val)) * round(log1p(min(ops.general.replace_inf, abs(node.val))), sigdigits=2))]
    end
end

count_terminals(trie::Trie) = sum(
    isempty(s.children) ? 1 : count_terminals(s) for s in values(trie.children)
        ; init = 0
)

