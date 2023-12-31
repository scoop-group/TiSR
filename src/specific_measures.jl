
""" Recursive complexity according to the dissertation of Kommenda 2018 except for the rule
    for ^, which was not addressed.
    Open TODO/caveats:
    - personally (VM), I think variables deserve less complexity than parameters
    - which complexity to use: { (lef * rig) + 1.0 } vs. { 2.0^(lef + rig) }
"""
global const recursive_compl_dict = Dict(
        :+       => (lef, rig) -> lef + rig,
        :-       => (lef, rig) -> lef + rig,
        :*       => (lef, rig) -> (lef * rig) + 1.0,
        :/       => (lef, rig) -> (lef * rig) + 1.0,
        :^       => (lef, rig) -> (lef * rig) + 1.0, # 2.0^(lef + rig) ?
        :pow_abs => (lef, rig) -> (lef * rig) + 1.0, # 2.0^(lef + rig) ?
    )

function recursive_compl(node::Node, ops)
    if node.ari == 2
        compl_lef = recursive_compl(node.lef, ops)
        compl_rig = recursive_compl(node.rig, ops)
        cur_op = Symbol(ops.binops[node.ind])

        compl = recursive_compl_dict[cur_op](compl_lef, compl_rig)

    elseif node.ari == 1
        compl_lef = recursive_compl(node.lef, ops)
        compl = 2.0^compl_lef

    elseif node.ari == 0
        compl = 2.0

    elseif node.ari == -1
        compl = 1.0
    end

    return min(floatmax(), compl)
end

function r_squared(pred, orig)
    total_sum_of_squares = sum(abs2, orig .- mean(orig))
    sum_of_squares_pred = sum(abs2, pred .- orig)
    1 - sum_of_squares_pred / total_sum_of_squares
end
