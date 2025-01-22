
""" Various pre-implemented fit quality measures.
"""
function get_measure_max_ae(residual, residual_relative, prediction, data, node, ops)
    maximum(abs, residual)
end

function get_measure_mae(residual, residual_relative, prediction, data, node, ops)
    mean(abs, residual)
end

function get_measure_mse(residual, residual_relative, prediction, data, node, ops)
    mean(abs2, residual)
end

function get_measure_minus_r2(residual, residual_relative, prediction, data, node, ops)
    get_minus_r2(prediction, data[end])
end

function get_measure_minus_abs_spearman(residual, residual_relative, prediction, data, node, ops)
    get_minus_abs_spearman(prediction, data[end])
end

function get_measure_mare(residual, residual_relative, prediction, data, node, ops)
    mean(abs, residual_relative)
end

function get_measure_q75_are(residual, residual_relative, prediction, data, node, ops)
    quantile(abs.(residual_relative), 0.75)
end

function get_measure_max_are(residual, residual_relative, prediction, data, node, ops)
    maximum(abs, residual_relative)
end

function get_measure_ms_processed_e(residual, residual_relative, prediction, data, node, ops)
    mean(abs2, ops.fitting.residual_processing(residual, eachindex(residual), ops) .* ops.data_descript.fit_weights)
end

""" Various pre-implemented complexity measures.
"""
function get_measure_compl(residual, residual_relative, prediction, data, node, ops)
    count_nodes(node)
end

function get_measure_weighted_compl(residual, residual_relative, prediction, data, node, ops)
    get_weighted_compl(node, ops)
end

function get_measure_recursive_compl(residual, residual_relative, prediction, data, node, ops)
    recursive_compl(node, ops)
end

function get_measure_n_params(residual, residual_relative, prediction, data, node, ops)
    length(list_of_param_nodes(node))
end

""" Recursive complexity according to the dissertation of Kommenda 2018 except for the rule
    for ^, which was not addressed.
    Open TODO/caveats:
    - personally (VM), I think variables deserve less complexity than parameters
    - which complexity to use: { (lef * rig) + 1.0 } vs. { 2.0^(lef + rig) }
"""
function recursive_compl(node::Node, ops)
    if node.ari == 2
        compl_lef = recursive_compl(node.lef, ops)
        compl_rig = recursive_compl(node.rig, ops)

        if ops.binops[node.ind] in (+, -)
            compl = compl_lef + compl_rig + 1
        else
            compl = (compl_lef * compl_rig) + 1
        end

    elseif node.ari == 1
        compl_lef = recursive_compl(node.lef, ops)
        compl = compl_lef^2 + 1

    elseif node.ari == 0
        compl = 1.5

    elseif node.ari == -1
        compl = 2.0
    end

    return min(floatmax(), compl)
end

# # more nuanced recursive_compl
# global const recursive_compl_dict = Dict(
#     :+           => (lef, rig) -> lef + rig,
#     :-           => (lef, rig) -> lef + rig,
#     :*           => (lef, rig) -> (lef + rig)           + 1.0,
#     :/           => (lef, rig) -> (1.5 * lef + 2 * rig) + 1.0,
#     :^           => (lef, rig) -> (1.5 * lef + 3 * rig) + 1.0,
#     :pow_abs     => (lef, rig) -> (1.5 * lef + 3 * rig) + 1.0,
#     # :pow_integer => (lef, rig) -> (1.5 * lef          ) + 1.0,
#     # :pow_float   => (lef, rig) -> (1.5 * lef + 1.0    ) + 1.0,
# )
#
# function recursive_compl(node::Node, ops)
#     if node.ari == 2
#         compl_lef = recursive_compl(node.lef, ops)
#         compl_rig = recursive_compl(node.rig, ops)
#         cur_op = Symbol(ops.binops[node.ind])
#
#         compl = recursive_compl_dict[cur_op](compl_lef, compl_rig)
#
#     elseif node.ari == 1
#         compl_lef = recursive_compl(node.lef, ops)
#         compl = 5.0 * compl_lef
#
#     elseif node.ari == 0
#         compl = 1.0
#
#     elseif node.ari == -1
#         compl = 1.0
#     end
#
#     return min(floatmax(), compl)
# end

""" calculate the weighted_coml of a node. The weights for the functions and terminals are provided
    by ops.grammar.weighted_compl_dict. Weights for variables and parameters can be set using "VAR"
    and "PARAM", respectively. For any funciton of terminal not specified, 1.0 is assumed.
"""
function get_weighted_compl(node, ops) # TODO: test
    if node.ari == 2
        cur_fun = string(ops.binops[node.ind])
        compl   = get(ops.grammar.weighted_compl_dict, cur_fun, 1.0)
        return compl + get_weighted_compl(node.lef, ops) + get_weighted_compl(node.rig, ops)
    elseif node.ari == 1
        cur_fun = string(ops.unaops[node.ind])
        compl   = get(ops.grammar.weighted_compl_dict, cur_fun, 1.0)
        return compl + get_weighted_compl(node.lef, ops)
    elseif node.ari == 0
        return get(ops.grammar.weighted_compl_dict, "VAR", 1.0)
    elseif node.ari == -1
        return get(ops.grammar.weighted_compl_dict, "PARAM", 1.0)
    end
end

function get_minus_r2(pred, orig)
    total_sum_of_squares = sum(abs2, orig .- mean(orig))
    sum_of_squares_pred = sum(abs2, pred .- orig)
    return -(1 - sum_of_squares_pred / total_sum_of_squares)
end

function get_minus_abs_spearman(pred, orig)
    minus_abs_spear = -abs(corspearman(pred, orig))
    return isfinite(minus_abs_spear) ? minus_abs_spear : 0.0
end

