
""" Various pre-implemented fit quality measures.
"""
function get_measure_max_ae(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[end]
    residual = target .- prediction
    @views maximum(abs, residual[inds])
end

function get_measure_mae(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[end]
    residual = target .- prediction
    @views mean(abs, residual[inds])
end

function get_measure_mse(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[end]
    residual = target .- prediction
    @views mean(abs2, residual[inds])
end

function get_measure_one_minus_r2(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[end]
    @views get_one_minus_r2(prediction[ind], target[ind])
end

function get_measure_one_minus_abs_spearman(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[end]
    @views get_one_minus_abs_spearman(prediction[inds], target[inds])
end

function get_measure_mare(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[end]
    rel_residual = @. (prediction - target) / (abs(target) + 1e-12)
    @views mean(abs, rel_residual[inds])
end

function get_measure_q75_are(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[end]
    rel_residual = @. (prediction - target) / (abs(target) + 1e-12)
    @views quantile(abs.(rel_residual[inds]), 0.75)
end

function get_measure_max_are(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[end]
    rel_residual = @. (prediction - target) / (abs(target) + 1e-12)
    maximum(abs, rel_residual[inds])
end

function get_measure_ms_processed_e(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[end]
    residual = target .- prediction
    @views mean(abs2, ops.fitting.residual_processing(residual, eachindex(residual), ops)[inds]
                            .* ops.data_descript.fit_weights[inds])
end

""" Various pre-implemented measures -> test versions.
"""
function get_measure_max_ae_test(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[3]
    residual = target .- prediction
    @views maximum(abs, residual[inds])
end

function get_measure_mae_test(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[3]
    residual = target .- prediction
    @views mean(abs, residual[inds])
end

function get_measure_mse_test(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[3]
    residual = target .- prediction
    @views mean(abs2, residual[inds])
end

function get_measure_one_minus_r2_test(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[3]
    @views get_one_minus_r2(prediction[ind], target[ind])
end

function get_measure_one_minus_abs_spearman_test(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[3]
    @views get_one_minus_abs_spearman(prediction[inds], target[inds])
end

function get_measure_mare_test(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[3]
    rel_residual = @. (prediction - target) / (abs(target) + 1e-12)
    @views mean(abs, rel_residual[inds])
end

function get_measure_q75_are_test(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[3]
    rel_residual = @. (prediction - target) / (abs(target) + 1e-12)
    @views quantile(abs.(rel_residual[inds]), 0.75)
end

function get_measure_max_are_test(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[3]
    rel_residual = @. (prediction - target) / (abs(target) + 1e-12)
    maximum(abs, rel_residual[inds])
end

function get_measure_ms_processed_e_test(prediction, target, node, ops)
    inds = ops.data_descript.split_inds[3]
    residual = target .- prediction
    @views mean(abs2, ops.fitting.residual_processing(residual, eachindex(residual), ops)[inds]
                            .* ops.data_descript.fit_weights[inds])
end


""" Various pre-implemented complexity measures.
"""
get_measure_compl(prediction, target, node, ops)           = count_nodes(node)
get_measure_weighted_compl(prediction, target, node, ops)  = get_weighted_compl(node, ops)
get_measure_recursive_compl(prediction, target, node, ops) = recursive_compl(node, ops)
get_measure_n_params(prediction, target, node, ops)        = length(list_of_param_nodes(node))

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

function get_one_minus_r2(pred, orig)
    total_sum_of_squares = sum(abs2, orig .- mean(orig))
    sum_of_squares_pred = sum(abs2, pred .- orig)
    return 1.0-(1.0 - sum_of_squares_pred / total_sum_of_squares)
end

function get_one_minus_abs_spearman(pred, orig)
    one_minus_abs_spear = 1.0-abs(corspearman(pred, orig))
    return isfinite(one_minus_abs_spear) ? one_minus_abs_spear : 1.0
end

