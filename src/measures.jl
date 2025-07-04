
""" Various pre-implemented fit quality measures.
"""
function get_measure_max_ae(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[end]
    residual = target .- prediction
    @views maximum(abs, residual[inds])
end

function get_measure_mae(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[end]
    residual = target .- prediction
    @views mean(abs, residual[inds])
end

function get_measure_mse(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[end]
    residual = target .- prediction
    @views mean(abs2, residual[inds])
end

function get_measure_one_minus_r2(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[end]
    @views get_one_minus_r2(prediction[ind], target[ind])
end

function get_measure_one_minus_abs_spearman(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[end]
    @views get_one_minus_abs_spearman(prediction[inds], target[inds])
end

function get_measure_mare(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    guard = inv(ops.general.replace_inf)
    inds = ops.data_descript.split_inds[end]
    rel_residual = @. (prediction - target) / (abs(target) + guard)
    @views mean(abs, rel_residual[inds])
end

function get_measure_q75_are(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    guard = inv(ops.general.replace_inf)
    inds = ops.data_descript.split_inds[end]
    rel_residual = @. (prediction - target) / (abs(target) + guard)
    @views quantile(abs.(rel_residual[inds]), 0.75)
end

function get_measure_max_are(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    guard = inv(ops.general.replace_inf)
    inds = ops.data_descript.split_inds[end]
    rel_residual = @. (prediction - target) / (abs(target) + guard)
    maximum(abs, rel_residual[inds])
end

function get_measure_ms_processed_e(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[end]
    residual = target .- prediction
    @views mean(abs2, ops.fitting.residual_processing(residual, eachindex(residual), ops)[inds] # TODO: are inds and eachindex(residual) right?
                            .* ops.data_descript.fit_weights[inds])
end

""" Various pre-implemented measures -> test versions.
"""
function get_measure_max_ae_test(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[3]
    residual = target .- prediction
    @views maximum(abs, residual[inds])
end

function get_measure_mae_test(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[3]
    residual = target .- prediction
    @views mean(abs, residual[inds])
end

function get_measure_mse_test(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[3]
    residual = target .- prediction
    @views mean(abs2, residual[inds])
end

function get_measure_one_minus_r2_test(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[3]
    @views get_one_minus_r2(prediction[ind], target[ind])
end

function get_measure_one_minus_abs_spearman_test(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[3]
    @views get_one_minus_abs_spearman(prediction[inds], target[inds])
end

function get_measure_mare_test(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    guard = inv(ops.general.replace_inf)
    inds = ops.data_descript.split_inds[3]
    rel_residual = @. (prediction - target) / (abs(target) + guard)
    @views mean(abs, rel_residual[inds])
end

function get_measure_q75_are_test(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    guard = inv(ops.general.replace_inf)
    inds = ops.data_descript.split_inds[3]
    rel_residual = @. (prediction - target) / (abs(target) + guard)
    @views quantile(abs.(rel_residual[inds]), 0.75)
end

function get_measure_max_are_test(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    guard = inv(ops.general.replace_inf)
    inds = ops.data_descript.split_inds[3]
    rel_residual = @. (prediction - target) / (abs(target) + guard)
    maximum(abs, rel_residual[inds])
end

function get_measure_ms_processed_e_test(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    inds = ops.data_descript.split_inds[3]
    residual = target .- prediction
    @views mean(abs2, ops.fitting.residual_processing(residual, eachindex(residual), ops)[inds]
                            .* ops.data_descript.fit_weights[inds])
end

""" Various pre-implemented complexity measures.
"""
function get_measure_compl(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    count_nodes(node)
end
function get_measure_weighted_compl(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    get_weighted_compl(node, ops)
end
function get_measure_n_params(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    length(list_of_param_nodes(node))
end
function get_measure_recursive_compl(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    recursive_compl(node, ops)
end

function get_measure_max_nodes_per_term(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    get_max_nodes_per_term(node, ops)
end

function get_measure_square_compl(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    get_max_nodes_per_term(node, ops) * get_weighted_compl(node, ops)
end

function get_measure_cross_compl(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    get_max_nodes_per_term(node, ops) + get_weighted_compl(node, ops)
end

""" calculate the weighted_coml of a node. The weights for the functions and terminals are provided
    by ops.grammar.weighted_compl_dict. Weights for variables and parameters can be set using "VAR"
    and "PARAM", respectively. For any funciton or terminal not specified, 3.0 is assumed.
"""
function get_weighted_compl(node, ops)::Float64
    if node.ari == 2
        cur_fun = string(ops.binops[node.ind])
        compl   = get(ops.grammar.weighted_compl_dict, cur_fun, 3.0)
        return compl + get_weighted_compl(node.lef, ops) + get_weighted_compl(node.rig, ops)
    elseif node.ari == 1
        cur_fun = string(ops.unaops[node.ind])
        compl   = get(ops.grammar.weighted_compl_dict, cur_fun, 3.0)
        return compl + get_weighted_compl(node.lef, ops)
    elseif node.ari == 0
        return get(ops.grammar.weighted_compl_dict, "VAR", 3.0)
    elseif node.ari == -1
        return get(ops.grammar.weighted_compl_dict, "PARAM", 3.0)
    end
end

function get_one_minus_r2(pred::Vector{T}, orig::Vector{T})::T where {T}
    total_sum_of_squares = sum(abs2, orig .- mean(orig))
    sum_of_squares_pred = sum(abs2, pred .- orig)
    return 1.0-(1.0 - sum_of_squares_pred / total_sum_of_squares)
end

function get_one_minus_abs_spearman(pred::AbstractArray{T}, orig::AbstractArray{T})::T where {T}
    one_minus_abs_spear = 1.0-abs(corspearman(pred, orig))
    return isfinite(one_minus_abs_spear) ? one_minus_abs_spear : 1.0
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

function recursive_compl(
    node::Node,
    ops;
    recursive_compl_dict = Dict(
        :+   => (lef, rig) -> lef + rig,
        :-   => (lef, rig) -> lef + rig,
        :*   => (lef, rig) -> (lef + rig)           + 1.0,
        :/   => (lef, rig) -> (1.5 * lef + 2 * rig) + 1.0,
        :^   => (lef, rig) -> (1.5 * lef + 3 * rig) + 1.0,
        :pow => (lef, rig) -> (1.5 * lef + 3 * rig) + 1.0,
    ),
    una_nest_pen = 3,
)
    if node.ari == 2
        compl_lef = recursive_compl(
            node.lef,
            ops,
            recursive_compl_dict = recursive_compl_dict,
            una_nest_pen = una_nest_pen
        )
        compl_rig = recursive_compl(
            node.rig,
            ops,
            recursive_compl_dict = recursive_compl_dict,
            una_nest_pen = una_nest_pen
        )
        cur_op = Symbol(ops.binops[node.ind])

        cur_fun = get(recursive_compl_dict, cur_op, (lef, rig) -> 3 * lef + 3 * rig + 1.0)
        compl = cur_fun(compl_lef, compl_rig)

    elseif node.ari == 1
        compl_lef = recursive_compl(
            node.lef,
            ops,
            recursive_compl_dict = recursive_compl_dict,
            una_nest_pen = una_nest_pen
        )
        compl = una_nest_pen * compl_lef

    elseif node.ari == 0
        compl = 1.0

    elseif node.ari == -1
        compl = 1.0 # + abs(round(node.val) - node.val) # non-integer penalty
    end

    return min(floatmax(), compl)
end

function get_measure_constr_vios(prediction::Vector{T}, target::Vector{T}, node, ops)::T where {T}
    isempty(ops.fitting.all_constr_f_select) && return 0.0

    function node_func(params::Vector{T1}, constr_data::Vector{Vector{T2}}) where {T1, T2} # evaluate a node at given parameters and given data
        TT = promote_type(T1, T2)
        node_ = convert(TT, node)
        pred, _ = eval_equation(node_, constr_data, ops)
        pred = ops.fitting.pre_residual_processing(pred, eachindex(constr_data[1]), ops)
        return pred
    end

    params = Float64[n.val for n in list_of_param_nodes(node)]
    return sum(abs, sum(abs2, f(node_func, params)) for f in ops.fitting.all_constr_f_select; init = 0.0)
end

