
""" Set the values of the parameter in a node.
"""
function set_params!(list_of_param, x)
    for i in eachindex(x)
        list_of_param[i].val = x[i]
    end
end

""" Returns the evalualted function for a given parameter vector.
"""
function fitting_eval(x, node, list_of_param, data, inds, ops)
    set_params!(list_of_param, x)

    dat = [view(d, inds) for d in data]
    pred, _ = eval_equation(node, dat, ops)
    pred = ops.fitting.pre_residual_processing(pred, inds, ops)
    return pred
end

""" Return the fitting residual. Pre-residual_processing! funciton is applied in-place
"""
function fitting_residual(x, node, list_of_param, data, inds, ops)
    pred = fitting_eval(x, node, list_of_param, data, inds, ops)
    view(data[end], inds) .- pred
end

""" Calculate the residual and apply pre- and post-processing.
"""
function fitting_objective(x::Vector{T}, node, data, ops)::Vector{T} where {T <: Number}
    data_ = data#[convert(Vector{T}, d) for d in data]
    node_ = convert(T, node)
    list_of_param_ = list_of_param_nodes(node_)

    res   = fitting_residual(x, node_, list_of_param_, data_, ops.data_descript.split_inds[1], ops)
    res  .= ops.fitting.residual_processing(res, ops.data_descript.split_inds[1], ops)
    res .*= view(ops.data_descript.fit_weights, ops.data_descript.split_inds[1])
    res  .= abs.(res)
    return res
end

""" Entrance/switch for parameter fitting. Depending on whether the equation has parameters
    or not, the correspoinding functions are called, and the parameters adapted in-place.
"""
function fit_n_eval!(node, data, ops; do_fit = true, do_constr = true)

    list_of_param = list_of_param_nodes(node)

    if length(list_of_param) > length(ops.data_descript.split_inds[1])
        return data[end], false
    end

    if do_fit && !isempty(list_of_param)
        if rand() < ops.fitting.NM_prob
            fitting_Optim!(
                node, data, ops, list_of_param, ops.fitting.NM_iter,
                Optim.NelderMead()
            )
        else
            fitting_LM!(node, data, ops, list_of_param, ops.fitting.max_iter)
        end

        if do_constr && (!isempty(ops.fitting.eq_constr) || !isempty(ops.fitting.ineq_constr))
            pred, valid = eval_equation(node, data, ops)
            valid || return pred, valid

            pred = ops.fitting.pre_residual_processing(pred, eachindex(pred), ops)
            mare  = mean(abs, (pred .- data[end]) ./ (data[end] .+ inv(ops.general.replace_inf)))

            if mare < ops.fitting.max_mare_for_constr_fit
                fitting_w_constr!(node, data, ops, list_of_param, ops.fitting.additional_constr_fit_iter)
            end
        end
    end

    # final evaluation
    pred, valid = eval_equation(node, data, ops)
    pred = ops.fitting.pre_residual_processing(pred, eachindex(pred), ops)
    return pred, valid
end

""" Create the cost function for the fitting and for the early stopping evaluation, fit the
    parameters in-place, and evaluate again using the optiized parameters. In case early
    stopping is used, the data is fit onto split_inds[1] and validated only split_inds[2].
    After fitting with early stopping, the parameters are selected for the best fit over all
    data."""
function fitting_LM!(node, data, ops, list_of_param, fit_iter)

    minim = x -> fitting_objective(x, node, data, ops)

    if ops.fitting.early_stop_iter > 0
        callback = trace -> early_stop_check(trace, node, list_of_param, data, ops,
            n=ops.fitting.early_stop_iter)
    else
        callback = _ -> false
    end

    x0 = [n.val for n in list_of_param]

    x_best, trace = lmfit(minim, x0, autodiff = :forward, max_iter = fit_iter,
        t_lim = ops.fitting.t_lim, rel_f_tol_5_iter = ops.fitting.rel_f_tol_5_iter,
        x_tol=0.0, g_tol=0.0,
        callback = callback
    )

    if ops.fitting.early_stop_iter > 0 && length(trace) > 2
        best_ind = argmin(
            x -> trace[x].value + trace[x].metadata["test_residual_norm"],
            2:length(trace)
        )
        x_best = trace[best_ind].metadata["x"]
    end

    set_params!(list_of_param, x_best)
    return trace[end].metadata["f_calls"]
end

""" Is called during each iteration of the fitting, if early_stop_iter > 0. The
    residual calculated for the current parameter and only for split_inds[2]. If this residual
    is monotonically increasing for early_stop_iter-iterations, the fitting is stopped.
"""
function early_stop_check(trace, node, list_of_param, data, ops; n=5)
    x = trace[end].metadata["x"]
    res_test = fitting_residual(x, node, list_of_param, data, ops.data_descript.split_inds[2], ops)
    res_test .= ops.fitting.residual_processing(res_test, ops.data_descript.split_inds[2], ops)
    res_test .*= view(ops.data_descript.fit_weights, ops.data_descript.split_inds[2])
    res_test .= abs.(res_test)

    trace[end].metadata["test_residual_norm"] = sum(abs2, res_test)

    length(trace) <= n + 1 && return false
    return issorted(trace[i].metadata["test_residual_norm"] for i in length(trace)-n:length(trace))
end

""" Fitting with Optim fitters.
"""
function fitting_Optim!(node, data, ops, list_of_param, fit_iter, algo)

    minim = x -> mean(abs2, fitting_objective(x, node, data, ops))

    x0 = [n.val for n in list_of_param]

    # callback # ----------------------------------------------------------------------------------
    function callback(tr)
        length(tr) > 5 || return false
        cur  = tr[end].value
        prev = tr[end-5].value
        return abs((prev - cur) / prev) < ops.fitting.rel_f_tol_5_iter
    end
    # ----------------------------------------------------------------------------------------------

    res = Optim.optimize(
        minim, x0, algo,
        Optim.Options(;
            show_warnings  = false, iterations       = fit_iter, time_limit     = ops.fitting.t_lim,
            show_trace     = false, store_trace      = true,     callback       = callback,
            x_abstol = 0.0, x_reltol = 0.0, f_abstol = 0.0, f_reltol = 0.0, g_abstol = 0.0,
        ), autodiff = :forward
    )
    x_best = res.trace[end].value < res.trace[1].value ? Optim.minimizer(res) : x0

    set_params!(list_of_param, x_best)
    return res.f_calls
end

""" Fitting with constraints using the penalty method.
"""
function fitting_w_constr!(node, data, ops, list_of_param, max_iter)

    minim = x -> fitting_objective(x, node, data, ops)

    function node_func(params::Vector{T1}, constr_data::Vector{Vector{T2}}) where {T1, T2} # evaluate a node at given parameters and given data
        T = promote_type(T1, T2)
        node_ = convert(T, node)
        list_of_param_ = list_of_param_nodes(node_)
        fitting_eval(params, node_, list_of_param_, constr_data, eachindex(constr_data[1]), ops)
    end

    function eq_constr(x::Vector{T})::Vector{T} where {T}
        reduce(vcat, f(node_func, x) for f in ops.fitting.eq_constr; init=eltype(x)[])
    end

    function ineq_constr(x::Vector{T})::Vector{T} where {T}
        reduce(vcat, f(node_func, x) for f in ops.fitting.ineq_constr; init=eltype(x)[])
    end

    x0 = [n.val for n in list_of_param]

    x_best, trace, stop_msg, valid = penalty_method(
        minim, eq_constr, ineq_constr, x0,
        max_iter                   = max_iter,
        b_track_iter               = 3,
        constr_tol                 = ops.fitting.constr_tol,
        obj_mare_tol               = ops.fitting.max_mare_for_constr_fit,
        init_constr_penalty_factor = 1e-1,
        incr_constr_penalty_factor = 1.5,
    )

    sort!(trace, lt=(t1, t2) -> lessthan_trace_w_constr(
        t1,
        t2,
        ops.fitting.constr_tol,
        ops.fitting.max_mare_for_constr_fit
    ))
    x_best .= trace[1].x

    set_params!(list_of_param, x_best)

    return
end

function lessthan_trace_w_constr(t1, t2, constr_tol, max_mare)
    if t1.obj_mare > max_mare || t2.obj_mare > max_mare
        return t1.obj_mare < t2.obj_mare
    elseif t1.constr_val > constr_tol || t2.constr_val > constr_tol
        return t1.constr_val < t2.constr_val
    else
        return t1.obj_mare < t2.obj_mare
    end
end

# # ==================================================================================================
# # DynamicExpressions
# # ==================================================================================================
# function convert_to_dynamic_expression(node, ops)
#     if node.ari == 2
#         return ops.binops[node.ind](
#             convert_to_dynamic_expression(node.lef, ops),
#             convert_to_dynamic_expression(node.rig, ops),
#         )
#     elseif node.ari == 1
#         return ops.unaops[node.ind](
#             convert_to_dynamic_expression(node.lef, ops)
#         )
#     elseif node.ari == 0
#         return DynamicExpressions.Expression(
#             DynamicExpressions.Node{Float64}(feature=node.ind); ops.dynam_expr.operators,
#             ops.dynam_expr.variable_names
#         )
#     elseif node.ari == -1
#         # return DynamicExpressions.Node{Float64}(node.val)
#         return node.val
#     end
# end

# function fitting_NM!(node, data, ops, list_of_param, fit_iter)
#
#     mnode = convert_to_dynamic_expression(node, ops)
#
#     x0, refs = DynamicExpressions.get_constants(mnode)
#     X        = reduce(hcat, data)'[:, ops.data_descript.split_inds[1]]
#
#     minim = x -> begin
#         DynamicExpressions.set_constants!(mnode, x, refs)
#         res, valid = TiSR.DynamicExpressions.eval_tree_array(mnode, X, ops.dynam_expr.operators)
#         # TODO: add pre_residual_processing
#         res .= X[end, :] .- res
#         # TODO: add residual_processing
#         res .*= view(ops.data_descript.fit_weights, ops.data_descript.split_inds[1])
#         return mean(abs2, res)
#     end
#
#     res = Optim.optimize(
#         minim, x0,
#         Optim.NelderMead(),
#         Optim.Options(;
#             show_warnings  = false, iterations     = fit_iter, time_limit = ops.fitting.t_lim,
#             show_trace     = false, store_trace    = true,
#             g_abstol       = 0.0,   g_reltol       = 0.0,
#             outer_g_abstol = 0.0,   outer_g_reltol = 0.0,
#         ),
#     )
#     x_best = res.trace[end].value < res.trace[1].value ? Optim.minimizer(res) : x0
#
#     set_params!(list_of_param, x_best)
#     return res.f_calls
# end
