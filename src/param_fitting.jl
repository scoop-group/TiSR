
# ==================================================================================================
# main fitting functions
# ==================================================================================================
""" Entrance/switch for parameter fitting. Depending on whether the equation has parameters
    or not, the correspoinding functions are called, and the residual statistics calculated,
    if the fitting was successful (valid).
"""
function fit_n_eval!(node, data, ops, fit_iter)
    list_of_param = list_of_param_nodes(node)

    if length(list_of_param) > length(ops.data_descript.split_inds[1])
        return data[end], false
    end

    if fit_iter > 0 && !isempty(list_of_param)

        if !iszero(ops.fitting.lasso_factor)
            max_iter_NW = 3
        else
            max_iter_NW = 0
        end

        max_iter_LM = fit_iter - max_iter_NW

        fitting_LM!(node, data, ops, list_of_param, max_iter_LM)

        if !iszero(max_iter_NW)
            fitting_NW!(node, data, ops, list_of_param, max_iter_NW)
        end
    end

    # final evaluation
    prediction, valid = eval_equation(node, data, ops)
    return prediction, valid
end

""" Create the cost function for the fitting and for the early stopping evaluation, fit the
    parameters in-place, and evaluate again using the optiized parameters. In case early
    stopping is used, the data is fit onto split_inds[1] and validated only split_inds[2].
    After fitting with early stopping, the parameters are selected for the best fit over all
    data."""
function fitting_LM!(node, data, ops, list_of_param, max_iter)

    minim = x -> fitting_objective(x, node, data, ops)

    if ops.fitting.early_stop_iter > 0
        callback = trace -> early_stop_check(trace, node, list_of_param, data, ops,
            n=ops.fitting.early_stop_iter)
    else
        callback = trace -> false
    end

    x0 = [n.val for n in list_of_param]

    x_best, trace = lmfit(minim, x0, autodiff = :forward, max_iter = ops.fitting.max_iter,
        t_lim = ops.fitting.t_lim, rel_f_tol_5_iter = ops.fitting.rel_f_tol_5_iter,
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
end

""" Is called during each iteration of the fitting, if early_stop_iter > 0. The
    residual calculated for the current parameter and only for split_inds[2]. If this residual
    is monotonically increasing for early_stop_iter-iterations, the fitting is stopped.
"""
function early_stop_check(trace, node, list_of_param, data, ops; n=5)
    x = trace[end].metadata["x"]
    res_test = fitting_residual(x, node, list_of_param, data, ops.data_descript.split_inds[2], ops)[1]
    res_test .= ops.fitting.residual_processing(res_test, ops.data_descript.split_inds[2], ops)
    res_test .*= view(ops.data_descript.fit_weights, ops.data_descript.split_inds[2])
    res_test .= abs.(res_test)

    trace[end].metadata["test_residual_norm"] = sum(abs2, res_test)

    length(trace) <= n + 1 && return false
    return issorted(trace[i].metadata["test_residual_norm"] for i in length(trace)-n:length(trace))
end

# ==================================================================================================
# helper functions
# ==================================================================================================

""" Set the values of the parameter in a node.
"""
function set_params!(list_of_param, x)
    for i in eachindex(x)
        list_of_param[i].val = x[i]
    end
end

""" Return the fitting residual. Pre-residual_processing! funciton is applied in-place
"""
function fitting_residual(x, node, list_of_param, data, inds, ops)
    set_params!(list_of_param, x)

    dat = [view(d, inds) for d in data]
    pred, valid = eval_equation(node, dat, ops)
    ops.fitting.pre_residual_processing!(pred, inds, ops)
    dat[end] .- pred, valid # TODO: maybe remove the valid here?
end

""" Calculate the residual and apply pre- and post-processing.
"""
function fitting_objective(x::Vector{T}, node, data, ops)::Vector{T} where {T <: Number}
    data_ = data#[convert(Vector{T}, d) for d in data]
    node_ = convert(T, node)
    list_of_param_ = list_of_param_nodes(node_)

    res   = fitting_residual(x, node_, list_of_param_, data_, ops.data_descript.split_inds[1], ops)[1]
    res  .= ops.fitting.residual_processing(res, ops.data_descript.split_inds[1], ops)
    res .*= view(ops.data_descript.fit_weights, ops.data_descript.split_inds[1])
    res  .= abs.(res)
    return res
end

# ==================================================================================================
# new one
# ==================================================================================================

function fitting_NW!(node, data, ops, list_of_param, max_iter)

    minim = x -> mean(abs2, fitting_objective(x, node, data, ops)) + ops.fitting.lasso_factor * sum(abs, x)

    x0 = [n.val for n in list_of_param]

    # # callback # ----------------------------------------------------------------------------------
    # function callback(tr)
    #     length(tr) > 5 || return false
    #     cur  = tr[end].value
    #     prev = tr[end-5].value
    #     return (prev - cur) / prev < ops.fitting.rel_f_tol_5_iter
    # end
    # # ----------------------------------------------------------------------------------------------

    res = Optim.optimize(
        minim, x0,
        # Optim.Newton(;linesearch=LineSearches.BackTracking()),
        Optim.LBFGS(;linesearch=LineSearches.BackTracking()),
        # Optim.GradientDescent(;linesearch=LineSearches.BackTracking()),
        Optim.Options(;
            show_warnings=false, iterations=max_iter, time_limit=ops.fitting.t_lim,
            show_trace=false,    store_trace=true,# callback = callback
        ), autodiff=:forward
    )
    x_best = res.trace[end].value < res.trace[1].value ? Optim.minimizer(res) : x0

    set_params!(list_of_param, x_best)
end

