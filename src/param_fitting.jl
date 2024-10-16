
# ==================================================================================================
# main fitting functions
# ==================================================================================================
""" Entrance/switch for parameter fitting. Depending on whether equation has parameters or not,
    the correspoinding functions are called, and the residual statistics calculated, if the
    fitting was successful (valid).
"""
function fit_n_eval!(indiv, data, ops, fit_iter)
    node = indiv.node
    list_of_param = list_of_param_nodes(node)

    if length(list_of_param) > length(data[1])
        indiv.valid = false
        return
    end

    if fit_iter > 0 && length(list_of_param) > 0

        if !iszero(ops.fitting.lasso_factor)
            max_iter_NW = 3
        else
            max_iter_NW = 0
        end

        max_iter_LM = fit_iter - max_iter_NW

        residual, valid = residual_after_fitting_LM(node, data, ops, list_of_param, max_iter_LM)

        if !iszero(max_iter_NW)
            residual, valid = residual_after_fitting_newton(node, data, ops, list_of_param, max_iter_NW)
        end
    else
        pred, valid = eval_equation(node, data, ops)
        residual = data[end] .- pred
    end

    # calculate statistics for the residual # ------------------------------------------------------
    if valid && isfinite(sum(residual))
        indiv.mae = mean(abs, residual)
        indiv.max_ae = maximum(abs, residual)
        indiv.mse = mean(abs2, residual)

        indiv.minus_r2 = get_minus_r2(data[end] .- residual, data[end])
        indiv.minus_abs_spearman = get_minus_abs_spearman(data[end] .- residual, data[end])

        relative_ref = any(d == 0 for d in data[end]) ? max.(abs.(data[end]), 0.1) : data[end]

        indiv.mare = mean(abs, residual ./ relative_ref)
        indiv.q75_are = quantile(abs.(residual ./ relative_ref), 0.75)
        indiv.max_are = maximum(abs, residual ./ relative_ref)

        indiv.ms_processed_e = mean(abs2, ops.fitting.residual_processing(residual, eachindex(residual), ops) .* ops.data_descript.fit_weights)
        indiv.valid = true
    else
        indiv.valid = false
    end
end

""" Create the cost function for the fitting and for the early stopping evaluation, fit the
    parameters in-place, and evaluate again using the optiized parameters. In case early
    stopping is used, the data is fit onto split_inds[1] and validated only split_inds[2].
    After fitting with early stopping, the parameters are selected for the best fit over all
    data and the final residual and the whether the fitting was successful is return for all
    data.
"""
function residual_after_fitting_LM(node, data, ops, list_of_param, max_iter)

    minim = x -> fitting_objective(x, node, data, ops)

    if ops.fitting.early_stop_iter > 0
        callback = trace -> early_stop_check(trace, node, list_of_param, data, ops,
            n=ops.fitting.early_stop_iter)
    else
        callback = trace -> false
    end

    x0 = getfield.(list_of_param, :val)

    x_best, trace = lmfit(minim, x0, autodiff = :forward, max_iter = ops.fitting.max_iter,
        t_lim = ops.fitting.t_lim, rel_f_tol_5_iter = ops.fitting.rel_f_tol_5_iter,
        callback = callback,
    )

    if ops.fitting.early_stop_iter > 0 && length(trace) > 2
        best_ind = argmin(
            x -> trace[x].value + trace[x].metadata["test_residual_norm"],
            2:length(trace)
        )
        x_best = trace[best_ind].metadata["x"]
    end

    for i in eachindex(list_of_param)
        list_of_param[i].val = x_best[i]
    end

    residual, valid = fitting_residual(x_best, node, list_of_param, data, eachindex(data[1]), ops)
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
""" Return the fitting residual. Pre-residual_processing! funciton is applied in-place
"""
function fitting_residual(x, node, list_of_param, data, inds, ops)
    for i in eachindex(x)
        list_of_param[i].val = x[i]
    end

    dat = [view(d, inds) for d in data]

    arr, valid = eval_equation(node, dat, ops)
    ops.fitting.pre_residual_processing!(arr, inds, ops)
    dat[end] .- arr, valid
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

    # if ops.fitting.lasso_factor > 0.0
    #     push!(res, sum(abs, x) * ops.fitting.lasso_factor)
    # end
    return res
end

# ==================================================================================================
# new one
# ==================================================================================================

function residual_after_fitting_newton(node, data, ops, list_of_param, max_iter)

    # create the closure
    minim = x -> mean(abs2, fitting_objective(x, node, data, ops)) + ops.fitting.lasso_factor * sum(abs, x)

    x0 = getfield.(list_of_param, :val)

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

    for i in eachindex(list_of_param)
        list_of_param[i].val = x_best[i]
    end

    residual, valid = fitting_residual(x_best, node, list_of_param, data, eachindex(data[1]), ops)

    return residual, valid
end

