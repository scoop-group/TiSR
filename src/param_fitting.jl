
# ==================================================================================================
# helper functions
# ==================================================================================================
""" Return the fitting residual. Pre-residual_processing! funciton is applied in-place
"""
function fitting_residual(x, node, list_of_param, data, inds, ops)
    for i in eachindex(x)
        list_of_param[i].val = x[i]
    end

    dat = [view(d, inds) for d in data] # dat = [d[inds] for d in data]

    arr, valid = eval_equation(node, dat[1:end-1], ops)
    ops.fitting.pre_residual_processing!(arr, inds)
    dat[end] .- arr, valid
end

# ==================================================================================================
# main fitting functions
# ==================================================================================================
""" Entrance/switch for parameter fitting. Depending on whether equation has parameters or not,
    the correspoinding functions are called, and the residual statistics calculated, if the
    fitting was successful (valid).
"""
function fit_n_eval!(indiv, data, ops)
    node = indiv.node
    list_of_param = list_of_param_nodes(node)

    if ops.fitting.max_iter > 0 && length(list_of_param) > 0
        residual, valid = residual_after_fitting(node, data, ops, list_of_param)

    else
        pred, valid = eval_equation(node, data[1:end-1], ops)
        residual = data[end] .- pred
    end

    # calculate statistics for the residual # ------------------------------------------------------
    if (valid = valid && isfinite(sum(residual)))
        indiv.mae = mean(abs, residual)
        indiv.max_ae = maximum(abs, residual)
        indiv.mse = mean(abs2, residual)

        indiv.minus_r2 = -r_squared(data[end] .- residual, data[end])

        relative_ref = any(d == 0 for d in data[end]) ? max.(abs.(data[end]), 0.1) : data[end]

        indiv.mare = mean(abs, residual ./ relative_ref)
        indiv.q75_are = quantile(abs.(residual ./ relative_ref), 0.75)
        indiv.max_are = maximum(abs, residual ./ relative_ref)

        indiv.ms_processed_e = mean(abs2, ops.fitting.residual_processing(residual, eachindex(residual)) .* ops.data_descript.fit_weights)
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
function residual_after_fitting(node, data, ops, list_of_param)

    function minim(x::Vector{T})::Vector{T} where {T <: Number}
        data_ = data#[convert(Vector{T}, d) for d in data]
        node_ = convert(T, node)
        list_of_param_ = list_of_param_nodes(node_)

        res   = fitting_residual(x, node_, list_of_param_, data_, ops.data_descript.split_inds[1], ops)[1]
        res  .= ops.fitting.residual_processing(res, ops.data_descript.split_inds[1])
        res .*= view(ops.data_descript.fit_weights, ops.data_descript.split_inds[1])
        res  .= abs.(res)
        return res
    end

    if ops.fitting.early_stop_iter > 0
        callback = trace -> early_stop_check(trace, node, list_of_param, data, ops,
            n=ops.fitting.early_stop_iter)
    else
        callback = trace -> false
    end

    x0 = getfield.(list_of_param, :val)

    x_best, trace = lmfit(
        minim,
        x0,
        autodiff           = :forward,
        x_tol              = 1e-16,
        g_tol              = 1e-16,
        max_iter           = ops.fitting.max_iter,
        t_lim              = ops.fitting.t_lim,
        rel_f_tol_5_iter   = ops.fitting.rel_f_tol_5_iter,
        callback           = callback,
        proceed_cautiously = false,
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
    test_res = abs.(ops.fitting.residual_processing(
        fitting_residual(
            trace[end].metadata["x"],
            node,
            list_of_param,
            data,
            ops.data_descript.split_inds[2],
            ops
        )[1],
        ops.data_descript.split_inds[2]
    ) .* view(ops.data_descript.fit_weights, ops.data_descript.split_inds[2]))

    trace[end].metadata["test_residual_norm"] = sum(abs2, test_res)

    length(trace) <= n + 1 && return false
    return issorted(trace[i].metadata["test_residual_norm"] for i in length(trace)-n:length(trace))
end

