
""" Recursive functions using multiple dispatch to extract a value from arbitrarily
    nested dual numbers.
"""
extract_from_dual(v::AbstractFloat) = v
extract_from_dual(v::Number) = extract_from_dual(v.value)

"""
    penalty_method(
        obj_func::Function,
        eq_constr::Function,
        ineq_constr::Function,
        x0::Vector{Float64};
        max_iter::Int64=1000,
        b_track_iter::Int64=3,
        b_track_factor::Float64=0.5,
        x_atol::Float64=1e-6,
        constr_tol::Float64=1e-8,
        g_tol::Float64=1e-8,
        obj_mare_tol::Float64=Inf,
        init_constr_penalty_factor::Float64=1e-1,
        incr_constr_penalty_factor::Float64=1.5
    )

Implements a penalty method for constrained optimization, minimizing the objective function `obj_func` subject to equality constraints `eq_constr` and inequality constraints `ineq_constr`, starting from an initial guess `x0`.

# Arguments
- `obj_func::Function`: The objective function to minimize, taking a vector `x::Vector{Float64}` as input and returning a `Vector{Float64}`.
- `eq_constr::Function`: A vector-valued function representing equality constraints, where `eq_constr(x::Vector{Float64})::Vector{Float64}` returns a vector that equals zero for feasible points.
- `ineq_constr::Function`: A vector-valued function representing inequality constraints, where `ineq_constr(x::Vector{Float64})::Vector{Float64}` returns a vector with non-positive elements for feasible points.
- `x0::Vector{Float64}`: Initial guess for the optimization variable.

# Keyword Arguments
- `max_iter::Int64=1000`: Maximum number of iterations for the optimization process.
- `b_track_iter::Int64=3`: Number of iterations to track for backtracking line search.
- `b_track_factor::Float64=0.5`: Backtracking factor for line search (0 < factor < 1).
- `x_atol::Float64=1e-6`: Absolute tolerance for convergence in the optimization variable `x`.
- `constr_tol::Float64=1e-8`: Tolerance for constraint satisfaction.
- `g_tol::Float64=1e-8`: Tolerance for the gradient norm of the Lagrangian.
- `obj_mare_tol::Float64=Inf`: Tolerance for the maximum absolute relative error in the objective function.
- `init_constr_penalty_factor::Float64=1e-1`: Initial penalty factor for constraints.
- `incr_constr_penalty_factor::Float64=1.5`: Factor by which the penalty parameter is increased in each iteration.

# Returns
- `x::Vector{Float64}`: result point
- optimization trace
- stop message
- whether the optimization exited gracefully

# Notes
- The penalty method converts the constrained optimization problem into a sequence of unconstrained problems by adding penalty terms for constraint violations.
- The algorithm iteratively solves these unconstrained problems, increasing the penalty factor until constraints are satisfied within `constr_tol`.
- Convergence is determined based on `x_atol`, `g_tol`, and `obj_mare_tol`.
"""
function penalty_method(
    obj_func::Function,
    eq_constr::Function,
    ineq_constr::Function,
    x0::Vector{Float64};
    max_iter::Int64                     = 1000,
    b_track_iter::Int64                 = 3,
    b_track_factor::Float64             = 0.5,
    x_atol::Float64                     = 1e-6,
    constr_tol::Float64                 = 1e-8,
    g_tol::Float64                      = 1e-8,
    obj_mare_tol::Float64               = Inf,
    init_constr_penalty_factor::Float64 = 1e-1,
    incr_constr_penalty_factor::Float64 = 1.5,
)
    x      = copy(x0)
    p      = init_constr_penalty_factor
    result = DiffResults.HessianResult(x)

    dx         = similar(x)
    f_evals    = 0

    constr_func = x -> sum(abs2, eq_constr(x)) + sum(abs2, max.(0.0, ineq_constr(x)))

    trace = [(
        iter = 0, x = copy(x), p = p,
        obj_val    = mean(abs2, obj_func(x)),
        obj_mare   = mean(abs, obj_func(x)),
        constr_val = constr_func(x),
        alpha=NaN, f_evals = 0
    )]

    obj_n_constr_val = zeros(3)

    pen_obj_func! = (x, p, obj_n_constr_val, obj_func, constr_func) -> begin
        obj_residual = obj_func(x)
        obj_val      = mean(abs2, obj_residual)
        obj_mare     = mean(abs,  obj_residual)
        constr_val   = constr_func(x)

        # save components of objective function for use outside
        obj_n_constr_val[1] = extract_from_dual(obj_val)
        obj_n_constr_val[2] = extract_from_dual(obj_mare)
        obj_n_constr_val[3] = extract_from_dual(constr_val)

        return obj_val + p * constr_val
    end

    stop_msg = ""

    for iter in 1:max_iter
        ForwardDiff.hessian!(result, x -> pen_obj_func!(x, p, obj_n_constr_val, obj_func, constr_func), x)

        H = DiffResults.hessian(result)
        g = DiffResults.gradient(result)
        v = DiffResults.value(result)

        if !(isfinite(sum(H)) && isfinite(sum(g)))
            return x, trace, "non finite", false
        end

        # lev_mar hack # ---------------------------------------------------------------------------
        lev_mar = 1e-12
        while rank(H) < size(H, 1) && lev_mar < 1e12
            H .+= (I(size(H, 1)) .* lev_mar)
            lev_mar *= 10.0
        end
        if lev_mar >= 1e12
            return x, trace, "signular exception", false
        end

        # step # -----------------------------------------------------------------------------------
        dx .= -H \ g

        # backtracking # ---------------------------------------------------------------------------
        alpha = 1.0
        suff_decr_factor = 1e-5

        for _ in 1:b_track_iter
            x_trial     = x .+ alpha .* dx
            trial_merit = pen_obj_func!(x_trial, p, obj_n_constr_val, obj_func, constr_func)

            if trial_merit < v + suff_decr_factor * alpha * (g' * dx)
                break
            end
            alpha *= b_track_factor
        end
        #-------------------------------------------------------------------------------------------
        x .+= (alpha * dx)

        f_evals += trunc(Int64, log(0.5, alpha) + 1.0)

        # save state # -----------------------------------------------------------------------------
        push!(trace, (
            iter = iter, x = copy(x), p = p,
            obj_val    = obj_n_constr_val[1],
            obj_mare   = obj_n_constr_val[2],
            constr_val = obj_n_constr_val[3],
            alpha = alpha, f_evals = f_evals
        ))

        # update penalty
        if abs(obj_n_constr_val[3]) > constr_tol
            p *= incr_constr_penalty_factor
        end

        # check stopping criteria # ----------------------------------------------------------------
        if norm(g, 2) < g_tol
            stop_msg = "g_tol reached"
        elseif norm(alpha * dx, 2) < x_atol
            stop_msg = "x_atol reached"
        elseif trace[end].obj_mare > obj_mare_tol
            stop_msg = "exceed obj_mare_tol"
        end
        stop_msg == "" || break
    end
    return x, trace, stop_msg, true
end


# ==================================================================================================
# archive
# ==================================================================================================

# -> an oler version, but as augmented lagrangian method
# function penalty_method(
#     obj_func,
#     eq_constr,
#     ineq_constr,
#     x0;
#     max_iter                   = 1000,
#     b_track_iter               = 3,
#     b_track_factor             = 0.5,
#     constr_tol                 = 1e-8,
#     x_atol                     = 1e-8,
#     g_tol                      = 1e-8,
#     obj_mare_tol               = Inf,
#     init_constr_penalty_factor = 1e-1,
#     incr_constr_penalty_factor = 1.5,
# )
#     x      = copy(x0)
#     p      = init_constr_penalty_factor
#     result = DiffResults.HessianResult(x)
#
#     dx         = similar(x)
#     alpha      = 0.0
#     f_evals    = 0
#
#     h_vios = eq_constr(x)
#     g_vios = max.(0.0, ineq_constr(x))
#
#     h_lagr = zeros(length(h_vios))
#     g_lagr = zeros(length(g_vios))
#
#     trace = [(
#         iter = 0, x = copy(x), p = p,
#         obj_val    = mean(abs2, obj_func(x)),
#         obj_mare   = mean(abs, obj_func(x)),
#         constr_val = sum(abs2, h_vios) + sum(abs2, g_vios),
#         alpha=NaN, f_evals = 0
#     )]
#
#     obj_n_constr_val = zeros(3)
#
#     pen_obj_func! = (x, p, h_vios, g_vios, h_lagr, g_lagr, obj_n_constr_val, obj_func, eq_constr, ineq_constr) -> begin
#         obj_residual = obj_func(x)
#         obj_val      = mean(abs2, obj_residual)
#         obj_mare     = mean(abs,  obj_residual) # TODO: not actual mare?
#
#         h_vios_ = eq_constr(x)
#         g_vios_ = max.(0.0, ineq_constr(x))
#
#         constr_val = sum(abs2, h_vios_) + sum(abs2, g_vios_)
#
#         # save components of objective function for use outside
#         obj_n_constr_val[1] = extract_from_dual(obj_val)
#         obj_n_constr_val[2] = extract_from_dual(obj_mare)
#         obj_n_constr_val[3] = extract_from_dual(constr_val)
#         h_vios             .= extract_from_dual.(h_vios_)
#         g_vios             .= extract_from_dual.(g_vios_)
#
#         return obj_val + p * constr_val + dot(h_lagr, h_vios) + dot(g_lagr, g_vios)
#     end
#
#     stop_msg = ""
#
#     for iter in 1:max_iter
#         ForwardDiff.hessian!(result, x -> pen_obj_func!(x, p, h_vios, g_vios, h_lagr, g_lagr, obj_n_constr_val, obj_func, eq_constr, ineq_constr), x)
#
#         if !(isfinite(sum(DiffResults.hessian(result))) && isfinite(sum(DiffResults.gradient(result))))
#             return x, trace, "non finite", false
#         end
#
#         dx .= -DiffResults.hessian(result) \ DiffResults.gradient(result)
#
#         # backtracking # -----------------------------------------------------------------------
#         alpha = 1.0
#         suff_decr_factor = 1e-5
#
#         for _ in 1:b_track_iter
#             x_trial     = x .+ alpha .* dx
#             trial_merit = pen_obj_func!(x_trial, p, h_vios, g_vios, h_lagr, g_lagr, obj_n_constr_val, obj_func, eq_constr, ineq_constr)
#
#             if trial_merit < DiffResults.value(result) + suff_decr_factor * alpha * (DiffResults.gradient(result)' * dx)
#                 break
#             end
#         end
#         #---------------------------------------------------------------------------------------
#         x .+= (alpha * dx)
#
#         f_evals += trunc(Int64, log(0.5, alpha) + 1.0)
#
#         # save state # -----------------------------------------------------------------------------
#         push!(trace, (
#             iter = iter, x = copy(x), p = p,
#             obj_val    = obj_n_constr_val[1],
#             obj_mare   = obj_n_constr_val[2],
#             constr_val = obj_n_constr_val[3],
#             alpha = alpha, f_evals = f_evals
#         ))
#
#         # update penalty
#         if abs(obj_n_constr_val[3]) > constr_tol
#             h_lagr .+= h_vios * p
#             g_lagr .+= g_vios * p
#             p *= incr_constr_penalty_factor
#         end
#
#         # check stopping criteria # ----------------------------------------------------------------
#         if norm(DiffResults.gradient(result), 2) < g_tol
#             stop_msg = "g_tol reached"
#             break
#         elseif norm(alpha * dx, 2) < x_atol
#             stop_msg = "x_atol reached"
#             break
#         elseif trace[end].obj_mare > obj_mare_tol
#             stop_msg = "exceed obj_mare_tol"
#             break
#         end
#     end
#     return x, trace, stop_msg, true
# end

# function backtracking(pen_obj_func!, x, dx, p, obj_n_constr_val, result;
#     alpha = 1.0, b_track_factor = 0.5, b_track_iter = 5, suff_decr_factor = 1e-5
# )
#     init_merit = DiffResults.value(result)
#     init_grad  = DiffResults.gradient(result)
#
#     iter = 0
#     trial_merit = Inf
#     while iter < b_track_iter || (isnan(trial_merit) && iter < 10)
#         iter += 1
#         x_trial = x .+ alpha .* dx
#
#         trial_merit = pen_obj_func!(x_trial, p, obj_n_constr_val)
#
#         if trial_merit < init_merit + suff_decr_factor * alpha * (init_grad' * dx)
#             return alpha
#         end
#
#         alpha *= b_track_factor
#     end
#     return alpha
# end

# --------------------------------------------------------------------------------------
# lev_mar = 1e-12
# while rank(H) < size(H, 1) && lev_mar < 1e12
#     H .+= (I(size(H, 1)) .* lev_mar)
#     lev_mar *= 10.0
#     println("padding")
# end
# if lev_mar >= 1e12
#     return x, trace, "signular exception", false
# end
# --------------------------------------------------------------------------------------

# for _ in 1:max_inner_iter
# end
