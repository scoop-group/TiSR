
using OptimBase
using LinearAlgebra
using ForwardDiff
import NLSolversBase: value, jacobian

function lmfit(f, x0::AbstractArray; autodiff=:finite, kwargs...)
    r = f(x0)
    R = OnceDifferentiable(f, x0, copy(r); inplace=false, autodiff=autodiff)
    levenberg_marquardt(R, x0; kwargs...)
end

# ==============================================================================
# core algorithm
# ==============================================================================
struct LevenbergMarquardt <: Optimizer end
Base.summary(::LevenbergMarquardt) = "Levenberg-Marquardt"

"""
    `levenberg_marquardt(f, g, initial_x; <keyword arguments>`

Returns the argmin over x of `sum(f(x).^2)` using the Levenberg-Marquardt
algorithm, and an estimate of the Jacobian of `f` at x.

The function `f` should take an input vector of length n and return an output
vector of length m. The function `g` is the Jacobian of f, and should return an m x
n matrix. `initial_x` is an initial guess for the solution.

Implements box constraints as described in Kanzow, Yamashita, Fukushima (2004; J
Comp & Applied Math).

# Keyword arguments
* `x_tol::Real=1e-8`: search tolerance in x
* `g_tol::Real=1e-12`: search tolerance in gradient
* `maxIter::Integer=1000`: maximum number of iterations
* `min_step_quality=1e-3`: for steps below this quality, the trust region is shrinked
* `good_step_quality=0.75`: for steps above this quality, the trust region is expanded
* `lambda::Real=10`: (inverse of) initial trust region radius
* `lambda_increase=10.0`: `lambda` is multiplied by this factor after step below min quality
* `lambda_decrease=0.1`: `lambda` is multiplied by this factor after good quality steps
* `lower,upper=[]`: bound solution to these limits
"""

function levenberg_marquardt(df::OnceDifferentiable, initial_x::AbstractVector{T};
    x_tol::Real=1e-8, g_tol::Real=1e-12,
    rel_f_tol_5_iter::AbstractFloat=0.0,
    max_iter::Integer=1000,
    lambda=T(10), lambda_increase::Real=3.0, lambda_decrease::Real=0.4,
    min_step_quality::Real=1e-3, good_step_quality::Real=0.75,
    callback=x -> false, t_lim::Float64=floatmax(1.0), proceed_cautiously=false
    ) where T

    # First evaluation
    value_jacobian!!(df, initial_x)

    @assert (0 <= min_step_quality < 1) " 0 <= min_step_quality < 1 must hold."
    @assert (0 < good_step_quality <= 1) " 0 < good_step_quality <= 1 must hold."
    @assert (min_step_quality < good_step_quality) "min_step_quality < good_step_quality must hold."

    # other constants
    MAX_LAMBDA = 1e16 # minimum trust region radius
    MIN_LAMBDA = 1e-16 # maximum trust region radius
    MIN_DIAGONAL = 1e-6 # lower bound on values of diagonal matrix used to regularize the trust region step

    x = copy(initial_x)
    delta_x = copy(initial_x)

    converged = false
    iteration = 0

    res_vect = similar(value(df))
    residual = sum(abs2, value(df))

    # Create buffers
    n = length(x)
    JᵀJ = Matrix{T}(undef, n, n)
    Jᵀr = Vector{T}(undef, n)
    x_trial = Vector{T}(undef, n)
    Jdelta_buffer = similar(value(df))

    # and an alias for the jacobian
    J = jacobian(df)
    delta_x = Array{T}(undef, n)

    # Maintain a trace of the system # -----------------------------------------
    tr = OptimBase.OptimizationTrace{LevenbergMarquardt}()
    d = Dict("lambda" => lambda, "x" => initial_x)
    os = OptimBase.OptimizationState{LevenbergMarquardt}(iteration, sum(abs2, value(df)), NaN, d)
    push!(tr, os)

    # start iterations # -------------------------------------------------------
    t0 = time()

    while (~converged && iteration < max_iter) && (time() - t0 < t_lim)
        iteration += 1

        # calculate the jacobian, but only if new x # ------------------------------
        jacobian!(df, x) # has alias J

        DtD = vec(sum(abs2, J, dims=1))

        clamp!(DtD, MIN_DIAGONAL, Inf)

        # calc JᵀJ + λ ⋅ diag(DtD) (left side)
        mul!(JᵀJ, transpose(J), J)
        @simd for i in 1:n
            @inbounds JᵀJ[i, i] += lambda * DtD[i]
        end

        # calc -Jᵀ⋅r (right side)
        mul!(Jᵀr, transpose(J), value(df))
        rmul!(Jᵀr, -1)

        # solve LGS to get delta_x # -----------------------------------------------
        if isfinite(sum(JᵀJ)) && rank(JᵀJ) == size(JᵀJ, 2)
            delta_x .= JᵀJ \ Jᵀr

        elseif isfinite(sum(JᵀJ)) && proceed_cautiously
            try
                delta_x .= JᵀJ \ Jᵀr
            catch
                lambda == MAX_LAMBDA && break
                lambda = min(lambda_increase * lambda, MAX_LAMBDA)
                continue
            end
        else
            break
        end

        # calculate step quality and determine dampening for next step # -----------
        # predicted residual according to model
        mul!(Jdelta_buffer, J, delta_x)
        Jdelta_buffer .= Jdelta_buffer .+ value(df)
        predicted_residual = sum(abs2, Jdelta_buffer)

        # compute res_vect for new potential x inplace according to NLSolversBase
        # value(obj, cache, state) -> doesn't update df besides function calls
        x_trial .= x .+ delta_x
        value(df, res_vect, x_trial)

        # update the sum of squares
        trial_residual = sum(abs2, res_vect)

        # step quality = residual change / predicted residual change
        rho = (trial_residual - residual) / (predicted_residual - residual)
        if trial_residual < residual && rho > min_step_quality
            # apply the step to x
            x .= x_trial
            copyto!(df.x_f, x)
            copyto!(value(df), res_vect)
            residual = trial_residual

            if rho > good_step_quality
                # increase trust region radius
                lambda = max(lambda_decrease * lambda, MIN_LAMBDA)
            end

        else
            # decrease trust region radius
            lambda = min(lambda_increase * lambda, MAX_LAMBDA)
        end

        # update the trace # -------------------------------------------------------
        g_norm = norm(J' * value(df), Inf)
        d = Dict("norm(g(x))" => g_norm,
                 "dx" => delta_x,
                 "step_norm" => norm(delta_x, 2),
                 "lambda" => lambda,
                 "x" => x,
                 )
        os = OptimBase.OptimizationState{LevenbergMarquardt}(iteration, sum(abs2, value(df)), g_norm, d)
        push!(tr, os)

        # check convergence criteria # ---------------------------------------------
        rel_tol_conv = (iteration > 5 &&
        (tr[end-min(5, length(tr) - 1)].value - os.value) / tr[end-min(5, length(tr) - 1)].value < rel_f_tol_5_iter)

        if (g_norm < g_tol
            || norm(delta_x) < x_tol * (x_tol + norm(x))
            || rel_tol_conv)
            converged = true
        end

        # do the callback # --------------------------------------------------------
        if callback(tr)
            break
        end
    end
    x, tr
end


