
# ==============================================================================
# from other files
# ==============================================================================
using OptimBase
using LinearAlgebra
using ForwardDiff
import NLSolversBase: value, jacobian

using LinearSolve

function lmfit_manual(func_val!::Function, func_jac!::Function, val, jac, x0::AbstractArray; autodiff=:finite, kwargs...)
    levenberg_marquardt_manual(func_val!, func_jac!, val, jac, x0; kwargs...)
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

# I think a smarter way to do this *might* be to create a type similar to `OnceDifferentiable`
# and the like. This way we could not only merge the two functions, but also have a convenient
# way to provide an autodiff-made acceleration when someone doesn't provide an `avv`.
# it would probably be very inefficient performace-wise for most cases, but it wouldn't hurt to have it somewhere
function levenberg_marquardt_manual(func_val!, func_jac!, val, jac, initial_x::AbstractVector{T};
    x_tol::Real=1e-8, g_tol::Real=1e-12, 
    max_iter::Integer=1000,
    rel_f_tol_5_iter::AbstractFloat=0.0,
    lambda=T(10), lambda_increase::Real=3.0, lambda_decrease::Real=0.4,
    min_step_quality::Real=1e-3, good_step_quality::Real=0.75, 
    callback=x -> false, t_lim::Float64=floatmax(1.0), proceed_cautiously=false
    ) where T

    # T = Float64
    # x_tol=1e-8
    # g_tol=1e-12
    # max_iter=1000
    # lambda=T(10)
    # lambda_increase=3.0
    # lambda_decrease=0.4
    # min_step_quality=1e-3
    # good_step_quality=0.75
    # callback=x -> false
    # t_lim=floatmax(1.0)
    # proceed_cautiously = true

    # initial_x = x0
    # func_val! = f_res_pro_val!
    # func_jac! = f_jac!

    # First evaluation
    valid = func_val!(val, initial_x)
    func_jac!(jac, initial_x)

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

    res_vect = similar(val)
    residual = sum(abs2, val)

    # Create buffers
    n = length(x)
    JᵀJ = Matrix{T}(undef, n, n)
    Jᵀr = Vector{T}(undef, n)
    x_trial = Vector{T}(undef, n)
    Jdelta_buffer = similar(val)

    # and an alias for the jacobian
    J = jac
    delta_x = Array{T}(undef, n)

    # Maintain a trace of the system # -----------------------------------------
    tr = OptimBase.OptimizationTrace{LevenbergMarquardt}()
    d = Dict("status" => "", "lambda" => lambda, "x" => initial_x)
    os = OptimBase.OptimizationState{LevenbergMarquardt}(iteration, sum(abs2, val), NaN, d)
    push!(tr, os)

    # start iterations # -------------------------------------------------------
    t0 = time()

    calc_jac = false

    while (~converged && iteration < max_iter) && (time() - t0 < t_lim)
        iteration += 1

    # calculate the jacobian, but only if new x # ------------------------------
        if calc_jac
            func_jac!(J, x)
        end

        DtD = vec(sum(abs2, J, dims=1))

        clamp!(DtD, MIN_DIAGONAL, Inf)
        
        # calc JᵀJ + λ ⋅ diag(DtD) (left side)
        mul!(JᵀJ, transpose(J), J)
        @simd for i in 1:n
            @inbounds JᵀJ[i, i] += lambda * DtD[i]
        end

        # calc -Jᵀ⋅r (right side)
        mul!(Jᵀr, transpose(J), val)
        rmul!(Jᵀr, -1)

    # solve LGS to get delta_x # -----------------------------------------------
        if isfinite(sum(JᵀJ)) && rank(JᵀJ) == size(JᵀJ, 2)
            tr[end].metadata["status"] = tr[end].metadata["status"] * "fullrank_"
            delta_x .= JᵀJ \ Jᵀr
            
        elseif isfinite(sum(JᵀJ)) && proceed_cautiously
            try
                delta_x .= JᵀJ \ Jᵀr
                tr[end].metadata["status"] = tr[end].metadata["status"] * "proceededcautiously_"
                
            catch
                lambda == MAX_LAMBDA && break
                lambda = min(lambda_increase * lambda, MAX_LAMBDA)
                tr[end].metadata["status"] = tr[end].metadata["status"] * "failedcautiously_"
                continue
            end
        else
            tr[end].metadata["status"] = tr[end].metadata["status"] * "notfinite"
            break
        end

    # calculate step quality and determine dampening for next step # -----------
        # predicted residual according to model
        mul!(Jdelta_buffer, J, delta_x)
        Jdelta_buffer .= Jdelta_buffer .+ val
        predicted_residual = sum(abs2, Jdelta_buffer)

        # compute res_vect for new potential x inplace according to NLSolversBase 
        # value(obj, cache, state) -> doesn't update df besides function calls
        x_trial .= x .+ delta_x

        valid = func_val!(res_vect, x_trial)

        if !valid 
            tr[end].metadata["status"] = tr[end].metadata["status"] * "nonvalid"
            lambda == MAX_LAMBDA && break
            lambda = min(lambda_increase * lambda, MAX_LAMBDA)
            continue
        end

        # update the sum of squares
        trial_residual = sum(abs2, res_vect)

        # step quality = residual change / predicted residual change
        rho = (trial_residual - residual) / (predicted_residual - residual)
        if trial_residual < residual && rho > min_step_quality
            # apply the step to x
            x .= x_trial
            val .= res_vect
            residual = trial_residual

            calc_jac = true

            # # increase trust region radius (from Ceres solver)
            # lambda = max(lambda * max(1 / 3, 1.0 - (2.0 * rho - 1.0)^3), MIN_LAMBDA)
            
            if rho > good_step_quality
                # increase trust region radius
                lambda = max(lambda_decrease * lambda, MIN_LAMBDA)
            end

        else
            # decrease trust region radius
            lambda = min(lambda_increase * lambda, MAX_LAMBDA)
            calc_jac = false
        end

    # update the trace # -------------------------------------------------------
        g_norm = norm(J' * val, Inf)
        d = Dict("status" => "",
                 "norm(g(x))" => g_norm, 
                 "dx" => delta_x, 
                 "step_norm" => norm(delta_x, 2),
                 "lambda" => lambda, 
                 "x" => x)            

        os = OptimBase.OptimizationState{LevenbergMarquardt}(iteration, sum(abs2, val), g_norm, d)
        push!(tr, os)


    # check convergence criteria # ---------------------------------------------
        rel_tol_conv = rel_f_tol_5_iter > 0 && (iteration > 5 && 
        (tr[end-min(5, length(tr) - 1)].value - os.value) / tr[end-min(5, length(tr) - 1)].value < rel_f_tol_5_iter)

        if (g_norm < g_tol 
            || norm(delta_x) < x_tol * (x_tol + norm(x))
            || rel_tol_conv)
            converged = true

            tr[end].metadata["status"] = tr[end].metadata["status"] * "converged"
        end

    # do the callback # --------------------------------------------------------
        if callback(tr)
            tr[end].metadata["status"] = tr[end].metadata["status"] * "callbackbreak"
            break
        end

    end
    x, tr
end



