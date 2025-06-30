
data = rand(100, 5)
ops, data_vect = Options(data, fit_weights=ones(size(data,1)))

# TODO: test fitting_residual(x, node, list_of_param, data, inds, ops)
# TODO: test fitting_objective(x::Vector{T}, node, data, ops)::Vector{T} where {T <: Number}
# TODO: test fit_n_eval!node, data, ops, fit_iter)
#   - early stopping
#   - pre_residual_processing
#   - residual_processing
#   - t_lim
#   - rel_f_tol_5_iter
# TODO: test fit_weights


@testset "fitting_LM!" begin

    improvements = Float64[]

    i = 0
    while i < 100
        node = TiSR.grow_equation(rand(3:7), ops)
        data_vect[end], valid = TiSR.eval_equation(node, data_vect, ops)
        valid || continue

        list_of_param = TiSR.list_of_param_nodes(node)
        !isempty(list_of_param) || continue

        for n in list_of_param
            n.val = n.val * (1 + (randn() * 0.1))
        end

        prediction, valid = TiSR.eval_equation(node, data_vect, ops)
        valid || continue
        mse_before = sum(abs2, prediction .- data_vect[end])
        0 < mse_before < Inf || continue

        TiSR.fitting_LM!(node, data_vect, ops, list_of_param, ops.fitting.max_iter)

        prediction, valid = TiSR.eval_equation(node, data_vect, ops)
        valid || continue

        mse_after = sum(abs2, (prediction .- data_vect[end]))

        impr = mse_after / mse_before

        push!(improvements, impr)
        i += 1
    end

    @test maximum(improvements) <= 1.0
    @test count(<(1e-1), improvements) > 70
end

