
# TODO: test lasso regression
# TODO: test early stopping
# TODO: test pre_residual_processing! = (x, ind) -> x,
# TODO: test residual_processing      = (x, ind) -> x,
# TODO: text t_lim
# TODO: text rel_f_tol_5_iter

# make some preparations # ------------------------------------------------------------------------

data = rand(100, 5)
ops, data_vect = Options(data)

@testset "residual after fitting" begin  # integration test

    improvements = []

    i = 0
    while i < 100
        node = TiSR.grow_equation(rand(3:7), ops)
        data_vect[end], valid = TiSR.eval_equation(node, data_vect, ops)

        valid || continue

        list_of_param = list_of_param_nodes(node)

        length(list_of_param) > 0 || continue

        orig_params = deepcopy(getfield.(list_of_param, :val))
        x0 = rand(length(orig_params))
        setfield!.(list_of_param, :val, x0)

        mse_before = mean(abs2, TiSR.fitting_residual(x0, node, list_of_param, data_vect, eachindex(data_vect[1]), ops)[1])

        isnan(mse_before) && continue

        residual, valid = residual_after_fitting(node, data_vect, ops, list_of_param)

        mse_after = mean(abs2, TiSR.eval_equation(node, data_vect, ops)[1] .- data_vect[end])

        impr = (mse_before - mse_after) / mse_before

        isnan(impr) && continue

        push!(improvements, impr)
        i += 1
    end

    @test median(imp[1] for imp in improvements) > 0.9
end

@testset "lasso_regression" begin




end


