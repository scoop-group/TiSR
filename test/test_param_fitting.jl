
data = rand(100, 5)
ops, data_vect = Options(data)

# TODO: test set_params!(list_of_param, x)
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

        length(list_of_param) > 0 || continue

        orig_params = deepcopy(getfield.(list_of_param, :val))
        x0 = rand(length(orig_params))
        setfield!.(list_of_param, :val, x0)

        mse_before = mean(abs2, TiSR.fitting_residual(x0, node, list_of_param, data_vect, eachindex(data_vect[1]), ops))

        isfinite(mse_before) || continue

        TiSR.fitting_LM!(node, data_vect, ops, list_of_param, ops.fitting.max_iter)
        prediction, valid = TiSR.eval_equation(node, data_vect, ops)

        valid || continue

        mse_after = mean(abs2, prediction .- data_vect[end])

        impr = mse_after / mse_before

        isfinite(impr) || continue

        push!(improvements, impr)
        i += 1
    end
    @test count(<(1e-3), [imp for imp in improvements]) > 70
end

data = rand(100, 5)
ops, data_vect = Options(
    data,
    fitting = fitting_params(;
        max_iter                = 15,
        early_stop_iter         = 0,
        t_lim                   = Inf,
        rel_f_tol_5_iter        = 0.0,
        lasso_factor            = 1e-1,
        pre_residual_processing = (x, ind, ops) -> x,
        residual_processing     = (x, ind, ops) -> x,
    )
)

@testset "fittig_NW!" begin

    # first check if no lasso does no minimize the redunded ones
    trash_param  = Float64[]

    i = 0
    while i < 100
        node = TiSR.grow_equation(rand(3:7), ops)

        # remove all bloat from the equation
        str1 = TiSR.node_to_string(node, ops)
        while true
            TiSR.simplify_unary_of_param!(node)
            TiSR.simplify_binary_of_param!(node)
            TiSR.simplify_binary_across_1_level!(node, ops)
            TiSR.replace_same_subst_n_div!(node, ops)
            str2 = TiSR.node_to_string(node, ops)
            str1 == str2 && break
            str1 = str2
        end

        data_vect[end], valid = TiSR.eval_equation(node, data_vect, ops)

        valid                             || continue
        (TiSR.maxim_tree_depth(node) > 2) || continue

        # add part to the equation
        node_elect = TiSR.random_node(node; mode=2)

        node_elect.lef = Node(
            findfirst(isequal(+), ops.binops),
            node_elect.lef,
            Node(0.0)
        )

        list_of_param = TiSR.list_of_param_nodes(node)

        ind_of_added = findfirst(n -> iszero(n.val), list_of_param)

        for n in list_of_param
            n.val = rand()
        end

        # add noise to data
        data_vect[end] .*= (1.0 .+ randn(length(data_vect[end])) .* 0.01)

        TiSR.fitting_LM!(node, data_vect, ops, list_of_param, 10)
        prediction, valid = TiSR.eval_equation(node, data_vect, ops)
        before_lasso = abs(list_of_param[ind_of_added].val)

        valid || continue

        TiSR.fitting_NW!(node, data_vect, ops, list_of_param, 10)
        prediction, valid = TiSR.eval_equation(node, data_vect, ops)
        after_lasso = abs(list_of_param[ind_of_added].val)

        valid || continue

        push!(trash_param, after_lasso/before_lasso)
        i += 1
    end

    # histogram(clamp.(trash_param, 0.0, 2.0), nbins=20)
    # sort(trash_param)

    @test count(<(1), trash_param) > 50
end

