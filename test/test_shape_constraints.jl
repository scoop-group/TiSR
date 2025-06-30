
pow(x, y) = abs(x)^y

@testset "calc constr vios" begin

    # data
    data_matr = rand(100, 2)
    data_matr[:, 2] .= @. 2.0 * data_matr[:, 1] + 1.0

    # define constraint manually
    monotonicity1 = (func, params) -> begin
        constr_data = [collect(-2.0:0.1:2.0)]
        pred = func(params, constr_data)
        sum(x -> max(0.0, x), diff(pred) / 0.1)
    end

    # prepare options
    eq_constr_vec   = Function[]
    ineq_constr_vec = Function[]
    all_constr_f_select = Function[monotonicity1]

    ops, data_vect = Options(data_matr,
        fitting = fitting_params(;
            max_iter                   = 10,
            rel_f_tol_5_iter           = 1e-2 * 0.01,
            all_constr_f_select        = all_constr_f_select,
            eq_constr                  = eq_constr_vec,
            ineq_constr                = ineq_constr_vec,
            max_mare_for_constr_fit    = 0.1,
            additional_constr_fit_iter = 10,
        ),
    )

    # build the node
    eq = "2.0 * v1 + 1.0"
    node = TiSR.string_to_node(eq, ops)

    manually = sum(x -> max(0.0, x), diff([2.0 * x + 1.0 for x in -2.0:0.1:2.0]) / 0.1)^2

    # calculate the violations
    @test isapprox(TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops), manually)

    # # again with obtain monotonicity functions # TODO: test also the "obtain"-functions
    # all_constr_f_select = [obtain_monotonicity_constr_func_wo_params([collect(-1.9:0.1:2.0)], 1, "smaller", 0.0)] # start inverval 1 later, as not finite diff
    # ops, data = prepare_ops(data_matr, eq_constr_vec, ineq_constr_vec, all_constr_f_select, 0.1);
    #
    # @test isapprox(TiSR.eval_constr_vios(node, ops), sum(i^2 for i in diff([2.0 * x for x in -2.0:0.1:2.0])/0.1))

    # ==================================================================================================
    #
    # ==================================================================================================
    data_matr = rand(100, 2)
    data_matr[:, 2] .= @. 2.0 * data_matr[:, 1] + 1.0

    monotonicity1 = (func, params) -> begin
        constr_data = [collect(-2.0:0.1:2.0)]
        pred = func(params, constr_data)
        cons = diff(pred) / 0.1
        sum(x -> max(0.0, x)^2, -cons)
    end

    eq_constr_vec   = Function[]
    ineq_constr_vec = Function[]
    all_constr_f_select = Function[monotonicity1]

    ops, data_vect = Options(data_matr,
        fitting = fitting_params(;
            max_iter                   = 10,
            rel_f_tol_5_iter           = 1e-2 * 0.01,
            all_constr_f_select        = all_constr_f_select,
            eq_constr                  = eq_constr_vec,
            ineq_constr                = ineq_constr_vec,
            max_mare_for_constr_fit    = 0.1,
            additional_constr_fit_iter = 10,
        ),
    )

    # buils the node # ---------------------------------------------------------------------------------
    eq = "2.0 * v1 + 1.0"
    node = TiSR.string_to_node(eq, ops)

    # calculate the violations # -----------------------------------------------------------------------
    @test isapprox(TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops), 0.0)

    # # again with obtain monotonicity functions # TODO: test also "obtain"-functions
    # all_constr_f_select = [obtain_monotonicity_constr_func_wo_params([collect(-1.9:0.1:2.0)], 1, "larger", 0.0)] # start inverval 1 later, as not finite diff
    # ops, data = prepare_ops(data_matr, eq_constr_vec, ineq_constr_vec, all_constr_f_select, 0.1);
    #
    # @test isapprox(TiSR.eval_constr_vios(node, ops), 0.0)

    # ==================================================================================================
    #
    # ==================================================================================================
    data_matr = rand(100, 2)
    data_matr[:, 2] .= @. 2.0 * sin(data_matr[:, 1]) + 1.0

    # constraints # ------------------------------------------------------------------------------------
    monotonicity1 = (func, params) -> begin
        constr_data = [collect(-2.0:0.1:2.0)]
        pred = func(params, constr_data)
        cons = diff(pred) / 0.1
        sum(x -> max(0.0, x), -cons)
    end

    eq_constr_vec   = Function[]
    ineq_constr_vec = Function[]
    all_constr_f_select = Function[monotonicity1]

    ops, data_vect = Options(data_matr,
        fitting = fitting_params(;
            max_iter                   = 10,
            rel_f_tol_5_iter           = 1e-2 * 0.01,
            all_constr_f_select        = all_constr_f_select,
            eq_constr                  = eq_constr_vec,
            ineq_constr                = ineq_constr_vec,
            max_mare_for_constr_fit    = 0.1,
            additional_constr_fit_iter = 10,
        ),
    )

    # buils the node # ---------------------------------------------------------------------------------
    eq = "2.0 * sin(v1) + 1.0"
    node = TiSR.string_to_node(eq, ops)

    # calculate the violations # -----------------------------------------------------------------------
    manually = sum(x -> max(0.0, x / 0.1), -diff([2.0 * sin(x) for x in -2.0:0.1:2.0]))^2

    @test isapprox(TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops), manually)

    # # again with obtain monotonicity functions # TODO: test also "obtain"-functions
    # all_constr_f_select = [obtain_monotonicity_constr_func_wo_params([collect(-1.9:0.1:2.0)], 1, "larger", 0.0)] # start inverval 1 later, as not finite diff
    # ops, data = prepare_ops(data_matr, eq_constr_vec, ineq_constr_vec, all_constr_f_select, 0.1);
    #
    # @test isapprox(TiSR.eval_constr_vios(node, ops), sum(x -> max(0.0, x / 0.1)^2, -diff([2.0 * sin(x) for x in -2.0:0.1:2.0])), rtol=0.1) # high rtol, as finite diff is inferior

    # ==================================================================================================
    #
    # ==================================================================================================

    data_matr = rand(100, 2)
    data_matr[:, 2] .= @. 2.0 * sin(data_matr[:, 1]) + 1.0 + data_matr[:, 1]^2

    # constraints # ------------------------------------------------------------------------------------
    monotonicity1 = (func, params) -> begin
        constr_data = [collect(-2.0:0.1:2.0)]
        pred = func(params, constr_data)
        cons = diff(pred) / 0.1
        sum(x -> max(0.0, x), -cons)
    end

    curvature = (func, params) -> begin
        constr_data = [collect(-2.0:0.1:2.0)]
        pred = func(params, constr_data)
        curv = diff(diff(pred) / 0.1) / 0.1
        sum(abs, 2.0 .- curv)
    end

    eq_constr_vec   = Function[]
    ineq_constr_vec = Function[]
    all_constr_f_select = Function[monotonicity1, curvature]

    # prepare options # --------------------------------------------------------------------------------
    ops, data_vect = Options(data_matr,
        binops = (+,   -,   *,   /,   ^, pow),
        fitting = fitting_params(;
            max_iter                   = 10,
            rel_f_tol_5_iter           = 1e-2 * 0.01,
            all_constr_f_select        = all_constr_f_select,
            eq_constr                  = eq_constr_vec,
            ineq_constr                = ineq_constr_vec,
            max_mare_for_constr_fit    = 0.1,
            additional_constr_fit_iter = 10,
        ),
    )

    # buils the node # ---------------------------------------------------------------------------------
    eq = "2.0 * sin(v1) + 1.0 + pow(v1, 2)"
    node = TiSR.string_to_node(eq, ops)

    # calculate the violations # -----------------------------------------------------------------------
    eq_constr_vio   = sum(abs, 2.0 .- diff(diff([2.0 * sin(x) + 1.0 + x^2 for x in -2.0:0.1:2.0])/0.1)/0.1)
    ineq_constr_vio = sum(xx -> max(0.0, xx), -diff([2.0 * sin(x) + 1.0 + x^2 for x in -2.0:0.1:2.0])/0.1)

    @test isapprox(eq_constr_vio^2 + ineq_constr_vio^2, TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops))

    # ==================================================================================================
    #
    # ==================================================================================================

    data_matr = rand(100, 2)
    data_matr[:, 2] .= @. 2.0 * data_matr[:, 1] + 1.0 + data_matr[:, 1]^2

    monotonicity1 = (func, params) -> begin
        constr_data = [collect(0.1:0.3:2.0)]
        pred = func(params, constr_data)
        cons = diff(pred) / 0.1
        sum(x -> max(0.0, x)^2, -cons)
    end

    curvature = (func, params) -> begin
        constr_data = [collect(0.1:0.3:2.0)]
        curv = ForwardDiff.derivative(x -> ForwardDiff.derivative(xx -> func(params, [xx]), x), constr_data[1])
        sum(abs2, 2.0 .- curv)
    end

    eq_constr_vec   = Function[]
    ineq_constr_vec = Function[]
    all_constr_f_select = Function[curvature, monotonicity1]

    # prepare options # --------------------------------------------------------------------------------
    ops, data_vect = Options(data_matr,
        binops = (+,   -,   *,   /,   ^, pow),
        fitting = fitting_params(;
            max_iter                   = 10,
            rel_f_tol_5_iter           = 1e-2 * 0.01,
            all_constr_f_select        = all_constr_f_select,
            eq_constr                  = eq_constr_vec,
            ineq_constr                = ineq_constr_vec,
            max_mare_for_constr_fit    = 0.1,
            additional_constr_fit_iter = 10,
        ),
    )

    eq = "2.0 * v1 + 1.0 + pow(v1, 2)"
    node = TiSR.string_to_node(eq, ops)

    @test TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops) == 0.0
end

@testset "fitting with constraints" begin

    data_matr = rand(100, 2) .* 10
    data_matr[:, 2] .= @. 0.1 * data_matr[:, 1]^3 +  12.2 * data_matr[:, 1]^2 + -3.2 * data_matr[:, 1] + -5

    monotonicity1 = obtain_monotonicity_constr_func([collect(0.1:0.5:2.0)], 1, ">", 0.0)

    eq_constr_vec = Function[]
    ineq_constr_vec = Function[monotonicity1]
    all_constr_f_select = Function[monotonicity1]

    # prepare options # --------------------------------------------------------------------------------
    ops, data_vect = Options(data_matr,
        binops = (+,   -,   *,   /,   ^, pow),
        fitting = fitting_params(;
            max_iter                   = 20,
            rel_f_tol_5_iter           = 1e-2 * 0.01,
            all_constr_f_select        = all_constr_f_select,
            eq_constr                  = eq_constr_vec,
            ineq_constr                = ineq_constr_vec,
            max_mare_for_constr_fit    = 1.0,
            additional_constr_fit_iter = 10,
        ),
    )

    # builds the node # --------------------------------------------------------------------------------
    eq = "0.1 * v1^3 +  12.2 * v1^2 + -3.2 * v1 + -5"
    node = TiSR.string_to_node(eq, ops)
    list_of_param = TiSR.list_of_param_nodes(node)
    p_orig = [n.val for n in list_of_param]

    ii = 0
    results_vio = Float64[]
    results_obj = Float64[]
    # println("1000 fits with monotonicity constraint..")
    @time while ii < 1000
        for i in eachindex(p_orig)
            list_of_param[i].val = p_orig[i] + 1.0 * randn()
        end

        pred, valid = TiSR.eval_equation(node, data_vect, ops)
        mare1 = sum(abs2, (pred .- data_vect[end]) ./ data_vect[end])

        valid || continue

        vio1 = TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops)
        vio1 == 0 && continue

        ii += 1

        pred, valid, coutner = TiSR.fit_n_eval!(node, data_vect, ops, 10)

        mare2 = sum(abs2, (pred .- data_vect[end]) ./ data_vect[end])

        vio2 = TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops)

        push!(results_vio, vio2 / vio1)
        push!(results_obj, mare2 / mare1)
    end

    @test count(i -> results_obj[i] > 1.0 && results_vio[i] > 1.0, eachindex(results_obj)) == 0
    @test count(v -> v < 0.1, results_obj) > 800
    @test count(v -> v < 0.1, results_vio) > 800
    @test count(i -> results_obj[i] < 0.1 && results_vio[i] < 0.1, eachindex(results_obj)) > 800

    # reference speed without constraints #  -------------------------------------------------------
    # eq_constr_vec = Function[]
    # ineq_constr_vec = Function[]
    # all_constr_f_select = Function[monotonicity_wo_params]
    #
    # ops, data = prepare_ops(data_matr, eq_constr_vec, ineq_constr_vec, all_constr_f_select, 0.1)
    #
    # # builds the node # --------------------------------------------------------------------------------
    # eq = "0.1 * v1^3 +  12.2 * v1^2 + -3.2 * v1 + -5"
    # node = TiSR.string_to_node(eq, ops)
    # list_of_param = TiSR.list_of_param_nodes(node)
    # p_orig = getfield.(list_of_param, :val)
    #
    # ii = 0
    # results = Float64[]
    # println("1000 fits without monotonicity constraint..")
    # @time while ii < 1000
    #     for i in eachindex(p_orig)
    #         list_of_param[i].val = p_orig[i] + 2.0 * randn()
    #     end
    #     vio1 = TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops)
    #     vio1 == 0 && continue
    #     ii += 1
    #
    #     pred, valid, coutner = TiSR.fit_n_eval!(node, data_vect, ops, 10)
    #     vio2 = TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops)
    #     push!(results, vio2 / vio1)
    # end

    # ==================================================================================================
    # constraints in the way
    # ==================================================================================================
    data_matr = rand(100, 2) .* 10
    data_matr[:, 2] .= @. 0.1 * data_matr[:, 1]^3 +  12.2 * data_matr[:, 1]^2 + -3.2 * data_matr[:, 1] + -5

    monotonicity1 = (func, params) -> begin
        constr_data = typeof(params)[collect(0.1:0.5:2.0)]
        pred = func(params, constr_data)
        cons = diff(pred) / 0.1
        sum(x -> max(0.0, x), cons)
    end

    curvature = (func, params) -> begin
        constr_data = typeof(params)[collect(0.1:0.1:2.0)]
        curv = ForwardDiff.derivative(x -> ForwardDiff.derivative(xx -> func(params, [xx]), x), constr_data[1])
        mean(x -> max(0.0, x), 1.0 .- curv) #.- 1.0
    end

    eq_constr_vec = Function[]
    ineq_constr_vec = Function[monotonicity1, curvature]
    all_constr_f_select = Function[monotonicity1, curvature]

    # prepare options # --------------------------------------------------------------------------------
    ops, data_vect = Options(data_matr,
        binops = (+,   -,   *,   /,   ^, pow),
        fitting = fitting_params(;
            max_iter                   = 20,
            rel_f_tol_5_iter           = 1e-2 * 0.01,
            all_constr_f_select        = all_constr_f_select,
            eq_constr                  = eq_constr_vec,
            ineq_constr                = ineq_constr_vec,
            max_mare_for_constr_fit    = 1.0,
            additional_constr_fit_iter = 5,
        ),
    )

    # builds the node # --------------------------------------------------------------------------------
    eq = "0.1 * v1^3 +  12.2 * v1^2 + -3.2 * v1 + -5"
    node = TiSR.string_to_node(eq, ops)
    list_of_param = TiSR.list_of_param_nodes(node)
    p_orig = getfield.(list_of_param, :val)

    ii = 0
    results = Float64[]
    @time while ii < 1000
        p = p_orig .+ 2.0 * randn(6)
        for i in eachindex(p_orig)
            list_of_param[i].val = p[i]
        end
        vio1 = TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops)
        vio1 == 0 && continue
        ii += 1

        pred, valid, coutner = TiSR.fit_n_eval!(node, data_vect, ops, 10)

        vio2 = TiSR.get_measure_constr_vios(rand(5), rand(5), node, ops)
        vio2 / vio1

        node
        push!(results, vio2 / vio1)
    end

    @test count(<(1e-3), results) > 500
end
