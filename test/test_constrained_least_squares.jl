
# set parameters

n_points                   = 100
b_track_iter               = 3
g_tol                      = 1e-4
max_iter                   = 1000
init_constr_penalty_factor = 1e-1
incr_constr_penalty_factor = 1.5
obj_mare_tol               = Inf
constr_tol                 = 1e-6
x_atol                     = 1e-3

@testset "fitting a straight line" begin

    # prepare data for first set of test # ---------------------------------------------------------
    X      = collect(-5:0.1:5)
    p_orig = [2.0, -0.2]
    func   = p -> @. p[1] * X[:, 1] + p[2]
    y      = func(p_orig)

    minim = p -> abs.(func(p) .- y)

    # test no 1 # ----------------------------------------------------------------------------------
    # constraints not in the way to fitting
    # y-offset must be lower than 0
    eq_constr_vec   = Function[]
    ineq_constr_vec = Function[p -> p[2]]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i

        x0 = rand(2) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = constr_tol,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([2.0, -0.2], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 2 # ----------------------------------------------------------------------------------
    eq_constr_vec   = Function[]
    ineq_constr_vec = Function[p -> -p[2]]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(2) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = constr_tol,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([2.0, 0.0], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 3 # ----------------------------------------------------------------------------------
    eq_constr_vec   = Function[p -> p[1] + p[2]]
    ineq_constr_vec = Function[]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(2) .* 100.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = constr_tol,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([1.8105263157894735, -1.810526315789474], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 4 # ----------------------------------------------------------------------------------
    eq_constr_vec   = Function[p -> p[1] + p[2]]
    ineq_constr_vec = Function[p -> -p[2]]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(2) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = constr_tol,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([0.0, 0.0], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 5 # ---------------------------------------------------------------------------------- # TODO: numeric problems
    eq_constr_vec   = Function[p -> p[1]^2 + p[2]^2 - 3.0]
    ineq_constr_vec = Function[p -> p[1]]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(2) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = constr_tol,
            x_atol                     = 1e-12,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )

        trace

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([5.419597569105593e-13, 1.7320508075687382], x, atol=1e-1) || isapprox([0.0, -1.732], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 6 # ----------------------------------------------------------------------------------
    eq_constr_vec   = Function[p -> p[1]^2 + p[2]^2 - 3.0]
    ineq_constr_vec = Function[p -> -p[1], p -> p[2]]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(2) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = 0.0,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([1.7299171646732987, -0.08594535108067346], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8
end

@testset "fitting cubic function" begin

    # prepare data for first set of test # ---------------------------------------------------------
    X      = collect(-5:0.1:5.0)
    p_orig = [1.0, -2.2, 3.2, -5.0]
    func   = p -> @. p[1] * X[:, 1]^3 +  p[2] * X[:, 1]^2 + p[3] * X[:, 1] + p[4]
    y      = func(p_orig)

    minim = p -> abs.(func(p) .- y)

    # test no 1 # ----------------------------------------------------------------------------------
    # constraints not in the way to fitting
    eq_constr_vec   = Function[]
    ineq_constr_vec = Function[p -> p[2], p -> p[4]]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(4) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = 0.0,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox(p_orig, x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 2 # ----------------------------------------------------------------------------------
    # some larger 0, which originally are not
    eq_constr_vec   = Function[]
    ineq_constr_vec = Function[p -> -p[2], p -> -p[4]]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(4) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = 0.0,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([1.0000000000000002, 4.0986427954591977e-16, 3.199999999999997, 3.896590271698892e-15], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 3 # ----------------------------------------------------------------------------------
    # constraints slightly in the way
    # parameters must opposites
    eq_constr_vec   = Function[p -> p[1] + p[2]]
    ineq_constr_vec = Function[]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(4) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = 0.0,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([1.1588983050847461, -1.1588983050847461, 0.7691737288135527, -13.849364406779658], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 4 # ----------------------------------------------------------------------------------
    eq_constr_vec   = Function[p -> p[1] + p[2]]
    ineq_constr_vec = Function[p -> -p[2]]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(4) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = 0.0,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([4.1977519351613885e-14, -1.9014277143923402e-14, 18.49799999999919, -23.699999999999214], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 5 # ----------------------------------------------------------------------------------
    eq_constr_vec   = Function[p -> p[1]^2 + p[2]^2 - 3.0]
    ineq_constr_vec = Function[p -> -p[1]]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(4) .* 100.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = 0.0,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor, # TODO: incr_constr lower than 2 might be better?
        )

        push!(mean_iters, trace[end].iter)

        (
           isapprox([0.9285, -1.462, 4.294, -11.27], x, atol=1e-1)
        || isapprox([1.424, -0.9865, -3.281, -15.31], x, atol=1e-1)
        || isapprox([1.614, 0.6275, -6.199, -29.03], x, atol=1e-1)
        || isapprox([1.38, 1.047, -2.607, -32.6], x, atol=1e-1)
        || isapprox([1.143, -1.302, 1.018, -12.64], x, atol=1e-1)

        )
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 6 # ----------------------------------------------------------------------------------
    eq_constr_vec   = Function[p -> p[1] + p[2]^2 + p[3]^2 + p[4]^2 - 10.0]
    ineq_constr_vec = Function[]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    # currently very slow to reach the final point
    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(4) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = 0.0,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor, # TODO: incr_constr lower than 2 might be better?
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([1.1005226505638641, -2.420918922553285, 1.362336184502689, -1.0875058805610434], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 7 # ----------------------------------------------------------------------------------
    eq_constr_vec = Function[]

    monotonicity1 = x -> begin
        pred = func(x)
        cons = pred[1:end-1] .- pred[2:end]
        -cons
    end

    ineq_constr_vec = Function[monotonicity1]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(4) .* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = 0.0,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor, # TODO: incr_constr lower than 2 might be better?
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox([-2.3960339765925e-5, -0.0001332847686767352, 0.001731037601924186, -23.69125424650505], x, atol=1e-1)
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 8

    # test no 8 # ----------------------------------------------------------------------------------
    eq_constr_vec = Function[]

    monotonicity1 = x -> begin
        pred = func(x)
        cons = pred[1:end-1] .- pred[2:end]
        cons
    end

    ineq_constr_vec = Function[monotonicity1]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:10) do i
        x0 = rand(4) #.* 10.0

        x, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(x0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = 0.0,
            x_atol                     = x_atol,
            g_tol                      = 1e-8,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor, # TODO: incr_constr lower than 2 might be better?
        )

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || (
            isapprox([1.0, -2.2, 3.2000000000000006, -5.000000000000001], x, atol=1e-1)
            || isapprox([3.5769606938293808, -5.114418373171924, 1.485214606483339, 22.009465330410876], x, atol=1e-1)
            || isapprox([3.241084088125545, -4.614077573743877, 2.1765168088174494, 15.495641309407556], x, atol=1e-1)
            || isapprox([4.116947940620407, -6.093320476145818, 2.9537365068500696, 28.088160168328447], x, atol=1e-1)
            || isapprox([1.5332621952965348, -2.616853047041113, 1.4448552873842777, -1.4846424037383907], x, atol=1e-1)
            || isapprox([1.8760020130095463, -3.252020211077348, 1.9550401193051212, 3.9202247209475063], x, atol=1e-1)
        )
    end
    # println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 1 # there are so many solution

    # test no 10 # ----------------------------------------------------------------------------------
    # curvature = x -> begin
    #     pred = func(x)
    #     slope = pred[1:end-1] .- pred[2:end]
    #     curvature = slope[1:end-1] .- slope[2:end]
    #     curvature .+ 0.2
    # end
    #
    # eq_constr_vec = Function[]
    #
    # ineq_constr_vec = Function[curvature]
    # eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    # ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])
    #
    # mean_iters = Int64[]
    # results = map(1:10) do i
    #     x0 = rand(4) .* 10.0
    #     x, trace, stop_msg = ip_gauss_newton(minim, eq_constr, ineq_constr, copy(x0), max_iter=max_iter, g_tol=g_tol, b_track_iter=b_track_iter, constr_penalty_factor=constr_penalty_factor)
    #     push!(mean_iters, trace[end].iter)
    #     stop_msg == "converged - KKT point" || isapprox([9.269080469790372e-12, -9.999999999995882, 18.498002256625984, 61.299999999979434], x, atol=1e-1)
    # end
    # println("mean_iters -> ", mean(mean_iters))
    #
    # @test count(results) > 8
end

@testset "fitting cubic function - nonlinear" begin


    p_orig = [0.1, 3.0, 12.2, 2.0, -3.2, 1.0, -5.0]

    X = rand(100) .* 10

    func = p -> @. p[1] * X^p[2] + p[3] * X^p[4] + p[5] * X^p[6] + p[7]

    y = func(p_orig)

    minim = p -> abs.(func(p) .- y)

    # test no 1 # ----------------------------------------------------------------------------------
    # constraints not in the way to fitting
    eq_constr_vec   = Function[]
    ineq_constr_vec = Function[p -> -p[2], p -> p[7]]
    eq_constr   = x -> reduce(vcat, c(x) for c in eq_constr_vec;   init = eltype(x)[])
    ineq_constr = x -> reduce(vcat, c(x) for c in ineq_constr_vec; init = eltype(x)[])

    mean_iters = Int64[]
    results = map(1:100) do i
        p0 = p_orig .* (1.0 + (randn() * 0.1))

        p, trace, stop_msg, valid = TiSR.penalty_method(
            minim,
            eq_constr,
            ineq_constr,
            copy(p0),
            max_iter                   = max_iter,
            b_track_iter               = b_track_iter,
            constr_tol                 = constr_tol,
            x_atol                     = x_atol,
            g_tol                      = g_tol,
            obj_mare_tol               = obj_mare_tol,
            init_constr_penalty_factor = init_constr_penalty_factor,
            incr_constr_penalty_factor = incr_constr_penalty_factor,
        )
        trace
        stop_msg

        push!(mean_iters, trace[end].iter)
        stop_msg == "converged - KKT point" || isapprox(p_orig, p, atol=2e-1)
    end
    println("mean_iters -> ", mean(mean_iters))

    @test count(results) > 50 # many local minima -> hard to check better

end


