
using Test
using OrderedCollections
using Statistics
using DataFrames

include("hardcoded_equations.jl")
include("../src/options.jl")
include("../src/node_n_eval_n_utilities.jl")
include("../src/param_fitting.jl")
include("../src/genetic_ops.jl")
include("../src/ParameterEstimation/levenberg_marquardt/levenberg_marquardt.jl")

# make some preparations # ------------------------------------------------------------------------
ln_abs(x) = log(abs(x))
pow_abs(x, y) = abs(x)^y

@testset "residual after fitting" begin  # integration test

    ops, data = Options(
        general=general_params(),
        unaops=(sin, cos, exp, ln_abs),
        binops=(+, -, *, /, pow_abs, ^),
        fitting=fitting_params(
            max_iter=200,
            early_stop_iter=0,
            t_lim=Inf,
            rel_f_tol_5_iter=0.,
            residual_processing=(x, ind) -> x,
        ),
        p_unaops=(0.1, 0.1, 0.1, 0.1),
        p_binops=(1.0, 1.0, 1.0, 0.1, 0.5, 0.5),
        data_descript=data_descript(
            rand(100, 5)
        ),
    );


    improvements = []

    i = 0
    while i < 50
        node = grow_equation(rand(3:7), ops)
        data[end], valid = eval_equation(node, data, ops)

        valid || continue

        list_of_param = list_of_param_nodes(node)

        length(list_of_param) > 0 || continue

        orig_params = deepcopy(getfield.(list_of_param, :val))
        x0 = rand(length(orig_params))
        setfield!.(list_of_param, :val, x0)

        mse_before = mean(abs2, fitting_residual(x0, node, list_of_param, data, eachindex(data[1]), ops)[1])

        isnan(mse_before) && continue

        residual, valid = residual_after_fitting(node, data, ops, list_of_param)

        mse_after = mean(abs2, eval_equation(node, data, ops)[1] .- data[end])

        impr = (mse_before - mse_after) / mse_before

        isnan(impr) && continue

        after_params = deepcopy(getfield.(list_of_param, :val))

        push!(improvements, (impr, norm(after_params .- orig_params, 2)))
        i += 1
    end

    @test median(imp[1] for imp in improvements) > 0.9
    # @test median(imp[2] for imp in improvements) < 1e-3  # norm2 of parameter estimate
end

@testset "lasso_regression" begin




end


