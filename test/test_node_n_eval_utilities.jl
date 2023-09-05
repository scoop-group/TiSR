

using DataFrames
using Test
using OrderedCollections
using Random
using ForwardDiff

include("hardcoded_equations.jl")
include("../options.jl")
include("../node_n_eval_n_utilities.jl")

# make some preparations # ------------------------------------------------------------------------
ln_abs(x) = log(abs(x))
pow_abs(x, y) = abs(x)^y

ops, data = Options(
    general=general_params(),
    unaops=(sin, cos, exp, ln_abs),
    binops=(+, -, *, /, pow_abs),
    p_unaops=(0.1, 0.1, 0.1, 0.1),
    p_binops=(1.0, 1.0, 1.0, 0.1, 0.5), data_descript=data_descript(
        rand(100, 10)
    )
);

@testset "node_to_string" begin
    eqs_dict = hardcoded_equations(ops)

    # foreach(display, node_to_string.(eqs_dict.vals, Ref(ops), sigdigits=6))
    @test all(node_to_string.(eqs_dict.vals, Ref(ops), sigdigits=6) .== [
        "((1.45455*-42.0)+sin(1.0))",
        "((1.45455*-42.0)-exp(pow_abs(1.0, v1)))",
        "((1.45455*exp(pow_abs(1.0, v1)))-exp(pow_abs(1.0, v1)))",
        "exp((((v1*0.282828)+pow_abs(v8, 0.998436))*sin(((v4*19.9833)-(v9/9.42395e11)))))",
        "exp(((sin(((v4*19.9833)-(v9/9.42395e11)))+pow_abs(v8, 0.998436))*sin(((v4*19.9833)-(v9/9.42395e11)))))",
    ])

    # foreach(display, node_to_string.(eqs_dict.vals, Ref(ops), params=true, sigdigits=15))
    @test all(node_to_string.(eqs_dict.vals, Ref(ops), sigdigits=15) .== [
        "((1.45454545*-42.0)+sin(1.0))",
        "((1.45454545*-42.0)-exp(pow_abs(1.0, v1)))",
        "((1.45454545*exp(pow_abs(1.0, v1)))-exp(pow_abs(1.0, v1)))",
        "exp((((v1*0.28282828)+pow_abs(v8, 0.9984357489357))*sin(((v4*19.9832754918)-(v9/9.42394578e11)))))",
        "exp(((sin(((v4*19.9832754918)-(v9/9.42394578e11)))+pow_abs(v8, 0.9984357489357))*sin(((v4*19.9832754918)-(v9/9.42394578e11)))))",
    ])

    @test all(node_to_string.(eqs_dict.vals, Ref(ops), sigdigits=3, unique_nodes=true) .== [
        "((1.45*-42.0)+sin(1.0))",
        "((1.45*-42.0)-exp(pow_abs(1.0, v1)))",
        "((1.45*exp(pow_abs(1.0, v1)))-{exp({pow_abs(1.0, v1)})})",
        "exp((((v1*0.283)+pow_abs(v8, 0.998))*sin(((v4*20.0)-(v9/9.42e11)))))",
        "exp(((sin(((v4*20.0)-(v9/9.42e11)))+pow_abs(v8, 0.998))*{sin({({(v4*20.0)}-{(v9/9.42e11)})})}))",
    ])

    @test all(node_to_string.(eqs_dict.vals, Ref(ops), sigdigits=3) .== node_to_string.(deepcopy.(eqs_dict.vals), Ref(ops), sigdigits=3))
end

@testset "isapprox" begin

    eqs_dict1 = hardcoded_equations(ops)
    eqs_dict2 = hardcoded_equations(ops)

    node_simple = deepcopy(eqs_dict1["1_simple"])
    node_complex = deepcopy(eqs_dict1["5_complex_directed_acyclic_graph"])

    node_simple.lef.lef.val += 1e-3
    node_complex.lef.lef.rig.rig.val += 1e-5

    eqs_dict2["copy_modif_1_simple"] = node_simple
    eqs_dict2["copy_modif_5_complex_directed_acyclic_graph"] = node_complex


    for (i, eq1) in enumerate(eqs_dict1.vals)                                                        # same same?
        for (j, eq2) in enumerate(eqs_dict2.vals)
            @test !(isapprox(eq1, eq2, rtol=0.0) ⊻ i == j)
        end
    end

    for (i, eq1) in enumerate(eqs_dict1.vals)                                                        # 1e-9 - same
        for (j, eq2) in enumerate(eqs_dict2.vals)
            @test !(isapprox(eq1, eq2, rtol=1e-9) ⊻ i == j)
        end
    end

    for (i, eq1) in enumerate(eqs_dict1.vals)                                                        # 1e-4
        for (j, eq2) in enumerate(eqs_dict2.vals)
            @test !(isapprox(eq1, eq2, rtol=1e-4) ⊻ (
                i == j || (i == 5 && j == 7)
            ))
        end
    end

    for (i, eq1) in enumerate(eqs_dict1.vals)                                                        # 1e-2
        for (j, eq2) in enumerate(eqs_dict2.vals)
            @test !(isapprox(eq1, eq2, rtol=1e-2) ⊻ (
                i == j || (i == 5 && j == 7) || (i == 1 && j == 6)
            ))
        end
    end

    node = deepcopy(eqs_dict1.vals[1])                                                                # and test with one of the simple nodes
    node.lef.rig.val = 1e100

    @test isapprox(node, eqs_dict1.vals[1], rtol=Inf)
end

@testset "copy_node and friends" begin
    eqs_dict = hardcoded_equations(ops)
    @test all(isapprox.(copy_node.(eqs_dict.vals), eqs_dict.vals))
    @test all(isapprox.(eqs_dict.vals, deepcopy.(eqs_dict.vals)))
end

@testset "count_nodes and friends" begin
    eqs_dict = hardcoded_equations(ops)
    @test all(count_nodes.(eqs_dict.vals) .== [6.0, 8.0, 11.0, 17.0, 22.0])
    @test all(getindex.(count_nodes_unique.(eqs_dict.vals), 1) .== [6.0, 8.0, 8.0, 17.0, 15.0])
end

@testset "maxim_tree_depth" begin
    eqs_dict = hardcoded_equations(ops)
    @test all(maxim_tree_depth.(eqs_dict.vals) .== [3, 4, 5, 6, 7])
end

@testset "eval_equation with the hardcoded_equations" begin
    eqs_dict = hardcoded_equations(ops)

    data = [rand(100) for _ in 1:ops.data_descript.n_vars]

    eq1(data) = @. ((1.45454545 * -42.0) + sin(1.0)) + data[1] * 0.0
    eq2(data) = @. ((1.45454545 * -42.0) - exp(pow_abs(1.0, data[1])))
    eq3(data) = @. ((1.45454545 * exp(pow_abs(1.0, data[1]))) - exp(pow_abs(1.0, data[1])))
    eq4(data) = @. exp((((data[1] * 0.28282828) + pow_abs(data[8], 0.9984357489357)) * sin(((data[4] * 19.9832754918) - (data[9] / 9.42394578e11)))))
    eq5(data) = @. exp(((sin(((data[4] * 19.9832754918) - (data[9] / 9.42394578e11))) + pow_abs(data[8], 0.9984357489357)) * sin(((data[4] * 19.9832754918) - (data[9] / 9.42394578e11)))))

    @test eq1(data) == eval_equation(eqs_dict.vals[1], data, ops, data[1][1])[1]
    @test eq2(data) == eval_equation(eqs_dict.vals[2], data, ops, data[1][1])[1]
    @test eq3(data) == eval_equation(eqs_dict.vals[3], data, ops, data[1][1])[1]
    @test eq4(data) == eval_equation(eqs_dict.vals[4], data, ops, data[1][1])[1]
    @test eq5(data) == eval_equation(eqs_dict.vals[5], data, ops, data[1][1])[1]

    # @btime eq1(data);
    # @btime eq2(data);
    # @btime eq3(data);

    # t1 = @benchmark eval_equation(eqs_dict_.vals[1], data, ops, data[1][1])[1];
    # t2 = @benchmark eval_equation(eqs_dict_.vals[2], data, ops, data[1][1])[1];
    # t3 = @benchmark eval_equation(eqs_dict_.vals[3], data, ops, data[1][1])[1];
    # Ts = [t1, t2, t3]
    # println("Execution times of the three hardcoded equations -> ", [round(mean(t.times) * 1e-3, sigdigits=3) for t in Ts], "μs")
end


