
# TODO: test copy_wo_copy
# TODO: does this copy with the same type?
# TODO: test list of param nodes
# TODO: test encode_single_char(node::Node, ops)

# TODO:  test deepcopy and copy

include("hardcoded_equations.jl")

data = rand(100, 10)
ops, data_vect = Options(data)

@testset "convert_node" begin
    node = nothing
    while true
        node = TiSR.grow_equation(5, ops)
        !isempty(TiSR.list_of_param_nodes(node)) && break
    end

    @test node isa Node{Float64}
    node2 = TiSR.convert_node(node, 1f0)
    @test node2 isa Node{Float32}
    node3 = TiSR.convert_node(node2, ForwardDiff.Dual(1.0))
    @test node3 isa Node{ForwardDiff.Dual{Nothing, Float64, 0}}
end

@testset "bad_array" begin
    @test !TiSR.bad_array(rand(5))
    @test !TiSR.bad_array(ones(5))
    @test !TiSR.bad_array(zeros(5))

    a = [ForwardDiff.Dual(1.0) for _ in 1:5]
    @test !TiSR.bad_array(a)

    a = [rand(5)..., Inf]
    @test TiSR.bad_array(a)

    a = [rand(5)..., NaN]
    @test TiSR.bad_array(a)
end

@testset "eval_equation" begin
    eqs_dict = hardcoded_equations(ops)

    data = [rand(100) for _ in 1:ops.data_descript.n_vars]

    eq1(data) = @. ((1.45454545 * -42.0) + sin(1.0)) + data[1] * 0.0
    eq2(data) = @. ((1.45454545 * -42.0) - exp((1.0 ^ data[1])))
    eq3(data) = @. ((1.45454545 * exp((1.0 ^ data[1]))) - exp((1.0 ^ data[1])))
    eq4(data) = @. exp((((data[1] * 0.28282828) + (data[8] ^ 0.9984357489357)) * sin(((data[4] * 19.9832754918) - (data[9] / 9.42394578e11)))))
    eq5(data) = @. exp(((sin(((data[4] * 19.9832754918) - (data[9] / 9.42394578e11))) + (data[8] ^ 0.9984357489357)) * sin(((data[4] * 19.9832754918) - (data[9] / 9.42394578e11)))))

    @test eq1(data) == TiSR.eval_equation(eqs_dict.vals[1], data, ops)[1]
    @test eq2(data) == TiSR.eval_equation(eqs_dict.vals[2], data, ops)[1]
    @test eq3(data) == TiSR.eval_equation(eqs_dict.vals[3], data, ops)[1]
    @test eq4(data) == TiSR.eval_equation(eqs_dict.vals[4], data, ops)[1]
    @test eq5(data) == TiSR.eval_equation(eqs_dict.vals[5], data, ops)[1]

    # TODO: test eval with bad values -> -5^1.5
    # TODO: test ForwardDiff with eval_equation
    # TODO: test log(-1) etc

end

@testset "node_to_string" begin
    eqs_dict = hardcoded_equations(ops)

    # foreach(display, TiSR.node_to_string.(eqs_dict.vals, Ref(ops), sigdigits=6))
    @test all(TiSR.node_to_string.(eqs_dict.vals, Ref(ops), sigdigits=6) .== [
        "((1.45455*-42.0)+sin(1.0))",
        "((1.45455*-42.0)-exp((1.0^v1)))",
        "((1.45455*exp((1.0^v1)))-exp((1.0^v1)))",
        "exp((((v1*0.282828)+(v8^0.998436))*sin(((v4*19.9833)-(v9/9.42395e11)))))",
        "exp(((sin(((v4*19.9833)-(v9/9.42395e11)))+(v8^0.998436))*sin(((v4*19.9833)-(v9/9.42395e11)))))",
    ])

    # foreach(display, TiSR.node_to_string.(eqs_dict.vals, Ref(ops), sigdigits=15))
    @test all(TiSR.node_to_string.(eqs_dict.vals, Ref(ops), sigdigits=15) .== [
        "((1.45454545*-42.0)+sin(1.0))",
        "((1.45454545*-42.0)-exp((1.0^v1)))",
        "((1.45454545*exp((1.0^v1)))-exp((1.0^v1)))",
        "exp((((v1*0.28282828)+(v8^0.9984357489357))*sin(((v4*19.9832754918)-(v9/9.42394578e11)))))",
        "exp(((sin(((v4*19.9832754918)-(v9/9.42394578e11)))+(v8^0.9984357489357))*sin(((v4*19.9832754918)-(v9/9.42394578e11)))))",
    ])

    @test all(TiSR.node_to_string.(eqs_dict.vals, Ref(ops), sigdigits=3) .== TiSR.node_to_string.(deepcopy.(eqs_dict.vals), Ref(ops), sigdigits=3))
end

@testset "isapprox" begin

    eqs_dict1 = hardcoded_equations(ops)
    eqs_dict2 = hardcoded_equations(ops)

    node_simple = deepcopy(eqs_dict1["1"])
    node_complex = deepcopy(eqs_dict1["5"])

    node_simple.lef.lef.val += 1e-3
    node_complex.lef.lef.rig.rig.val += 1e-5

    eqs_dict2["copy_modif_1"] = node_simple
    eqs_dict2["copy_modif_5"] = node_complex


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

@testset "deepcopy and friends" begin
    eqs_dict = hardcoded_equations(ops)
    @test all(isapprox.(deepcopy.(eqs_dict.vals), eqs_dict.vals))
    @test all(isapprox.(eqs_dict.vals, deepcopy.(eqs_dict.vals)))
end

@testset "count_nodes and friends" begin
    eqs_dict = hardcoded_equations(ops)
    @test all(TiSR.count_nodes.(eqs_dict.vals) .== [6.0, 8.0, 11.0, 17.0, 22.0])
end

@testset "maxim_tree_depth" begin
    eqs_dict = hardcoded_equations(ops)
    @test all(TiSR.maxim_tree_depth.(eqs_dict.vals) .== [3, 4, 5, 6, 7])
end

@testset "eval_equation with the hardcoded_equations" begin
    eqs_dict = hardcoded_equations(ops)

    data = [rand(100) for _ in 1:ops.data_descript.n_vars]

    eq1(data) = @. ((1.45454545 * -42.0) + sin(1.0)) + data[1] * 0.0
    eq2(data) = @. ((1.45454545 * -42.0) - exp((1.0 ^ data[1])))
    eq3(data) = @. ((1.45454545 * exp((1.0 ^ data[1]))) - exp((1.0 ^ data[1])))
    eq4(data) = @. exp((((data[1] * 0.28282828) + (data[8] ^ 0.9984357489357)) * sin(((data[4] * 19.9832754918) - (data[9] / 9.42394578e11)))))
    eq5(data) = @. exp(((sin(((data[4] * 19.9832754918) - (data[9] / 9.42394578e11))) + (data[8] ^ 0.9984357489357)) * sin(((data[4] * 19.9832754918) - (data[9] / 9.42394578e11)))))

    @test eq1(data) == TiSR.eval_equation(eqs_dict.vals[1], data, ops)[1]
    @test eq2(data) == TiSR.eval_equation(eqs_dict.vals[2], data, ops)[1]
    @test eq3(data) == TiSR.eval_equation(eqs_dict.vals[3], data, ops)[1]
    @test eq4(data) == TiSR.eval_equation(eqs_dict.vals[4], data, ops)[1]
    @test eq5(data) == TiSR.eval_equation(eqs_dict.vals[5], data, ops)[1]
end


