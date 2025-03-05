
data = rand(100, 10)
ops, data_vect = Options(data)

include("hardcoded_equations.jl") # tests also the Node constructors

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

    # test with many many # ------------------------------------------------------------------------
    data = [rand(100) .- 0.5 for _ in 1:ops.data_descript.n_vars] # also negative values; should not error
    global v1 = data[1]
    global v2 = data[2]
    global v3 = data[3]
    global v4 = data[4]
    global v5 = data[5]
    global v6 = data[6]
    global v7 = data[7]
    global v8 = data[8]
    global v9 = data[9]

    function gen_valid_node(ops, data, compl)
        node = TiSR.grow_equation(5, ops)
        pred, valid = TiSR.eval_equation(node, data, ops)
        if valid
            return node, pred
        else
            return gen_valid_node(ops, data, compl)
        end
    end

    max_diffs = map(1:500) do _
        node, pred = gen_valid_node(ops, data, rand(1:6))
        str = TiSR.node_to_string(node, ops) # TODO: continue here

        julia_pred = eval(Meta.parse("@. " * str))

        a_diff = abs.((julia_pred .- pred))

        if maximum(a_diff) < 1e-14 || any(iszero, julia_pred)
            return maximum(a_diff)
        end
        r_diff = abs.(a_diff ./ julia_pred)

        return maximum(min.(a_diff, r_diff))
    end

    @test count(>(1e-10), max_diffs) < 10

    # TODO: test ForwardDiff with eval_equation
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

@testset "deepcopy and copy" begin
    eqs_dict = hardcoded_equations(ops)

    # test whether node is the same
    @test all(isapprox.(deepcopy.(eqs_dict.vals), eqs_dict.vals))
    @test all(isapprox.(eqs_dict.vals, deepcopy.(eqs_dict.vals)))

    # test whether node is not the same instance -> perform mutation on copyied
    copied = deepcopy.(eqs_dict.vals)
    for node in copied
        TiSR.hoist_mutation!(node, ops)
    end
    @test all(.!isapprox.(copied, eqs_dict.vals))

    # TODO: test whether type stays the same
    # copied = map(copied) do node
    #     TiSR.convert_node(node, 1f0)
    # end
    # copied isa Vector{Node{Float32}}
    #
    # deepcopy.(copied)
    #
    # deepcopy.(eqs_dict.vals)
    #
    # @test all(.!isapprox.(copied, eqs_dict.vals))
end

@testset "copy_node_wo_copy!" begin
    eqs_dict = hardcoded_equations(ops)

    copied = deepcopy.(eqs_dict.vals)
    for node in copied
        TiSR.hoist_mutation!(node, ops)
    end
    @assert all(.!isapprox.(copied, eqs_dict.vals))

    TiSR.copy_node_wo_copy!.(copied, eqs_dict.vals)
    @test all(isapprox.(copied, eqs_dict.vals))
end

@testset "count_nodes" begin
    eqs_dict = hardcoded_equations(ops)
    @test all(TiSR.count_nodes.(eqs_dict.vals) .== [6.0, 8.0, 11.0, 17.0, 22.0])
end

@testset "maxim_tree_depth" begin
    eqs_dict = hardcoded_equations(ops)
    @test all(TiSR.maxim_tree_depth.(eqs_dict.vals) .== [3, 4, 5, 6, 7])
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


    for (i, eq1) in enumerate(eqs_dict1.vals)
        for (j, eq2) in enumerate(eqs_dict2.vals)
            @test !(isapprox(eq1, eq2, rtol=0.0) ⊻ i == j)
        end
    end

    for (i, eq1) in enumerate(eqs_dict1.vals)
        for (j, eq2) in enumerate(eqs_dict2.vals)
            @test !(isapprox(eq1, eq2, rtol=1e-9) ⊻ i == j)
        end
    end

    for (i, eq1) in enumerate(eqs_dict1.vals)
        for (j, eq2) in enumerate(eqs_dict2.vals)
            @test !(isapprox(eq1, eq2, rtol=1e-4) ⊻ (
                i == j || (i == 5 && j == 7)
            ))
        end
    end

    for (i, eq1) in enumerate(eqs_dict1.vals)
        for (j, eq2) in enumerate(eqs_dict2.vals)
            @test !(isapprox(eq1, eq2, rtol=1e-2) ⊻ (
                i == j || (i == 5 && j == 7) || (i == 1 && j == 6)
            ))
        end
    end

    node = deepcopy(eqs_dict1.vals[1])
    node.lef.rig.val = 1e100

    @test isapprox(node, eqs_dict1.vals[1], rtol=Inf)
end

@testset "Base.:(==)" begin
    eqs_dict = hardcoded_equations(ops).vals
    copied = deepcopy.(eqs_dict)

    copied[1].lef.lef.val *= 1.1
    copied[2].lef.lef.val *= 1.1
    copied[3].lef.lef.val *= 1.1
    copied[4].lef.lef.lef.rig.val *= 1.1
    copied[5].lef.lef.lef.lef.rig.rig.val *= 1.1

    @assert all(isapprox.(eqs_dict, copied, rtol=Inf))
    @assert all(.!isapprox.(eqs_dict, copied, rtol=0))
    @test all(.!.==(eqs_dict, copied))
end

@testset "list_of_param_nodes" begin
    eqs_dict = hardcoded_equations(ops).vals

    list_of_params = TiSR.list_of_param_nodes.(eqs_dict)
    list_of_params isa Vector{Vector{Node}}

    vals = [[round(n.val, sigdigits=3) for n in list_of_param] for list_of_param in list_of_params]

    @test all(vals .== [
        [1.45, -42.0, 1.0],
        [1.45, -42.0, 1.0],
        [1.45, 1.0, 1.0],
        [0.283, 0.998, 20.0, 9.42e11],
        [20.0, 9.42e11, 0.998, 20.0, 9.42e11],
    ])
end

@testset "list_of_param_nodes" begin
    eqs_dict = hardcoded_equations(ops).vals

    node = eqs_dict[1]
    node.ari = 0

    @assert TiSR.count_nodes(node) == 1
    @assert TiSR.count_nodes(node.lef) > 1

    null_node = Node(0.0)
    TiSR.clean_trash_nodes!(node, null_node)
    @test TiSR.count_nodes(node.lef) == 1
end


