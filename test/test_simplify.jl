
# TODO: test drastic_simplify!(node, ops; threshold=1e-1, full=false) # TODO: can go in infinite loop
# TODO: test drastic_simplify_!(node, ops; threshold=1e-1)
# TODO: test get_drastic_simplify_nodes(node, ops; threshold=1e-1)
# TODO: test is_drastic_simplifyable(node, ops; threshold=1e-1)

pow2(x) = x^2
pow_abs(x, y) = abs(x)^y

data = rand(100, 10)
ops, data_vect = Options(
    data,
    binops = (+,   -,   *,   /,   ^, pow_abs),
    unaops = (exp, log, sin, cos, abs, pow2, sqrt),
)

@testset "simplify_unary_of_param!" begin
    unary_of_param_regex = r"\(\(-?\d+(\.\d+)?\)\)"

    i = 0
    while i < 1000
        node = TiSR.grow_equation(rand(3:5), ops)
        node_str = TiSR.node_to_string(node, ops)
        occursin(unary_of_param_regex, node_str) || continue
        i += 1

        @test TiSR.simplify_unary_of_param!(node)

        node_str = TiSR.node_to_string(node, ops)
        @test !occursin(unary_of_param_regex, node_str)
        @test !TiSR.simplify_unary_of_param!(node)
    end

end

@testset "simplify_binary_of_param!" begin
    binary_of_param1_regex = r"\(\(-?\d+(\.\d+)\)\,\(-?\d+(\.\d+)?\)\)"
    binary_of_param2_regex = r"\(\(-?\d+(\.\d+)?\)[\+\-\*/\^]\(-?\d+(\.\d+)?\)"

    i = 0
    while i < 1000
        node = TiSR.grow_equation(rand(3:5), ops)
        node_str = TiSR.node_to_string(node, ops)
        (occursin(binary_of_param1_regex, node_str) || occursin(binary_of_param2_regex, node_str)) || continue
        i += 1

        @test TiSR.simplify_binary_of_param!(node)

        node_str = TiSR.node_to_string(node, ops)
        @test !(occursin(binary_of_param1_regex, node_str) || occursin(binary_of_param2_regex, node_str))
        @test !TiSR.simplify_binary_of_param!(node)
    end
end

@testset "reorder_add_n_mul!" begin

    node_strs = [
        "abs(log((v2+(0.343))))"
        "sqrt(pow2((exp(v3)^exp((0.848)))))"
        "(sin((pow_abs(v3,v3)-sqrt(v4)))/(abs(log(v4))-(sqrt((0.394))^abs(v3))))"
        "abs(sqrt((exp(v1)/abs(v2))))"
        "exp((pow_abs((0.611),v2)*abs((0.375))))"
        "log(((v1*(0.627))+sin((0.92))))"
        "pow2((sin(pow_abs(v2,(0.031)))+(((0.884)-(0.206))-(v2*(0.303)))))"
        "pow_abs(((0.278)+(0.696)),((0.264)/(0.613)))"
        "pow_abs(abs((0.0186)),(v1+(0.671)))"
        "sin(cos(exp(pow_abs(v4,(0.592)))))+(v4+v1)"
    ]

    nodes = [TiSR.string_to_node(n, ops) for n in node_strs]

    for n in nodes
        TiSR.reorder_add_n_mul!(n, ops)
    end

    reordered = [TiSR.node_to_string(n, ops) for n in nodes]

    reordered_right=[
        "abs(log(((0.343)+v2)))"
        "sqrt(pow2((exp(v3)^exp((0.848)))))"
        "(sin((pow_abs(v3,v3)-sqrt(v4)))/(abs(log(v4))-(sqrt((0.394))^abs(v3))))"
        "abs(sqrt((exp(v1)/abs(v2))))"
        "exp((abs((0.375))*pow_abs((0.611),v2)))"
        "log((sin((0.92))+((0.627)*v1)))"
        "pow2((sin(pow_abs(v2,(0.031)))+(((0.884)-(0.206))-((0.303)*v2))))"
        "pow_abs(((0.278)+(0.696)),((0.264)/(0.613)))"
        "pow_abs(abs((0.0186)),((0.671)+v1))"
        "(sin(cos(exp(pow_abs(v4,(0.592)))))+(v1+v4))"
    ]

    reordered .== reordered_right

    @test reordered == reordered_right

    for n in nodes
        @test !TiSR.reorder_add_n_mul!(n, ops)
    end
end

@testset "replace_same_subst_n_div!!" begin
    for _ in 1:100
        node = TiSR.grow_equation(rand(3:5), ops)
        TiSR.count_nodes(node) > 1 || continue
        node_elect = TiSR.random_node(node, mode=1)
        lefrig = TiSR.mutate_left(node_elect, 1) ? :lef : :rig

        which_op = rand((-, /))
        new_node = TiSR.Node(2, findfirst(==(which_op), ops.binops));
        new_node.lef = TiSR.grow_equation(rand(1:5), ops, method=:full)
        new_node.rig = copy(new_node.lef)
        setfield!(node_elect, lefrig, new_node)

        compl_bef = TiSR.count_nodes(node)

        @test TiSR.replace_same_subst_n_div!(node, ops)

        while TiSR.replace_same_subst_n_div!(node, ops)
        end

        @test TiSR.count_nodes(node) < compl_bef
        @test !TiSR.replace_same_subst_n_div!(node, ops)
    end
end

@testset "simplify_binary_across_1_level!" begin

    node_strs = [
        "((abs(((((v2+(1.0))/(1.0))*(1.0))+(1.0)))^log(v1))^pow_abs((v3/v4),pow_abs(v2,v3)))"
        "pow2(exp(cos(((((pow2(v1)/(1.0))*(1.0))*(1.0))))))"
        "sqrt(((((0.799)+v2)+(1.0))-(1.0)))"
    ]

    nodes = [TiSR.string_to_node(n, ops) for n in node_strs]

    for n in nodes
        @test TiSR.simplify_binary_across_1_level!(n, ops)
    end

    simplified = [TiSR.node_to_string(n, ops) for n in nodes]

    right_ones = [
        "((abs((((1.0)*(v2+(1.0)))+(1.0)))^log(v1))^pow_abs((v3/v4),pow_abs(v2,v3)))"
        "pow2(exp(cos(((1.0)*pow2(v1)))))"
        "sqrt(((0.799)-v2))"
    ]

    @test right_ones == simplified
    for n in nodes
        @test !TiSR.simplify_binary_across_1_level!(n, ops)
    end
end

@testset "div_to_mul_param!" begin
    for _ in 1:100
        node = TiSR.grow_equation(rand(3:5), ops)
        TiSR.count_nodes(node) > 1 || continue

        node_elect = TiSR.random_node(node, mode=1)
        lefrig = TiSR.mutate_left(node_elect, 1) ? :lef : :rig
        div_param_node = TiSR.Node(2, findfirst(==(/), ops.binops));
        div_param_node.lef = getfield(node_elect, lefrig)
        div_param_node.rig = TiSR.Node(2.0)
        setfield!(node_elect, lefrig, div_param_node)

        str = TiSR.node_to_string(node, ops)
        @test occursin("/(2.0)", str)

        @test TiSR.div_to_mul_param!(node, ops)
        str = TiSR.node_to_string(node, ops)

        @test !occursin("/(2.0)", str)
        @test occursin("*(0.5)", str)
        @test !TiSR.div_to_mul_param!(node, ops)
    end
end

@testset "apply_simple_simplifications" begin
    # just a test whether the equaitons get shorter on average. All components have are tested
    shorter = map(1:1000) do _
        node = TiSR.grow_equation(rand(4:7), ops)
        bef = TiSR.count_nodes(node)
        TiSR.apply_simple_simplifications!(node, ops)
        bef - TiSR.count_nodes(node)
    end
    @test all(shorter .>= 0)
    @test count(shorter .> 0) > 100
end

data = rand(100, 5)
ops, data_vect = Options(
    data,
    binops   = (+,   -,   *,   /,   ^,),
    unaops   = (exp, log, sin, cos, abs, pow2, sqrt, inv),
)

@testset "drastic_simplify!" begin

    node_strs = [
        "(abs((0.38))-(v4-v1))"
        "(abs((0.38))-1e-7-(v4-v1+1e-7))"
        "sqrt(((((0.98)-v1)*((0.24)*(0.89)))/cos(cos(v3))))"
        "sqrt(((((0.98)-v1)*1e-7*(v2*(0.89)))/cos(cos(v3))))"
        "sqrt(((((0.98)-v1)*((0.24)*(0.89)))/cos(cos(v3))^1e-7))"
        "sin(exp(sqrt((0.25))))"
        "sqrt((cos((v1-v4))+((v4/(0.24))+cos(v3))))"
        "sqrt((cos((v1-v4))+1e-7/((v4/(0.24))+cos(v3))))"
        "sqrt((cos((v1-v4)*(1.0))+((v4/(0.24))+cos(v3))))"
        "sin(((cos((0.49))^(v3+(0.92)))^((v2*v3)^((0.044)+v2))))"
        "sin(((cos((0.49))^(v3+(0.92)))^((v2*v3/(1.0))^((0.044)+v2))))"
        "sin(sqrt(((0.27)-(0.87)))^(1.0))"
        "sin(sqrt(((0.27)-(0.87))))"
        "sin(sqrt((cos(v2)/(v3-v4))))"
        "cos(pow2(pow2((0.56))))"
        "sin(sin((0.72)))"
        "log(log(sqrt((0.3))))"
        "log(log(sqrt((0.3))))+1e-7^log(log(sqrt((0.3))))"
        "log(log(sqrt((0.3))))+(1.0)^log(log(sqrt((0.3))))"
    ]

    nodes = [TiSR.string_to_node(n, ops) for n in node_strs]

    for n in nodes
        TiSR.drastic_simplify!(n, ops; threshold=1e-6, full=true)
    end

    simplified = [TiSR.node_to_string(n, ops) for n in nodes]

    right_ones = [
        "(abs((0.38))-(v4-v1))"
        "(abs((0.38))-(v4-v1))"
        "sqrt(((((0.98)-v1)*((0.24)*(0.89)))/cos(cos(v3))))"
        "sqrt((1.0e-7))"
        "sqrt((((0.98)-v1)*((0.24)*(0.89))))"
        "sin(exp(sqrt((0.25))))"
        "sqrt((cos((v1-v4))+((v4/(0.24))+cos(v3))))"
        "sqrt(cos((v1-v4)))"
        "sqrt((cos((v1-v4))+((v4/(0.24))+cos(v3))))"
        "sin(((cos((0.49))^(v3+(0.92)))^((v2*v3)^((0.044)+v2))))"
        "sin(((cos((0.49))^(v3+(0.92)))^((v2*v3)^((0.044)+v2))))"
        "sin(sqrt(((0.27)-(0.87))))"
        "sin(sqrt(((0.27)-(0.87))))"
        "sin(sqrt((cos(v2)/(v3-v4))))"
        "cos(pow2(pow2((0.56))))"
        "sin(sin((0.72)))"
        "log(log(sqrt((0.3))))"
        "log(log(sqrt((0.3))))"
        "(log(log(sqrt((0.3))))+(1.0))"
    ]

    @test simplified == right_ones
end

@testset "node_to_symbolic" begin

    mse_diffs = Float64[]
    compl_diffs = Float64[]

    while length(mse_diffs) < 10000
        node = TiSR.grow_equation(rand(3:5), ops, method=:full)
        TiSR.apply_simple_simplifications!(node, ops)

        compl = TiSR.count_nodes(node)
        compl > 5 || continue

        pred, valid = TiSR.eval_equation(node, data_vect, ops)
        valid     || continue

        mse = mean(abs2, pred .- data_vect[end])

        symbolic = nothing
        try
            symbolic = TiSR.node_to_symbolic(node, ops)
        catch
            continue
        end
        re_noded = TiSR.string_to_node(symbolic, ops)

        re_noded_compl = TiSR.count_nodes(re_noded)
        compl_diff     = compl - re_noded_compl
        push!(compl_diffs, compl_diff)

        re_noded_mse   = mean(abs2, TiSR.eval_equation(re_noded, data_vect, ops)[1] .- data_vect[end])
        mse_diff       = mse - re_noded_mse
        push!(mse_diffs, abs(mse_diff))
    end

    @test count(<=(1e-3), mse_diffs) > 8000
end

