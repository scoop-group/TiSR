
# TODO: test simplify_to_string(node::Node, ops::Options; sigdigits=15, use_simplify=false) # TODO: make nicer, this was quick and dirty
# TODO: test apply_simple_simplifications!(node, ops)
# TODO: test div_to_mul_param!(node, ops) # TODO: test
# TODO: test replace_same_subst_n_div!(node, ops)


pow2(x) = x^2
pow_abs(x, y) = abs(x)^y

data = rand(100, 10)
ops, data_vect = Options(
    data,
    binops                             = (+,   -,   *,   /,   ^, pow_abs),
    unaops                             = (exp, log, sin, cos, abs, pow2, sqrt),
)

@testset "simplify_unary_of_param!" begin
    unary_of_param_regex = r"\(-?\d+(\.\d+)?\)"

    i = 0
    while i < 1000
        node = TiSR.grow_equation(rand(3:5), ops)
        node_str = TiSR.node_to_string(node, ops)
        occursin(unary_of_param_regex, node_str) || continue
        i += 1

        TiSR.simplify_unary_of_param!(node)

        node_str = TiSR.node_to_string(node, ops)
        @test !occursin(unary_of_param_regex, node_str)
    end
end

@testset "simplify_binary_of_param!" begin
    binary_of_param1_regex = r"\(-?\d+(\.\d+)\,-?\d+(\.\d+)?\)"
    binary_of_param2_regex = r"\(-?\d+(\.\d+)?[\+\-\*/\^]-?\d+(\.\d+)?\)"

    i = 0
    while i < 1000
        node = TiSR.grow_equation(rand(3:5), ops)
        node_str = TiSR.node_to_string(node, ops)
        (occursin(binary_of_param1_regex, node_str) || occursin(binary_of_param2_regex, node_str)) || continue
        i += 1

        TiSR.simplify_binary_of_param!(node)

        node_str = TiSR.node_to_string(node, ops)
        @test !(occursin(binary_of_param1_regex, node_str) || occursin(binary_of_param2_regex, node_str))
    end
end

@testset "reorder_add_n_mul!" begin

    # nodes = [TiSR.grow_equation(rand(3:5), ops, method=:full) for _ in 1:10]
    # for n in nodes
    #     println(n)
    # end

    node_strs = [
        "(cos(0.156)/(v4-0.714))"
        "(exp(pow_abs(v4,v4))*pow2(log(0.416)))"
        "(log(v3)+pow_abs(0.357,v1))"
        "(pow_abs(0.544,v4)-(0.431/0.526))"
        "abs(log((v2+0.343)))"
        "abs(sin(cos(0.48)))"
        "exp(pow_abs(sin(sin(v3)),pow_abs(log(0.352),cos(v1))))"
        "pow_abs((exp(sqrt(v1))^cos(exp(0.879))),exp(((v4/v4)-(0.0202/0.313))))"
        "sin((exp(0.327)+pow2(v4)))"
        "sqrt(pow2((exp(v3)^exp(0.848))))"
        "(sin((pow_abs(v3,v3)-sqrt(v4)))/(abs(log(v4))-(sqrt(0.394)^abs(v3))))"
        "abs(sqrt((exp(v1)/abs(v2))))"
        "exp((pow_abs(0.611,v2)*abs(0.375)))"
        "log(((v1*0.627)+sin(0.92)))"
        "log(log((v3/0.52)))"
        "pow2((sin(pow_abs(v2,0.031))+((0.884-0.206)-(v2*0.303))))"
        "pow_abs((0.278+0.696),(0.264/0.613))"
        "pow_abs(abs(0.0186),(v1+0.671))"
        "pow_abs(sqrt((pow_abs(0.293,v1)-(v2-0.938))),pow_abs((exp(v4)^(0.664-0.28)),(exp(v1)*(0.326+0.313))))"
        "sin(cos(exp(pow_abs(v4,0.592))))+(v4+v1)"
    ]

    nodes = [TiSR.string_to_node(n, ops) for n in node_strs]

    for n in nodes
        TiSR.reorder_add_n_mul!(n, ops)
    end

    reordered = [TiSR.node_to_string(n, ops) for n in nodes]

    reordered_right=[
        "(cos(0.156)/(v4-0.714))"
        "(exp(pow_abs(v4,v4))*pow2(log(0.416)))"
        "(log(v3)+pow_abs(0.357,v1))"
        "(pow_abs(0.544,v4)-(0.431/0.526))"
        "abs(log((0.343+v2)))"
        "abs(sin(cos(0.48)))"
        "exp(pow_abs(sin(sin(v3)),pow_abs(log(0.352),cos(v1))))"
        "pow_abs((exp(sqrt(v1))^cos(exp(0.879))),exp(((v4/v4)-(0.0202/0.313))))"
        "sin((exp(0.327)+pow2(v4)))"
        "sqrt(pow2((exp(v3)^exp(0.848))))"
        "(sin((pow_abs(v3,v3)-sqrt(v4)))/(abs(log(v4))-(sqrt(0.394)^abs(v3))))"
        "abs(sqrt((exp(v1)/abs(v2))))"
        "exp((abs(0.375)*pow_abs(0.611,v2)))"
        "log((sin(0.92)+(0.627*v1)))"
        "log(log((v3/0.52)))"
        "pow2((sin(pow_abs(v2,0.031))+((0.884-0.206)-(0.303*v2))))"
        "pow_abs((0.278+0.696),(0.264/0.613))"
        "pow_abs(abs(0.0186),(0.671+v1))"
        "pow_abs(sqrt((pow_abs(0.293,v1)-(v2-0.938))),pow_abs((exp(v4)^(0.664-0.28)),(exp(v1)*(0.326+0.313))))"
        "(sin(cos(exp(pow_abs(v4,0.592))))+(v1+v4))"
    ]

    reordered .== reordered_right

    @test reordered == reordered_right
end

@testset "simplify_binary_across_1_level!" begin

    # nodes = [TiSR.grow_equation(rand(3:5), ops, method=:full) for _ in 1:10]
    # for n in nodes
    #     println(n)
    # end

    node_strs = [
        "((0.113*v1)+sin(v3))"
        "((abs(v2)^log(v1))^pow_abs((v3/v4),pow_abs(v2,v3)))"
        "(exp(abs(v3))-sin(pow_abs(v1,v1)))"
        "(exp(v1)+(0.977^v3))"
        "(pow2((v2+v1))*(sin(0.151)^(0.236*0.0362)))"
        "log((cos(0.443)^log(0.23)))"
        "pow2(exp(cos(pow2(v1))))"
        "pow_abs((v4+v2),(v2/v4))"
        "sqrt((0.799+v2))"
        "((abs(((((v2+1.0)/1.0)*1.0)+1.0))^log(v1))^pow_abs((v3/v4),pow_abs(v2,v3)))"
        "pow2(exp(cos(((((pow2(v1)/1.0)*1.0)*1.0)))))"
        "sqrt((((0.799+v2)+1.0)-1.0))"
    ]

    nodes = [TiSR.string_to_node(n, ops) for n in node_strs]

    for n in nodes
        TiSR.simplify_binary_across_1_level!(n, ops)
    end

    simplified = [TiSR.node_to_string(n, ops) for n in nodes]

    right_ones = [
        "((0.113*v1)+sin(v3))"
        "((abs(v2)^log(v1))^pow_abs((v3/v4),pow_abs(v2,v3)))"
        "(exp(abs(v3))-sin(pow_abs(v1,v1)))"
        "(exp(v1)+(0.977^v3))"
        "(pow2((v2+v1))*(sin(0.151)^(0.236*0.0362)))"
        "log((cos(0.443)^log(0.23)))"
        "pow2(exp(cos(pow2(v1))))"
        "pow_abs((v4+v2),(v2/v4))"
        "sqrt((0.799+v2))"
        "((abs(((1.0*(v2+1.0))+1.0))^log(v1))^pow_abs((v3/v4),pow_abs(v2,v3)))"
        "pow2(exp(cos((1.0*pow2(v1)))))"
        "sqrt((0.799-v2))"
    ]

    @test right_ones == simplified
end

data = rand(100, 5)
ops, data_vect = Options(
    data,
    binops   = (+,   -,   *,   /,   ^,),
    unaops   = (exp, log, sin, cos, abs, pow2, sqrt),
)

@testset "drastic_simplify!" begin

    node_strs = [
        "(abs(0.38)-(v4-v1))"
        "(abs(0.38)-1e-7-(v4-v1+1e-7))"
        "sqrt((((0.98-v1)*(0.24*0.89))/cos(cos(v3))))"
        "sqrt((((0.98-v1)*1e-7*(v2*0.89))/cos(cos(v3))))"
        "sqrt((((0.98-v1)*(0.24*0.89))/cos(cos(v3))^1e-7))"
        "sin(exp(sqrt(0.25)))"
        "sqrt((cos((v1-v4))+((v4/0.24)+cos(v3))))"
        "sqrt((cos((v1-v4))+1e-7/((v4/0.24)+cos(v3))))"
        "sqrt((cos((v1-v4)*1.0)+((v4/0.24)+cos(v3))))"
        "sin(((cos(0.49)^(v3+0.92))^((v2*v3)^(0.044+v2))))"
        "sin(((cos(0.49)^(v3+0.92))^((v2*v3/1.0)^(0.044+v2))))"
        "sin(sqrt((0.27-0.87))^1.0)"
        "sin(sqrt((0.27-0.87)))"
        "sin(sqrt((cos(v2)/(v3-v4))))"
        "cos(pow2(pow2(0.56)))"
        "sin(sin(0.72))"
        "log(log(sqrt(0.3)))"
        "log(log(sqrt(0.3)))+1e-7^log(log(sqrt(0.3)))"
        "log(log(sqrt(0.3)))+1.0^log(log(sqrt(0.3)))"
    ]

    nodes = [TiSR.string_to_node(n, ops) for n in node_strs]

    for n in nodes
        TiSR.drastic_simplify!(n, ops; threshold=1e-6, full=true)
    end

    simplified = [TiSR.node_to_string(n, ops) for n in nodes]

    right_ones = [
        "(abs(0.38)-(v4-v1))"
        "(abs(0.38)-(v4-v1))"
        "sqrt((((0.98-v1)*(0.24*0.89))/cos(cos(v3))))"
        "sqrt(1.0e-7)"
        "sqrt(((0.98-v1)*(0.24*0.89)))"
        "sin(exp(sqrt(0.25)))"
        "sqrt((cos((v1-v4))+((v4/0.24)+cos(v3))))"
        "sqrt(cos((v1-v4)))"
        "sqrt((cos((v1-v4))+((v4/0.24)+cos(v3))))"
        "sin(((cos(0.49)^(v3+0.92))^((v2*v3)^(0.044+v2))))"
        "sin(((cos(0.49)^(v3+0.92))^((v2*v3)^(0.044+v2))))"
        "sin(sqrt((0.27-0.87)))"
        "sin(sqrt((0.27-0.87)))"
        "sin(sqrt((cos(v2)/(v3-v4))))"
        "cos(pow2(pow2(0.56)))"
        "sin(sin(0.72))"
        "log(log(sqrt(0.3)))"
        "log(log(sqrt(0.3)))"
        "(log(log(sqrt(0.3)))+1.0)"
    ]

    @test simplified == right_ones
end

@testset "node_to_symbolic" begin

    # TODO: following symbolic do not convert to nodes -> 1//0 ; Inf ;

    mse_diffs = Float64[]
    compl_diffs = Float64[]

    while length(mse_diffs) < 10000

        node = TiSR.grow_equation(rand(3:5), ops, method=:full)

        TiSR.apply_simple_simplifications!(node, ops)

        compl = TiSR.count_nodes(node)

        compl > 5 || continue

        mse       = mean(abs2, TiSR.eval_equation(node, data_vect, ops)[1])

        symbolic = nothing
        try
            symbolic = TiSR.node_to_symbolic(node, ops)
        catch
            # println("----symbolics")
            # @show node
            continue
        end

        re_noded = nothing
        try
            re_noded       = TiSR.string_to_node(symbolic, ops)
        catch
            # println("----re_noded")
            # @show symbolic
            # @show re_noded
            continue
        end

        re_noded_compl = TiSR.count_nodes(re_noded)
        compl_diff     = compl - re_noded_compl

        re_noded_mse   = mean(abs2, TiSR.eval_equation(re_noded, data_vect, ops)[1])
        mse_diff       = mse - re_noded_mse

        push!(compl_diffs, compl_diff)
        push!(mse_diffs, abs(mse_diff))
    end

    @test count(<=(1e-3), mse_diffs) > 8000
end

@testset "simplify_w_symbolic_utils!" begin

    mse_diffs = Float64[]
    compl_diffs = Float64[]

    while length(mse_diffs) < 10000

        node = TiSR.grow_equation(rand(3:5), ops, method=:full)

        TiSR.apply_simple_simplifications!(node, ops)

        compl = TiSR.count_nodes(node)

        compl > 5 || continue

        mse = mean(abs2, TiSR.eval_equation(node, data_vect, ops)[1])
        prev_node = TiSR.node_to_string(node, ops, sigdigits=2)

        try
            TiSR.simplify_w_symbolic_utils!(node, ops; use_simplify=false, though_polyform=false)
        catch
            # @show node
            continue
        end

        re_noded_compl = TiSR.count_nodes(node)
        compl_diff     = compl - re_noded_compl

        # if compl_diff > 5
        #     println("----")
        #     println(prev_node)
        #     println(TiSR.node_to_string(node, ops, sigdigits=2))
        # end

        re_noded_mse   = mean(abs2, TiSR.eval_equation(node, data_vect, ops)[1])
        mse_diff       = mse - re_noded_mse

        push!(compl_diffs, compl_diff)
        push!(mse_diffs, abs(mse_diff))
    end

    @test count(<=(1e-5), mse_diffs) > 8000
end




