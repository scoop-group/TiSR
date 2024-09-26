
@testset "get_max_nodes_per_term" begin

    # make some preparations # ---------------------------------------------------------------------
    data = rand(100, 10)
    ops, data_vect = Options(
        data,
        binops   = (+, -, *, /),
        p_binops = (0.0, 0.0, 1.0, 1.0),
        unaops   = (),
        p_unaops = ()
    )

    for _ in 1:1000
        node = TiSR.grow_equation(3, ops, method=:full)

        max_nodes_before = TiSR.get_max_nodes_per_term(node, ops)

        @assert TiSR.count_nodes(node) == max_nodes_before

        new_node = TiSR.addterm_mutation!(node, ops, subtree_depth=4)

        # this may happen when grow_equation in asym-mode makes a smaller node
        if TiSR.count_nodes(new_node) <= 7 
            continue
        end

        max_nodes_after = TiSR.get_max_nodes_per_term(node, ops)
    end

    # some more manual tests # ---------------------------------------------------------------------
    ops, data_vect = Options(data)

    nodes = TiSR.string_to_node.(
        [
            "(((sin(0.0945) * (0.119 + 0.473)) + 0.103) + exp(0.912))",                                                                     
            "((log(log((v5 - 0.323))) + cos(((0.914 / v7) ^ abs(v2)))) + abs((v3 + 0.422)))",                                               
            "(((abs(0.696) ^ abs(0.583)) + log((v1 / v5))) + abs((abs(v4) / 0.62)))",                                                       
            "((((v1 ^ v7) / (v3 - v3)) + (log(0.688) + (0.675 - v8))) + (cos(v2) / sin(0.253)))",                                           
            "(((abs(log(0.465)) * abs(exp(v5))) + (exp(0.6) - (0.233 + v1))) + 0.243)",                                                     
            "((exp((v9 * v5)) + sin(0.809)) + (0.817 - 0.72))",                                                                             
            "((cos(abs(v7)) + (0.197 ^ sin(0.931))) + log(cos(v2)))",                                                                       
            "(((v6 / 0.688) + abs(v1)) + (0.736 + 0.682))",                                                                                 
            "((((0.132 + sin(v1)) - (exp(0.358) / (v6 ^ 0.17))) + ((exp(0.448) - exp(v5)) + ((v1 / v8) + (0.865 + 0.151)))) + exp(0.195))", 
            "((log(exp(0.808)) + (0.39 / v7)) + abs(cos((v9 - 0.152))))",                                                                   
        ],
        Ref(ops)
    )

    @assert TiSR.get_max_nodes_per_term.(nodes, Ref(ops)) == [6, 7, 5, 7, 7, 4, 4, 3, 6 ,5]
end

# make some preparations # ------------------------------------------------------------------------
@testset "trim_to_max_nodes_per_term!" begin
    """ sufficient tests, but depends on get_max_nodes_per_term (which it is not sufficient)
    """

    data = rand(100, 10)

    ops, data_vect = Options(
        data,
        grammar=TiSR.grammar_params(;max_nodes_per_term=5)
    )

    for _ in 1:1000
        node = TiSR.grow_equation(rand(4:7), ops, method=:full)
        max_nodes_before = TiSR.get_max_nodes_per_term(node, ops)
        TiSR.trim_to_max_nodes_per_term!(node, ops)
        max_nodes_after = TiSR.get_max_nodes_per_term(node, ops)
        @assert max_nodes_after <= 5
    end
end

@testset "check_legal_function_nesting" begin


    # check if any invalid, inspite of empty illegal_dict
    data = rand(100, 10)
    ops, data_vect = Options(data)

    for _ in 1:1000
        node = TiSR.grow_equation(rand(4:7), ops, method=:full)
        @assert TiSR.check_legal_function_nesting(node, ops)
    end

    # check if any invalid, inspite of illegal operation not present in function set
    ops, data_vect = Options(
        data,
        binops=(+, -, *, /),
        p_binops=(1.0, 1.0, 1.0, 1.0),
        grammar=TiSR.grammar_params(;
            illegal_dict = Dict(
                :^ => (lef = (), rig = (+, -, *, /, sin, cos)),
            )
        )
    )

    for _ in 1:1000
        node = TiSR.grow_equation(rand(4:7), ops, method=:full)
        @assert TiSR.check_legal_function_nesting(node, ops)
    end

    # make some manual checks with examples
    ops, data_vect = Options(
        data,
        grammar=TiSR.grammar_params(;
            illegal_dict = Dict(
                :^ => (lef = (), rig = (+, -, *, /, sin, cos)),
                :/ => (lef = (), rig = (+, -)),
            )
        )
    )

    nodes = TiSR.string_to_node.(
        [
            "(((log(0.0149) - (v1 ^ 0.725)) * log((0.899 - v3))) ^ sin(sin((v6 / 0.127))))",
            "log((cos(abs(v6)) - (abs(v6) + exp(v5))))",
            "exp((cos(abs(0.881)) - exp(cos(0.395))))",
            "((log((0.556 / 0.303)) * ((0.629 ^ v1) + (v7 ^ v3))) / ((abs(0.576) + (v8 + 0.329)) / (log(0.0768) - (0.143 / 0.34))))",
            "((sin(log(v7)) / ((v1 / v5) * sin(v6))) + abs((abs(0.779) / (0.16 - v7))))",
            "exp(((abs(v8) / (v2 ^ v6)) - ((0.701 * v5) ^ abs(0.0541))))",
            "((cos(cos(0.931)) * (exp(v9) * log(v3))) / (log((v6 ^ v5)) / abs(abs(v4))))",
            "abs(abs(exp(exp(0.874)))) / ((v1 / v5) * sin(v6))",
            "cos(cos(cos(cos(v2))))",
            "((((v5 * v4) * cos(v5)) / log(sin(0.22))) + log(((v3 * 0.667) / (v1 * v7))))",
        ],
        Ref(ops)
    )

    @assert all([TiSR.check_legal_function_nesting(n, ops) for n in nodes] .== [0, 1, 1, 0, 0, 1, 1, 1, 1, 1])
end                                                                       
