
data = rand(100, 10)
ops, data_vect = Options(data)

@testset "is_legal_nesting" begin

    # check if any invalid, inspite of empty illegal_dict
    data = rand(100, 10)
    ops, data_vect = Options(data)

    for _ in 1:1000
        node = TiSR.grow_equation(rand(4:7), ops, method=:full)
        @test TiSR.is_legal_nesting(node, ops)
    end

    # check if any invalid, inspite of illegal operation not present in function set
    ops, data_vect = Options(
        data,
        binops=(+, -, *, /),
        grammar=TiSR.grammar_params(;
            illegal_dict = Dict(
                "^" => (lef = (), rig = ("+", "-", "*", "/", "sin", "cos")),
            )
        )
    )

    for _ in 1:1000
        node = TiSR.grow_equation(rand(4:7), ops, method=:full)
        @test TiSR.is_legal_nesting(node, ops)
    end

    # make some manual checks with examples
    ops, data_vect = Options(
        data,
        grammar=TiSR.grammar_params(;
            illegal_dict = Dict(
                "^" => (lef = (), rig = ("+", "-", "*", "/", "sin", "cos")),
                "/" => (lef = (), rig = ("+", "-")),
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

    @test all([TiSR.is_legal_nesting(n, ops) for n in nodes] .== [0, 1, 1, 0, 0, 1, 1, 1, 1, 1])

    # some more manual tests
    nodes = TiSR.string_to_node.(
        [
          "exp(sin(0.533))",
          "abs(abs(abs(v6)))",
          "abs(cos(exp(v6)))",
          "abs(log(exp(v4)))",
          "log(exp(abs(v1)))",
          "abs(sin((v1 ^ v2)))",
          "cos(abs((v4 / v2)))",
          "log(abs((v8 * v4)))",
          "(v8 - sin((v2 - v1)))",
           "sin(log(log(0.0358)))",
           "abs(log(exp(abs(v7))))",
           "exp(log( (cos(v9) * v9)))",
           "log((exp(v3) - cos(v1)))",
           "sin(cos(cos((v6 * v9))))",
           "abs(((v9 - v3) ^ sin(v3)))",
           "abs(abs(abs((v8 + 0.794 ))))",
           "sin(((v9 - v3) + (v2 - v4)))",
           "(cos(cos(v9)) - cos(sin(v7)))",
           "cos(abs(abs((0.13 - 0.658))))",
           "sin((log(v8) - abs(exp( v7))))",
           "log(sin((cos(v6) ^ (v1 ^ v5))))",
           "(log((v3 * v2)) + abs((0.318 * v1)))",
           "exp((abs(v9) / (v5 + (0.493 - v8))))",
           "log( (exp(exp(v3)) ^ log((v8 * v7))))",
           "(log(sin(exp(v4))) / cos(abs(exp(v4))))",
           "(exp(((0.982 * v5) * (v7 / 0.137))) - v2)",
           "sin(( cos((v9 * v4)) / abs((0.221 / v3))))",
           "log(((sin(v5) ^ v1) + (abs(v2) + (v4 * v8))))",
           "(log(sin(0.445)) / ((v8 - v1) / (v2 - 0.932)))",
           "sin((log(log(v5)) + ((v3 - 0.829) + log(v4))))",
           "log((sin(cos(0.143)) * (abs(0.832) + log(v5))))",
           "(cos(((v4 / v6) * v1)) ^ log(exp((v4 + 0.502))))",
           "sin(((abs(0.473) - (v9 / 0.265)) * log(abs(v6))))",
           "(((v2 * v1) * sin(0.429)) - ((v5 + v4) + ( v7 / v7)))",
           "(sin(((v5 * v6) ^ sin(v5))) * cos((sin(0.0645) + v7)))",
           "((cos(0.0623) / (v1 - v3)) * ((v7 * 0.177) * (v2 ^ v1)))" ,
           "((exp(log(0.67)) / cos(abs(v3))) - abs((log(v2) * exp(v7))))",
           "(exp((abs(0.562) ^ log(v5))) / cos(((v4 - v5) / (v7 ^ v4))))",
           "((((v8 - v1) / (v4 ^ v2)) + abs(log(v9))) - cos(sin(abs(v2))))",
           "((cos(exp(v7)) - (log(v5) + (0.618 * v8))) * cos(sin(exp(v8)) ))",
           "(exp((cos(v9) ^ cos(v3))) ^ log(((0.119 * 0.404) + (v1 * v8))))",
           "((log(log(v6)) / ((v2 * v7) / sin(v5))) + sin((sin(v6) + abs(v5))))",
           "(abs(cos(log(v7))) / (exp((v3 * v7)) / ((v1 * 0.124) ^ (v3 + v3))))",
           "((log((v5 - v7)) * (cos(v6) ^ cos(v6))) * ( log(v5) * (v1 / (v3 * v9))))",
           "(log((cos(v2) + (v6 ^ v4))) * ((cos(v2) / (v6 / v5)) / sin((v1 * v4))))",
           "(((sin(0.0138) + v2) - (abs(0.347) - abs(v9))) * exp(log((v8 - 0.594))))",
           "(cos((v3 + (v4 / v1))) ^ ((log(v8) + (v2 / v3)) / ((v5 - v7) ^ exp(v6))))",
           "((((v1 * v6) / cos(0.602)) + (cos(0.837) + sin(v2))) - log(log((v7 + v2))))",
           "(((cos(v3) - log(0.779)) ^ (sin(0.389) / v1)) * (abs((v4 / 0.0874)) * abs((0.88 / v2))))",
           "((cos(log(v2)) ^ ((v1 / 0.297) / log(v8))) / (((v6 ^ v4) + (v5 * 0.131)) + log((v1 /v1))))"
        ],
        Ref(ops)
    )

    # make some manual checks with examples
    ops, data_vect = Options(
        data,
        grammar=TiSR.grammar_params(;
            illegal_dict = Dict(
                "cos" => (lef = ("sin", "cos"), rig = ())
            )
        )
    )

    @test findall([!TiSR.is_legal_nesting(n, ops) for n in nodes]) == [14, 18, 35, 39, 40]

    ops, data_vect = Options(
        data,
        grammar=TiSR.grammar_params(;
            illegal_dict = Dict(
                "cos" => (lef = ("sin", "cos"), rig = ()),
                "sin" => (lef = ("sin", "cos", "PARAM"), rig = ())
            )
        )
    )

    @test findall([!TiSR.is_legal_nesting(n, ops) for n in nodes]) == [1, 10, 14, 18, 21, 27, 29, 30, 31, 33, 34, 35, 39, 40, 42, 46, 49]

    ops, data_vect = Options(
        data,
        grammar=TiSR.grammar_params(;
            illegal_dict = Dict(
                "cos" => (lef = ("sin", "cos"), rig = ()),
                "sin" => (lef = ("sin", "cos", "PARAM"), rig = ()),
                "^" => (lef = (), rig = ("+", "-", "*", "/", "sin", "cos", "VAR")),
            )
        )
    )

    @test findall([!TiSR.is_legal_nesting(n, ops) for n in nodes]) == [1, 6, 10, 14, 15, 18, 21, 24, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 50]
end

@testset "get_max_nodes_per_term" begin

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

    @test TiSR.get_max_nodes_per_term.(nodes, Ref(ops)) == [6, 7, 5, 7, 7, 4, 4, 3, 6 ,5]
end

