

using DataFrames
using Test
using OrderedCollections
using Random
using ForwardDiff
using SymbolicUtils

include("hardcoded_equations.jl")
include("../src/options.jl")
include("../src/node_n_eval_n_utilities.jl")
include("../src/genetic_ops.jl")
include("../src/string_to_node.jl")
include("../src/simplify.jl")
include("../src/write_to_excel.jl")

# make some preparations # ------------------------------------------------------------------------
ops, data = Options(
    general=general_params(),
    unaops=(sin, cos, exp, log),
    binops=(+, -, *, /, ^),
    p_unaops=(0.1, 0.1, 0.1, 0.1),
    p_binops=(1.0, 1.0, 1.0, 0.1, 0.5), data_descript=data_descript(
        rand(100, 10)
    )
);

@testset "string_to_node" begin

    for _ in 1:100
        node = grow_equation(rand(3:7), ops)
        node_str = node_to_string(node, ops, sigdigits=3)
        rebuild_node = string_to_node(node_str, ops)
        rebuild_node_str = node_to_string(rebuild_node, ops, sigdigits=3)
        @test node_str == rebuild_node_str
    end

    for _ in 1:100
        node = grow_equation(rand(3:5), ops)
        sym_node = node_to_symbolic(node, ops)
        rebuild_node = string_to_node(sym_node, ops)
        @test all(isapprox.(eval_equation(node, data, ops)[1], eval_equation(rebuild_node, data, ops)[1], rtol=1e-10))
    end
end



