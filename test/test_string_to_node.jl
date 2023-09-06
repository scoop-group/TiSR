

using Statistics
using Test
using DataFrames
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
    p_binops=(1.0, 1.0, 1.0, 0.1, 0.5), 
    data_descript=data_descript(
        rand(100, 10)
    )
);

@testset "string_to_node" begin

    for _ in 1:1000
        node = grow_equation(rand(3:7), ops)
        node_str = node_to_string(node, ops, sigdigits=3)
        rebuild_node = string_to_node(node_str, ops)
        rebuild_node_str = node_to_string(rebuild_node, ops, sigdigits=3)
        @test node_str == rebuild_node_str
    end

    counter = 0
    while counter < 1000
        node = grow_equation(rand(3:6), ops)

        # skip this node, if invalid -> in the actual algorithm, nodes are simplified after their 
        # evaluation as well # ---------------------------------------------------------------------
        res1, valid1 =  eval_equation(node, data, ops)
        valid1 || continue                               

        sym_node = node_to_symbolic(node, ops)

        rebuild_node = string_to_node(sym_node, ops)
        res2, valid2 =  eval_equation(rebuild_node, data, ops)

        # in seldom cases, the node may become invalid after simplification 
        # -> log(v2) * log(v2) is valid but log.(v2).^2.0 may be invalid, because we filter 
        # negative bases for the ^-operator by choice
        valid2 || continue 

        @test (all(isapprox.(res1, res2, rtol=1e-3)) || quantile(abs.(res1 .- res2), 0.80) < 1e-10)
        counter += 1
    end
end
    


