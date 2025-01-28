
data = rand(100, 10)
ops, data_vect = Options(data)

@testset "string_to_node rediscover" begin
    for _ in 1:1000
        node = TiSR.grow_equation(rand(3:7), ops)
        node_str = TiSR.node_to_string(node, ops, sigdigits=3)
        rebuild_node = TiSR.string_to_node(node_str, ops)
        rebuild_node_str = TiSR.node_to_string(rebuild_node, ops, sigdigits=3)
        @test node_str == rebuild_node_str
    end
end

@testset "node_to_symbolic" begin
    counter = 0
    while counter < 1000
        node = TiSR.grow_equation(rand(3:6), ops)
        TiSR.apply_simple_simplifications!(node, ops)

        # skip this node, if invalid -> in the actual algorithm, nodes are simplified after their
        # evaluation as well # ---------------------------------------------------------------------
        res1, valid1 =  TiSR.eval_equation(node, data_vect, ops)
        valid1 || continue

        sym_node = nothing
        try
            sym_node = TiSR.node_to_symbolic(node, ops)
        catch
            @show node
            @assert false
        end

        rebuild_node = TiSR.string_to_node(sym_node, ops)
        res2, valid2 = TiSR.eval_equation(rebuild_node, data_vect, ops)

        # in seldom cases, the node may become invalid after simplification
        # -> log(v2) * log(v2) is valid but log.(v2).^2.0 may be invalid, because we filter
        # negative bases for the ^-operator by choice
        valid2 || continue

        @test (all(isapprox.(res1, res2, rtol=1e-3)) || quantile(abs.(res1 .- res2), 0.80) < 1e-10)
        counter += 1
    end
end




