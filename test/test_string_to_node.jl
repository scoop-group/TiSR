
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

