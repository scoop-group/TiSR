
# TODO: subtree_mutation
# TODO: apply_genetic_operations

data = rand(100, 5)
ops, data_vect = Options(data)

@testset "grow_equation" begin

    # test full method
    for _ in 1:1000
        max_depth = rand(3:7)
        node = TiSR.grow_equation(max_depth, ops, method=:full)
        @assert TiSR.maxim_tree_depth(node) == max_depth
    end

    # test asym method -> test only whether tree does not reach depth of 5 sometimes -> should suffice
    max_depth = 5
    nodes = map(1:10000) do _
        node = TiSR.grow_equation(max_depth, ops, method=:asym)
    end

    @assert all([TiSR.maxim_tree_depth(node) <= max_depth for node in nodes])
    @assert any([TiSR.maxim_tree_depth(node) < max_depth for node in nodes])
end

@testset "random_node" begin

    distribs = map(1:1000) do _
        max_depth = rand(1:7)
        node = TiSR.grow_equation(max_depth, ops, method=:asym)

        obj_ids = map(1:10_000) do _
            rand_node = TiSR.random_node(node)
            obj_id = objectid(rand_node)
        end

        count_map_ = countmap(obj_ids)

        @assert length(count_map_) == TiSR.count_nodes(node)

        max_frequency_diff = maximum(abs, diff(collect(values(count_map_))))

        mean_freq = mean(collect(values(count_map_)))
        outest_lier = maximum(abs, mean_freq - freq for freq in values(count_map_))

        outest_lier / mean_freq
    end

    @assert quantile(distribs, 0.95) < 0.2

    for _ in 1:1000
        max_depth = rand(4:7)
        node = TiSR.grow_equation(max_depth, ops, method=:full)

        obj_ids_mode1 = map(1:1_000) do _
            rand_node = TiSR.random_node(node, mode=1)
            @assert TiSR.maxim_tree_depth(rand_node) > 1
            obj_id = objectid(rand_node)
        end

        obj_ids_mode2 = map(1:1_000) do _
            rand_node = TiSR.random_node(node, mode=2)
            @assert TiSR.maxim_tree_depth(rand_node) > 2
            obj_id = objectid(rand_node)
        end

        @assert TiSR.count_nodes(node) > length(unique(obj_ids_mode1)) > length(unique(obj_ids_mode2))
    end
end

@testset "mutate_left" begin
    """ this test is kind of pointless
    """

    for _ in 1:1000
        node = TiSR.grow_equation(rand(3:7), ops, method = :full)
        rand_node = TiSR.random_node(node, mode=1)

        min_depth = 1
        lefrig = TiSR.mutate_left(rand_node, min_depth)

        if rand_node.ari == 1
            @assert lefrig == true
        end 
    end

end

@testset "point_mutation" begin

    distances = map(1:1000) do _
        node = TiSR.grow_equation(rand(4:7), ops, method=:asym)
        num_nodes_before = TiSR.count_nodes(node)
        str1 = TiSR.encode_single_char(node, ops)

        TiSR.point_mutation!(node, ops)

        num_nodes_after = TiSR.count_nodes(node)
        str2 = TiSR.encode_single_char(node, ops)

        @assert num_nodes_after == num_nodes_after

        Levenshtein()(str1, str2)
    end

    @assert all(d in (0, 1) for d in distances)
    @assert quantile(distances, 0.7) == 1
end

@testset "insert_mutation!" begin

    distances = map(1:1000) do _
        node = TiSR.grow_equation(rand(4:7), ops, method=:asym)

        num_nodes_before = TiSR.count_nodes(node)
        str1 = TiSR.encode_single_char(node, ops)

        TiSR.insert_mutation!(node, ops)

        num_nodes_after = TiSR.count_nodes(node)
        str2 = TiSR.encode_single_char(node, ops)

        cmap1 = countmap(str1)
        cmap2 = countmap(str2)

        # make sure no node is lost during mutation
        for char in keys(cmap1)
            @assert cmap1[char] <= cmap2[char]
        end

        # check if more nodes after mutation
        @assert num_nodes_before < num_nodes_after

        # string distance?
        dist = Levenshtein()(str1, str2)
    end
    @assert all(d in (3, 5, 6) for d in distances)
end

@testset "hoist_mutation!" begin

    distances = map(1:1000) do _
        node = TiSR.grow_equation(rand(4:7), ops, method=:asym)
        TiSR.maxim_tree_depth(node, minim=3) > 2 || return 0

        num_nodes_before = TiSR.count_nodes(node)
        str1 = TiSR.encode_single_char(node, ops)

        TiSR.hoist_mutation!(node, ops)

        num_nodes_after = TiSR.count_nodes(node)
        str2 = TiSR.encode_single_char(node, ops)

        cmap1 = countmap(str1)
        cmap2 = countmap(str2)

        # make sure no node is lost during mutation
        for char in keys(cmap2)
            @assert cmap2[char] <= cmap1[char]
        end

        # check if more nodes after mutation
        @assert num_nodes_before > num_nodes_after

        # string distance?
        dist = Levenshtein()(str1, str2)
    end

    @assert quantile(distances, 0.5) < 8
end

@testset "crossover_mutation!" begin
    for _ in 1:1000
        node1 = TiSR.grow_equation(rand(4:7), ops, method=:asym)
        TiSR.maxim_tree_depth(node1, minim=3) > 2 || return 0
        node2 = TiSR.grow_equation(rand(4:7), ops, method=:asym)
        TiSR.maxim_tree_depth(node2, minim=3) > 2 || return 0

        str1_before = TiSR.encode_single_char(node1, ops)
        str2_before = TiSR.encode_single_char(node2, ops)

        num_nodes1_before = TiSR.count_nodes(node1)
        num_nodes2_before = TiSR.count_nodes(node2)

        TiSR.crossover_mutation!(node1, node2, ops)

        num_nodes1_after = TiSR.count_nodes(node1)
        num_nodes2_after = TiSR.count_nodes(node2)

        @assert num_nodes1_after + num_nodes2_after == num_nodes1_before + num_nodes2_before

        str1_after = TiSR.encode_single_char(node1, ops)
        str2_after = TiSR.encode_single_char(node2, ops)

        str_before = str1_before * str2_before
        str_after = str1_after * str2_after

        cmap1 = countmap(str_before)
        cmap2 = countmap(str_after)

        @assert length(cmap1) == length(cmap2)

        # make sure no node is lost during mutation
        for char in keys(cmap1)
            @assert cmap2[char] == cmap1[char]
        end
    end
end

@testset "subtree_mutation!" begin
    """ no idea how to test that. num_nodes can be higher or lower & edit distance can be 
        arbitrarily high.
    """
end

@testset "addterm_mutation!" begin
    distances = map(1:1000) do _
        node = TiSR.grow_equation(rand(4:7), ops, method=:asym)

        num_nodes_before = TiSR.count_nodes(node)
        str1 = TiSR.encode_single_char(node, ops)

        TiSR.addterm_mutation!(node, ops)

        num_nodes_after = TiSR.count_nodes(node)
        str2 = TiSR.encode_single_char(node, ops)

        cmap1 = countmap(str1)
        cmap2 = countmap(str2)

        # make sure no node is lost during mutation
        for char in keys(cmap1)
            @assert cmap2[char] >= cmap1[char]
        end

        # check if more nodes after mutation
        @assert num_nodes_before < num_nodes_after

        # string distance?
        dist = Levenshtein()(str1, str2)
    end
    @assert all(d in (5, 6, 8, 9, 10, 11, 12) for d in distances)
end






















