
# TODO: test add_from_bank_of_terms_mutation!(node, ops, bank_of_terms)
# TODO: test apply_genetic_operations!(indivs, ops, bank_of_terms;
# TODO: test rand_mult(;minn=2.0, maxx=10.0, sign_flip_prob=0.05)


data = rand(100, 5)
ops, data_vect = Options(data)

@testset "grow_equation" begin

    # test full method
    for _ in 1:1000
        max_depth = rand(3:7)
        node = TiSR.grow_equation(max_depth, ops, method=:full)
        @test TiSR.maxim_tree_depth(node) == max_depth
    end

    # test asym method -> test only whether tree does not reach depth of 5 sometimes -> should suffice
    max_depth = 5
    nodes = map(1:10000) do _
        node = TiSR.grow_equation(max_depth, ops, method=:asym)
    end

    @test all([TiSR.maxim_tree_depth(node) <= max_depth for node in nodes])
    @test any([TiSR.maxim_tree_depth(node) < max_depth for node in nodes])
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

        # test whether all nodes have been sampled at least ones
        @test length(count_map_) == TiSR.count_nodes(node)

        mean_freq = mean(collect(values(count_map_)))
        outest_lier = maximum(abs, mean_freq - freq for freq in values(count_map_))
        outest_lier / mean_freq
    end

    # test whether the sample distribution accross the nodes is approximately uniform
    @test quantile(distribs, 0.95) < 0.2

    for _ in 1:1000
        max_depth = rand(4:7)
        node = TiSR.grow_equation(max_depth, ops, method=:full)

        obj_ids_mode1 = map(1:1_000) do _
            rand_node = TiSR.random_node(node, mode=1)
            # test the mode parameter
            @test TiSR.maxim_tree_depth(rand_node) > 1
            obj_id = objectid(rand_node)
        end

        obj_ids_mode2 = map(1:1_000) do _
            rand_node = TiSR.random_node(node, mode=2)
            # test the mode parameter
            @test TiSR.maxim_tree_depth(rand_node) > 2
            obj_id = objectid(rand_node)
        end

        # make sure some nodes are never sampled
        @test TiSR.count_nodes(node) > length(unique(obj_ids_mode1)) > length(unique(obj_ids_mode2))
    end
end

@testset "mutate_left" begin
    for _ in 1:1000
        node = TiSR.grow_equation(rand(3:7), ops, method = :full)
        rand_node = TiSR.random_node(node, mode=1)

        min_depth = 1
        lefrig = TiSR.mutate_left(rand_node, min_depth)

        if rand_node.ari == 1
            @test lefrig == true
        end
    end

end

@testset "point_mutation1!" begin
    distances = map(1:1000) do _
        node = TiSR.grow_equation(rand(4:7), ops, method=:asym)
        num_nodes_before = TiSR.count_nodes(node)
        str1 = TiSR.encode_single_char(node, ops)

        TiSR.point_mutation1!(node, ops)

        num_nodes_after = TiSR.count_nodes(node)
        str2 = TiSR.encode_single_char(node, ops)

        @test num_nodes_after == num_nodes_before

        dist = Levenshtein()(str1, str2)
        @test dist in (0, 1)
        dist
    end
    @test quantile(distances, 0.7) == 1
end

@testset "point_mutation2!" begin
    distances = map(1:1000) do _
        node = TiSR.grow_equation(rand(4:7), ops, method=:asym)
        num_nodes_before = TiSR.count_nodes(node)
        str1 = TiSR.encode_single_char(node, ops)
        TiSR.point_mutation2!(node, ops)
        num_nodes_after = TiSR.count_nodes(node)
        str2 = TiSR.encode_single_char(node, ops)
        dist = num_nodes_after - num_nodes_before
        lev = Levenshtein()(str1, str2)
        dist, lev/num_nodes_before
    end

    @test length(unique(first.(distances))) > 5
    @test any(>(0), first.(distances))
    @test any(<(0), first.(distances))
    @test quantile(last.(distances), 0.5) < 0.5
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
            @test cmap1[char] <= cmap2[char]
        end

        # check if more nodes after mutation
        @test num_nodes_before < num_nodes_after

        # string distance?
        dist = Levenshtein()(str1, str2)
    end
    @test all(d in (3, 5, 6) for d in distances)
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
            @test cmap2[char] <= cmap1[char]
        end

        # check if more nodes after mutation
        @test num_nodes_before > num_nodes_after

        # string distance?
        dist = Levenshtein()(str1, str2)
    end

    @test quantile(distances, 0.5) < 8

    # test if smaller than 2-3
    node = TiSR.string_to_node("exp(v1)", ops)
    TiSR.hoist_mutation!(node, ops)
    @test TiSR.node_to_string(node, ops) == "v1"

    node = TiSR.string_to_node("v1 + v2", ops)
    TiSR.hoist_mutation!(node, ops)
    @test TiSR.node_to_string(node, ops) in ("v1", "v2")
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

        @test num_nodes1_after + num_nodes2_after == num_nodes1_before + num_nodes2_before

        str1_after = TiSR.encode_single_char(node1, ops)
        str2_after = TiSR.encode_single_char(node2, ops)

        str_before = str1_before * str2_before
        str_after = str1_after * str2_after

        cmap1 = countmap(str_before)
        cmap2 = countmap(str_after)

        @test length(cmap1) == length(cmap2)

        # make sure no node is lost during mutation
        for char in keys(cmap1)
            @test cmap2[char] == cmap1[char]
        end
    end
end

@testset "subtree_mutation!" begin

    # not the bests test for point_mutation2, but not the worst either
    distances = map(1:1000) do _
        node = nothing
        num_nodes_before = nothing
        while true
            node = TiSR.grow_equation(rand(4:7), ops, method=:asym)
            num_nodes_before = TiSR.count_nodes(node)
            num_nodes_before > 1 && break
        end
        str1 = TiSR.encode_single_char(node, ops)
        TiSR.subtree_mutation!(node, ops)
        num_nodes_after = TiSR.count_nodes(node)
        str2 = TiSR.encode_single_char(node, ops)
        dist = num_nodes_after - num_nodes_before
        lev = Levenshtein()(str1, str2)
        dist, lev/num_nodes_before
    end

    @test length(unique(first.(distances))) > 5
    @test any(>(0), first.(distances))
    @test any(<(0), first.(distances))
    @test quantile(last.(distances), 0.5) > 0.5
end

@testset "insert_times_param_mutation!" begin
    distances = map(1:1000) do _
        node = TiSR.grow_equation(rand(4:7), ops, method=:full)

        num_nodes_before = TiSR.count_nodes(node)
        str1 = TiSR.encode_single_char(node, ops)

        TiSR.insert_times_param_mutation!(node, ops)

        num_nodes_after = TiSR.count_nodes(node)
        str2 = TiSR.encode_single_char(node, ops)

        cmap1 = countmap(str1)
        cmap2 = countmap(str2)

        # make sure no node is lost during mutation
        for char in keys(cmap1)
            @test cmap2[char] >= cmap1[char]
        end

        # check if more nodes after mutation
        @test num_nodes_before < num_nodes_after

        # string distance?
        dist = Levenshtein()(str1, str2)
    end
    @test all(d == 5 for d in distances)
end

@testset "addterm_mutation!" begin
    distances = map(1:1000) do _
        node = TiSR.grow_equation(rand(4:7), ops, method=:asym)

        num_nodes_before = TiSR.count_nodes(node)
        str1 = TiSR.encode_single_char(node, ops)

        TiSR.addterm_mutation!(node, ops, subtree_depth=1)

        num_nodes_after = TiSR.count_nodes(node)
        str2 = TiSR.encode_single_char(node, ops)

        cmap1 = countmap(str1)
        cmap2 = countmap(str2)

        # make sure no node is lost during mutation
        for char in keys(cmap1)
            @test cmap2[char] >= cmap1[char]
        end

        # check if more nodes after mutation
        @test num_nodes_before < num_nodes_after

        # string distance?
        dist = Levenshtein()(str1, str2)
    end
    @test all(d in (5, 6, 8, 9, 10, 11, 12) for d in distances)
end























