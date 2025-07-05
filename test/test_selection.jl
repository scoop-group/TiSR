
# TODO: test parent_selection(pop)
# TODO: test normalize_objectives(indiv_obj_vals) -> median, offset
# TODO: test get_relative_fitness(indiv_obj_vals)

# TODO: test perform_population_selection!(pop, ops) -> niching, age, last front sampling, pareto, rank, crowding
# TODO: test perform_hall_of_fame_selection!(hall_of_fame, population, ops) -> niching, deepcopy, pareto

@testset "non_dominated_sort" begin
    x = 1.0:10.0
    x_data = repeat(x, 10)                                                                           # generate shifted 1/x points
    y_data = reduce(vcat, [1 ./ x .+ offs for offs in 1:10])
    data = [[xx, yy] for (xx, yy) in zip(x_data, y_data)]
    @assert length(data) == 100 && unique(length.(data)) == [2]

    inds = TiSR.non_dominated_sort(data)#; n_select=100, first_front=false)
    inds = sortperm(inds)[1:100]
    @test length(unique(inds)) == 100
    @test issorted(inds)                                                                             # shold not have changed the ind positions, as 10 fronts with 10 points each are created with increasing offset -> a little quick and dirty, as within each set of 10 points, they have no particular order
    @test all(i in eachindex(x_data) for i in inds)                                                  # return only valid inds

    inds = TiSR.first_pareto_front(data)
    @test length(unique(inds)) == 10                                                                 # first front consists of 10 points
    @test all(i < 11 for i in inds)
    @test all(i in eachindex(x_data) for i in inds)

    inds = TiSR.non_dominated_sort(data)#; n_select=40, first_front=false)
    inds = sortperm(inds)[1:40]
    @test length(unique(inds)) == 40
    @test all(i < 41 for i in inds)
    @test all(i in eachindex(x_data) for i in inds)

    inds = TiSR.non_dominated_sort(data)#; n_select=35, first_front=false)
    inds = sortperm(inds)[1:35]
    @test length(unique(inds)) == 35
    @test all(i in inds for i in 1:30)
    @test all(i in eachindex(x_data) for i in inds)

    # ----------------------------------------------------------------------------------------------
    x_data = repeat(1.0:10., 10)                                                                      # generate a grid of points
    y_data = sort(x_data)
    data = [[xx, yy] for (xx, yy) in zip(x_data, y_data)]
    @assert length(data) == 100 && unique(length.(data)) == [2]

    inds = TiSR.first_pareto_front(data)
    @test length(unique(inds)) == 1                                                                  # 1 point dominates all
    @test all(i in eachindex(x_data) for i in inds)

    inds = TiSR.non_dominated_sort(data)#; n_select=100, first_front=false)
    inds = sortperm(inds)[1:100]
    @test length(unique(inds)) == 100
    @test all(i in eachindex(x_data) for i in inds)

    manh_dist(p) = reduce(+, p)

    s_data = data[inds]
    @test issorted([manh_dist(p) for p in s_data])                                                   # the manhattan distance should increase in for each front, as they are diagonal lines
    @test all(i in eachindex(x_data) for i in inds)

    num = 37                                                                                         # separete a front to check whether that works as well -> using the crowding_distance_selection
    inds = TiSR.non_dominated_sort(data)#; n_select=num, first_front=false)
    inds = sortperm(inds)[1:num]
    @test length(unique(inds)) == num
    @test all(i in eachindex(x_data) for i in inds)

    s_data = data[inds]
    @test issorted([manh_dist(p) for p in s_data])
    @test all(i in eachindex(x_data) for i in inds)

    # ----------------------------------------------------------------------------------------------
    x_data = repeat(1.0:10., 10)                                                                     # generate a grid of points
    y_data = sort(x_data)
    z_data = ones(size(x_data))
    data = [[xx, yy, zz] for (xx, yy, zz) in zip(x_data, y_data, z_data)]
    @assert length(data) == 100 && unique(length.(data)) == [3]

    inds = TiSR.non_dominated_sort(data)#; n_select=25, first_front=false)
    inds = sortperm(inds)[1:25]
    @test length(inds) == 25                                                                         # this tests whether the function can handle if one variable is the same for all points -> there was a division by 0 error in the crowing distance before
end

@testset "crowding_distance" begin
    # Test case 1: Basic functionality with 2D points
    rows = [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]]
    expected = [Inf, 2.0, Inf]  # Boundary points have Inf, middle point has finite crowding distance
    @test TiSR.crowding_distance(rows) ≈ expected

    # Test case 2: All points have the same value in one dimension
    rows = [[1.0, 2.0], [2.0, 2.0], [3.0, 2.0]]
    expected = [Inf, 1.0, Inf]  # No variability in the second dimension
    @test TiSR.crowding_distance(rows) ≈ expected

    # Test case 3: Single point
    rows = [[1.0, 2.0]]
    expected = [0.0]  # Only one point, no neighbors, no distance
    @test TiSR.crowding_distance(rows) == expected

    # Test case 4: Two points
    rows = [[1.0, 2.0], [3.0, 4.0]]
    expected = [Inf, Inf]  # Boundary points are always Inf
    @test TiSR.crowding_distance(rows) == expected

    # Test case 5: 3D points
    rows = [[1.0, 2.0, 1.0], [2.0, 4.0, 2.0], [3.0, 6.0, 3.0]]
    expected = [Inf, 3.0, Inf]  # Distances calculated considering all dimensions
    @test TiSR.crowding_distance(rows) ≈ expected

    # Test case 6: Uniform points
    rows = [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]
    expected = [0.0, 0.0, 0.0]  # All values are the same, no distance
    @test TiSR.crowding_distance(rows) == expected
end

@testset "niching tests" begin
    x = 11.0:20.0
    x_data = repeat(x, 10)                                                                           # generate shifted 1/x points
    y_data = reduce(vcat, [1 ./ x .+ offs for offs in 11:20])
    data = [[xx, yy] for (xx, yy) in zip(x_data, y_data)]
    @assert length(data) == 100 && unique(length.(data)) == [2]

    num = 100                                                                                         # separete a front to check whether that works as well -> using the crowding_distance_selection
    inds = TiSR.non_dominated_sort(data)#; n_select=num, first_front=false)
    inds = sortperm(inds)[1:num]
    @test length(unique(inds)) == 100
    @test issorted(inds)                                                                             # shold not have changed the ind positions, as 10 fronts with 10 points each are created with increasing offset -> a little quick and dirty, as within each set of 10 points, they have no particular order
    @test all(i in eachindex(x_data) for i in inds)                                                  # return only valid inds

    # add some slighty shifted data -> should be removed by roundeding and niching
    for i in eachindex(data)
        push!(data, data[i] .* 1.01)
    end

    # round to 2 significant digits
    for i in eachindex(data)
        data[i] = round.(data[i], sigdigits=2)
    end

    sort!(data)

    @assert length(unique(data)) == 100
    num = 200                                                                                         # separete a front to check whether that works as well -> using the crowding_distance_selection
    inds = TiSR.non_dominated_sort(data)#; n_select=num, first_front=false)
    inds = sortperm(inds)[1:100]
    @test length(unique(inds)) == 100 # only 100, because only 100 are unique

    # now acutally the niching # -------------------------------------------------------------------
    x = 11.0:20.0
    x_data = repeat(x, 1)                                                                           # generate shifted 1/x points
    y_data = reduce(vcat, [100.0./x .+ offs for offs in 0:0])
    data = [[xx, yy] for (xx, yy) in zip(x_data, y_data)]
    @assert length(data) == 10 && unique(length.(data)) == [2]

    # add data, which is slightly better in one, but significantly worse in the other objective
    for i in 1:5
        data_ = copy(data[1])
        data_ .*= [0.999, 2.0].^i
        push!(data, data_)
    end

    @assert length(data) == 15
    @assert length(unique(data)) == 15

    sort!(data)
    orig_data = deepcopy(data)

    avg = 0.0
    for i in 1:10000
        shuffle!(data)
        fronts = TiSR.non_dominated_sort(data)
        inds = sortperm(fronts)[1:10]
        orig_inds = [findfirst(==(data[i]), orig_data) for i in inds]
        c = count(>(5), orig_inds)
        avg += (c - avg) / i
    end
    @test isapprox(avg, 6.666, rtol=0.1)

    # round to 2 significant digits
    data = [round.(d, sigdigits=2) for d in data]
    sort!(data)
    orig_data = deepcopy(data)

    @assert length(unique(data)) == 15

    # in this test, we make sure that, not all of the additionally added points are added, which
    # are slightly better on one, but significantly worse on the other
    avg = 0.0
    for i in 1:10000
        shuffle!(data)
        fronts = TiSR.non_dominated_sort(data)
        inds = sortperm(fronts)[1:10]
        orig_inds = [findfirst(==(data[i]), orig_data) for i in inds]
        c = count(>(5), orig_inds)
        avg += (c - avg) / i
    end
    @test avg == 9
end

