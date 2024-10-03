
@testset "non_dominated_sort" begin
    x = 1.0:10.0
    x_data = repeat(x, 10)                                                                           # generate shifted 1/x points
    y_data = reduce(vcat, [1 ./ x .+ offs for offs in 1:10])
    data = [[xx, yy] for (xx, yy) in zip(x_data, y_data)]
    @assert length(data) == 100 && unique(length.(data)) == [2]

    inds = TiSR.non_dominated_sort(data; n_select=100, first_front=false)
    @test length(unique(inds)) == 100
    @test issorted(inds)                                                                             # shold not have changed the ind positions, as 10 fronts with 10 points each are created with increasing offset -> a little quick and dirty, as within each set of 10 points, they have no particular order
    @test all(i in eachindex(x_data) for i in inds)                                                  # return only valid inds

    inds = TiSR.non_dominated_sort(data; n_select=100, first_front=true)
    @test length(unique(inds)) == 10                                                                 # first front consists of 10 points
    @test all(i < 11 for i in inds)
    @test all(i in eachindex(x_data) for i in inds)

    inds = TiSR.non_dominated_sort(data; n_select=40, first_front=false)
    @test length(unique(inds)) == 40
    @test all(i < 41 for i in inds)
    @test all(i in eachindex(x_data) for i in inds)

    inds = TiSR.non_dominated_sort(data; n_select=35, first_front=false)
    @test length(unique(inds)) == 35
    @test all(i in inds for i in 1:30)
    @test all(i in eachindex(x_data) for i in inds)

    # ----------------------------------------------------------------------------------------------
    x_data = repeat(1.0:10., 10)                                                                      # generate a grid of points
    y_data = sort(x_data)
    data = [[xx, yy] for (xx, yy) in zip(x_data, y_data)]
    @assert length(data) == 100 && unique(length.(data)) == [2]

    inds = TiSR.non_dominated_sort(data; n_select=100, first_front=true)
    @test length(unique(inds)) == 1                                                                  # 1 point dominates all
    @test all(i in eachindex(x_data) for i in inds)

    inds = TiSR.non_dominated_sort(data; n_select=100, first_front=false)
    @test length(unique(inds)) == 100
    @test all(i in eachindex(x_data) for i in inds)

    manh_dist(p) = reduce(+, p)

    s_data = data[inds]
    @test issorted([manh_dist(p) for p in s_data])                                                   # the manhattan distance should increase in for each front, as they are diagonal lines
    @test all(i in eachindex(x_data) for i in inds)

    num = 37                                                                                         # separete a front to check whether that works as well -> using the crowding_distance_selection
    inds = TiSR.non_dominated_sort(data; n_select=num, first_front=false)
    @test length(unique(inds)) == num
    @test all(i in eachindex(x_data) for i in inds)

    s_data = data[inds]
    @test issorted([manh_dist(p) for p in s_data])
    @test all(i in eachindex(x_data) for i in inds)

    # ----------------------------------------------------------------------------------------------
    x_data = repeat(1.0:10., 10)                                                                      # generate a grid of points
    y_data = sort(x_data)
    z_data = ones(size(x_data))
    data = [[xx, yy, zz] for (xx, yy, zz) in zip(x_data, y_data, z_data)]
    @assert length(data) == 100 && unique(length.(data)) == [3]

    @test length(TiSR.non_dominated_sort(data; n_select=25, first_front=false)) == 25                      # this tests whether the function can handle if one variable is the same for all points -> there was a division by 0 error in the crowing distance before
end

@testset "crowding_distance_selection" begin
    x = 1.0:10.0
    x_data = repeat(x, 10)                                                                           # generate a grid of points
    y_data = sort(x_data)
    data = [x_data, y_data]
    @assert length(data) == 2 && length.(data) == [100, 100]
    arr = reduce(hcat, data)
    all_inds = collect(1:size(arr, 1))

    num_select = 20
    inds = TiSR.crowding_distance_selection(arr, all_inds, num_select)
    @test length(unique(inds)) == num_select
    @test all(i in all_inds for i in inds)

    dist = TiSR.crowding_distance(arr)
    @test 2 <= count(!isfinite, dist) <= 4                                                           # as we have 2 dimensions, depending on the row sort, there may be 2...4 Inf

    num_select = 40 - 4                                                                              # check whether all the edge points of the 2d grid are selected
    inds = TiSR.crowding_distance_selection(arr, all_inds, num_select)
    @test all(length(intersect([1. 10.], arr[i, :])) >= 1 for i in inds)                             # boundary points should have either a 1 or a 10 as one of their coordinates
    @test all(i in all_inds for i in inds)

    num_select = 60                                                                                  # if more points are requested, are still all boundary points selected
    inds = TiSR.crowding_distance_selection(arr, all_inds, num_select)
    @test count(length(intersect([1. 10.], arr[i, :])) >= 1 for i in inds) == 36
    @test length(unique(inds)) == num_select
    @test all(i in all_inds for i in inds)
end

@testset "tournament_selection" begin
    x = 1.0:10.0
    x_data = repeat(x, 10)                                                                           # generate shifted 1/x points
    y_data = sort(x_data)
    data = [x_data, y_data]
    @assert length(data) == 2 && length.(data) == [100, 100]
    all_inds = collect(eachindex(x_data))

    compl_coef      = 0.1
    n_select        = 10
    tournament_size = 100

    fitness = x_data .+ compl_coef .* y_data
    fitness = 1 ./ fitness

    inds = TiSR.tournament_selection(
        fitness,
        all_inds;
        tournament_size = tournament_size,
        n_select        = n_select,
        modify          = false,
    )

    @test length(unique(inds)) == n_select
    @test all(i in all_inds for i in inds)

    fitness_select = sum(fitness[inds]) ./ length(inds)
    fitness_not_selected = sum(fitness[setdiff(all_inds, inds)]) ./ length(setdiff(all_inds, inds))
    @test fitness_select > fitness_not_selected

    tournament_size = 20

    selection_better = map(1:100) do _
        inds = TiSR.tournament_selection(
            fitness,
            all_inds;
            tournament_size = tournament_size,
            n_select        = n_select,
            modify          = false,
        )

        @test length(unique(inds)) == n_select
        @test all(i in all_inds for i in inds)

        fitness_select = sum(fitness[inds]) / length(inds)
        fitness_not_selected = sum(fitness[setdiff(all_inds, inds)]) / length(setdiff(all_inds, inds))
        fitness_select / fitness_not_selected
    end

    @test count(>(1), selection_better) > 90
end

@testset "niching tests" begin
    x = 11.0:20.0
    x_data = repeat(x, 10)                                                                           # generate shifted 1/x points
    y_data = reduce(vcat, [1 ./ x .+ offs for offs in 11:20])
    data = [[xx, yy] for (xx, yy) in zip(x_data, y_data)]
    @assert length(data) == 100 && unique(length.(data)) == [2]

    inds = TiSR.non_dominated_sort(data; n_select=100, first_front=false)
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
    inds = TiSR.non_dominated_sort(data; n_select=200, first_front=false)
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

    avg = 0.0
    for i in 1:1000
        inds = TiSR.non_dominated_sort(data; n_select=10, first_front=false)
        avg += (count(>(10), inds) - avg) / i
    end
    @test isapprox(avg, 3.333, rtol=0.1) # is slightly higher than 5/15, because the 5 are less crowded

    # round to 2 significant digits
    for i in eachindex(data)
        data[i] = round.(data[i], sigdigits=2)
    end

    @assert length(unique(data)) == 15

    inds = TiSR.non_dominated_sort(data; n_select=10, first_front=false)

    avg = 0.0
    for i in 1:1000
        inds = TiSR.non_dominated_sort(data; n_select=10, first_front=false)
        avg += (count(>(10), inds) - avg) / i
    end
    @test isapprox(avg, 0.0, rtol=0.1)
end

