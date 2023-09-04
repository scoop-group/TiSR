
# non-dominated sort # ---------------------------------------------------------
dominates(x, y) = all(i -> x[i] <= y[i], eachindex(x)) && any(i -> x[i] < y[i], eachindex(x))

""" Slightly adapted of a version found on Julia Discourse:
    https://discourse.julialang.org/t/fast-non-dominated-sorting/24164
    Non-dominated sorting, which selects n_select individuals. The selection in the last front is
    made using crowding_distance selection. The fronts are then concatenated to one indices
    vector. If the keysword first_front=true, only the unfiltered first front is returned.
"""
function non_dominated_sort(rows; n_select=200, first_front=false)
    arr = reduce(hcat, rows)

    fronts = Vector{Int64}[]
    individs = SVector{size(arr, 2)}.(eachrow(arr))
    indices = collect(eachindex(individs))

    if first_front
        n_select = 1
    end

    cur_num = 0
    while !isempty(individs) && cur_num < n_select
        red = [all(x -> !dominates(x, y), individs) for y in individs]
        next_front = indices[red]
        cur_num += length(next_front)
        push!(fronts, next_front)
        deleteat!(indices, red)
        deleteat!(individs, red)
    end

    if !first_front && (too_many = cur_num - n_select) > 0
        to_select = length(fronts[end]) - too_many
        inds = crowding_distance_selection(arr, fronts[end], to_select)
        deleteat!(fronts, length(fronts))
        push!(fronts, inds)
    end
    reduce(vcat, fronts)
end

""" Adapted a ChatGPT generated answer (21.04.2023) after I found a fundamental bug in my own
    implementation, which did not account for higher dimensions.
    Calculates the crowing distance between the points.
"""
function crowding_distance(points)
    n_points, n_dims = size(points)
    distances = zeros(n_points)

    for i in 1:n_dims                                                                                # compute the crowding distance for each point and each dimension
        sorted_indices = sortperm(points[:, i])                                                      # sort the points by their i-th objective value
        sorted_points = points[sorted_indices, :]
        min_value, max_value = sorted_points[1, i], sorted_points[end, i]
        min_value == max_value && continue

        distances[sorted_indices[1]] = Inf                                                           # set the crowding distance for the boundary points to infinity
        distances[sorted_indices[end]] = Inf

        for j in 2:n_points-1                                                                        # compute the crowding distance for the remaining points
            distances[sorted_indices[j]] += (sorted_points[j+1, i] - sorted_points[j-1, i]) / (max_value - min_value)
        end
    end
    distances
end

""" Calls the crowing distance function, sets the distances to the minimal value of 1e-9
    because wsample cannot handle 0, and returns the selected indices.
"""
function crowding_distance_selection(arr, inds, n_select)
    dists = crowding_distance(arr[inds, :])
    dists .= max.(dists, 1e-20)                                                                        # wsmample cannot work with 0 or nan
    replace!(dists, NaN => 1e-20)
    wsample(inds, dists, n_select, replace=false)
end

# tournament selection # -------------------------------------------------------
""" Tournament selection. The proprocessing of the fitness may be adapted. The inds passed into
    this function are modified and should not be used afterwards.
"""
function tournament_selection(fitness, inds; tournament_size=5, n_select=10, modify_inds=true, best=true)
    n_select > length(inds) && return inds

    fitness = 1.0 ./ fitness
    # fitness = -1.0 .* fitness

    if !modify_inds
        inds = deepcopy(inds)
    end

    selected = Int64[]

    for _ in 1:n_select
        if best
            shuffle!(inds)
            tournament = view(inds, 1:min(tournament_size, length(inds)))
            winner = argmax(i -> fitness[i], tournament)
        else
            tournament = rand(inds, tournament_size)
            winner = wsample(tournament, fitness[tournament])
        end

        push!(selected, winner)
        filter!(!isequal(winner), inds)
    end

    sort!(selected)
    return selected
end

