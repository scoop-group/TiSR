

# non-dominated sort # -----------------------------------------------------------------------------
dominates(x, y) = all(i -> x[i] <= y[i], eachindex(x)) && any(i -> x[i] < y[i], eachindex(x))

""" Slightly adapted of a version found on Julia Discourse:
    https://discourse.julialang.org/t/fast-non-dominated-sorting/24164
    Calculate the domination scores, i.e., number of pareto frontier, but
    disregard doubles.
"""
function non_dominated_sort(rows)
    fronts   = Vector{Int64}[]
    individs = SVector{length(rows[1])}.(rows)

    domination_score = fill(typemax(Int64), length(rows))
    indices = unique(i -> individs[i], 1:length(individs))
    keepat!(individs, indices)

    while !isempty(individs)
        red = [all(x -> !dominates(x, y), individs) for y in individs]
        next_front = indices[red]
        push!(fronts, next_front)
        deleteat!(indices, red)
        deleteat!(individs, red)
    end

    for (front, inds) in enumerate(fronts)
        domination_score[inds] .= front
    end

    return domination_score
end

""" Calculate the crowding distances.
"""
function crowding_distance(rows)
    n_pnts, n_dims = length(rows), length(rows[1])
    dists = zeros(n_pnts)
    s_inds = zeros(Int64, n_pnts)

    for i in 1:n_dims
        sortperm!(s_inds, rows, by=x->x[i])

        min_val, max_val = rows[s_inds[1]][i], rows[s_inds[end]][i]
        min_val == max_val && continue

        # set the crowding dist for the boundary pnts to infinity
        dists[s_inds[1]] = dists[s_inds[end]] = Inf

        # compute the crowding dist for the remaining points
        for j in 2:n_pnts-1
            dists[s_inds[j]] += (rows[s_inds[j+1]][i] - rows[s_inds[j-1]][i]) / (max_val - min_val)
        end
    end
    return dists
end

""" perfrom binary tournament selection for parent selection
"""
parent_selection(pop) = min(rand(pop), rand(pop))

""" Slightly adapted of a version found on Julia Discourse:
    https://discourse.julialang.org/t/fast-non-dominated-sorting/24164
    Find the first pareto frontier but disregards doubles.
"""
function first_pareto_front(rows)
    individs = SVector{length(rows[1])}.(rows)
    indices = unique(i -> individs[i], 1:length(individs))
    keepat!(individs, indices)
    red = [all(x -> !dominates(x, y), individs) for y in individs]
    return indices[red]
end

# tournament selection # ---------------------------------------------------------------------------
""" Tournament selection. The proprocessing of the fitness may be adapted. The inds passed into
    this function are modified and should not be used afterwards.
"""
function tournament_selection(fitness, inds; tournament_size=5, n_select=10, modify=true)
    n_select >= length(inds) && return inds

    if !modify
        inds = deepcopy(inds)
    end

    selected = Int64[]

    for _ in 1:n_select
        shuffle!(inds)
        _, win_ind = findmax(i -> fitness[i], view(inds, 1:min(tournament_size, length(inds))))
        push!(selected, popat!(inds, win_ind))
    end

    sort!(selected)
    return selected
end

""" Normalize vector of vectors using the median instead of the maximum.
"""
function normalize_objectives(indiv_obj_vals)
    indiv_obj_vals   = reduce(hcat, indiv_obj_vals)'
    indiv_obj_vals .-= minimum(indiv_obj_vals, dims=1)
    indiv_obj_vals ./= median(indiv_obj_vals, dims=1)
    return eachrow(indiv_obj_vals)
end

""" Calculate the relative fitness based on normalized objectives using the
    median instead of the maximum.
"""
function get_relative_fitness(indiv_obj_vals)
    normed  = normalize_objectives(indiv_obj_vals)
    fitness = sum.(normed)
    return -fitness
end

