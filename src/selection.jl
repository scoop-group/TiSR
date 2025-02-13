
""" Whether x dominates y.
"""
dominates(x, y) = all(i -> x[i] <= y[i], eachindex(x)) && any(i -> x[i] < y[i], eachindex(x))

""" Slightly adapted of a version found on Julia Discourse:
    [source](https://discourse.julialang.org/t/fast-non-dominated-sorting/24164)
    For each entry, the pareto front number is determined.
"""
function non_dominated_sort(individs::Vector{Vector{Float64}})
    fronts   = Vector{Int64}[]
    indices = collect(eachindex(individs))
    individs = deepcopy(individs)
    domination_score = fill(typemax(Int64), length(individs))

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

""" Find the first pareto front.
"""
first_pareto_front(individs::Vector{Vector{Float64}}) = findall(y -> all(x -> !dominates(x, y), individs), individs)

""" Calculate the crowding distances.
"""
function crowding_distance(individs)
    n_pnts, n_dims = length(individs), length(individs[1])
    dists = zeros(n_pnts)
    s_inds = zeros(Int64, n_pnts)
    for i in 1:n_dims
        sortperm!(s_inds, individs, by=x->x[i])
        min_val, max_val = individs[s_inds[1]][i], individs[s_inds[end]][i]
        min_val == max_val && continue
        dists[s_inds[1]] = dists[s_inds[end]] = Inf
        length(s_inds) == 2 && continue # would lead to NaN
        for j in 2:n_pnts-1
            dists[s_inds[j]] += (individs[s_inds[j+1]][i] - individs[s_inds[j-1]][i]) / (max_val - min_val)
        end
    end
    return dists
end

""" perfrom binary tournament selection for parent selection
"""
parent_selection(pop) = min(rand(pop), rand(pop))

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
    indiv_obj_vals ./= (median(indiv_obj_vals, dims=1) .+ 1e-1)
    return eachrow(indiv_obj_vals)
end

""" Calculate the relative fitness based on normalized objectives using the
    median instead of the maximum.
"""
get_relative_fitness(indiv_obj_vals) = -sum.(indiv_obj_vals)

