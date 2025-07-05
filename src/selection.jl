
""" Whether x dominates y.
"""
dominates(x, y) = all(i -> x[i] <= y[i], eachindex(x)) && any(i -> x[i] < y[i], eachindex(x))

""" Slightly adapted of a version found on Julia Discourse:
    [source](https://discourse.julialang.org/t/fast-non-dominated-sorting/24164)
    For each entry, the pareto front number is determined.
"""
function non_dominated_sort(individs::Vector{Vector{Float64}})
    #@timeit to "non_dominated_sort" begin
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
    #end # @timeit

    return domination_score
end

""" Find the first pareto front.
"""
first_pareto_front(individs::Vector{Vector{Float64}}) = findall(y -> all(x -> !dominates(x, y), individs), individs)

""" Calculate the crowding distances.
"""
function crowding_distance(individs)
    #@timeit to "crowding_distance" begin
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
                dists[s_inds[j]] += (individs[s_inds[j+1]][i] - individs[s_inds[j-1]][i]) / (max_val - min_val) # TODO: maybe use * instead of +
            end
        end
    #end # @timeit
    return dists
end

""" Perfrom binary tournament selection for parent selection.
"""
parent_selection(pop; tournament_size=2) = minimum(rand(pop) for _ in 1:tournament_size)

""" Perform the population selection.
"""
function perform_population_selection!(pop, ops, isle)
    sort!(pop, by=i->i.age)

    if ops.general.constraint_island_function(isle)
        selection_objectives = ops.selection.selection_objectives
    else
        selection_objectives = filter(!=(:constr_vio), ops.selection.selection_objectives)
    end

    # extract objectives # -------------------------------------------------------------
    indiv_obj_vals = [
        Float64[
            indiv.measures[obj] for obj in ops.selection.selection_objectives
        ]
        for indiv in pop
    ]

    # apply niching
    if ops.selection.population_niching_sigdigits > 0
        indiv_obj_vals = [
            round.(indiv, sigdigits=ops.selection.population_niching_sigdigits)
            for indiv in indiv_obj_vals
        ]
        unique_inds = unique(i -> indiv_obj_vals[i], 1:length(indiv_obj_vals))
        keepat!(indiv_obj_vals, unique_inds)
        keepat!(pop, unique_inds)
    end

    # determine rank and crowding for all individuals
    ranks    = non_dominated_sort(indiv_obj_vals)
    crowding = crowding_distance(indiv_obj_vals)

    for i in eachindex(pop)
        pop[i].rank     = ranks[i]
        pop[i].crowding = crowding[i]
    end

    length(pop) > ops.general.pop_per_isle || return

    # pareto selection # ---------------------------------------------------------------------------
    sort!(pop)
    keepat!(pop, 1:ops.general.pop_per_isle)
end

""" Perfroms the hall_of_fame selection.
"""
function perform_hall_of_fame_selection!(hall_of_fame, population, ops)
    for isle in 1:ops.general.num_islands
        prepend!(hall_of_fame, deepcopy.(population[isle])) # TODO: maybe only add, and deepcopy only hall_of_fame of itself in the end?
    end

    indiv_obj_vals = [
        Float64[
            indiv.measures[obj] for obj in ops.selection.hall_of_fame_objectives
        ]
        for indiv in hall_of_fame
    ]

    # apply niching
    if ops.selection.hall_of_fame_niching_sigdigits > 0
        indiv_obj_vals = [
            round.(indiv, sigdigits=ops.selection.hall_of_fame_niching_sigdigits)
            for indiv in indiv_obj_vals
        ]
        unique_inds = unique(i -> indiv_obj_vals[i], 1:length(indiv_obj_vals))
        keepat!(indiv_obj_vals, unique_inds)
        keepat!(hall_of_fame, unique_inds)
    end

    selection_inds = first_pareto_front(indiv_obj_vals)

    keepat!(hall_of_fame, selection_inds)
end

