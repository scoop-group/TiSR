
function generational_loop(data, ops ;start_pop=Node[])
# ==================================================================================================
# some preparation
# ==================================================================================================
    data_type = eltype(data[1])

    # create dict for the simulation progression # -------------------------------------------------
    keep_track_of = [
        "time", "generation", "mean age hall_of_fame",
        "mean compl", "mean recursive_compl", "mean n_params",
        "min mae", "min mse", "min max_ae", "min minus_r2", "min mare", "min q75_are",
        "min max_are", "min ms_processed_e"
    ]

    prog_dict = OrderedDict(key => data_type[] for key in keep_track_of)
    cur_prog_dict = OrderedDict(key => convert(data_type, 0.0) for key in keep_track_of)

    # initialize data structures # -----------------------------------------------------------------
    population =   [Individual[] for _ in 1:ops.general.num_islands]
    children =     [Individual[] for _ in 1:ops.general.num_islands]
    hall_of_fame = Individual[]

    # generate initial random children # -----------------------------------------------------------
    new_nodes = [
        start_pop[isle:ops.general.num_islands:end]
        for isle in 1:ops.general.num_islands
    ]

# ==================================================================================================
# start generational loop
# ==================================================================================================
    t_start = time()
    gen = 0.0

    while (gen <= ops.general.n_gens && time() - t_start < ops.general.t_lim)
        gen += 1.0

        foreach(indiv -> indiv.age += 1, hall_of_fame)

        for isle in 1:ops.general.num_islands

            foreach(indiv -> indiv.age += 1, population[isle])

# ==================================================================================================
# genetic operations
# ==================================================================================================
            # create new children # ----------------------------------------------------------------
            if gen > 1.0
                shuffle!(population[isle])

                for i in 1:max(length(population[isle]), (2 * ops.general.pop_per_isle  - length(population[isle])))
                    push!(new_nodes[isle], copy_node(getfield(population[isle][mod1(i, length(population[isle]))], :node)))
                end

                apply_genetic_operations!(new_nodes[isle], ops)
            else
                while length(new_nodes[isle]) + length(population[isle]) < 2 * ops.general.pop_per_isle
                    push!(new_nodes[isle], grow_equation(ops.general.init_tree_depth, ops))
                end
            end

            # grow them up # -----------------------------------------------------------------------
            children[isle] = Array{Individual}(undef, length(new_nodes[isle]))
            if ops.general.multihreadding
                Threads.@threads for ii in eachindex(new_nodes[isle])
                    children[isle][ii] = Individual(new_nodes[isle][ii], data, ops)
                end
            else
                for ii in eachindex(new_nodes[isle])
                    children[isle][ii] = Individual(new_nodes[isle][ii], data, ops)
                end
            end

            filter!(indiv -> indiv.valid, children[isle])
            new_nodes[isle] = Individual[]

# ==================================================================================================
# add children to population
# ==================================================================================================
            append!(population[isle], children[isle])

        end # for isle in 1:ops.general.num_islands

# ==================================================================================================
# remove apperently same by MAE and MSE
# ==================================================================================================
        if ops.general.prevent_doubles > 0
            if ops.general.prevent_doubles_across_islands
                remove_doubles_across_islands!(population, ops)
            else
                for isle in 1:ops.general.num_islands
                    remove_doubles_by_structure!(population[isle])
                    remove_doubles!(population[isle], ops)
                end
            end
        end

# ==================================================================================================
# selection
# ==================================================================================================
        for isle in 1:ops.general.num_islands

            selection_inds = Int64[]

            if ops.selection.n_pareto_select_per_isle > 0
                append!(selection_inds, non_dominated_sort(
                [
                        getfield.(population[isle], obj)
                        for obj in ops.selection.selection_objectives
                    ],
                    n_select=ops.selection.n_pareto_select_per_isle,
                    first_front=false)
                )
            end

            if ops.general.pop_per_isle - ops.selection.n_pareto_select_per_isle > 0
                append!(selection_inds, tournament_selection(
                (
                    reduce(+, scale .* getfield.(population[isle], field)
                        for (scale, field) in ops.selection.tournament_selection_fitness)
                ),
                    setdiff(eachindex(population[isle]), selection_inds),
                    tournament_size=ops.selection.tournament_size,
                    n_select = ops.general.pop_per_isle - ops.selection.n_pareto_select_per_isle
                )
                )
            end

            sort!(selection_inds)
            keepat!(population[isle], selection_inds)
        end

# ==================================================================================================
# migration
# ==================================================================================================
        if gen % ops.general.migration_interval == 0
            for _ in 1:ops.general.n_migrations
                emmigrate_island = rand(1:ops.general.num_islands)
                immigrate_island = mod1(emmigrate_island + rand((1, -1)), ops.general.num_islands)

                push!(
                    population[immigrate_island],
                    deepcopy(
                        rand(population[emmigrate_island])
                    )
                )
            end
        end

# ==================================================================================================
# hall of fame
# ==================================================================================================
        for isle in 1:ops.general.num_islands
            append!(hall_of_fame, copy.(population[isle])) # deepcopy quite expansive
        end

        remove_doubles!(hall_of_fame, ops)
        remove_doubles_by_structure!(hall_of_fame)

        selection_inds = non_dominated_sort(
            [getfield.(hall_of_fame, obj) for obj in ops.selection.hall_of_fame_objectives]
            ; first_front=true)

        keepat!(hall_of_fame, selection_inds)

# ==================================================================================================
# every couple of generations
# ==================================================================================================
        t_since = time() - t_start

        if length(prog_dict["time"]) == 0 || t_since - prog_dict["time"][end] > 5.0

            # current KPIs # -----------------------------------------------------------------------
            get_for_prog = [t_since, gen,
                mean(getfield.(hall_of_fame, :age)),
                mean(getfield.(hall_of_fame, :compl)),
                mean(getfield.(hall_of_fame, :recursive_compl)),
                mean(getfield.(hall_of_fame, :n_params)),
                minimum(getfield.(hall_of_fame, :mae)),
                minimum(getfield.(hall_of_fame, :mse)),
                minimum(getfield.(hall_of_fame, :max_ae)),
                minimum(getfield.(hall_of_fame, :minus_r2)),
                minimum(getfield.(hall_of_fame, :mare)),
                minimum(getfield.(hall_of_fame, :q75_are)),
                minimum(getfield.(hall_of_fame, :max_are)),
                minimum(getfield.(hall_of_fame, :ms_processed_e))]

            for i in eachindex(keep_track_of)
                push!(prog_dict[keep_track_of[i]], get_for_prog[i])
                cur_prog_dict[keep_track_of[i]] = get_for_prog[i]
            end

            display(cur_prog_dict)
            println("\n", round(Int64, t_since ÷ 60), " min  ", round(Int64, t_since % 60), " sec")
        end
    end

# ==================================================================================================
# Post-pare for return
# ==================================================================================================
    population = reduce(vcat, population)

    population = Dict(string(f) => getfield.(population, f)
        for f in fieldnames(typeof(population[1]))
    )

    hall_of_fame = Dict(string(f) => getfield.(hall_of_fame, f)
        for f in fieldnames(typeof(hall_of_fame[1]))
    )

    return hall_of_fame, population, prog_dict
end


