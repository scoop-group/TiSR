
function generational_loop(data, ops ;start_pop=Node[])

    @assert length(data) == ops.data_descript.n_vars "please use the data variable which was returned by the Options constructor."

# ==================================================================================================
# some preparation
# ==================================================================================================
    data_type = eltype(data[1])

    # create dict for the simulation progression # -------------------------------------------------
    keep_track_of = [
        "time", "generation", "mean age hall_of_fame",
        "cur_max_compl", "mean compl", "mean recursive_compl", "mean n_params",
        "min mae", "min mse", "min max_ae", "min minus_r2", "min minus_abs_spearman", "min mare", "min q75_are",
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
    stop_msg = ""

    while true
        gen += 1.0

        cur_max_compl = maximum(indiv.compl for indiv in hall_of_fame; init=5)

        foreach(indiv -> indiv.age += 1, hall_of_fame)

        for isle in 1:ops.general.num_islands

            foreach(indiv -> indiv.age += 1, population[isle])

# ==================================================================================================
# genetic operations
# ==================================================================================================
            # create new children # ----------------------------------------------------------------
            while length(new_nodes[isle]) + length(population[isle]) < 0.6 * ops.general.pop_per_isle
                push!(new_nodes[isle], grow_equation(ops.grammar.init_tree_depth, ops, method=:asym))
            end

            # perform mutations # ------------------------------------------------------------------
            if length(population[isle]) > 0.4 * ops.general.pop_per_isle
                shuffle!(population[isle])

                for i in 1:max(length(population[isle]), (2 * ops.general.pop_per_isle  - length(population[isle])))
                    push!(new_nodes[isle], deepcopy(getfield(population[isle][mod1(i, length(population[isle]))], :node)))
                end

                apply_genetic_operations!(new_nodes[isle], ops)
            end

            # grow them up # -----------------------------------------------------------------------
            children[isle] = Array{Individual}(undef, length(new_nodes[isle]))
            if ops.general.multihreadding
                Threads.@threads for ii in eachindex(new_nodes[isle])
                    children[isle][ii] = Individual(new_nodes[isle][ii], data, ops, cur_max_compl)
                end
            else
                for ii in eachindex(new_nodes[isle])
                    children[isle][ii] = Individual(new_nodes[isle][ii], data, ops, cur_max_compl)
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
# remove apperently same by rounded MAE and MSE
# ==================================================================================================
        if ops.general.remove_doubles_sigdigits > 0
            if ops.general.remove_doubles_across_islands
                remove_doubles_across_islands!(population, ops)
            else
                for isle in 1:ops.general.num_islands
                    remove_doubles!(population[isle], ops)
                end
            end
        end

# ==================================================================================================
# selection
# ==================================================================================================
        for isle in 1:ops.general.num_islands
            isempty(population[isle]) && continue

            selection_inds = Int64[]

            indiv_obj_vals = [
                Float64[
                    round(getfield(indiv, obj), sigdigits=ops.selection.population_niching_sigdigits)
                    for obj in ops.selection.selection_objectives
                ]
                for indiv in population[isle]
            ]

            # Pareto selection # -------------------------------------------------------------------
            if ops.selection.n_pareto_select_per_isle > 0
                selected = non_dominated_sort(
                    indiv_obj_vals,
                    n_select=ops.selection.n_pareto_select_per_isle,
                    first_front=false
                )

                append!(selection_inds, selected)
            end

            # tournament selection # ---------------------------------------------------------------
            if ops.general.pop_per_isle - ops.selection.n_pareto_select_per_isle > 0
                selected = tournament_selection(
                (
                    reduce(+, scale .* getfield.(population[isle], field)
                        for (scale, field) in ops.selection.tournament_selection_fitness)
                ),
                    setdiff(eachindex(population[isle]), selection_inds),
                    tournament_size=ops.selection.tournament_size,
                    n_select = ops.general.pop_per_isle - ops.selection.n_pareto_select_per_isle
                )

                append!(selection_inds, selected)
            end

            sort!(selection_inds)
            keepat!(population[isle], selection_inds)
        end

# ==================================================================================================
# migration
# ==================================================================================================
        if gen % ops.general.migration_interval == 0
            emmigrate_island = rand(1:ops.general.num_islands)
            immigrate_island = mod1(emmigrate_island + rand((1, -1)), ops.general.num_islands)

            !isempty(population[emmigrate_island]) || continue

            push!(
                population[immigrate_island],
                popat!(
                    population[emmigrate_island],
                    rand(1:length(population[emmigrate_island]))
                )
            )
        end

# ==================================================================================================
# hall of fame
# ==================================================================================================
        for isle in 1:ops.general.num_islands
            append!(hall_of_fame, deepcopy.(population[isle]))
        end

        indiv_obj_vals = [
            Float64[
                round(getfield(indiv, obj), sigdigits=ops.selection.hall_of_fame_niching_sigdigits)
                for obj in ops.selection.hall_of_fame_objectives
            ]
            for indiv in hall_of_fame
        ]

        selection_inds = non_dominated_sort(indiv_obj_vals; first_front=true)

        keepat!(hall_of_fame, selection_inds)

# ==================================================================================================
# every couple of generations
# ==================================================================================================
        t_since = time() - t_start

        if length(prog_dict["time"]) == 0 || t_since - prog_dict["time"][end] > 5.0

            # current KPIs # -----------------------------------------------------------------------
            get_for_prog = [t_since, gen,
                mean(getfield.(hall_of_fame, :age)),
                cur_max_compl,
                mean(getfield.(hall_of_fame, :compl)),
                mean(getfield.(hall_of_fame, :recursive_compl)),
                mean(getfield.(hall_of_fame, :n_params)),
                minimum(getfield.(hall_of_fame, :mae)),
                minimum(getfield.(hall_of_fame, :mse)),
                minimum(getfield.(hall_of_fame, :max_ae)),
                minimum(getfield.(hall_of_fame, :minus_r2)),
                minimum(getfield.(hall_of_fame, :minus_abs_spearman)),
                minimum(getfield.(hall_of_fame, :mare)),
                minimum(getfield.(hall_of_fame, :q75_are)),
                minimum(getfield.(hall_of_fame, :max_are)),
                minimum(getfield.(hall_of_fame, :ms_processed_e))]

            for i in eachindex(keep_track_of)
                push!(prog_dict[keep_track_of[i]], get_for_prog[i])
                cur_prog_dict[keep_track_of[i]] = get_for_prog[i]
            end

            if ops.general.print_progress
                display(cur_prog_dict)
                println("\n", round(Int64, t_since ÷ 60), " min  ", round(Int64, t_since % 60), " sec")
            end
        end

        if gen >= ops.general.n_gens
            stop_msg = "reached maximum number of generations"
            break
        elseif time() - t_start >= ops.general.t_lim
            stop_msg = "reached time limit"
            break
        elseif ops.general.callback(hall_of_fame, population, ops)
            stop_msg = "callback returned true"
            break
        end
    end

    # final display of current KPIs # --------------------------------------------------------------
    t_since = time() - t_start
    cur_max_compl = maximum(indiv.compl for indiv in hall_of_fame; init=5)

    get_for_prog = [t_since, gen,
        mean(getfield.(hall_of_fame, :age)),
        cur_max_compl,
        mean(getfield.(hall_of_fame, :compl)),
        mean(getfield.(hall_of_fame, :recursive_compl)),
        mean(getfield.(hall_of_fame, :n_params)),
        minimum(getfield.(hall_of_fame, :mae)),
        minimum(getfield.(hall_of_fame, :mse)),
        minimum(getfield.(hall_of_fame, :max_ae)),
        minimum(getfield.(hall_of_fame, :minus_r2)),
        minimum(getfield.(hall_of_fame, :minus_abs_spearman)),
        minimum(getfield.(hall_of_fame, :mare)),
        minimum(getfield.(hall_of_fame, :q75_are)),
        minimum(getfield.(hall_of_fame, :max_are)),
        minimum(getfield.(hall_of_fame, :ms_processed_e))]

    for i in eachindex(keep_track_of)
        push!(prog_dict[keep_track_of[i]], get_for_prog[i])
        cur_prog_dict[keep_track_of[i]] = get_for_prog[i]
    end

    if ops.general.print_progress
        display(cur_prog_dict)
        println("\n", round(Int64, t_since ÷ 60), " min  ", round(Int64, t_since % 60), " sec")
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

    return hall_of_fame, population, prog_dict, stop_msg
end


