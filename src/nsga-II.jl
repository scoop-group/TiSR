
function generational_loop(data, ops)
    population = [Individual[] for _ in 1:ops.general.num_islands]
    new_nodes  = [Node[] for _ in 1:ops.general.num_islands]
    return generational_loop(data, ops, population, new_nodes)
end

function generational_loop(data, ops, start_pop::Vector{Individual})
    population = [
        start_pop[isle:ops.general.num_islands:end]
        for isle in 1:ops.general.num_islands
    ]
    new_nodes  = [Node[] for _ in 1:ops.general.num_islands]
    return generational_loop(data, ops, population, new_nodes)
end

function generational_loop(data, ops, start_pop::Vector{String})
    population = [Individual[] for _ in 1:ops.general.num_islands]
    start_pop  = Node[string_to_node(eq, ops) for eq in start_pop]
    new_nodes  = [
        start_pop[isle:ops.general.num_islands:end]
        for isle in 1:ops.general.num_islands
    ]
    return generational_loop(data, ops, population, new_nodes)
end

function generational_loop(
    data,
    ops,
    population::Vector{Vector{Individual}},
    new_nodes::Vector{Vector{Node}}
)

    children     = [Individual[] for _ in 1:ops.general.num_islands]
    hall_of_fame = Individual[]

# ==================================================================================================
# prepar user interrupt
# ==================================================================================================
    # Create a channel for communication between tasks
    input_channel = Channel{String}(1)

    # Start listening for user input asynchronously
    @async listen_for_input(input_channel)

# ==================================================================================================
# some preparation
# ==================================================================================================
    @assert length(data) == ops.data_descript.n_vars + 1 "please use the data variable which was returned by the Options constructor."

    # create dict for the simulation progression # -------------------------------------------------
    keep_track_of = [
        "time", "generation", "mean age hall_of_fame",
        "cur_max_compl", "mean compl", "mean recursive_compl", "mean n_params",
        "min mae", "min mse", "min max_ae", "min minus_r2", "min minus_abs_spearman", "min mare", "min q75_are",
        "min max_are", "min ms_processed_e"
    ]

    prog_dict     = OrderedDict(key => eltype(data[1])[]             for key in keep_track_of)
    cur_prog_dict = OrderedDict(key => convert(eltype(data[1]), 0.0) for key in keep_track_of)

# ==================================================================================================
# start generational loop
# ==================================================================================================
    t_start  = time()
    gen      = 0.0
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
            if ops.general.multithreading
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
            length(population[isle]) > ops.general.pop_per_isle || continue

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

            if ops.general.plot_hall_of_fame
                compl          = [indiv.compl          for indiv in hall_of_fame]
                ms_processed_e = [indiv.ms_processed_e for indiv in hall_of_fame]

                display(scatterplot(
                    compl,
                    ms_processed_e,
                    yscale           = :log10,
                    title            = "hall of fame",
                    xlabel           = "complexity",
                    ylabel           = "log10 of ms_processed_e",
                    marker           = :circle,
                    unicode_exponent = false,
                    xlim             = (0, ops.grammar.max_compl),
                    ylim             = (
                        10^(floor(log10(minimum(ms_processed_e)))),
                        10^(ceil(log10(maximum(ms_processed_e))))
                    ),
                ))
            end
        end


# ==================================================================================================
# population shuffle
# ==================================================================================================
        if gen % ops.general.population_shuffle_interval == 0
            split_list(list, n) = length(list) <= n ? [list] : [list[1:n], split_list(list[n+1:end], n)...]
            population = reduce(vcat, population)
            shuffle!(population)
            population = split_list(population, ceil(Int, length(population) / ops.general.num_islands))
        end

# ==================================================================================================
# termination criteria
# ==================================================================================================
        if gen >= ops.general.n_gens
            stop_msg = "reached maximum number of generations"
            break
        elseif time() - t_start >= ops.general.t_lim
            stop_msg = "reached time limit"
            break
        elseif ops.general.callback(hall_of_fame, population, ops)
            stop_msg = "callback returned true"
            break
        elseif isready(input_channel)
            user_input = take!(input_channel)  # Get input from the channel
            if user_input == "qq"
                println("Detected 'qq'. Exiting...")
                stop_msg = "user interrupt"
                break
            elseif length(user_input) > 0
                println("type qq and enter to finish the equation search early")
            end
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

    return hall_of_fame, population, prog_dict, stop_msg
end

"""
    listen_for_input(input_channel::Channel)

Listens for user input in a blocking fashion and sends the input to the given channel. 
This function runs asynchronously and continuously reads from `stdin`. If the input 
"qq" is detected, it will send the input to the channel and terminate the listener.

### Arguments
- `input_channel::Channel`: A communication channel where user input is sent.

### Behavior
- The function reads input from the terminal using `readline(stdin)`.
- The input is passed to the `input_channel` for further processing in a separate task.
- If the user types "qq", the function stops, indicating the program should terminate.

### Source
Generated using OpenAI's ChatGPT. 
"""
function listen_for_input(input_channel::Channel)
    while true
        input = readline(stdin)  # Blocking read, but in its own task
        put!(input_channel, input)  # Send the input to the channel
        if input == "qq"
            break  # Exit the listener task if 'qq' is entered
        end
    end
end
