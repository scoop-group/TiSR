
function generational_loop(data, ops)
    population = [Individual[] for _ in 1:ops.general.num_islands]
    children   = [Individual[] for _ in 1:ops.general.num_islands]
    return generational_loop(data, ops, population, children)
end

function generational_loop(data, ops, start_pop::Vector{Individual})
    population = [
        start_pop[isle:ops.general.num_islands:end]
        for isle in 1:ops.general.num_islands
    ]
    children   = [Individual[] for _ in 1:ops.general.num_islands]
    return generational_loop(data, ops, population, children)
end

function generational_loop(data, ops, start_pop::Vector{String})
    population = [Individual[] for _ in 1:ops.general.num_islands]
    start_pop  = Node[string_to_node(eq, ops) for eq in start_pop]
    children   = [
        [
            Individual(node, ops)
            for node in start_pop[isle:ops.general.num_islands:end]
        ]
        for isle in 1:ops.general.num_islands
    ]
    return generational_loop(data, ops, population, children)
end

function generational_loop(
    data,
    ops,
    population::Vector{Vector{Individual}},
    children::Vector{Vector{Individual}},
)
    timer = TimerOutput()
# ==================================================================================================
# prepare bank_of_terms
# ==================================================================================================
    bank_of_terms  = Node[string_to_node(eq, ops) for eq in ops.grammar.bank_of_terms]

# ==================================================================================================
# check if bank_of_terms and start_pop legal
# ==================================================================================================
    if !isempty(ops.grammar.illegal_dict)
        for node in bank_of_terms
            @assert is_legal_nesting(node, ops) "'$node' from bank_of_terms has some illegal nestings"
        end
        for new_node in population
            for node in new_node
                @assert is_legal_nesting(node, ops) "'$node' from start_pop has some illegal nestings"
            end
        end
    end

# ==================================================================================================
# initialize data structures
# ==================================================================================================
    hall_of_fame = Individual[]

# ==================================================================================================
# prepare user interrupt
# ==================================================================================================
    stdin_reader = watch_stream(stdin)

# ==================================================================================================
# misc
# ==================================================================================================
    null_node = Node(0.0)

    @assert length(data) == ops.data_descript.n_vars + 1 "please use the data variable which was returned by the Options constructor."

    prog_dict = OrderedDict(
        k => Float64[]
        for k in ["time", "generation", "mean age", "cur_max_compl", "mean compl",
                  ["min " * string(m) for m in keys(ops.measures)]...]
    )

# ==================================================================================================
# start generational loop
# ==================================================================================================
    t_start  = time()
    gen      = 0.0
    stop_msg = ""

    while true
        gen += 1.0

        cur_max_compl = maximum(indiv.measures[:compl] for indiv in hall_of_fame; init=ops.grammar.min_compl)

        foreach(indiv -> indiv.age += 1, hall_of_fame)

        for isle in 1:ops.general.num_islands

            foreach(indiv -> indiv.age += 1, population[isle])

# ==================================================================================================
# genetic operations
# ==================================================================================================
            # create new children # ----------------------------------------------------------------
            while length(children[isle]) + length(population[isle]) < 0.6 * ops.general.pop_per_isle
                push!(children[isle],
                    Individual(grow_equation(ops.grammar.init_tree_depth, ops), ops)
                )
            end

            @timeit timer "mutations" begin
                if length(population[isle]) > 0.4 * ops.general.pop_per_isle

                    # select parents
                    shuffle!(population[isle])
                    for i in 1:ops.general.n_children
                        if ops.general.parent_selection
                            push!(children[isle], copy(parent_selection(population[isle])))
                        else
                            push!(children[isle], copy(population[isle][mod1(i, length(population[isle]))]))
                        end
                    end

                    apply_genetic_operations!(children[isle], ops, bank_of_terms)

                    for _ in 1:ops.general.n_refitting
                        push!(children[isle], copy(rand(population[isle])))
                    end
                end
            end

            if ops.general.fitting_island_function(isle)
                fit_iter = ops.fitting.max_iter
            else
                fit_iter = 0
            end

            # grow them up # -----------------------------------------------------------------------
            if ops.general.multithreading
                Threads.@threads :greedy for ii in eachindex(children[isle])
                    fit_individual!(children[isle][ii], data, ops, cur_max_compl, fit_iter, timer)
                end
            else
                @timeit timer "Individual" begin
                    for ii in eachindex(children[isle])
                        fit_individual!(children[isle][ii], data, ops, cur_max_compl, fit_iter, timer)
                    end
                end
            end

            filter!(indiv -> indiv.valid, children[isle])

# ==================================================================================================
# add children to population
# ==================================================================================================
            append!(population[isle], children[isle])
            children[isle] = Individual[]

        end # for isle in 1:ops.general.num_islands

# ==================================================================================================
# remove individuals
# ==================================================================================================
        @timeit timer "remove doubles and old individuals" begin
            if ops.general.remove_doubles_sigdigits > 0
                for isle in 1:ops.general.num_islands
                    remove_doubles!(population[isle], ops)
                end
            end

            for isle in 1:ops.general.num_islands
                filter!(i -> i.age <= ops.general.max_age, population[isle])
            end
        end

# ==================================================================================================
# selection
# ==================================================================================================
        @timeit timer "selection" begin
            for isle in 1:ops.general.num_islands
                length(population[isle]) > ops.general.pop_per_isle || continue

                selection_inds = Int64[]

                indiv_obj_vals = [
                    Float64[
                        round(indiv.measures[obj], sigdigits=ops.selection.population_niching_sigdigits)
                        for obj in ops.selection.selection_objectives
                    ]
                    for indiv in population[isle]
                ]

                # determine rank and crowding for all individuals -> overkill, if no parent selection
                ranks    = non_dominated_sort_clean(indiv_obj_vals)
                crowding = crowding_distance_clean(indiv_obj_vals)

                for i in eachindex(population[isle])
                    population[isle][i].rank     = ranks[i]
                    population[isle][i].crowding = crowding[i]
                end

                # Pareto selection # -------------------------------------------------------------------
                if ops.selection.n_pareto_select_per_isle > 0
                    sort!(population[isle]) # apply non_dominated_sort
                    append!(selection_inds, 1:ops.selection.n_pareto_select_per_isle)
                end

                # tournament selection # ---------------------------------------------------------------
                if ops.general.pop_per_isle - ops.selection.n_pareto_select_per_isle > 0
                    remaining_inds = setdiff(eachindex(population[isle]), selection_inds)

                    if ops.general.pop_per_isle - length(selection_inds) < length(remaining_inds)
                        fitness  = get_relative_fitness(indiv_obj_vals)
                        selected = tournament_selection(fitness, remaining_inds,
                            tournament_size = ops.selection.tournament_size,
                            n_select        = ops.general.pop_per_isle - length(selection_inds)
                        )
                    else
                        selected = remaining_inds
                    end

                    append!(selection_inds, selected)
                end

                # apply selection
                sort!(selection_inds)
                keepat!(population[isle], selection_inds)
            end

            # remove doubles across islands
            if ops.general.remove_doubles_across_islands && ops.general.remove_doubles_sigdigits > 0
                remove_doubles_across_islands!(population, ops)
            end
        end # TODO: put the for isle together

# ==================================================================================================
# migration
# ==================================================================================================
        @timeit timer "migration" begin
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

            # hall of fame migration # -----------------------------------------------------------------
            if gen % ops.general.hall_of_fame_migration_interval == 0
                indiv = copy(rand(hall_of_fame))
                indiv.age = 0
                push!(population[rand(1:ops.general.num_islands)], indiv)
            end
        end

# ==================================================================================================
# hall of fame
# ==================================================================================================
        @timeit timer "update hall of fame" begin
            for isle in 1:ops.general.num_islands
                append!(hall_of_fame, copy.(population[isle]))
            end

            indiv_obj_vals = [
                Float64[
                    round(indiv.measures[obj], sigdigits=ops.selection.hall_of_fame_niching_sigdigits)
                    for obj in ops.selection.hall_of_fame_objectives
                ]
                for indiv in hall_of_fame
            ]

            selection_inds = non_dominated_sort(indiv_obj_vals; first_front=true) # TODO: cleanup

            keepat!(hall_of_fame, selection_inds)
        end

# ==================================================================================================
# every couple of generations
# ==================================================================================================
        t_since = time() - t_start

        @timeit timer "plotting, logging, and garbage prevention" begin
            if length(prog_dict["time"]) == 0 || t_since - prog_dict["time"][end] > 5.0

                # overwrite references to trash nodes with the null node
                for isle in population
                    for indiv in isle
                        clean_trash_nodes!(indiv.node, null_node)
                    end
                end
                for indiv in hall_of_fame
                    clean_trash_nodes!(indiv.node, null_node)
                end

                # current KPIs # -----------------------------------------------------------------------
                get_for_prog = OrderedDict([
                    "time"          => t_since,
                    "generation"    => gen,
                    "mean age"      => mean(i.age for i in hall_of_fame),
                    "cur_max_compl" => cur_max_compl,
                    "mean compl"    => mean(i.measures[:compl] for i in hall_of_fame),
                    ["min " * string(m) => minimum([i.measures[m] for i in hall_of_fame])
                     for m in keys(ops.measures)]...
               ])

                if ops.general.print_progress
                    display(get_for_prog)
                    println("\n", round(Int64, t_since ÷ 60), " min  ", round(Int64, t_since % 60), " sec | type q and enter to finish early")
                end

                for k in keys(get_for_prog)
                    push!(prog_dict[k], get_for_prog[k])
                end

                if ops.general.plot_hall_of_fame
                    compl          = [indiv.measures[:compl]          for indiv in hall_of_fame]
                    ms_processed_e = [indiv.measures[:ms_processed_e] for indiv in hall_of_fame]

                    plt = scatterplot(
                        compl,
                        clamp.(ms_processed_e, 1e-30, 1e30),
                        yscale           = :log10,
                        title            = "hall of fame",
                        xlabel           = "compl",
                        ylabel           = "log10 of ms_processed_e",
                        marker           = :circle,
                        unicode_exponent = false,
                        xlim             = (0, ops.grammar.max_compl),
                        ylim             = (
                            10^(floor(log10(minimum(ms_processed_e)))),
                            10^(ceil(log10(maximum(ms_processed_e))))
                        ),
                        compact=true
                    )

                    if ops.general.print_hall_of_fame
                        sort!(hall_of_fame, by=i->i.measures[:compl])

                        inds_to_show = round.(Int64, collect(range(1, length(hall_of_fame), length=16)))
                        unique!(inds_to_show)

                        eq_strs = [simplify_to_string(hall_of_fame[i].node, ops, sigdigits=2) for i in inds_to_show]

                        for (ii, i) in enumerate(inds_to_show)
                            label!(plt, :r, ii, replace(
                                eq_strs[ii],
                                " " => "", r"(\d)\.0\b" => s"\1"
                            ))
                        end
                    end

                    display(plt)
                end
            end
        end

# ==================================================================================================
# island extinction
# ==================================================================================================
        @timeit timer "island extinction" begin
            if gen % ops.general.island_extinction_interval == 0
                emmigrate_island = rand(1:ops.general.num_islands)

                # both directions, decreasing with distance
                offsets = -trunc(Int64, 0.5 * ops.general.num_islands):trunc(Int64, 0.5 * ops.general.num_islands)
                offsets = filter(!=(0), offsets)
                probs   = (1 ./ abs.(offsets)).^2

                while !isempty(population[emmigrate_island])
                    indiv = popat!(
                        population[emmigrate_island],
                        rand(1:length(population[emmigrate_island]))
                    )

                    if rand() < ops.general.migrate_after_extinction_prob
                        immigrate_island = mod1(
                            emmigrate_island + wsample(offsets, probs),
                            ops.general.num_islands
                        )
                        push!(population[immigrate_island], indiv)
                    end
                end
            end
        end

# ==================================================================================================
# termination criteria
# ==================================================================================================
    @timeit timer "termination criteria (inlc. callback)" begin
            if gen >= ops.general.n_gens
                stop_msg = "reached maximum number of generations"
                break
            elseif time() - t_start >= ops.general.t_lim
                stop_msg = "reached time limit"
                break
            elseif ops.general.callback(hall_of_fame, population, ops)
                stop_msg = "callback returned true"
                break
            elseif check_for_user_quit(stdin_reader)
                stop_msg = "graceful user interrupt"
                break
            end
        end
    end

    # final display of current KPIs # --------------------------------------------------------------
    t_since = time() - t_start
    cur_max_compl = maximum(indiv.measures[:compl] for indiv in hall_of_fame; init=ops.grammar.min_compl)

    get_for_prog = OrderedDict([
        "time"          => t_since,
        "generation"    => gen,
        "mean age"      => mean(i.age for i in hall_of_fame),
        "cur_max_compl" => cur_max_compl,
        "mean compl"    => mean(i.measures[:compl] for i in hall_of_fame),
        ["min " * string(m) => minimum([i.measures[m] for i in hall_of_fame])
         for m in keys(ops.measures)]...
    ])

    if ops.general.print_progress
        display(get_for_prog)
        println("\n", round(Int64, t_since ÷ 60), " min  ", round(Int64, t_since % 60), " sec")
    end

    for k in keys(get_for_prog)
        push!(prog_dict[k], get_for_prog[k])
    end

    close_reader!(stdin_reader)
# ==================================================================================================
# Post-pare for return
# ==================================================================================================
    population = reduce(vcat, population)

    return hall_of_fame, population, prog_dict, stop_msg, timer
end

""" Copied from SymbolicRegression.jl, after own failed attempts.
"""
struct StdinReader{ST}
    can_read_user_input::Bool
    stream::ST
end

""" Start watching stream (like stdin) for user input.
    Copied from SymbolicRegression.jl, after own failed attempts.
"""
function watch_stream(stream)
    can_read_user_input = isreadable(stream)

    can_read_user_input && try
        Base.start_reading(stream)
        bytes = bytesavailable(stream)
        if bytes > 0
            # Clear out initial data
            read(stream, bytes)
        end
    catch err
        if isa(err, MethodError)
            can_read_user_input = false
        else
            throw(err)
        end
    end
    return StdinReader(can_read_user_input, stream)
end

"""Close the stdin reader and stop reading.
   Copied from SymbolicRegression.jl, after own failed attempts.
"""
function close_reader!(reader::StdinReader)
    if reader.can_read_user_input
        Base.stop_reading(reader.stream)
    end
end

"""Check if the user typed 'q' and <enter> or <ctl-c>.
   Copied from SymbolicRegression.jl, after own failed attempts.
"""
function check_for_user_quit(reader::StdinReader)::Bool
    if reader.can_read_user_input
        bytes = bytesavailable(reader.stream)
        if bytes > 0
            # Read:
            data = read(reader.stream, bytes)
            control_c = 0x03
            quit = 0x71
            if length(data) > 1 && (data[end] == control_c || data[end - 1] == quit)
                return true
            end
        end
    end
    return false
end
