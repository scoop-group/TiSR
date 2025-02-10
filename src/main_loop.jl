
function generational_loop(data::Vector{Vector{Float64}}, ops)
    population = [Individual[] for _ in 1:ops.general.num_islands]
    children   = [Individual[] for _ in 1:ops.general.num_islands]
    return generational_loop(data, ops, population, children)
end

function generational_loop(data::Vector{Vector{Float64}}, ops, start_pop::Vector{Individual})
    population = [
        start_pop[isle:ops.general.num_islands:end]
        for isle in 1:ops.general.num_islands
    ]
    children   = [Individual[] for _ in 1:ops.general.num_islands]
    return generational_loop(data, ops, population, children)
end

function generational_loop(data::Vector{Vector{Float64}}, ops, start_pop::Vector{String})
    population = [Individual[] for _ in 1:ops.general.num_islands]
    start_pop  = Node[string_to_node(eq, ops) for eq in start_pop]
    children   = [
        [
            Individual(node)
            for node in start_pop[isle:ops.general.num_islands:end]
        ]
        for isle in 1:ops.general.num_islands
    ]
    return generational_loop(data, ops, population, children)
end

function generational_loop(data::Vector{Vector{Float64}}, ops,
    population::Vector{Vector{Individual}},
    children::Vector{Vector{Individual}},
)
    @assert length(data) == ops.data_descript.n_vars + 1 "please use the data variable which was returned by the Options constructor."

    # prepare and check bank_of_terms and start_pop # ----------------------------------------------
    bank_of_terms  = Node[string_to_node(eq, ops) for eq in ops.grammar.bank_of_terms]

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

    # misc # ---------------------------------------------------------------------------------------
    hall_of_fame = Individual[]                      # initialize data structures

    prog_dict = OrderedDict(k => Float64[]           # initialize progress dict
        for k in ["time", "generation", "mean age", "cur_max_compl", "mean compl",
                  ["min " * string(m) for m in keys(ops.measures)]...]
    )

    stdin_reader = watch_stream(stdin)               # prepare user interrupt

    null_node = Node(0.0)                            # create the null node for some manual garbification

    t_start  = time()
    gen      = 0.0
    stop_msg = ""

# ==================================================================================================
# start generational loop
# ==================================================================================================
    while true
        gen += 1.0

        cur_max_compl = maximum(indiv.measures[:compl] for indiv in hall_of_fame; init=ops.grammar.min_compl)

        foreach(indiv -> indiv.age += 1, hall_of_fame)

        for isle in 1:ops.general.num_islands
            one_isle_one_generation!(
                population[isle],
                children[isle],
                bank_of_terms,
                data,
                ops,
                ops.general.fitting_island_function(isle) ? ops.fitting.max_iter : 0,
                cur_max_compl
            )
        end

# ==================================================================================================
# inter-isle
# ==================================================================================================
        # remove doubles across islands # ----------------------------------------------------------
        if ops.general.remove_doubles_across_islands && ops.general.remove_doubles_sigdigits > 0
            remove_doubles_across_islands!(population, ops)
        end

        # migration # ------------------------------------------------------------------------------
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

        # hall of fame # ---------------------------------------------------------------------------
        for isle in 1:ops.general.num_islands
            prepend!(hall_of_fame, copy.(population[isle]))
        end

        # extract objectives # -------------------------------------------------------------
        indiv_obj_vals = [
            Float64[
                indiv.measures[obj] for obj in ops.selection.hall_of_fame_objectives
            ]
            for indiv in hall_of_fame
        ]

        # apply niching
        indiv_obj_vals = [
            round.(indiv, sigdigits=ops.selection.hall_of_fame_niching_sigdigits)
            for indiv in indiv_obj_vals
        ]
        unique_inds = unique(i -> indiv_obj_vals[i], 1:length(indiv_obj_vals))
        keepat!(indiv_obj_vals, unique_inds)
        keepat!(hall_of_fame, unique_inds)

        selection_inds = first_pareto_front(indiv_obj_vals)

        keepat!(hall_of_fame, selection_inds)

# ==================================================================================================
# every couple of generations
# ==================================================================================================
        t_since = time() - t_start

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

                    inds_to_show = round.(Int64, collect(range(1, length(hall_of_fame), length=15)))
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

# ==================================================================================================
# island extinction
# ==================================================================================================
        if gen % ops.general.island_extinction_interval == 0
            emmigrate_island = rand(1:ops.general.num_islands)

            # both directions, decreasing with distance
            offsets = -trunc(Int64, 0.25 * ops.general.num_islands):trunc(Int64, 0.25 * ops.general.num_islands)
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

# ==================================================================================================
# termination criteria
# ==================================================================================================
        if gen >= ops.general.n_gens
            stop_msg = "reached maximum number of generations"
            break
        elseif time() - t_start >= ops.general.t_lim
            stop_msg = "reached time limit"
            break
        elseif ops.general.callback(hall_of_fame, population, gen, prog_dict, ops)
            stop_msg = "callback returned true"
            break
        elseif check_for_user_quit(stdin_reader)
            stop_msg = "graceful user interrupt"
            break
        end
    end # while true

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

    # Post-pare for return
    population = reduce(vcat, population)

    return hall_of_fame, population, prog_dict, stop_msg
end

function one_isle_one_generation!(pop, chil, bank_of_terms, data, ops, fit_iter, cur_max_compl)

    foreach(indiv -> indiv.age += 1, pop)

# ==================================================================================================
# genetic operations
# ==================================================================================================
    # create new children # ----------------------------------------------------------------
    while length(chil) + length(pop) < 0.6 * ops.general.pop_per_isle
        push!(chil,
            Individual(grow_equation(ops.grammar.init_tree_depth, ops))
        )
    end

    if length(pop) > 0.4 * ops.general.pop_per_isle

        # select parents
        shuffle!(pop)
        for i in 1:ops.general.n_children
            if ops.general.parent_selection
                push!(chil, copy(parent_selection(pop)))
            else
                push!(chil, copy(pop[mod1(i, length(pop))]))
            end
        end

        apply_genetic_operations!(chil, ops, bank_of_terms)

        for _ in 1:ops.general.n_refitting
            push!(chil, copy(rand(pop)))
        end
    end

    # grow them up # -----------------------------------------------------------------------
    if ops.general.multithreading
        Threads.@threads :greedy for ii in eachindex(chil)
            fit_individual!(chil[ii], data, ops, cur_max_compl, fit_iter)
        end
    else
        for ii in eachindex(chil)
            fit_individual!(chil[ii], data, ops, cur_max_compl, fit_iter)
        end
    end

    filter!(indiv -> indiv.valid, chil)

    # add children to population
    prepend!(pop, chil) # prepend, so that children are preferred during niching
    empty!(chil)

    # remove individuals
    if ops.general.remove_doubles_sigdigits > 0
        remove_doubles!(pop, ops)
    end

    filter!(i -> i.age <= ops.general.max_age, pop)

# ==================================================================================================
# selection
# ==================================================================================================
    if length(pop) > ops.general.pop_per_isle
        selection_inds = Int64[]

        # extract objectives # -------------------------------------------------------------
        indiv_obj_vals = [
            Float64[
                indiv.measures[obj] for obj in ops.selection.selection_objectives
            ]
            for indiv in pop
        ]

        # apply niching
        indiv_obj_vals = [
            round.(indiv, sigdigits=ops.selection.population_niching_sigdigits)
            for indiv in indiv_obj_vals
        ]
        unique_inds = unique(i -> indiv_obj_vals[i], 1:length(indiv_obj_vals))
        keepat!(indiv_obj_vals, unique_inds)
        keepat!(pop, unique_inds)

        # determine rank and crowding for all individuals
        ranks    = non_dominated_sort(indiv_obj_vals)
        crowding = crowding_distance(indiv_obj_vals)

        for i in eachindex(pop)
            pop[i].rank     = ranks[i]
            pop[i].crowding = crowding[i]
        end

        # Pareto selection # -------------------------------------------------------------------
        if ops.selection.n_pareto_select_per_isle > 0
            # sort!(pop)
            # append!(selection_inds, 1:ops.selection.n_pareto_select_per_isle)

            n_front = 1
            while true
                n_required = ops.selection.n_pareto_select_per_isle - length(selection_inds)
                n_required > 0 || break
                front = findall(==(n_front), ranks)
                if n_required < length(front)
                    front = wsample(front, replace([crowding[f] for f in front], Inf => 1e100, 0.0 => 1e-100), n_required, replace=false)
                end
                append!(selection_inds, front)
                n_front += 1
            end

        end

        # tournament selection # ---------------------------------------------------------------
        if ops.general.pop_per_isle - ops.selection.n_pareto_select_per_isle > 0
            remaining_inds = setdiff(eachindex(pop), selection_inds)

            if ops.general.pop_per_isle - length(selection_inds) < length(remaining_inds)
                indiv_obj_vals = normalize_objectives(indiv_obj_vals)
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
        keepat!(pop, selection_inds)
    end
end

# ==================================================================================================
# user interrupt
# ==================================================================================================

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
