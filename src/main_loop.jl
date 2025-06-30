
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

    t_start        = time()
    gen            = 0.0
    stop_msg       = ""
    eval_counters  = fill(0, ops.general.num_islands)
    expression_log = [Dict{UInt64, Int8}() for _ in 1:ops.general.num_islands]

    #reset_timer!(to) # @timeit

# ==================================================================================================
# start generational loop
# ==================================================================================================
    while true
        gen += 1.0

        cur_max_compl = maximum(indiv.measures[:compl] for indiv in hall_of_fame; init=ops.grammar.min_compl)

        foreach(indiv -> indiv.age += 1, hall_of_fame)

        #@timeit to "one generation" begin
            if ops.general.multithreading
                Threads.@threads :greedy for isle in 1:ops.general.num_islands
                    eval_counters[isle] += one_isle_one_generation!(
                        population[isle],
                        children[isle],
                        bank_of_terms,
                        data,
                        ops,
                        cur_max_compl,
                        expression_log[isle],
                        isle,
                    )
                end
            else
                for isle in 1:ops.general.num_islands
                    eval_counters[isle] += one_isle_one_generation!(
                        population[isle],
                        children[isle],
                        bank_of_terms,
                        data,
                        ops,
                        cur_max_compl,
                        expression_log[isle],
                        isle
                    )
                end
            end
        #end # @timeit

# ==================================================================================================
# inter-isle
# ==================================================================================================
        #@timeit to "inter isle" begin
            if gen % ops.general.migration_interval == 0
                perform_migration!(population, ops)
            end

            if gen % ops.general.hall_of_fame_migration_interval == 0
                indiv = deepcopy(rand(hall_of_fame))
                indiv.age = 0
                push!(population[rand(1:ops.general.num_islands)], indiv)
            end

            #@timeit to "garbage collect expression_log" begin
                if gen % ops.general.seen_merge_interval == 0
                   merged_log = mergewith(safe_plus, expression_log...)
                   #@show length(merged_log) # DEBUG expression_log
                   expression_log = [deepcopy(merged_log) for _ in 1:ops.general.num_islands]
                end
                if gen % ops.general.seen_forget_interval == 0
                    for exp_log in expression_log
                       while true
                           map!(v -> v ÷ Int8(2), values(exp_log))
                           filter!(kv -> kv[2] > 0, exp_log)
                           length(exp_log) < ops.general.seen_gc_length && break
                       end
                    end
                end
            #end # @timeit

            if gen % ops.general.island_extinction_interval == 0
                perform_island_extinction!(population, gen, ops)
            end

            #@timeit to "hall_of_fame_selection" begin
                perform_hall_of_fame_selection!(hall_of_fame, population, ops)
            #end # @timeit
        #end # @timeit
# ==================================================================================================
# every couple of generations
# ==================================================================================================
        t_since = time() - t_start

        #@timeit to "every couple of gens" begin
            if isempty(prog_dict["time"]) || t_since - prog_dict["time"][end] > 5.0
                #println("reject_rate -> ", reject_rate[1] / reject_rate[2]) # DEBUG expression_log
                #reject_rate .= 0                                            # DEBUG expression_log

                # GC.gc() # no idea why that is necessary
                clean_trash_nodes!(population, null_node)
                clean_trash_nodes!(hall_of_fame, null_node)

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
                    println("\n$(round(Int64, t_since ÷ 60)) min $(round(Int64, t_since % 60)) sec | type q and enter to finish early")
                    if !isempty(prog_dict["time"])
                        println("\n" * @sprintf("%.3e evals per second", sum(eval_counters) / (t_since - prog_dict["time"][end])))
                    end
                    eval_counters .= 0
                end

                # check for replace_inf behavior
                obj_replace_inf_error = findfirst(
                    all(indiv.measures[obj] == ops.general.replace_inf for indiv in hall_of_fame)
                    for obj in ops.selection.hall_of_fame_objectives
                )
                if !isnothing(obj_replace_inf_error)
                    @warn """all individual in the hall_of_fame have the ops.general.replace_inf for
                        objective '$(ops.selection.hall_of_fame_objectives[obj_replace_inf_error])'.
                        Please check for correctness or increase the ops.general.replace_inf value."""
                end

                for k in keys(get_for_prog)
                    push!(prog_dict[k], get_for_prog[k])
                end

                #@timeit to "plot and show equations" begin
                    if ops.general.show_hall_of_fame
                        show_hall_of_fame(hall_of_fame, ops)
                    end
                #end # @timeit
            end
        #end # @timeit

# ==================================================================================================
# termination criteria
# ==================================================================================================
        if gen >= ops.general.n_gens
            stop_msg = "reached maximum number of generations"
            break
        elseif time() - t_start >= ops.general.t_lim
            stop_msg = "reached time limit"
            break
        elseif ops.general.callback(hall_of_fame, population, gen, t_since, prog_dict, ops)
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

    if ops.general.show_hall_of_fame
        show_hall_of_fame(hall_of_fame, ops)
    end

    close_reader!(stdin_reader)

    # Post-pare for return
    population = reduce(vcat, population)

    #display(to) # @timeit

    return (
        hall_of_fame, population, prog_dict, stop_msg
    )
end

function one_isle_one_generation!(pop, chil, bank_of_terms, data, ops, cur_max_compl, expression_log, isle; trial=1)

    eval_counter = 0

    foreach(indiv -> indiv.age += 1, pop)

    # create new children # ----------------------------------------------------------------
    #@timeit to "new children" begin
        while length(chil) + length(pop) < ops.general.pop_per_isle
            push!(chil,
                Individual(grow_equation(
                    ops.grammar.init_tree_depth,
                    ops,
                    param_prob = ops.general.fitting_island_function(isle) ? 2.0 : 0.0
                ))
            )
        end
    #end # @timeit

    # genetic operations # -------------------------------------------------------------------------
    #@timeit to "genetic ops" begin
        if 0.8 * ops.general.pop_per_isle < length(pop) < 1.2 * ops.general.pop_per_isle
            for _ in 1:ops.general.n_children
                push!(chil, fastcopy(parent_selection(pop)))
            end

            apply_genetic_operations!(chil, ops, bank_of_terms)

            for _ in 1:ops.general.n_refitting
                push!(chil, fastcopy(rand(pop)))
            end
        end
    #end # @timeit

    # fitting and evaluation # ---------------------------------------------------------------------
    #@timeit to "Individual" begin
        for ii in eachindex(chil)
            eval_counter += fit_individual!(chil[ii], data, ops, cur_max_compl, expression_log, isle)
        end
    #end # @timeit

    filter!(indiv -> indiv.valid, chil)

    # add children to population
    prepend!(pop, chil) # prepend, so that children are preferred during niching
    empty!(chil)

    # remove individuals # -------------------------------------------------------------------------
    filter!(i -> i.age <= ops.general.max_age, pop)

    # selection # ----------------------------------------------------------------------------------
    if length(pop) > ops.general.pop_per_isle
        #@timeit to "selection" begin
            perform_population_selection!(pop, ops, isle)
        #end # @timeit
    elseif isempty(pop) && trial < 100
        println("all individuals filtered, redoing generation")
        one_isle_one_generation!(pop, chil, bank_of_terms, data, ops, cur_max_compl, expression_log, isle, trial=trial+1)
    elseif isempty(pop)
        throw("Failed redoing the generation 100 times. All individuals are filtered out. Possible filters: illegal_dict, custom_check_legal, nonfinite evaluation, some of the defined measues is nonfinite.")
    end
    return eval_counter
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

""" Performs the migration of individuals form one island to another
"""
function perform_migration!(population, ops)
    emmigrate_island = rand(1:ops.general.num_islands)
    immigrate_island = mod1(emmigrate_island + rand((1, -1)), ops.general.num_islands)

    if !isempty(population[emmigrate_island])
        push!(
            population[immigrate_island],
            popat!(
                population[emmigrate_island],
                rand(1:length(population[emmigrate_island]))
            )
        )
    end
end

""" Perform island extinction. Chooses a random island, and distributes the best
    indiviuals into the surrounding islands, while emptying the chosen island.
"""
function perform_island_extinction!(population, gen, ops)
    if ops.general.island_extinction_rotation
        emmigrate_island = Int(div(gen, ops.general.island_extinction_interval))
    else
        emmigrate_island = rand(1:ops.general.num_islands)
    end

    offsets = 1:ops.general.migrate_after_extinction_dist

    sort!(population[emmigrate_island])
    for _ in 1:ops.general.migrate_after_extinction_num
        indiv = popfirst!(population[emmigrate_island])
        immigrate_island = mod1(
            emmigrate_island + rand(offsets),
            ops.general.num_islands
        )
        push!(population[immigrate_island], indiv)
        isempty(population[emmigrate_island]) && break
    end
    empty!(population[emmigrate_island])
end

""" Plots the current hall_of_fame fame and prints a selection of the equations
    after simplifying them.
"""
function show_hall_of_fame(hall_of_fame, ops)
    #@timeit to "scatterplot" begin
        compl          = [indiv.measures[:compl]          for indiv in hall_of_fame]
        ms_processed_e = [indiv.measures[:ms_processed_e] for indiv in hall_of_fame]

        if any(!(1e-100 < m < 1e100) for m in ms_processed_e)
            println("ms_processed_e clamped to in-between 1e-100 and 1e100 for hall_of_fame plot")
            clamp!(ms_processed_e, 1e-100, 1e100)
        end

        if length(unique(ms_processed_e)) > 1
            ymin = 10^(floor(log10(minimum(ms_processed_e))))
            ymax = 10^(ceil(log10(maximum(ms_processed_e))))
        else
            ymin = 10^(floor(log10(minimum(ms_processed_e)))-1)
            ymax = 10^(ceil(log10(maximum(ms_processed_e)))+1)
        end

        plt = scatterplot(compl, ms_processed_e,
            yscale           = :log10,
            title            = "hall of fame",
            xlabel           = "compl",
            ylabel           = "log10 of ms_processed_e",
            marker           = :circle,
            unicode_exponent = false,
            xlim             = (0, ops.grammar.max_compl),
            ylim             = (ymin, ymax),
            compact          = true,
        )
    #end # @timeit

    #@timeit to "simplify equations to show" begin
        sort!(hall_of_fame, by=i->i.measures[:ms_processed_e], rev=true)

        inds_to_show = round.(Int64, collect(range(1, length(hall_of_fame), length=15)))
        unique!(inds_to_show)

        eq_strs = [simplify_to_string(hall_of_fame[i].node, ops, sigdigits=2) for i in inds_to_show]

        for i in eachindex(eq_strs)
            label!(plt, :r, i, replace(
                eq_strs[i],
                " " => "",
                r"(\d)\.0\b" => s"\1",
                r"(\d)\.0(\D)" => s"\1\2",
            ))
        end
    #end # @timeit
    display(plt)
end

""" Plus, which does not overflow integers, but replaces with typemax in that case.
"""
safe_plus(x, y) = typemax(x) - x > y ? x + y : typemax(x)
