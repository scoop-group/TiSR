
mutable struct Individual   #{T}
    node::Node

    compl::Float64
    recursive_compl::Float64
    age::Float64

    n_params::Float64

    ms_processed_e::Float64 # T

    mae::Float64            # T
    mse::Float64            # T
    max_ae::Float64         # T

    minus_r2::Float64       # T

    mare::Float64           # T
    q75_are::Float64        # T
    max_are::Float64        # T

    valid::Bool

    function Individual(node, data, ops)
        for _ in 1:3
            simplify_unary_of_param!(node)
            simplify_binary_of_param!(node)
            simplify_binary_across_1_level!(node, ops)
        end

        trim_to_max_compl!(node, ops)
        reorder_add_n_mul!(node, ops)

        indiv = new(node)

        if !isempty(ops.illegal_dict) && !check_legel(node, ops)
            indiv.valid = false
            return indiv
        end

        fit_n_eval!(indiv, data, ops)

        indiv.valid || return indiv

        indiv.age = 0.0
        indiv.compl = count_nodes_unique(node)
        indiv.recursive_compl = recursive_compl(node, ops)
        indiv.n_params = length(list_of_param_nodes(node))

        return indiv
    end
end

# ==================================================================================================
# some helpers
# ==================================================================================================
""" prints or displays an Individual
"""
Base.show(io::IO, indiv::Individual) = println(io, "Indiv($(node_to_string(indiv.node, __operators; sigdigits=3)))")

""" Compares two individuals and returns true if they are approx. same.
"""
function Base.isapprox(indiv1::Individual, indiv2::Individual; rtol=5)
    Base.isapprox(indiv1.mae, indiv2.mae, rtol=rtol) && Base.isapprox(indiv1.mse, indiv2.mse, rtol=rtol)
end

# ==================================================================================================
# different remove doubles -> maybe combine somehow with views?
# ==================================================================================================
""" Removes similar individual from vector of individuals. The individuals are considered
    similar, if they have a similar MAE and MSE. Of those, the one with the least complexity is
    chosen.
    CAVEAT: The current approach uses as the current best_of_cluster, which is replaced if a
    shorter similar one is found, rather than first searching for all similar ones and then
    choosing the one with the least complexity. This may lead to moving clusters, which should
    not be a problem, if the rtol is chosen small enough.
"""
function remove_doubles!(individs::Vector{Individual}, ops)
    eachind = collect(eachindex(individs))
    unique_inds = Int64[]

    while !isempty(eachind)
        best_of_cluster = popfirst!(eachind)

        i = 1
        while i <= length(eachind)
            if isapprox(
                individs[best_of_cluster],
                individs[eachind[i]],
                rtol=ops.general.prevent_doubles
            )
                next_ind = popat!(eachind, i)
                if individs[next_ind].compl < individs[best_of_cluster].compl
                    best_of_cluster = next_ind
                end
            else
                i += 1
            end
        end
        push!(unique_inds, best_of_cluster)
    end
    sort!(unique_inds)
    keepat!(individs, unique_inds)
end

""" Removes similar individual from vector of vector of individuals. The individuals are
    considered similar, if they have a similar MAE and MSE. Of those, the one with the least
    complexity is chosen.
    CAVEAT: The current approach uses as the current best_of_cluster, which is replaced if a
    shorter similar one is found, rather than first searching for all similar ones and then
    choosing the one with the least complexity. This may lead to moving clusters, which should
    not be a problem, if the rtol is chosen small enough.
"""
function remove_doubles_across_islands!(individs::Vector{Vector{Individual}}, ops)
    get_indiv = (individs, nest_ind) -> individs[nest_ind[1]][nest_ind[2]]

    eachind = [(isle, ind) for isle in 1:ops.general.num_islands for ind in eachindex(individs[isle])]
    unique_inds = Tuple{Int64, Int64}[]

    while !isempty(eachind)
        best_of_cluster = popfirst!(eachind)

        i = 1
        while i <= length(eachind)
            if isapprox(
                get_indiv(individs, best_of_cluster),
                get_indiv(individs, eachind[i]),
                rtol=ops.general.prevent_doubles
            )
                next_ind = popat!(eachind, i)
                if get_indiv(individs, next_ind).compl < get_indiv(individs, best_of_cluster).compl
                    best_of_cluster = next_ind
                end
            else
                i += 1
            end
        end
        push!(unique_inds, best_of_cluster)
    end

    unique_inds = [
        [
            unique_inds[i][2]
            for i in 1:length(unique_inds)
            if unique_inds[i][1] == isle
        ] for isle in 1:ops.general.num_islands
    ]

    foreach(isle -> sort!(unique_inds[isle]), 1:ops.general.num_islands)
    foreach(isle -> keepat!(individs[isle], unique_inds[isle]), 1:ops.general.num_islands)
end
