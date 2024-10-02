
# TODO: integration test for Individual()

data = rand(100, 10)
ops, data_vect = Options(
    data,
    general                            = general_params(
        remove_doubles_sigdigits       = 2,
        remove_doubles_across_islands  = false, # TODO: experiment
    ),
)


@testset "Base.copy(indiv::Individual)" begin

    indiv = nothing
    while isnothing(indiv)
        node  = TiSR.grow_equation(rand(4:7), ops, method=:full)
        list_of_param = TiSR.list_of_param_nodes(node)
        !isempty(list_of_param) || continue
        indiv = TiSR.Individual(node, data_vect, ops, Inf)
        indiv.valid || continue
    end

    indiv2 = copy(indiv)

    for field in fieldnames(typeof(indiv))
        @test getfield(indiv, field) == getfield(indiv, field)

        if field == :node
            list_of_param2 = TiSR.list_of_param_nodes(indiv2.node)
            list_of_param2[1].val += 1.0
            @test indiv.node != indiv2.node
        elseif field == :valid
        else
            setfield!(indiv2, field, getfield(indiv, field) + 1.0)
            @test getfield(indiv, field) != getfield(indiv2, field)
        end
    end
end

@testset "_remove_doubles_helper!" begin
    
    # create population island with mae and mse such that all are different
    indiv_obj_vals = [[i, 10.0-i] for i in 1:5]

    # copy individuals but with different complexities
    for i in 1:5
        for j in 1:5
            push!(indiv_obj_vals, [(indiv_obj_vals[i] .+ rand(2) .* 0.01)..., j])
        end
    end

    deleteat!(indiv_obj_vals, 1:5)

    unique_inds = collect(eachindex(indiv_obj_vals))
    for i in eachindex(indiv_obj_vals)
        indiv_obj_vals[i] .= round.(indiv_obj_vals[i], sigdigits=ops.general.remove_doubles_sigdigits)
    end

    TiSR._remove_doubles_helper!(indiv_obj_vals, unique_inds)

    @test all(isone(i[end]) for i in indiv_obj_vals)
end

@testset "remove_doubles!" begin
    
    avg = 0.0
    
    for iters in 1:100
        indivs = TiSR.Individual[]

        while length(indivs) < 5
            node  = TiSR.grow_equation(rand(4:7), ops)

            indiv = TiSR.Individual(node, data_vect, ops, Inf)
            indiv.valid || continue
            push!(indivs, indiv)

        end

        # copy individuals but with different complexities
        for i in 1:5
            for j in 1:5
                indiv = copy(indivs[i])
                indiv.mae   += rand() * 0.0001
                indiv.mse   += rand() * 0.0001
                indiv.compl += j
                push!(indivs, indiv)
            end
        end

        shuffle!(indivs)

        TiSR.remove_doubles!(indivs, ops)

        avg += (length(indivs) - avg)/iters
    end

    @test isapprox(avg, 5, rtol=0.1)
end

@testset "remove_doubles_across_islands!" begin
    
    avg = 0.0
    
    for iters in 1:100
        indivs = TiSR.Individual[]

        while length(indivs) < 5
            node  = TiSR.grow_equation(rand(4:7), ops)

            indiv = TiSR.Individual(node, data_vect, ops, Inf)
            indiv.valid || continue
            push!(indivs, indiv)

        end

        # copy individuals but with different complexities
        for i in 1:5
            for j in 1:5
                indiv = copy(indivs[i])
                indiv.mae   += rand() * 0.0001
                indiv.mse   += rand() * 0.0001
                indiv.compl += j
                push!(indivs, indiv)
            end
        end

        shuffle!(indivs)

        populations_with_islands = [indivs[i:i+4] for i in 1:5:length(indivs)]

        TiSR.remove_doubles_across_islands!(populations_with_islands, ops)

        pop = reduce(vcat, populations_with_islands)

        avg += (length(pop) - avg)/iters
    end

    @test isapprox(avg, 5, rtol=0.1)
end


