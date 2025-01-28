
# TODO: integration test for Individual()

data = rand(100, 10)
ops, data_vect = Options(
    data,
    general                            = general_params(
        remove_doubles_sigdigits       = 2,
    ),
)

@testset "Base.copy(indiv::Individual)" begin

    indiv = nothing
    while isnothing(indiv)
        node  = TiSR.grow_equation(rand(4:7), ops, method=:full)
        list_of_param = TiSR.list_of_param_nodes(node)
        !isempty(list_of_param) || continue
        indiv = TiSR.Individual(node)
        TiSR.fit_individual!(indiv, data_vect, ops, ops.grammar.max_compl, 0)
        indiv.valid || continue
    end

    indiv2 = copy(indiv)

    @test !(indiv === indiv2)
    for (k, v) in indiv
        @test string(v) == string(getfield(indiv2, k))
    end

    indiv.measures[:compl] += 1
    @test indiv.measures[:compl] != indiv2.measures[:compl]
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
            indiv = TiSR.Individual(node)
            TiSR.fit_individual!(indiv, data_vect, ops, ops.grammar.max_compl, 0)
            indiv.valid || continue
            push!(indivs, indiv)
        end

        # copy individuals but with different complexities
        for i in 1:5
            for j in 1:5
                indiv = copy(indivs[i])
                indiv.measures[:mae]   += rand() * 0.0001
                indiv.measures[:mse]   += rand() * 0.0001
                indiv.measures[:compl] += j
                push!(indivs, indiv)
            end
        end

        shuffle!(indivs)

        TiSR.remove_doubles!(indivs, ops)

        avg += (length(indivs) - avg)/iters
    end

    @test isapprox(avg, 5, rtol=0.1)
end


