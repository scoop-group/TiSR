
# TODO: integration test fit_individual
# TODO: integration test cur_max_compl
# TODO: test remove_doubles_across_islands!(individs::Vector{Vector{Individual}}, ops)

data = rand(100, 10)
ops, data_vect = Options(
    data,
    general                            = general_params(
        remove_doubles_sigdigits       = 2,
    )
)

@testset "Base.isless" begin
    indiv1 = TiSR.Individual(TiSR.grow_equation(rand(3:5), ops, method = :asym))
    indiv2 = TiSR.Individual(TiSR.grow_equation(rand(3:5), ops, method = :asym))

    indiv1.rank = 1
    indiv2.rank = 2
    @test indiv1 < indiv2

    indiv1.rank = 1
    indiv1.crowding = 2.0
    indiv2.rank = 1
    indiv2.crowding = 1.0
    @test indiv1 < indiv2

    vect = [indiv1, indiv2]
    @test vect == sort(vect)
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
                indiv = deepcopy(indivs[i])
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


