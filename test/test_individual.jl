
# TODO: integration test cur_max_compl
# TODO: integration test fit_iter
# TODO: integration test fit_individual!(indiv, data, ops, cur_max_compl, fit_iter)
# TODO: test fastcopy(indiv::Individual) = Individual(deepcopy(indiv.node))

data = rand(100, 10)
ops, data_vect = Options(
    data,
    general = general_params(
        num_islands              = 3,
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
