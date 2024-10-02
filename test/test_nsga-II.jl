

@testset "stop_msg and callback" begin

    # test n_gens # --------------------------------------------------------------------------------
    data = rand(100, 10)
    ops, data_vect = Options(
        data,
        general=general_params(
            n_gens=10,
        )
    );

    hall_of_fame, population, prog_dict, stop_msg = generational_loop(data_vect, ops);

    @test stop_msg == "reached maximum number of generations"
    @test maximum(prog_dict["generation"]) <= 10
    @test maximum(i.age for i in population) <= 10
    @test maximum(i.age for i in hall_of_fame) <= 10

    # test t_lim # ---------------------------------------------------------------------------------
    data = rand(100, 10)
    ops, data_vect = Options(
        data,
        general=general_params(
            t_lim=10,
        )
    );

    hall_of_fame, population, prog_dict, stop_msg = generational_loop(data_vect, ops);

    @test stop_msg == "reached time limit"
    @test 10 < prog_dict["time"][end] < 15

    # test callback # ------------------------------------------------------------------------------
    data_matr = rand(100, 3)
    data_matr[:, end] .= 3.0 .* (data_matr[:, 1] .* 5.0 .+ data_matr[:, 2])

    ops, data_vect = Options(
        data_matr,
        general=general_params(
            t_lim=Inf,
            callback = (hall_of_fame, population, ops) -> any(i.compl <= 7 && i.mare < 1e-5 for i in hall_of_fame)
        )
    );

    hall_of_fame, population, prog_dict, stop_msg = generational_loop(data_vect, ops);
    
    @assert stop_msg == "callback returned true"
    @assert any(i.compl <= 7 && i.mare <= 1e-5 for i in hall_of_fame)
end

