

# TODO: test start_pop with finished equation, invalids, strings, individuals
# TODO: test bank_of_terms with finished equation and with illegal nesting
# TODO: test cur_max_compl
# TODO: test hall_of_fame migration
# TODO: test prog_dict
# TODO: test t_since
# TODO: test one_isle_one_generation!(pop, chil, bank_of_terms, data, ops, fit_iter, cur_max_compl)
# TODO: test max_age
# TODO: test fitting_isle/fit_iter
# TODO: test perform_migration!(population, ops)
# TODO: test perform_island_extinction!(population, ops)
# TODO: test hall_of_fame
# TODO: many integration tests

@testset "stop_msg and callback" begin

    # test n_gens # --------------------------------------------------------------------------------
    data = rand(100, 10)
    ops, data_vect = Options(
        data,
        general=general_params(
            n_gens=10,
            print_progress     = false,
            plot_hall_of_fame  = false,
            print_hall_of_fame = false,
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
            t_lim=10.0,
            print_progress     = false,
            plot_hall_of_fame  = false,
            print_hall_of_fame = false,
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
            callback = (hall_of_fame, population, gen, prog_dict, ops) -> any(i.measures[:compl] <= 7 && i.measures[:mare] < 1e-5 for i in hall_of_fame),
            print_progress     = false,
            plot_hall_of_fame  = false,
            print_hall_of_fame = false,
        )
    );

    hall_of_fame, population, prog_dict, stop_msg = generational_loop(data_vect, ops);

    @test stop_msg == "callback returned true"
    @test any(i.measures[:compl] <= 7 && i.measures[:mare] <= 1e-5 for i in hall_of_fame)
end

