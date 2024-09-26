
function hardcoded_equations(ops)

    eqs_dict = OrderedDict{String, Node}()

    eqs_dict["1"] = Node(                                                                # manually create a simple node
        findfirst(isequal(+), ops.binops),
        Node(
            findfirst(isequal(*), ops.binops),
            Node(1.45454545),
            Node(-42.0)
        ),
        Node(
            findfirst(isequal(sin), ops.unaops),
            Node(1.0)
        ))

    eqs_dict["2"] = Node(                                                              # manually create a second node and use different rounding
        findfirst(isequal(-), ops.binops),
        Node(
            findfirst(isequal(*), ops.binops),
            Node(1.45454545),
            Node(-42.0)
        ),
        Node(findfirst(isequal(exp), ops.unaops),
            Node(
                findfirst(isequal(^), ops.binops),
                Node(1.0),
                Node(1)
            )
        ))

    eqs_dict["3"] = deepcopy(eqs_dict["2"])
    eqs_dict["3"].lef.rig = deepcopy(eqs_dict["3"].rig)

    eqs_dict["4"] = Node(
        findfirst(isequal(exp), ops.unaops),
        Node(
            findfirst(isequal(*), ops.binops),
            Node(
                findfirst(isequal(+), ops.binops),
                Node(
                    findfirst(isequal(*), ops.binops),
                    Node(1),
                    Node(0.28282828),
                ),
                Node(
                    findfirst(isequal(^), ops.binops),
                    Node(8),
                    Node(0.9984357489357),
                ),
            ),
            Node(
                findfirst(isequal(sin), ops.unaops),
                Node(
                    findfirst(isequal(-), ops.binops),
                    Node(
                        findfirst(isequal(*), ops.binops),
                        Node(4),
                        Node(19.9832754918),
                    ),
                    Node(
                        findfirst(isequal(/), ops.binops),
                        Node(9),
                        Node(94.2394578e10),
                    ),
                ),
            ),
        ),
    )

    eqs_dict["5"] = deepcopy(eqs_dict["4"])
    eqs_dict["5"].lef.lef.lef = deepcopy(eqs_dict["5"].lef.rig)

    return eqs_dict
end




# gaussian
# magnetic manipulator
# vanderwaals
# gravity
# relativity
# spring force
#

# BenchmarkProblem1 = (
#     name = "I.12.4 - magnitude of electric field",
#     str = "E = q1 / (4 * pi * epsilon * r^2)",
#     func = (q1, r) -> q1 / (4 * pi * 8.854e-12 * r^2),
#     node = Node(
#             findfirst(==(/), ops.binops),
#                   Node(1),
#                   Node(findfirst(*),
#                        Node(4.0),
#                        Node(findfirst(
#                  )
#     )
# )



