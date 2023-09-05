
using OrderedCollections


include("../src/options.jl")
include("../src/node_n_eval_n_utilities.jl")

function hardcoded_equations(ops)
    eqs_dict = OrderedDict{String, Node}()

    eqs_dict["1_simple"] = Node(                                                                # manually create a simple node
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

    eqs_dict["2_simple_w_var"] = Node(                                                              # manually create a second node and use different rounding
        findfirst(isequal(-), ops.binops),
        Node(
            findfirst(isequal(*), ops.binops),
            Node(1.45454545),
            Node(-42.0)
        ),
        Node(findfirst(isequal(exp), ops.unaops),
            Node(
                findfirst(isequal(pow_abs), ops.binops),
                Node(1.0),
                Node(1)
            )
        ))

    eqs_dict["3_simple_w_var_directed_acyclic_graph"] = deepcopy(eqs_dict["2_simple_w_var"])
    eqs_dict["3_simple_w_var_directed_acyclic_graph"].lef.rig = eqs_dict["3_simple_w_var_directed_acyclic_graph"].rig

    eqs_dict["4_complex"] = Node(
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
                    findfirst(isequal(pow_abs), ops.binops),
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

    eqs_dict["5_complex_directed_acyclic_graph"] = deepcopy(eqs_dict["4_complex"])
    eqs_dict["5_complex_directed_acyclic_graph"].lef.lef.lef = eqs_dict["5_complex_directed_acyclic_graph"].lef.rig

    return eqs_dict
end
