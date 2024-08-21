
""" Builds the equations from top down using "full"-method. Parameters of unary operators are
    prevented and power to parameters are enforced (if specified in ops). The probability of the
    operators are adapted in case a particular nesting is prohibited according to illegal_dict.
"""
grow_equation(
    rem_depth::Int, ops::Options;
    unaweights=Float64[ops.p_unaops...], binweights=Float64[ops.p_binops...],
    no_param=false, method = "not full"
)::Node = grow_equation_(rem_depth, ops, unaweights, binweights, no_param, method)

function grow_equation_(rem_depth::Int, ops::Options, unaweights, binweights, no_param, method)::Node

    if rem_depth <= 1 || (method != "full" && rand() < 0.3^rem_depth)
        if !no_param && rand(Bool)
            next_node = Node(rand())
        else
            next_node = Node(rand(1:ops.data_descript.n_vars))
        end
    else

        rand_ = rand() * (sum(unaweights) + sum(binweights))

        if length(ops.unaops) > 0 && rand_ <= sum(unaweights)
            op_ind = findfirst(x -> sum(unaweights[1:x]) >= rand_, 1:length(ops.unaops))
            next_node = Node(1, op_ind)

            if length(keys(ops.illegal_dict)) > 0
                cur_fun = Symbol(ops.unaops[op_ind])
                ill_nexts = get(ops.illegal_dict, cur_fun, ())

                adjust_ops_weights!(unaweights, binweights, ill_nexts, ops)
            end

            next_node.lef =     grow_equation_(rem_depth - 1, ops, copy(unaweights), copy(binweights), true, method)
        else
            op_ind = findfirst(x -> (sum(unaweights) + sum(binweights[1:x])) >= rand_, 1:length(ops.binops))
            next_node = Node(2, op_ind)

            cur_fun = Symbol(ops.binops[op_ind])

            if length(keys(ops.illegal_dict)) > 0
                ill_nexts = get(ops.illegal_dict, cur_fun, ())

                adjust_ops_weights!(unaweights, binweights, ill_nexts, ops)
            end

            next_node.lef =     grow_equation_(rem_depth - 1, ops, unaweights, binweights, false, method)
            if ops.general.pow_abs_param && cur_fun in (:pow_abs, :^)
                next_node.rig = Node(rand())
            elseif next_node.lef.ari == -1
                next_node.rig = grow_equation_(rem_depth - 1, ops, copy(unaweights), copy(binweights), true, method)
            else
                next_node.rig = grow_equation_(rem_depth - 1, ops, copy(unaweights), copy(binweights), false, method)
            end
        end
    end
    return next_node
end

""" Check for illegal nestings in existing nodes.
"""
check_legal(node, ops; unaweights=Float64[ops.p_unaops...], binweights=Float64[ops.p_binops...]) = check_legal_(node, ops, unaweights, binweights)

function check_legal_(node, ops, unaweights, binweights)
    isempty(ops.illegal_dict) && return true
    node.ari <= 0 && return true

    if node.ari == 2
        (binweights[node.ind] == 0.0) && return false

        cur_fun = Symbol(ops.binops[node.ind])
        ill_nexts = get(ops.illegal_dict, cur_fun, ())
        adjust_ops_weights!(unaweights, binweights, ill_nexts, ops)

        legal = check_legal_(node.lef, ops, copy(unaweights), copy(binweights))

        !legal && return false

        cur_fun in (:pow_abs, :^) && ops.general.pow_abs_param && node.rig.ari > -1 && return false

        return check_legal_(node.rig, ops, copy(unaweights), copy(binweights))

    elseif node.ari == 1
        (unaweights[node.ind] == 0.0) && return false

        cur_fun = Symbol(ops.unaops[node.ind])
        ill_nexts = get(ops.illegal_dict, cur_fun, ())
        adjust_ops_weights!(unaweights, binweights, ill_nexts, ops)

        return check_legal_(node.lef, ops, copy(unaweights), copy(binweights))
    end
end

""" Adjust the weight-vectors for the operator selection probability of the grow_equation and
    check_legal function to prevent illegal nestings according to illegal_dict.
"""
function adjust_ops_weights!(unaweights, binweights, ill_nexts, ops)
    for fun in ill_nexts
        ind = findfirst(isequal(fun), ops.unaops)
        if !isnothing(ind)
            unaweights[ind] = 0.0
        end

        ind = findfirst(isequal(fun), ops.binops)
        if !isnothing(ind)
            binweights[ind] = 0.0
        end
    end
end

# ==================================================================================================
# genetic operations helpers
# ==================================================================================================
""" Return a random sub-node of a given node. Has multiple modes:
    0 -> any node
    1 -> at least one level above terminal
    2 -> at least two levels above terminal
"""
function random_node(node::Node; mode=0)

    maxim_tree_depth(node, minim=mode + 2) <= mode && return node

    while true
        node_elect, valid = random_node_(node, mode)
        valid && return node_elect
    end
end

""" Traversal for the random_node function. Implemented along the lines of
    SymbolicRegression.jl
"""
function random_node_(node::Node, mode)
    if node.ari <= 0
        valid = maxim_tree_depth(node, minim=mode + 2) > mode
        return node, valid
    end

    b = 0
    c = 0
    if node.ari >= 1
        b = count_nodes(node.lef)
    end
    if node.ari == 2
        c = count_nodes(node.rig)
    end

    i = rand(1:(1+b+c))
    if i <= b
        node_elect, valid = random_node_(node.lef, mode)
        if valid
            return node_elect, valid
        end
    elseif i == b + 1
        valid = maxim_tree_depth(node, minim=mode + 2) > mode
        return node, valid
    end
    if c > 0
        node_elect, valid = random_node_(node.rig, mode)
        return node_elect, valid
    else
        return node, false
    end
end

""" Entrace/switch function for the genetic operations. If nodes are shallow, only insert-,
    point, and addterm are allowed. For crossover, if there are not any nodes deeper than 2
    left, point mutation is applied. Some of the mutations may not mutate by change. This is not
    tracked or avoided. In this case, the individual is just refit with it's new parameters as
    start point.
"""
function apply_genetic_operations!(
    nodes,
    ops;
    (one_node_muts!)=[
        insert_mutation!, point_mutation!, addterm_mutation!, # ___ below ar deeper than 2
        hoist_mutation!, subtree_mutation!, drastic_simplify!
    ]
)
    eachind = collect(eachindex(nodes))

    if ops.general.always_drastic_simplify
        drastic_inds = findall(drastic_simplify!(n, ops, potential=true, threshold=1e-3) for n in nodes)
        drastic_nodes = [copy_node(nodes[i]) for i in drastic_inds]
        foreach(n -> drastic_simplify!(n, ops, potential=false, threshold=1e-3), drastic_nodes)

        for _ in 1:3
            simplify_unary_of_param!.(drastic_nodes)
            simplify_binary_of_param!.(drastic_nodes)
            simplify_binary_across_1_level!.(drastic_nodes, Ref(ops))
        end

        append!(nodes, drastic_nodes)
    end

    while !isempty(eachind) # TODO: maybe invert -> sample mutation first and then find appropriate node(s) -> closer to selected mutation probabilities
        node = nodes[popfirst!(eachind)]

        until_mut = maxim_tree_depth(node, minim=3) > 2 ? length(ops.mutation) : 3

        rand_mutation = rand() * ops.mutation[until_mut]

        if rand_mutation < ops.mutation[6]
            one_node_muts![findfirst(i -> rand_mutation < ops.mutation[i], 1:6)](node, ops)

        elseif rand_mutation < ops.mutation[7]
            try
                simplify_w_symbolic_utils!(node, ops; through_polyform=rand() < 0.2)               # add as parameter?
            catch
            end

        elseif length(eachind) > 1 && rand_mutation <= ops.mutation[8]
            ind = findfirst(i -> maxim_tree_depth(nodes[i], minim=3) > 2, eachind)
            if isnothing(ind)
                point_mutation!(node, ops)
            else
                node2 = nodes[popat!(eachind, ind)]
                crossover_mutation!(node, node2, ops)
            end
        end
    end
end

""" Helper function, which decidies whether a mutation is applied to the left or right child node.
"""
function mutate_left(node, ops, min_depth)
    (node.ari == 1                                                                    # -> definetly lef, if ari == 1
     ||
     (node.ari == 2 &&
      ((ops.general.pow_abs_param && Symbol(ops.binops[node.ind]) in (:pow_abs, :^))  # -> definetly lef, if pow_abs_param
       ||
       !(maxim_tree_depth(node.rig, minim=min_depth + 2) >= min_depth))               # -> definetly lef, if rig not enough
      ||
      (rand(Bool) && maxim_tree_depth(node.lef, minim=min_depth + 2) >= min_depth)    # -> maybe lef, if lef enough
    )
    )
end

# ==================================================================================================
# genetic operations
# ==================================================================================================
""" Changes one node to another without modifying the structure of the node.
"""
function point_mutation!(node, ops)
    node_elect = random_node(node, mode=0)

    if node_elect.ari == 2
        node_elect.ind = wsample(1:length(ops.binops), collect(ops.p_binops)) # TODO: unnecessary collect

        if ops.general.pow_abs_param && Symbol(ops.binops[node_elect.ind]) in (:pow_abs, :^) # TODO: move this to grammar check
            node_elect.rig = Node(rand())
        end

    elseif node_elect.ari == 1
        node_elect.ind = wsample(1:length(ops.unaops), collect(ops.p_unaops)) # TODO: unnecessary collect

    elseif node_elect.ari == 0
        node_elect.ind = rand(1:ops.data_descript.n_vars)

    elseif node_elect.ari == -1
        node_elect.val *= rand_mult(;minn=0.5, maxx=2.0)
    end
end

""" Inserts or prepend an operation.
"""
function insert_mutation!(node, ops; subtree_depth=2)
    new_node = grow_equation(subtree_depth, ops, method = "full")
    lefrig1 = mutate_left(new_node, ops, 1) ? :lef : :rig

    if rand(1:count_nodes(node)) == 1
        orig_node = copy_node(node)
        copy_node_wo_copy!(node, new_node)

        setfield!(node, lefrig1, orig_node)
    else
        node_elect = random_node(node, mode=1)
        lefrig2 = mutate_left(node_elect, ops, 1) ? :lef : :rig

        setfield!(new_node, lefrig1, getfield(node_elect, lefrig2))
        setfield!(node_elect, lefrig2, new_node)
    end
end

""" Removes an operation.
"""
function hoist_mutation!(node, ops)
    node_elect = random_node(node, mode=2)

    lefrig1 = mutate_left(node_elect, ops, 2) ? :lef : :rig
    sub1 = getfield(node_elect, lefrig1)

    lefrig2 = mutate_left(sub1, ops, 1) ? :lef : :rig
    sub2 = getfield(sub1, lefrig2)
    setfield!(node_elect, lefrig1, sub2)
end

""" Combines two nodes to create two new ones.
"""
function crossover_mutation!(node1, node2, ops)
    node_elect1 = random_node(node1, mode=1)
    node_elect2 = random_node(node2, mode=1)

    lefrig1 = mutate_left(node_elect1, ops, 2) ? :lef : :rig
    lefrig2 = mutate_left(node_elect2, ops, 2) ? :lef : :rig

    temp = getfield(node_elect1, lefrig1)
    setfield!(node_elect1, lefrig1, copy_node(getfield(node_elect2, lefrig2)))
    setfield!(node_elect2, lefrig2, copy_node(temp))
end

# redundand mutations # ----------------------------------------------------------------------------
""" Replaces a subtree of a random subtree.
"""
function subtree_mutation!(node, ops; subtree_depth=3) # TODO: mabe rand depth
    node_elect = random_node(node, mode=1)
    lefrig = mutate_left(node_elect, ops, 1) ? :lef : :rig
    setfield!(node_elect, lefrig, grow_equation(subtree_depth, ops))
end

""" Adds a top-level term.
"""
function addterm_mutation!(node, ops; subtree_depth=2)
    orig_node = copy_node(node)
    copy_node_wo_copy!(node, Node(2, findfirst(isequal(+), ops.binops)))
    node.lef = orig_node
    node.rig = grow_equation(subtree_depth, ops)
end

