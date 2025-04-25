
""" Builds the equations using either the full-method (method = :full) or a
    probabilistic, assymetric grow method (method = :asym).
"""
grow_equation(
    rem_depth::Int, ops::Options;
    unaweights = ones(length(ops.unaops)), binweights = ones(length(ops.binops)),
    method = :asym,
    param_prob = 2.0,
)::Node = grow_equation_(rem_depth, ops, unaweights, binweights, method, param_prob)

function grow_equation_(rem_depth::Int, ops::Options, unaweights, binweights, method, param_prob)::Node
    if rem_depth <= 1 || (method == :asym && rand() < 0.3^rem_depth)

        if rand() < (param_prob / (ops.data_descript.n_vars + param_prob))       # parameter twice as likely as any one variable
            next_node = Node(randn() * 10.0)                           # parameter
        else
            next_node = Node(rand(1:ops.data_descript.n_vars)) # variable
        end
    else
        rand_ = rand() * (sum(unaweights; init=0.0) + sum(binweights; init=0.0))

        if length(ops.unaops) > 0 && rand_ <= sum(unaweights; init=0.0)
            op_ind = findfirst(x -> sum(unaweights[1:x]; init=0.0) >= rand_, 1:length(ops.unaops))
            next_node = Node(1, op_ind)

            next_node.lef = grow_equation_(rem_depth - 1, ops, unaweights, binweights, method, param_prob)
        else
            op_ind = findfirst(x -> (sum(unaweights; init=0.0) + sum(binweights[1:x]; init=0.0)) >= rand_, 1:length(ops.binops))
            next_node = Node(2, op_ind)

            next_node.lef = grow_equation_(rem_depth - 1, ops, unaweights, binweights, method, param_prob)
            next_node.rig = grow_equation_(rem_depth - 1, ops, unaweights, binweights, method, param_prob)
        end
    end
    return next_node
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
        node_elect = random_node_(node)
        if maxim_tree_depth(node_elect, minim=mode + 2) > mode
            return node_elect
        end
    end
end

""" Traversal for the random_node function. Implemented along the lines of
    SymbolicRegression.jl.
"""
function random_node_(node)
    node.ari <= 0 && return node

    n_all = count_nodes(node)
    rand_ = rand() * n_all
    rand_ <= 1.0 && return node

    n_lef = count_nodes(node.lef)
    rand_ <= n_lef + 1.0 && return random_node_(node.lef)

    return random_node_(node.rig)
end

""" Entrace/switch function for the genetic operations. If nodes are shallow, only insert-,
    point, and addterm are allowed. For crossover, if there are not any nodes deeper than 2
    left, point mutation is applied. Some of the mutations may not mutate by change. This is not
    tracked or avoided. In this case, the individual is just refit with it's new parameters as
    start point.
"""
function apply_genetic_operations!(indivs, ops, bank_of_terms;
    (one_node_muts!)=(
        insert_mutation!, point_mutation!, point_mutation2!, addterm_mutation!,  # ___ below ar deeper than 2
        insert_times_param_mutation!, hoist_mutation!, subtree_mutation!, drastic_simplify!
    ),
   (repeated_muts!)=(insert_mutation!, point_mutation!, point_mutation2!, hoist_mutation!),
)
    eachind = shuffle(eachindex(indivs))

    # mutations # ----------------------------------------------------------------------------------
    while !isempty(eachind)
        indiv = indivs[popfirst!(eachind)]

        until_mut     = maxim_tree_depth(indiv.node, minim = 3) > 2 ? length(ops.mutation.mut_probs) : 4
        rand_mutation = rand() * sum(ops.mutation.mut_probs[1:until_mut])
        mut_ind       = findfirst(i -> rand_mutation <= sum(ops.mutation.mut_probs[1:i]), eachindex(ops.mutation.mut_probs))

        if mut_ind == 11 && !isempty(eachind) && any(i -> maxim_tree_depth(indivs[i].node, minim=3) > 2, eachind)
            ind = findfirst(i -> maxim_tree_depth(indivs[i].node, minim=3) > 2, eachind)
            indiv2 = indivs[popat!(eachind, ind)]
            crossover_mutation!(indiv.node, indiv2.node, ops)
        elseif mut_ind == 10
            add_from_bank_of_terms_mutation!(indiv.node, ops, bank_of_terms)
        elseif false #mut_ind == 9
        elseif mut_ind == 8 && is_drastic_simplifyable(indiv.node, ops)
            drastic_simplify!(indiv.node, ops)
        elseif mut_ind == 7
            subtree_mutation!(indiv.node, ops)
        elseif mut_ind == 5
            insert_times_param_mutation!(indiv.node, ops)
        elseif mut_ind == 4
            addterm_mutation!(indiv.node, ops)
        else
            # repeated mutations -> insert, point, point2, hoist
            max_muts = max(1.0, count_nodes(indiv.node) * ops.mutation.max_muts_ratio)
            i_muts = 0
            while i_muts < max_muts
                i_muts += 1
                rand_mutation = rand() * sum(ops.mutation.multiple_mut_probs)
                mut_ind = findfirst(i -> rand_mutation <= sum(ops.mutation.multiple_mut_probs[1:i]), eachindex(ops.mutation.multiple_mut_probs))

                repeated_muts![mut_ind](indiv.node, ops)
                rand() < ops.mutation.p_multiple_mutations || break
            end
        end
    end
end

""" Helper function, which decidies whether a mutation is applied to the left or right child node.
"""
function mutate_left(node, min_depth)
    if node.ari == 1
        return true
    elseif node.ari == 2
        suff_depth_lef = maxim_tree_depth(node.lef, minim=min_depth + 2) >= min_depth
        suff_depth_rig = maxim_tree_depth(node.rig, minim=min_depth + 2) >= min_depth

        if suff_depth_lef && suff_depth_rig
            return rand(Bool)
        elseif suff_depth_lef
            return true
        else
            return false
        end
    else
        @assert false
    end
end

# ==================================================================================================
# genetic operations
# ==================================================================================================
""" Changes one node to another without modifying the structure of the node.
"""
function point_mutation!(node, ops)
    node_elect = random_node(node, mode=0)
    if node_elect.ari == 2
        node_elect.ind = rand(1:length(ops.binops))
    elseif node_elect.ari == 1
        node_elect.ind = rand(1:length(ops.unaops))
    elseif node_elect.ari == 0
        node_elect.ind = rand(1:ops.data_descript.n_vars)
    elseif node_elect.ari == -1
        node_elect.val *= rand_mult()
    end
end

""" return a positive or negative float with an absolute within the specified range.
"""
function rand_mult(;minn=2.0, maxx=10.0, sign_flip_prob=0.05)
    sign_flip = rand() < sign_flip_prob ? -1.0 : 1.0
    r = rand() * (maxx - minn) + minn
    r *= sign_flip
    return (rand(Bool) ? r : inv(r))
end

""" Changes arity 2 to 1 and vv, as well as arity 0 to -1 and vv.
"""
function point_mutation2!(node, ops)
    node_elect = random_node(node, mode=0)

    if node_elect.ari == 2
        if !isempty(ops.unaops)
            node_elect.ari = 1
            node_elect.ind = rand(1:length(ops.unaops))
            if rand(Bool)
                node_elect.lef = copy(node_elect.rig)
            end
        else
            point_mutation!(node, ops)
        end

    elseif node_elect.ari == 1
        if !isempty(ops.binops)
            node_elect.ari = 2
            node_elect.ind = rand(1:length(ops.binops))
            node_elect.rig = grow_equation(1, ops)
        else
            point_mutation!(node, ops)
        end

    elseif node_elect.ari == 0
        node_elect.ari = -1
        node_elect.val = rand()

    elseif node_elect.ari == -1
        node_elect.ari = 0
        node_elect.ind = rand(1:ops.data_descript.n_vars)
    end
end

""" Inserts or prepend an operation.
"""
function insert_mutation!(node, ops; subtree_depth=2)
    new_node = grow_equation(subtree_depth, ops, method = :full)
    lefrig1 = mutate_left(new_node, 1) ? :lef : :rig

    if rand(1:count_nodes(node)) == 1
        orig_node = copy(node)
        copy_node_wo_copy!(node, new_node)

        setfield!(node, lefrig1, orig_node)
    else
        node_elect = random_node(node, mode=1)
        lefrig2 = mutate_left(node_elect, 1) ? :lef : :rig

        setfield!(new_node, lefrig1, getfield(node_elect, lefrig2))
        setfield!(node_elect, lefrig2, new_node)
    end
end

""" Removes an operation.
"""
function hoist_mutation!(node, ops)
    node.ari <= 0 && return
    node_elect = random_node(node, mode=1)
    lefrig = mutate_left(node_elect, 1) ? :lef : :rig
    copy_node_wo_copy!(node_elect, getfield(node_elect, lefrig))
end

""" Combines two nodes to create two new ones.
"""
function crossover_mutation!(node1, node2, ops)
    node_elect1 = random_node(node1, mode=2)
    node_elect2 = random_node(node2, mode=2)

    lefrig1 = mutate_left(node_elect1, 2) ? :lef : :rig
    lefrig2 = mutate_left(node_elect2, 2) ? :lef : :rig

    temp = getfield(node_elect1, lefrig1)
    setfield!(node_elect1, lefrig1, copy(getfield(node_elect2, lefrig2)))
    setfield!(node_elect2, lefrig2, copy(temp))
end

# redundand mutations # ----------------------------------------------------------------------------
""" Replaces a subtree of a random subtree.
"""
function subtree_mutation!(node, ops; subtree_depth=0)
    subtree_depth = iszero(subtree_depth) ? rand(2:5) : subtree_depth
    node_elect = random_node(node, mode=1)
    lefrig = mutate_left(node_elect, 1) ? :lef : :rig
    setfield!(node_elect, lefrig, grow_equation(subtree_depth, ops))
end

""" Adds a *parameter to a random subnode.
"""
function insert_times_param_mutation!(node, ops)
    node_elect = random_node(node, mode=1)
    lefrig = mutate_left(node_elect, 1) ? :lef : :rig

    times_param_node = Node(2, findfirst(==(*), ops.binops))
    times_param_node.lef = Node(1.0)
    times_param_node.rig = getfield(node_elect, lefrig)

    setfield!(node_elect, lefrig, times_param_node)
end

""" Adds a top-level term.
"""
function addterm_mutation!(node, ops; subtree_depth=0)
    subtree_depth = iszero(subtree_depth) ? rand(2:5) : subtree_depth
    orig_node = copy(node)
    copy_node_wo_copy!(node, Node(2, findfirst(==(+), ops.binops)))
    node.lef = orig_node
    node.rig = grow_equation(subtree_depth, ops)
end

""" add from bank of terms
"""
function add_from_bank_of_terms_mutation!(node, ops, bank_of_terms)
    orig_node = copy(node)
    copy_node_wo_copy!(node, Node(2, findfirst(==(+), ops.binops)))
    node.lef = orig_node
    node.rig = copy(rand(bank_of_terms))
end
