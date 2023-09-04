
# SymbolicRegression.jl

There are several parts of this code that are inspired by the SymbolicRegression.jl package. Most
importantly, it was used as a starting point on how the basis, i.e., the Node that form the
equations and how they are evaluated. However, most of the code was completely rewritten without
reference to SymbolicRegression.jl, as I was curious to build it up myself and understand the
process completely, after I had understood the Node and eval_equation parts. The version of
SymbolicRegression.jl in questions is 51b7da4 of May 1, 2022. 

- Node -> Node
    - approximately same, but other names and used more values in the "degree" to save the
      "constant" file

- copyNode -> copy_node
    - build for different node structure

- nodeToSymbolic -> node_to_symbolic
    - adopted the way to create SymbolicUtils symbols without the @syms macro

- randomNode -> random_node
    - added a parameter and the feature that ensures that the funciton returns a node with a minimum
      maximal depth

- prependRandomOp
    - added this type of mutation, but added it to the insert_mutation!

- Options -> Options
    - implemented a very similar struct for to carry all the settings but broke it up into smaller
      pieces

- countNodes -> count_nodes
    - essentially the same

- evalTreeArray -> eval_equation
    - made a very naive, simplified version, which is less performant, but easier to read and reason
      through

- return_on_nonfinite -> bad_array
    - converted the macro to a inline function

- PopMember -> Individual
    - added several more fit-measures at the cost of performance
    - more of the individual-specific calculations inside the individual creation

- combineOperators -> simplify_binary_across_1_level!
    - adopted the idea to make this simplification, but implemented it differently

- optFunc & optimizeConstants
    - general structure and setup of optimizing the parameters inside the node

# DynamicExpressions.jl

I adopted the same typing structure of the Node struct. Version 4638851, August 28, 2023.

