
using ForwardDiff
using OrderedCollections
using Statistics
using StatsBase
using StringDistances
using Test

import Pkg
Pkg.develop(path=".")

using TiSR

include("test_node_n_eval_utilities.jl")
include("test_string_to_node.jl")
include("test_param_fitting.jl")
include("test_selection.jl")
include("test_genetic_ops.jl")
include("test_grammar.jl")
include("test_individual.jl")
include("test_nsga-II.jl")


