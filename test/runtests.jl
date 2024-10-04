
using UnicodePlots
using ForwardDiff
using OrderedCollections
using Statistics
using StatsBase
using StringDistances
using Test
using Random

import Pkg
Pkg.develop(path=".")

using TiSR

include("test_node_n_eval_n_utilities.jl")
include("test_individual.jl")
include("test_string_to_node.jl")
include("test_selection.jl")
include("test_genetic_ops.jl")
include("test_grammar.jl")
include("test_nsga-II.jl")
include("test_param_fitting.jl")
include("test_simplify.jl")

# include("test_misc_helpers.jl")
# include("test_save_results.jl")
# include("test_specific_measures.jl")
# include("test_ParameterEstimation")
# include("test_options.jl")
