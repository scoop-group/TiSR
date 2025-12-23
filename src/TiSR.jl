
module TiSR

export Options,
    Node,
    generational_loop,
    save_to_csv,
    save_to_fwf,
    save_to_excel,
    data_split_params,
    general_params,
    measure_params,
    selection_params,
    fitting_params,
    mutation_params,
    grammar_params,
    obtain_value_constr_func,
    obtain_monotonicity_constr_func

using Statistics
using Random
using Dates
import Base: isapprox, ==, isless, show, copy, deepcopy, iterate, length, hash

using ForwardDiff
using OrderedCollections
using StatsBase: wsample, sample, corspearman
using SymbolicUtils
using SymbolicUtils.Code
using StaticArrays
using LinearAlgebra
using Optim
using LineSearches
# import DynamicExpressions
# import Zygote

#using TimerOutputs # @timeit
#const to = TimerOutput() # @timeit

using DataFrames
using XLSX
using YAML
using UnicodePlots
using Printf

using DiffResults

include("options.jl")
include("node_n_eval_n_utilities.jl")
include("genetic_ops.jl")
include("simplify.jl")
include("selection.jl")
include("param_fitting.jl")
include("string_to_node.jl")
include("individual.jl")
include("measures.jl")
include("main_loop.jl")
include("save_results.jl")
include("ParameterEstimation/levenberg_marquardt.jl")
include("grammar.jl")

# for shape constraints
include("forwarddiff_vector_extension.jl")
include("obtain_constraint_funcs.jl")
include("ParameterEstimation/penalty_method.jl")

#global reject_rate = [0, 0] # DEBUG expression_log

end # module
