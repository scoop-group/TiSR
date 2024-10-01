
module TiSR

export Options,
    Node,
    generational_loop,
    save_to_csv,
    save_to_fwf,
    save_to_excel,
    general_params,
    selection_params,
    fitting_params,
    mutation_params,
    grammar_params

using Statistics
using Random
using Dates
import Base: isapprox, show, deepcopy

using ForwardDiff
using OrderedCollections
using StatsBase: wsample, sample, corspearman
using SymbolicUtils
using SymbolicUtils.Code
using StaticArrays
using LinearAlgebra

using DataFrames
using XLSX
using UnicodePlots

include("options.jl")
include("node_n_eval_n_utilities.jl")
include("genetic_ops.jl")
include("simplify.jl")
include("selection.jl")
include("param_fitting.jl")
include("string_to_node.jl")
include("individual.jl")
include("specific_measures.jl")
include("nsga-II.jl")
include("save_results.jl")
include("ParameterEstimation/levenberg_marquardt/levenberg_marquardt.jl")
include("misc_helpers.jl")
include("grammar.jl")

end # module
