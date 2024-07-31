__precompile__()
module Mobility_prediction

using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots
using LinearAlgebra, ScikitLearn, Random, MLJ, PyCall, Conda

function __init__()
    skl = pyimport("sklearn.ensemble")
    jl = pyimport("joblib")
    pd = pyimport("padelpy")
    global skl, jl, pd
end

include("All functions.jl")


export smiles_to_mobility

end # module Mobility_prediction
