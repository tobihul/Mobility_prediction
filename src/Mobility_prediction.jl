__precompile__()
module Mobility_prediction

using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots
using LinearAlgebra, ScikitLearn, Random, MLJ, PyCall, Conda

const skl = Ref{PyObject}()
const jl = Ref{PyObject}()
const pd = Ref{PyObject}()

function __init__()
    skl[] = pyimport("sklearn.ensemble")
    jl[] = pyimport("joblib")
    pd[] = pyimport("padelpy")
end

include("All functions.jl")


export smiles_to_mobility

end # module Mobility_prediction
