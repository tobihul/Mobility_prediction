module Mobility_prediction

__precompile__()

using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots
using LinearAlgebra, ScikitLearn, Random, MLJ, PyCall, Conda

# Declare module-level variables
const skl = Ref{PyObject}()
const jl = Ref{PyObject}()
const pd = Ref{PyObject}()
const rf_cl = Ref{Any}()
const Precompiled_smiles = Ref{DataFrame}()

function __init__()
    skl[] = pyimport("sklearn.ensemble")
    jl[] = pyimport("joblib")
    pd[] = pyimport("padelpy")

    # Initialize the random forest classifier and precompiled SMILES
    rf_cl[] = jl[].load("optimized_random_forest_classifier_RepoRT.joblib")
    Precompiled_smiles[] = CSV.read("Precompiled SMILES.csv", DataFrame)
end

# Include your functions file
include("All functions.jl")


export smiles_to_mobility

end # module Mobility_prediction
