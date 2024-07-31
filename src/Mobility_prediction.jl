__precompile__()
module Mobility_prediction



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

    pkg_dir = Base.dirname(@__DIR__)
    rf_path = joinpath(pkg_dir, "optimized_random_forest_classifier_RepoRT.joblib")
    csv_path = joinpath(pkg_dir,"Precompiled_SMILES.csv")

    # Load the model and data
    rf_cl[] = jl[].load(rf_path)
    Precompiled_smiles[] = CSV.read(csv_path, DataFrame)
end

# Include your functions file
include("All functions.jl")


export smiles_to_mobility

end # module Mobility_prediction
