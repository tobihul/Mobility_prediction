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

function setup_environment()
    ENV["PYTHON"] = ""  # Ensure PyCall uses Conda's Python
    Conda.pip_interop(true, Conda.ROOTENV)  # Enable pip interop in Conda root environment

    # List of required Python packages
    required_packages = ["joblib", "padelpy"]

    # Install each package if not already installed
    for pkg in required_packages
        try
            run(`$(Conda.ROOTENV)/bin/pip show $pkg`)  # Check if package is already installed
        catch
            Conda.pip("install", [pkg], Conda.ROOTENV)  # Install package if not found
        end
    end

    Pkg.build("PyCall")  # Rebuild PyCall to use the newly installed packages
end

function __init__()

    # Activate Conda environment
    Conda.add("scikit-learn=1.5.1")
   
   # Get the directory of the package source file
    src_dir = dirname(pathof(Mobility_prediction))
    
    # Go up one level to get the root directory of the package
    pkg_dir = dirname(src_dir)

    # Construct paths to files in the root directory
    rf_path = joinpath(pkg_dir, "optimized_random_forest_classifier_RepoRT.joblib")
    csv_path = joinpath(pkg_dir, "Precompiled SMILES.csv")

    # Print paths for debugging
    println("Random Forest path: $rf_path")
    println("CSV path: $csv_path")

    # Import necessary modules
    skl[] = pyimport("sklearn.ensemble")
    jl[] = pyimport("joblib")
    pd[] = pyimport("padelpy")

    # Load the model and data
    rf_cl[] = jl[].load(rf_path)
    Precompiled_smiles[] = CSV.read(csv_path, DataFrame)
end

# Include your functions file
include("smiles_to_mobility.jl")


export smiles_to_mobility

end # module Mobility_prediction
