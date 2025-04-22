__precompile__(true)
module Mobility_prediction

using CSV, Statistics, DataFrames, LinearAlgebra, ScikitLearn, Random
using PyCall, Conda


const jl = Ref{PyObject}()
const pd = Ref{PyObject}()
const rf_cl = Ref{Any}()
const Precompiled_smiles = Ref{DataFrame}()

using Conda

ENV["PYTHON"] = ""  

Conda.pip_interop(true, Conda.ROOTENV)  

Conda.pip("install", ["joblib"], Conda.ROOTENV)

Conda.pip("install", ["padelpy"], Conda.ROOTENV)

try
    run(`java -version`)  # Check if Java is already installed
catch
    println("Java not found. Installing OpenJDK...")
    Conda.add("openjdk")
end


function __init__()

    try
        readdir("$(Conda.ROOTENV)\\pkgs\\scikit-learn-1.6.1-py310hf2a6c47_0") 
    catch
        println("Scikitlearn version not found installing 1.6.1...")
        Conda.add("scikit-learn=1.6.1")
    end
    
   
   
    src_dir = dirname(pathof(Mobility_prediction))
    
    
    pkg_dir = dirname(src_dir)

    
    rf_path = joinpath(pkg_dir, "optimized_random_forest_classifier_RepoRT_improved pH3 and 2.6.joblib")
    csv_path = joinpath(pkg_dir, "Precompiled SMILES new.csv")

    
    println("Random Forest path: $rf_path")
    println("CSV path: $csv_path")

   
    jl[] = pyimport("joblib")
    pd[] = pyimport("padelpy")

    # Load the model and data
    rf_cl[] = jl[].load(rf_path)
    Precompiled_smiles[] = CSV.read(csv_path, DataFrame)
end


include("smiles_to_mobility.jl")


export smiles_to_mobility

end 
