read(sdf_file)
sdf_file = "C:\\Users\\uqthulle\\Downloads\\S1\\OPERA_KOC\\TR_KOC_545.txt"
sdf_file_2 = "C:\\Users\\uqthulle\\Downloads\\S1\\OPERA_KOC\\TST_KOC_184.txt"
using DataFrames

function extract_smiles_and_logkoc(file_path::String)
    smiles_list = String[]
    logkoc_list = Float64[]

    current_smiles = nothing
    current_logkoc = nothing

    open(file_path, "r") do io
        while !eof(io)
            line = strip(readline(io))

            if startswith(line, "> <SMILES>")
                current_smiles = strip(readline(io))
            elseif startswith(line, "> <LogKOC>")
                current_logkoc = parse(Float64, strip(readline(io)))
            elseif startswith(line, "\$\$\$\$")
                if current_smiles !== nothing && current_logkoc !== nothing
                    push!(smiles_list, current_smiles)
                    push!(logkoc_list, current_logkoc)
                end
                # Reset for next molecule
                current_smiles = nothing
                current_logkoc = nothing
            end
        end
    end

    return DataFrame(SMILES = smiles_list, LogKOC = logkoc_list)
end

# Run it
df = extract_smiles_and_logkoc(sdf_file)
df_2 = extract_smiles_and_logkoc(sdf_file_2)
df_final = vcat(df, df_2)
first(df, 5)

CSV.write("C:\\Users\\uqthulle\\Downloads\\OPERA experimental KOC used in model.csv", df_final)
SMILES_opera = df_final.SMILES


mobility_opera_exp = smiles_to_mobility("C:\\Users\\uqthulle\\Downloads", SMILES_opera )
using PyCall

@pyimport sklearn
println("scikit-learn version: ", sklearn.__version__)

RepoRT_FP = CSV.read("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Finerprints\\All RepoRT fingerprints.csv", DataFrame)
REACH_FP = CSV.read("C:\\Users\\uqthulle\\Downloads\\REACH fingerprints final.csv", DataFrame)
REACH_FP = select!(REACH_FP, [1, 4:ncol(REACH_FP)...])
rename!(REACH_FP, names(REACH_FP)[1] => :x1)
All_FP_data = vcat(RepoRT_FP, REACH_FP)

All_FPs = Matrix(All_FP_data[:,2:end])

rf_cl = joblib.load("optimized_random_forest_classifier_RepoRT_improved pH3 and 2.6.joblib")

predicted_class = ScikitLearn.predict(rf_cl,All_FPs)

predicted_probability = Int(round.(maximum(ScikitLearn.predict_proba(rf_cl,All_FPs))*100, digits = 0))

Precompiled_SMILES = DataFrame(SMILES = All_FP_data.x1, Predicted_class = predicted_class,
Probability = predicted_probability)

CSV.write("Precompiled SMILES new.csv", Precompiled_SMILES)

