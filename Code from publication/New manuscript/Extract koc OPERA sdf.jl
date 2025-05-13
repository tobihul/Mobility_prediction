read(sdf_file)
sdf_file = "C:\\Users\\uqthulle\\Downloads\\S1\\OPERA_KOC\\TR_KOC_545.txt"
sdf_file_2 = "C:\\Users\\uqthulle\\Downloads\\S1\\OPERA_KOC\\TST_KOC_184.txt"

sdf_file = "R:\\PHD2024TH-Q6813\\Research files\\Models and other documents\\TR_KOC_545.txt"
sdf_file_2 = "R:\\PHD2024TH-Q6813\\Research files\\Models and other documents\\TST_KOC_184.txt"
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

RepoRT_data = CSV.read("R:\\PHD2024TH-Q6813\\Research files\\Code\\Final Mobility model\\RPLC_data_pH3and2.6_updated.csv", DataFrame)
RepoRT_FP = CSV.read("R:\\PHD2024TH-Q6813\\Research files\\Code\\Final Mobility model\\pH 3 and 2.6 RepoRT fingerprints final.csv", DataFrame)
REACH_FP = CSV.read("R:\\PHD2024TH-Q6813\\Research files\\Code\\Final Mobility model\\REACH fingerprints final.csv", DataFrame)

SMILES_RepoRT = RepoRT_data.SMILES
SMILES_REACH = REACH_FP.canonicalsmiles
REACH_FP = REACH_FP[:,4:end]
#REACH_FP = select!(REACH_FP, [1, 4:ncol(REACH_FP)...])
#rename!(REACH_FP, names(REACH_FP)[1] => :x1)
All_FP_data = vcat(RepoRT_FP, REACH_FP)
All_FPs = Matrix(All_FP_data)

All_smiles = vcat(SMILES_RepoRT, SMILES_REACH)

rf_cl = joblib.load("optimized_random_forest_classifier_RepoRT_improved pH3 and 2.6 final.joblib")

predicted_class = ScikitLearn.predict(rf_cl,All_FPs)

probas = ScikitLearn.predict_proba(rf_cl,All_FPs)

predicted_probability = vec(maximum(probas, dims=2))

Precompiled_SMILES = DataFrame(SMILES = All_smiles, Predicted_class = predicted_class,
Probability = predicted_probability)

CSV.write("Precompiled SMILES new.csv", Precompiled_SMILES)

RepoRT_data = CSV.read("R:\\PHD2024TH-Q6813\\Research files\\Files for work PC\\Other csv and stuff\\RPLC_data_pH3and2.6_updated.csv", DataFrame)

overlapping = intersect(RepoRT_data.SMILES, df_final.SMILES)

indices_RepoRT = [findfirst(x -> x == smiles, RepoRT_data.SMILES) for smiles in overlapping]
indices_OPERA = [findfirst(x -> x == smiles, df_final.SMILES) for smiles in overlapping]

modifiers = RepoRT_data.Modifier[indices_RepoRT]
kocs = (df_final[:,2][indices_OPERA])

using GLM

scatter(kocs, modifers)

df = DataFrame(kocs = kocs, modifers = modifiers./100)

# Fit the linear model: modifers ~ kocs
model = lm(@formula(modifers ~ kocs), df)

# Get the predicted values
y_pred = GLM.predict(model, df)

# Calculate R^2
R2 = cor(kocs, modifiers./100)

# Plotting the scatter plot and the fitted line
scatter(kocs, modifiers./100, label=false, xlabel="log Koc", ylabel=" Fraction of organic Modifier (Φ)", dpi = 300,
 legend=:bottomright,
 ylims = (0,1))
plot!(kocs, y_pred, label="r = $(round(R2, digits=2))", color=:red)

savefig("R:\\PHD2024TH-Q6813\\Research files\\Plots\\Modifer koc relationship.png")



