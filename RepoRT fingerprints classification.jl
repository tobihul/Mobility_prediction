using CSV, Statistics, DataFrames, StatsPlots, MultivariateStats
using GLM, TypedTables, LinearAlgebra
folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], MACCS = String[], Pubchem_fps = String[],MW = Float64[], XlogP = Float64[], Retention_factors = Float64[], Modifier = Float64[] )
function interpolate_B_modifier(time::Float64, gradient::DataFrame)
    times = gradient[:,1]./run_time
    idx = searchsortedlast(times, time)
    if idx == 0
        return gradient[1,3]
    elseif idx == length(times)
        return gradient[end,3]
    elseif time > 1
        return gradient[end,3]
    else
        t1, t2 = times[idx], times[idx+1]
        B1, B2 = gradient[:,3][idx], gradient[:,3][idx+1]
        return B1 + (B2 - B1) * (time - t1) / (t2 - t1)
    end
end
for i = 1:374
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data", readdir(folder_path)[i])
    println(i)
    #Getting the SMILES and retention times
    rt_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "rtdata_canonical_success.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    Compound_Inchi = data_compounds[:,6]
    RT = data_compounds[:,4]

     #Getting the retention factors
     gradient_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "gradient.tsv"],"_"))
     gradient = CSV.read(gradient_path, DataFrame)
     local retention_factors
     Modifier = []
     try 
         run_time = gradient[end,1]
         retention_factors = RT./run_time
         for rtfs in retention_factors
            push!(Modifier, interpolate_B_modifier(rtfs, gradient))
         end
     catch
         retention_factors = repeat(["Not gradient info"], length(RT))
         Modifier = repeat(["Not gradient info"], length(RT))
     end

    #Getting the MACCS
    MACCS_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "fingerprints_maccs_canonical_success.tsv"],"_"))
    MACCS_info = CSV.read(MACCS_path, DataFrame)
    df_merged = innerjoin(data_compounds, MACCS_info, on = :id)
    MACCS = df_merged[:,end]
    #Getting the PubChemFingerprints
    Pubchem_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "fingerprints_pubchem_canonical_success.tsv"],"_"))
    Pubchem_info = CSV.read(Pubchem_path, DataFrame)
    df_merged_2 = innerjoin(data_compounds, Pubchem_info, on = :id)
    Pubchem_fps = df_merged_2[:,end]

    #Getting the LC mode info for the data set
    info_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    #Check if the column info is missing in the meta data
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"

        LC_mode = repeat(["No mode info"], length(MACCS))

    else

        LC_mode = repeat([LC_mode], length(MACCS))

    end

    #Getting MW and XLogP
    descriptor_path = joinpath(current_file, join(["$(readdir(folder_path)[1:end-1][i])", "descriptors_canonical_success.tsv"],"_"))
    data_descriptors = CSV.read(descriptor_path, DataFrame)
    MW = data_descriptors.MW
    XLogP = data_descriptors.XLogP

    
  
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(Inchikey = Compound_Inchi, LC_mode = LC_mode, MACCS = MACCS, Pubchem_fps = Pubchem_fps,MW = MW, XlogP = XLogP, Retention_factors = retention_factors, Modifier = Modifier)
    Final_table = vcat(Final_table, current_table)
end

# Separate entries with RP and HILIC modes
rp_Inchikey = Final_table[Final_table.LC_mode .== "RP", :].Inchikey
hilic_Inchikey = Final_table[Final_table.LC_mode .== "HILIC", :].Inchikey

# Find the intersection of CIDss for both modes
common_Inchikey = intersect(rp_Inchikey, hilic_Inchikey)

#Removing the ones that have been separated by both as we already have them
condition = in.(Final_table.Inchikey, Ref(common_Inchikey))
Final_table_unique = Final_table[.!condition, :]

Final_table_unique = dropmissing(Final_table_unique, :Retention_factors)
valuess = Final_table_unique.Retention_factors
indexes = Int[]
for i = 1:length(valuess)
    if typeof(valuess[i]) == Float64
        push!(indexes, i)
    end
end
indexes
Final_table_unique = Final_table_unique[indexes,:]

MACCS = Int.(zeros(length(Final_table_unique[:,2]),166))
for i = 1:length(Final_table_unique[:,3])
    string_bits = Final_table_unique[i,3]
    vectorized_bits = parse.(Int, split(string_bits, ","))
    for v in vectorized_bits
        MACCS[i,v] = 1
    end
end
MACCS
Pubchem_fps = Int.(zeros(length(Final_table_unique[:,4]),881))
for i = 1:length(Final_table_unique[:,4])
    string_bits = Final_table_unique[i,4]
    vectorized_bits = parse.(Int, split(string_bits, ","))
    for v in vectorized_bits
        Pubchem_fps[i,v] = 1
    end
end
Pubchem_fps

indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table_unique))
indices_HILIC = findall(row -> occursin("HILIC", row.LC_mode), eachrow(Final_table_unique))

X = [MACCS Pubchem_fps]
y = Final_table_unique[indices_RPLC,6]
X = X[indices_RPLC,:]


histogram(Final_table_unique[indices_HILIC,7], bins = 50, title = "Classification of HILIC retention indices",
dpi = 300, label = false)
separator = 0.18
vspan!([0, separator], color=:green, alpha=0.3, label = "Non-mobile")
vspan!([separator, separator*2], color=:yellow3, alpha=0.3, label = "Not very mobile")
vspan!([separator*2, separator*3], color=:darkorange, alpha=0.3, label = "Mobile")
vspan!([separator*3, maximum(Final_table_unique[indices_HILIC,7])], color=:red, alpha=0.3, label = "Very mobile")


savefig("C:\\Users\\uqthulle\\Documents\\Plots\\RPLC space with organic modifier.png")


histogram(Final_table_unique[indices_RPLC,7], title = "Classification of RPLC retention indices",
dpi = 300, label = false)
separator = 0.2
vspan!([0, separator], color=:red, alpha=0.3, label = "Very mobile")
vspan!([separator, separator*2], color=:darkorange, alpha=0.3, label = "Mobile")
vspan!([separator*2, separator*3], color=:yellow3, alpha=0.3, label = "Not very mobile")
vspan!([separator*3, maximum(Final_table_unique[indices_RPLC,7])], color=:green, alpha=0.3, label = "Non-mobile")


Table_RPLC = Final_table_unique[indices_RPLC, :]
label_RPLC = String[]  # Initialize an empty array to store strings
cutoff_RPLC = 0.2
for i = 1:length(Table_RPLC[:, 7])
    if Table_RPLC[i, 7] <= cutoff_RPLC
        push!(label_RPLC, "Very mobile")
    elseif Table_RPLC[i, 7] <= cutoff_RPLC * 2
        push!(label_RPLC, "Mobile")
    elseif Table_RPLC[i, 7] <= cutoff_RPLC * 3
        push!(label_RPLC, "Not very mobile")
    else
        push!(label_RPLC, "Non mobile")
    end
end
label_RPLC

X = [MACCS Pubchem_fps]
X = X[indices_RPLC]
y = label_RPLC

using ScikitLearn
using ScikitLearn.CrossValidation: cross_val_score
using ScikitLearn: @sk_import
@sk_import ensemble: RandomForestClassifier
@sk_import model_selection: train_test_split
using ScikitLearn.GridSearch: GridSearchCV

train, valid, test = partition(eachindex(y), 0.7, 0.2, rng=42)

rf_regressor = RandomForestClassifier()

@time fit!(rf_regressor, X_train, y_train)
    
y_pred_train = ScikitLearn.predict(rf_regressor, X_train)
y_pred_train_score = ScikitLearn.score(rf_regressor, X_train, y_train)

scatter(y_train, y_pred_train, xlabel = "test", ylabel = "pred", title = "RF model to predict retention factors from fingerprints RPLC", titlefont = font(10),
label = "train score = $(round(y_pred_train_score,digits = 2))")



######Plotting Modifier data
histogram(Final_table_unique[indices_RPLC,end-1], label = false)
histogram(Final_table_unique[indices_RPLC,end], yticks = false, xlims= (0,100))
scatter(Final_table_unique[indices_RPLC,end],Final_table_unique[indices_RPLC,end-1], dpi = 300,
title = "RPLC Modifier vs Retention factor", xlabel = "%B", ylabel = "Retention factor",
label = false)
scatter(Final_table_unique[indices_HILIC,end],Final_table_unique[indices_HILIC,end-1], dpi = 300,
title = "HLILC Modifier vs Retention factor", xlabel = "%B", ylabel = "Retention factor",
label = false)

scatter(Final_table_unique.MW[indices_HILIC], 
Final_table_unique.XlogP[indices_HILIC], 
zcolor = Final_table_unique.Modifier[indices_HILIC], c = :plasma, clims = (0,100)
,title = "HILIC Space with %B Modifier", dpi = 1000, xlabel = "MW", 
ylabel = "XlogP", cbar_title = "%Organic modifier", label = false,
markersize = 2, 
markerstrokewidth  = 0)
savefig("C:\\Users\\uqthulle\\Documents\\Plots\\2d histograpm rplc.png")



indices_0_2 = findall(x-> x <0.2, Final_table_unique[:,end-1])
histogram(Final_table_unique[:,end-1])
histogram(Final_table_unique[:,end]./100)
vline!([0.333, 0.666,1])
vline!([0.333,0.666,1])

histogram2d((Final_table_unique[indices_RPLC,end]./100), 
Final_table_unique[indices_RPLC,end-1],
ylabel = "Retention factor",
xlabel = "%B",
title = "2d histogram RF and organic modifier",
dpi = 300,
xlims = (0,1),
ylims = (0,1))

b_separator = 0.33333333
rf_separator = 0.33333333

plot!(rectangle(b_separator,rf_separator,0,0), opacity = 0.4, c = :red, label = "Very mobile")
plot!(rectangle(b_separator,rf_separator,b_separator,0.25), opacity = 0.4, c = :darkorange, label = "Moderately mobile")
plot!(rectangle(b_separator,rf_separator,b_separator*2,0.4), opacity = 0.4, c = :green, label = "Not mobile")

    
scatter((Final_table_unique[indices_RPLC,end]./100), Final_table_unique[indices_RPLC,6], zcolor =Final_table_unique[indices_RPLC,end-1] )


findall(Final_table_unique[indices_RPLC,end])

condition_1(x) = x < 12.5
condition_2(x) = x > 0.65

# Find indices where both conditions are satisfied
indices = findall(x -> condition_1(Final_table_unique[indices_RPLC,end][x]) && condition_2(Final_table_unique[indices_RPLC,end-1][x]), 1:min(length(Final_table_unique[indices_RPLC,end]), length(Final_table_unique[indices_RPLC,end-1])))

scatter(Final_table_unique[indices_RPLC,end][indices],Final_table_unique[indices_RPLC,end-1][indices], dpi = 300,
title = "RPLC Modifier vs Retention factor", xlabel = "%B", ylabel = "Retention factor",
label = false, xlims = (0,100), ylims = (0,1))

println(Final_table_unique[indices_RPLC,1][indices][10])
println(Final_table_unique[indices_RPLC,5][indices][10])
println(Final_table_unique[indices_RPLC,6][indices][10])

outlier_compounds = Final_table_unique[indices_RPLC,:][indices,:]

CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\Outlier compounds.csv", outlier_compounds)

indices_21 = findall(x-> x=="InChI=1S/C16H14O6/c1-21-13-3-2-8(4-10(13)18)14-7-12(20)16-11(19)5-9(17)6-15(16)22-14/h2-6,14,17-19H,7H2,1H3", Final_table_unique[indices_RPLC,1])

Final_table_unique[indices_RPLC,end-1:end][indices_21,:]




