using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats
using GLM, TypedTables, LinearAlgebra, ScikitLearn, Random, MLJ, MLDataUtils
using ScikitLearn: @sk_import
@sk_import ensemble: RandomForestRegressor
@sk_import model_selection: StratifiedKFold
@sk_import model_selection: train_test_split
cd("R:\\PHD2024TH-Q6813\\Code\\Regression")
folder_path = "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], MACCS = String[], Pubchem_fps = String[],MW = Float64[], XlogP = Float64[], Retention_factors = Float64[], Modifier = Float64[] )
MACCS_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\MACCS keys.txt")
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")
All_keys =  [MACCS_keys ;PubChem_keys]
All_keys_y = [MACCS_keys ;PubChem_keys]
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
function calculate_leverage(X,x)
    leverage = zeros(length(x[:,1]))
    b = pinv(X'*X)
    for i = 1:length(x[:,1])
        leverage[i] = (x[i,:])' * b * x[i,:] 
    end
    return leverage
end
function remove_overlapping(Results_table)

    # Separate entries with RP and HILIC modes
    rp_Inchikey = Results_table[Results_table.LC_mode .== "RP", :].Inchikey
    hilic_Inchikey = Results_table[Results_table.LC_mode .== "HILIC", :].Inchikey

    # Find the intersection of CIDss for both modes
    common_Inchikey = intersect(rp_Inchikey, hilic_Inchikey)

    #Removing the ones that have been separated by both as we already have them
    condition = in.(Results_table.Inchikey, Ref(common_Inchikey))
    Final_table_unique = Results_table[.!condition, :]

    return Final_table_unique
end
function remove_missing_and_outliers(Results_unique)
    #This removes rows that have missing retention factors due to no gradient data or other 
    valuess = Results_unique.Retention_factors
    indexes = Int[]
    for i = 1:length(valuess)
        if !ismissing(valuess[i]) && typeof(valuess[i]) == Float64
            push!(indexes, i)
        end
    end
    indexes
    Results_unique = Results_unique[indexes,:]

    #This removes outliers with high retention factor but low %B
    rtfs = Results_unique[:,end-1]
    pmbph = Results_unique[:,end]
    indexes = Int[]
    for i = 1:length(rtfs)
        if rtfs[i] > 0.65 && pmbph[i] < 12.5
            continue
        else
            push!(indexes, i)
        end
    end
    indexes
    Results_unique = Results_unique[indexes,:]
    return Results_unique
end
function format_fingerprints(Final_table_unique)
    #First we get MACCS
    MACCS = Int.(zeros(length(Final_table_unique[:,2]),166))
    for i = 1:length(Final_table_unique[:,3])
        string_bits = Final_table_unique[i,3]
        vectorized_bits = parse.(Int, split(string_bits, ","))
        for v in vectorized_bits
            MACCS[i,v] = 1
        end
    end
    #Now PubChem FPs
    Pubchem_fps = Int.(zeros(length(Final_table_unique[:,4]),881))
    for i = 1:length(Final_table_unique[:,4])
        string_bits = Final_table_unique[i,4]
        vectorized_bits = parse.(Int, split(string_bits, ","))
        for v in vectorized_bits
            Pubchem_fps[i,v] = 1
        end
    end
    return MACCS, Pubchem_fps
end
function custom_binning(y, num_bins)
    min_val = minimum(y)
    max_val = maximum(y)
    bin_width = (max_val - min_val) / num_bins
    bins = Array{Int}(undef, length(y))
    
    for i in 1:length(y)
        bin_index = ceil(Int, (y[i] - min_val) / bin_width)
        bins[i] = min(bin_index, num_bins)
    end
    
    return bins
end
function Stratifiedcv(X, groups, n_folds)
    # Initialize StratifiedKFold
    skf = StratifiedKFold(n_splits=n_folds)

    # Initialize variables to store indices for each fold
    train_indices = Vector{Vector{Int}}(undef, n_folds)
    test_indices = Vector{Vector{Int}}(undef, n_folds)

    # Iterate over the splits and store indices for each fold
    for (i, (train_idx, test_idx)) in enumerate(skf.split(X, groups))
        train_indices[i] = train_idx .+ 1  # Adjust for 1-based indexing
        test_indices[i] = test_idx .+ 1  # Adjust for 1-based indexing
    end
    return train_indices, test_indices
end
for i = 1:374
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\processed_data", readdir(folder_path)[i])
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

Final_table

#Removing compounds that occur both in RPLC and HILIC 
Final_table_unique = remove_overlapping(Final_table)

#Removing rows with missing retention factor and outliers
Final_table_unique = remove_missing_and_outliers(Final_table_unique)

#Formatting MACCS and PubChem Fingerprints for modelling
MACCS, PubChem_fps = format_fingerprints(Final_table_unique)

#Getting indices for RPLC and HILIC rows
indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table_unique))
indices_HILIC = findall(row -> occursin("HILIC", row.LC_mode), eachrow(Final_table_unique))



histogram(Final_table_unique[indices_HILIC,7], bins = 50, title = "Classification of HILIC retention indices",
dpi = 300, label = false)
separator = 0.18
vspan!([0, separator], color=:green, alpha=0.3, label = "Non-mobile")
vspan!([separator, separator*2], color=:yellow3, alpha=0.3, label = "Not very mobile")
vspan!([separator*2, separator*3], color=:darkorange, alpha=0.3, label = "Mobile")
vspan!([separator*3, maximum(Final_table_unique[indices_HILIC,7])], color=:red, alpha=0.3, label = "Very mobile")


savefig("C:\\Users\\uqthulle\\Documents\\Plots\\RPLC space with organic modifier.png")

Final_table_unique
histogram(filtered_RPLC_data[:,end]./100, title = "Classification of RPLC %B",
dpi = 300, label = false, xlims = (0,1))
separator = 0.333
vspan!([0, 0.2], color=:red, alpha=0.3, label = "Very mobile")
vspan!([0.2, 0.6], color=:darkorange, alpha=0.3, label = "Mobile")
vspan!([0.6, 1], color=:green, alpha=0.3, label = "Non-mobile")


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
histogram(Final_table_unique[:,end]./100, label = false, dpi = 300,
xlabel = "%B")
vline!([0.333, 0.666,1], label = false)
vline!([0.333,0.666,1], label = false)

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
using Plots
plot!(rectangle(b_separator,rf_separator,0,0), opacity = 0.4, c = :red, label = "Very mobile")
plot!(rectangle(b_separator,rf_separator,b_separator,0.25), opacity = 0.4, c = :darkorange, label = "Moderately mobile")
plot!(rectangle(b_separator,rf_separator,b_separator*2,0.4), opacity = 0.4, c = :green, label = "Not mobile")

    
scatter((Final_table_unique[indices_RPLC,end]./100), Final_table_unique[indices_RPLC,6], zcolor =Final_table_unique[indices_RPLC,end-1] )




condition_1(x) = x > 80
condition_2(x) = x > 0.8

# Find indices where both conditions are satisfied
indices = findall(x -> condition_1(Final_table_unique[indices_RPLC,end][x]) && condition_2(Final_table_unique[indices_RPLC,end-1][x]), 1:min(length(Final_table_unique[indices_RPLC,end]), length(Final_table_unique[indices_RPLC,end-1])))

scatter(Final_table_unique[indices_RPLC,end][indices],Final_table_unique[indices_RPLC,end-1][indices], dpi = 300,
title = "RPLC Modifier vs Retention factor", xlabel = "%B", ylabel = "Retention factor",
label = false, xlims = (0,100), ylims = (0,1))

println(Final_table_unique[indices_RPLC,1][indices][19])
println(Final_table_unique[indices_RPLC,5][indices][19])
println(Final_table_unique[indices_RPLC,6][indices][19])

outlier_compounds = Final_table_unique[indices_RPLC,:][indices,:]

CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\Outlier compounds.csv", outlier_compounds)

indices_21 = findall(x-> x=="InChI=1S/C16H14O6/c1-21-13-3-2-8(4-10(13)18)14-7-12(20)16-11(19)5-9(17)6-15(16)22-14/h2-6,14,17-19H,7H2,1H3", Final_table_unique[indices_RPLC,1])

