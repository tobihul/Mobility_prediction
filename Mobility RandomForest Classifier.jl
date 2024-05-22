using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats
using GLM, TypedTables, LinearAlgebra, ScikitLearn, Random, MLJ, MLDataUtils
using ScikitLearn: @sk_import
@sk_import ensemble: RandomForestClassifier
@sk_import model_selection: StratifiedKFold
folder_path = "C:\\Users\\uqthulle\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], MACCS = String[], Pubchem_fps = String[],MW = Float64[], XlogP = Float64[], Retention_factors = Float64[], Modifier = Float64[] )
MACCS_keys = readlines( "C:\\Users\\uqthulle\\Documents\\MACCS keys.txt")
PubChem_keys = readlines( "C:\\Users\\uqthulle\\Documents\\PubChem keys.txt")
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
        leverage[i] = pinv(x[i,:]) * b * x[i,:] 
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

"""RandomForestClassifier based on mobility"""

##Assigning the mobility classes First
Retention_factors = Final_table_unique[indices_RPLC,end-1]
p_b = Final_table_unique[indices_RPLC,end]./100

#Shuffling data and creating classes
Random.seed!(42)

X = [MACCS PubChem_fps]

X = (X[indices_RPLC,:])

y = []
for i in eachindex(Retention_factors)
    if Retention_factors[i] >= 0.666666 
        push!(y, "Non-mobile")
    elseif Retention_factors[i] <= 0.3333333
        push!(y, "Very mobile")
    else 
        push!(y, "Mobile")
    end
end

shuffle_indices = shuffle(collect(1:length(y)))

X = X[shuffle_indices,:]

y = y[shuffle_indices]

model,test_final = partition(eachindex(y), 0.9, rng=42)

#Manual Hyperparameter optimisation for ScikitLearn RF
n_estimators_t = 100

criterion = "gini"

min_samples_split_t = [2]

min_samples_leaf_t = [1]

random_state = 42

f_name = "Random Forest Classification RPLC"

df_results_3 = DataFrame(n_estimators = Int[], min_samples_split = Int[], min_samples_leaf = Int[], train_score = Float64[], test_score = Float64[])

df_results_3 = CSV.read("C:\\Users\\uqthulle\\Documents\\$f_name.csv ", DataFrame)

rf_cl = RandomForestClassifier(n_estimators = 10, random_state = 42)

ScikitLearn.fit!(rf_cl, X[model,:], y[model])

y_hat_train = ScikitLearn.predict(rf_cl,X[model,:])
y_hat_test = ScikitLearn.predict(rf_cl,X[test_final,:])

score_train = ScikitLearn.score(rf_cl, X[model,:], y[model])
score_test = ScikitLearn.score(rf_cl, X[test_final,:], y[test_final])

scatter(y[model], y_hat_train, label = "train = $(round(score_train, digits= 2))", xlims = (0,1), ylims = (0,1), dpi = 300)
scatter!(y[test_final],y_hat_test, label = "test = $(round(score_test, digits= 2))")

importance = rf_cl.feature_importances_

sorted_importance = sortperm(importance, rev = true)

labels = All_keys[sorted_importance]

bar(labels[1:15],sort(importance, rev = true)[1:15],
xrotation=30, 
dpi = 300,
title = "RPLC important variables Classification", 
bottom_margin = 8Plots.mm,
legend = false)

# Get predicted probabilities for test set
y_prob_test = ScikitLearn.predict_proba(rf_cl, X[test_final,:])

real = y[test_final]
pred = y_hat_test

sum(real .== "Mobile")
sum(real .== "Non-mobile")
sum(real .== "Very mobile")
threshold = collect(0:0.01:1)

num_classes = 3

TPRs = zeros(num_classes, length(threshold))
FPRs = zeros(num_classes, length(threshold))
y_hat_test
for c in 1:num_classes
    for j in eachindex(threshold)
        TP, FP, TN, FN = 0, 0, 0, 0
        for i in eachindex(y_hat_test)
            if y_prob_test[i, c] >= threshold[j]
                if real[i] == pred[i]
                    TP += 1
                else
                    FP += 1
                end
            else
                if real[i] != pred[i]
                    TN += 1
                else
                    FN += 1
                end
            end
        end
        TPRs[c, j] = TP / (TP + FN)
        FPRs[c, j] = FP / (FP + TN)
    end
end

# Plot ROC curves for each class
mean_TPR = mean(TPRs, dims = 1)
mean_FPR = mean(FPRs, dims = 1)
scatter(FPRs[1,:], TPRs[1,:],
title = "ROC curve even split %B mobile phase",
xlabel = "FPR", ylabel = "TPR", 
label = "Mobile", 
dpi = 300)
scatter!(FPRs[2,:], TPRs[2,:], label = "Non-mobile")
scatter!(FPRs[3,:], TPRs[3,:], label = "Very mobile")
scatter!(mean_FPR[1,:], mean_TPR[1,:], label = "Mean")
plot!(threshold, threshold, linestyle = :dash, label = false)

rf_cl.classes_

savefig("C:\\Users\\uqthulle\\Documents\\Plots\\ROC curve evel split B mobile phase.png")