using CSV, Statistics, DataFrames, PubChemCrawler, StatsPlots, MultivariateStats, Images
using GLM, TypedTables, LinearAlgebra, ScikitLearn, Random, MLDataUtils
using ScikitLearn: @sk_import
@sk_import ensemble: RandomForestClassifier
@sk_import model_selection: StratifiedKFold
@sk_import model_selection: train_test_split
folder_path = "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\processed_data"
Final_table = DataFrame(Inchikey = String[], LC_mode = String7[], MACCS = String[], Pubchem_fps = String[],MW = Float64[], XlogP = Float64[], Retention_factors = Float64[], Modifier = Float64[] )
MACCS_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\MACCS keys.txt")
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")
All_keys =  [MACCS_keys ;PubChem_keys]
All_keys_log_MW = [MACCS_keys ;PubChem_keys;"MW"]
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
        if rtfs[i] > 0.7 && pmbph[i] < 10
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
function TPR_FDR(c_matrix)
    TP_M = c_matrix[1,1]
    FN_M = c_matrix[2,1] + c_matrix[3,1]
    FP_M = c_matrix[1,2] + c_matrix[1,3]
    TPR_M = TP_M/(TP_M + FN_M)
    FDR_M = FP_M/(TP_M + FP_M)
    F1_M = (2*TP_M)/(2*TP_M + FP_M + FN_M)

    TP_NM = c_matrix[2,2]
    FN_NM = c_matrix[1,2] + c_matrix[3,2]
    FP_NM = c_matrix[2,1] + c_matrix[2,3]
    TPR_NM = TP_NM/(TP_NM + FN_NM)
    FDR_NM = FP_NM/(TP_NM + FP_NM)
    F1_NM = (2*TP_NM)/(2*TP_NM + FP_NM + FN_NM)

    TP_VM = c_matrix[3,3]
    FN_VM = c_matrix[1,3] + c_matrix[2,3]
    FP_VM = c_matrix[3,1] + c_matrix[3,2]
    TPR_VM = TP_VM/(TP_VM + FN_VM)
    FDR_VM = FP_VM/(TP_VM + FP_VM)
    F1_VM = (2*TP_VM)/(2*TP_VM + FP_VM + FN_VM)
    
    Classes = ["Mobile", "Non-Mobile", "Very mobile"]
    TPRs = [TPR_M, TPR_NM, TPR_VM]
    FDRs = [FDR_M, FDR_NM, FDR_VM]
    F1_scores = [F1_M, F1_NM, F1_VM]
    df = DataFrame(Class = Classes, TPR = TPRs, FDR = FDRs, F1_scores = F1_scores)

    return df
end
function remove_high_std(RPLC_data, std_threshold)
    unique_compounds = unique(RPLC_data[:,1])

    standard_devs = []

    for i in eachindex(unique_compounds)
        if i % 100 == 0
            println("$(round(i/length(unique_compounds)*100, digits = 2))%")
        end

        occurrences = findall(x-> x == unique_compounds[i],RPLC_data[:,1])

        push!(standard_devs,std(RPLC_data[occurrences,end]./100))
    
    end

    high_std_cmp = unique_compounds[findall(x-> x>=std_threshold, standard_devs)]

    filtered_RPLC_data = RPLC_data[.!in.(RPLC_data.Inchikey, Ref(high_std_cmp)), :]

    return filtered_RPLC_data

end
function train_test_split_no_leakage_classifier(filtered_RPLC_data, split)
    #removed_fp = 668
    #fingerprints = PubChem_fps[:, [1:removed_fp-1; removed_fp+1:end]]
    

    # Generate a random permutation of column indices
    #perm = randperm(size(PubChem_fps, 2))

    # Shuffle the columns of the matrix using the permutation
    #fingerprints = PubChem_fps[:, perm]
    fingerprints = PubChem_fps

    #fingerprints = fingerprints .- mean(fingerprints,dims = 1)

    #fingerprints = fingerprints ./ std(fingerprints,dims=1)

    #fingerprints = replace(fingerprints, NaN => 0.0)

    unique_compounds = unique(filtered_RPLC_data[:,1])

    p_b = filtered_RPLC_data[:,end]./100

    y = assign_labels(p_b)

    first_labels = strat_labels(unique_compounds,filtered_RPLC_data[:,1], y)
    
    train, test = train_test_split(unique_compounds, test_size = split, random_state = 42, stratify = first_labels, shuffle = true)
    
    # Initialize variables
    train_indices = []
    test_indices = []

    # Find indices for train and test data
    @time for i in eachindex(unique_compounds)
        if (i % 100) == 0
            println("$(round(i/length(unique_compounds)*100, digits = 2))%")
        end
        occurrences = findall(x -> x == unique_compounds[i], filtered_RPLC_data[:, 1])
        if unique_compounds[i] in train
            append!(train_indices, occurrences)
        else
            append!(test_indices, occurrences)
        end
    end
    Random.seed!(42)
     # Shuffle the indices
     shuffle!(train_indices)
     shuffle!(test_indices)
    # Extract train and test data
    X_train = fingerprints[train_indices, :]
    X_test = fingerprints[test_indices, :]
    y_train = y[train_indices]
    y_test = y[test_indices]

    return X_train, X_test, y_train, y_test, train_indices, test_indices, y
end
function assign_labels(Φ)
    y = String[]
    for i in eachindex(Φ)
        if Φ[i] >= 0.6
            push!(y, "Non-mobile")
        elseif Φ[i] <= 0.2
            push!(y, "Very mobile")
        else 
            push!(y, "Mobile")
        end
    end
    return y
end
function strat_labels(unique_vector, original_data, labels)

    index_first_occurrence = Vector{Int}(undef, length(unique_vector))

    for i in eachindex(unique_vector)

        if (i % 100) == 0
            println("$(round(i/length(unique_vector)*100, digits = 2))%")
        end
        index_first_occurrence[i] = findfirst(x-> x == unique_vector[i], original_data[:,1])

    end

    strat_labels = labels[index_first_occurrence]
    return strat_labels
end
function remove_outliers_IQR(RPLC_data)
    unique_compounds = unique(RPLC_data[:,1])

    non_outliers = []
    for i in eachindex(unique_compounds)
        if i % 100 == 0
            println("$(round(i/length(unique_compounds)*100, digits = 2))%")
        end

        occurrences = findall(x-> x == unique_compounds[i],RPLC_data[:,1])
        group = RPLC_data[occurrences,end]./100
        Q1 = quantile(group, 0.25)
        Q3 = quantile(group, 0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        indices_non_outliers = occurrences[(group .>= lower_bound) .& (group .<= upper_bound)]
        append!(non_outliers, indices_non_outliers)
    end

    filtered_RPLC_data = RPLC_data[non_outliers,:]

    return filtered_RPLC_data
end
for i = 1:374
    #Locating the file
    current_file = joinpath("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\processed_data", readdir(folder_path)[i])
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
         retention_factors = repeat(["No gradient info"], length(RT))
         Modifier = repeat(["No gradient info"], length(RT))
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

#Removing rows with missing retention factor and outliers
Final_table_unique = remove_missing_and_outliers(Final_table)

#Formatting MACCS and PubChem Fingerprints for modelling
MACCS, PubChem_fps = format_fingerprints(Final_table_unique)

#Getting indices for RPLC and HILIC rows
indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table_unique))

RPLC_data = Final_table_unique[indices_RPLC,:]

#############################################################################################################
#Removing outliers if IQR > 1.5 * threshold
@time filtered_RPLC_data = remove_outliers_IQR(RPLC_data)
filtered_RPLC_data
MACCS, PubChem_fps = format_fingerprints(filtered_RPLC_data)
PubChem_fps

#############################################################################################################
#Train/test split without data leakage

X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(filtered_RPLC_data, 0.1)

#################
#Cross-validation
X_train
X_test

unique_train = unique(filtered_RPLC_data[train_indices,1])

Φ = filtered_RPLC_data[train_indices,end]./100
y = assign_labels(Φ)
first_labels = strat_labels(unique_train,filtered_RPLC_data[train_indices,1],y)

skf = StratifiedKFold(n_splits=3, shuffle=true, random_state=42)

train_fold = []
test_fold = []

for (i, (train_index, test_index)) in enumerate(skf.split(unique_indices, first_labels))

    train_index_julia = Int64.(collect(train_index)).+1
    test_index_julia = Int64.(collect(test_index)).+1
    
    push!(train_fold, train_index_julia)
    push!(test_fold, test_index_julia)
end

function append_occurrences_folds(unique_train::Vector{String}, fold::Int)
    # Initialize variables
    train_indices_fold = Int[]
    test_indices_fold = Int[]

    # Find indices for train and test data
    @time for i in eachindex(unique_train)
        if (i % 10) == 0
            println("$(round(i/length(unique_train)*100, digits = 2))%")
        end
        occurrences = findall(x -> x == unique_train[i], filtered_RPLC_data[train_indices, 1])
        if unique_train[i] in unique_train[train_fold[fold]]
            append!(train_indices_fold, occurrences)
        else
            append!(test_indices_fold, occurrences)
        end
    end

    shuffle!(train_indices_fold)
    shuffle!(test_indices_fold)
    return train_indices_fold, test_indices_fold

end

train_fold_1,test_fold_1 = append_occurrences_folds(unique_train, 1) 

train_fold_2,test_fold_2 = append_occurrences_folds(unique_train, 2) 

train_fold_3,test_fold_3 = append_occurrences_folds(unique_train, 3) 

unique_train[train_fold_1]

sum(y_train[train_fold_1].=="Very mobile")
sum(y_train[train_fold_1].=="Non-mobile")
sum(y_train[train_fold_1].=="Mobile")

sum(y_train[train_fold_2].=="Very mobile")
sum(y_train[train_fold_2].=="Non-mobile")
sum(y_train[train_fold_2].=="Mobile")

sum(y_train[train_fold_3].=="Very mobile")
sum(y_train[train_fold_2].=="Non-mobile")
sum(y_train[train_fold_2].=="Mobile")
using MLJ
#Making the model CV
fingeprints = [MACCS PubChem_fps][train_indices,:]
function CV_3_RFC(n_estimatorss, max_featuress, max_depths, min_samples_splits, min_samples_leafs)

    rf_cl = RandomForestClassifier(n_estimators = n_estimatorss, max_features = max_featuress,
                                   max_depth = max_depths, min_samples_split = min_samples_splits,
                                   min_samples_leaf = min_samples_leafs,
                                   random_state = 42, class_weight = "balanced", n_jobs = -1)

    #Fold 1
    println("Fold1")
    ScikitLearn.fit!(rf_cl, X_train[train_fold_1,:], y_train[train_fold_1], )
    
    score_train_fold_1 = ScikitLearn.score(rf_cl, X_train[train_fold_1,:], y_train[train_fold_1])
    score_test_fold_1 = ScikitLearn.score(rf_cl, X_train[test_fold_1,:], y_train[test_fold_1])

    y_hat_train = ScikitLearn.predict(rf_cl,X_train[train_fold_1,:])
    y_hat_test = ScikitLearn.predict(rf_cl,X_train[test_fold_1,:])

    c_matrix_train_f1 = confusion_matrix(y_hat_train, y_train[train_fold_1])
    results_train_f1 = TPR_FDR(c_matrix_train_f1)
    F1_train_f1 = mean(results_train_f1[:,4])

    c_matrix_test_f1 = confusion_matrix(y_hat_test, y_train[test_fold_1])
    results_test_f1 = TPR_FDR(c_matrix_test_f1)
    F1_test_f1 = mean(results_test_f1[:,4])
    #Fold 2
    println("Fold2")
    ScikitLearn.fit!(rf_cl, X_train[train_fold_2,:], y_train[train_fold_2])
    
    score_train_fold_2 = ScikitLearn.score(rf_cl, X_train[train_fold_2,:], y_train[train_fold_2])
    score_test_fold_2 = ScikitLearn.score(rf_cl, X_train[test_fold_2,:], y_train[test_fold_2])
    
    y_hat_train = ScikitLearn.predict(rf_cl,X_train[train_fold_2,:])
    y_hat_test = ScikitLearn.predict(rf_cl,X_train[test_fold_2,:])

    c_matrix_train_f2 = confusion_matrix(y_hat_train, y_train[train_fold_2])
    results_train_f2 = TPR_FDR(c_matrix_train_f2)
    F1_train_f2 = mean(results_train_f2[:,4])

    c_matrix_test_f2 = confusion_matrix(y_hat_test, y_train[test_fold_2])
    results_test_f2 = TPR_FDR(c_matrix_test_f2)
    F1_test_f2 = mean(results_test_f2[:,4])
    #Fold 3
    println("Fold3")
    ScikitLearn.fit!(rf_cl, X_train[train_fold_3,:], y_train[train_fold_3])
    
    score_train_fold_3 = ScikitLearn.score(rf_cl, X_train[train_fold_3,:], y_train[train_fold_3])
    score_test_fold_3 = ScikitLearn.score(rf_cl, X_train[test_fold_3,:], y_train[test_fold_3])

    y_hat_train = ScikitLearn.predict(rf_cl,X_train[train_fold_3,:])
    y_hat_test = ScikitLearn.predict(rf_cl,X_train[test_fold_3,:])

    c_matrix_train_f3 = confusion_matrix(y_hat_train, y_train[train_fold_3])
    results_train_f3 = TPR_FDR(c_matrix_train_f3)
    F1_train_f3 = mean(results_train_f3[:,4])

    c_matrix_test_f3 = confusion_matrix(y_hat_test, y_train[test_fold_3])
    results_test_f3 = TPR_FDR(c_matrix_test_f3)
    F1_test_f3 = mean(results_test_f3[:,4])

    mean_train_score = mean([score_train_fold_1,score_train_fold_2,score_train_fold_3])
    mean_test_score = mean([score_test_fold_1,score_test_fold_2,score_test_fold_3])
    mean_train_F1 = mean([F1_train_f1, F1_train_f2, F1_train_f3])
    mean_test_F1 = mean([F1_test_f1, F1_test_f2, F1_test_f3])
    return mean_train_score, mean_test_score, mean_train_F1, mean_test_F1
end

df_results = DataFrame(n_estimators = Int[], max_features = Float64[], max_depth = Int[],
                       min_samples_split = Int[], min_samples_leaf = Int[], 
                       CV_3_train_score = Float64[],
                       CV_3_test_score = Float64[], CV_3_train_F1 = Float64[], CV_3_test_F1 = Float64[])


n_estimatorss = [10,20,50,100,200]
max_featuress = [1.0, 0.5, 0.2, 0.1]
max_depths = [10, 25, 100]
min_samples_splits = [2,4,6]
min_samples_leafs = [1,2,4]

for n_estimators in n_estimatorss
    for max_features in max_featuress
        for max_depth in max_depths
            for min_samples_split in min_samples_splits
                for min_samples_leaf in min_samples_leafs

                    @show n_estimators, max_features, max_depth, min_samples_split, min_samples_leaf 

                    train_score, test_score, F1_train, F1_test = CV_3_RFC(n_estimators, max_features, max_depth, min_samples_split, min_samples_leaf)

                    @show train_score, test_score

                    push!(df_results, [n_estimators, max_features, max_depth, min_samples_split, min_samples_leaf, train_score, test_score, F1_train, F1_test])

                end
            end
        end
    end
end

vscodedisplay(df_results)

CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\All CV scores mobility data full dataset.csv", df_results)

##############
#Model with best settings
using MLJ
rf_cl = RandomForestClassifier(n_estimators = 100, max_features = 0.5,
                                   max_depth = 100, min_samples_split = 2,
                                   min_samples_leaf = 2,
                                   random_state = 42, class_weight = "balanced", n_jobs = -1)


ScikitLearn.fit!(rf_cl, X_train, y_train)
depths = [maximum([tree.tree_.max_depth for tree in rf_cl.estimators_]) for _ in 1:length(rf_cl.estimators_)]
y_hat_train = ScikitLearn.predict(rf_cl,X_train)
y_hat_test = ScikitLearn.predict(rf_cl,X_test)

score_train = ScikitLearn.score(rf_cl, X_train, y_train)
score_test = ScikitLearn.score(rf_cl, X_test, y_test)

c_matrix = confusion_matrix(y_hat_train, y_train)

results = TPR_FDR(c_matrix)

scatter([1,2,3], results[:,2], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "TPR", shape = :diamond, xlims = (0.5,3.5), grid = :y, markersize = 4,
title = "Training data n = $(length(y_train))", titlefont = font(10),
legend = false)

scatter!([1,2,3], results[:,3], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "FDR", shape = :square, markersize = 4)

p_train = scatter!([1,2,3], results[:,4], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "F1 score", shape = :utriangle, dpi = 300, markersize = 4)

c_matrix = confusion_matrix(y_hat_test, y_test)

results = TPR_FDR(c_matrix)

scatter([1,2,3], results[:,2], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "TPR", shape = :diamond, xlims = (0.5,3.5), grid = :y, markersize = 4,
title = "Test data n = $(length(y_test))", titlefont = font(10))

scatter!([1,2,3], results[:,3], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "FDR", shape = :square, markersize = 4)

p_test = scatter!([1,2,3], results[:,4], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "F1 score", shape = :utriangle, dpi = 300, markersize = 4)


plot(p_train, p_test, size = (800,400), dpi = 300)

#Checking feature importances
importance = rf_cl.feature_importances_
sum(importance)
sorted_importance = sortperm(importance, rev = true)
sorted_importance_original
labels_original = PubChem_keys[sorted_importance_original]
labels = PubChem_keys[sorted_importance]

bar(labels[1:10],sort(importance, rev = true)[1:10],
xrotation=25, 
dpi = 300,
title = "RPLC important variables Classification", 
left_margin = 10Plots.mm,
bottom_margin = 10Plots.mm,
legend = false)

#Applicability domain

lev = calculate_leverage(X_train, X_test)

histogram(lev, dpi = 300,
xlabel = "hᵢ", label = "test", xlims = (0,0.2))

warning_h = (3 * length(X_train[1,:])+1)/length(X_train[:,1])

vline!([warning_h], label = "warning (h*)")

melamine = findall(x-> x =="InChI=1S/C3H6N6/c4-1-7-2(5)9-3(6)8-1/h(H6,4,5,6,7,8,9)",filtered_RPLC_data[:,1])

scatter(filtered_RPLC_data[melamine, end], ylims = (0,100))

sum(X_train[:,527])/length(X_train[:,1])*100



#Tracking misclassifications

indices = []
for i in eachindex(y_test)
    if y_test[i] == "Non-mobile" && y_hat_test[i] == "Non-mobile"
        push!(indices, i)
    end
end
indices

cmp_indices = test_indices[indices]


Inchi = filtered_RPLC_data[cmp_indices,1][7]

p_bs = filtered_RPLC_data[cmp_indices,end]

histogram(filtered_RPLC_data[cmp_indices,6])

filtered_RPLC_data[cmp_indices,5]

histogram(p_bs./100, bins = 10, dpi = 300, xlims = (0,1),
ylabel = "Frequency",
xlabel = "Φ",
label = "Non-mobile  misclassified as mobile",
legend = :topright, legendfont = font(12), xtickfont=font(12), 
ytickfont=font(12), 
guidefont=font(18),
ylims = (0,300))

vline!([0.6], linestyle = :dash, label = "Non-mobile label threshold",
c = :red)

cmps = findall(x-> x == Inchi, filtered_RPLC_data[:,1])
means = mean(filtered_RPLC_data[cmps,end]./100)
stds = std(filtered_RPLC_data[cmps,end]./100)

cmps = findall(x-> x == Inchi, filtered_RPLC_data[:,1])
means = mean(filtered_RPLC_data[cmps,end]./100)
stds = std(filtered_RPLC_data[cmps,end]./100)

histogram(filtered_RPLC_data[cmp_indices,end]./100, bins = 10)


scatter(filtered_RPLC_data[cmps,end]./100, ylims = (0,1),
dpi = 300, title = "All entries Bis(2-ethylhexyl) phthalate", 
xlabel = "Entry nr.", ylabel = "Φ", legend = false)




scatter(p_bs./100,  dpi = 300,
title = "Very mobile misclassified as mobile",
ylabel = "ϕ",
xlabel = "Compound nr",
legend = false)

scatter!(RTs, dpi = 300,
title = "Very mobile misclassified as mobile",
ylabel = "%B elution",
xlabel = "N compound",
label = "Retention factor")


MW = filtered_RPLC_data[:,5]
XlogP = filtered_RPLC_data[:,6]
Φ = filtered_RPLC_data[:,end]./100

scatter(MW, XlogP, alpha = 0.4, zcolor = Φ,
label = false, dpi = 300, color_label = "Φ")