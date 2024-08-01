include("All functions.jl")
folder_path = "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\raw_data"
Final_table = DataFrame(SMILES = String[], InChi = String[], LC_mode = String7[], Column_name = String[], CID = String[], Retention_factors = Float64[], Modifier = Float64[])
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")
for i = 1:374
    #Locating the file
    current_file_raw = joinpath("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\raw_data", readdir(folder_path)[i+1])
    current_file_processed = joinpath("C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\RepoRT-master\\processed_data", readdir(processed_folder_path)[1:end-1][i])
    println(i)
    #Getting the SMILES and retention times
    rt_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "rtdata.tsv"],"_"))
    data_compounds = CSV.read(rt_path, DataFrame)
    Compound_SMILES = data_compounds[:,7]
    InChI = data_compounds[:,8]
    RT = data_compounds[:,4]
    CID = data_compounds[:,5]
    #Getting the retention factors
    gradient_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "gradient.tsv"],"_"))
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
    #Getting the LC mode info for the data set
    info_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "info.tsv"],"_"))
    info = CSV.read(info_path, DataFrame)
    #Check if the column info is missing in the meta data
    LC_mode = info[1,3]
    if ismissing(LC_mode) || LC_mode == "NA"

        LC_mode = repeat(["No mode info"], length(RT))

    else

        LC_mode = repeat([LC_mode], length(RT))

    end

    #Metadata 
    meta_path = joinpath(current_file_raw, join(["$(readdir(folder_path)[1:end][i+1])", "metadata.tsv"],"_"))
    meta = CSV.read(meta_path, DataFrame)
    Column_name = meta[1,2]
    #Getting all the data in one table to append all the other datasets in this format
    current_table = DataFrame(SMILES = Compound_SMILES, InChi = InChI, LC_mode = LC_mode, Column_name = Column_name, CID = CID, Retention_factors = retention_factors, Modifier = Modifier)
    Final_table = vcat(Final_table, current_table)
end
Final_table

Final_table_unique = remove_missing_and_outliers(Final_table)

#Getting indices for RPLC and HILIC rows
indices_RPLC = findall(row -> occursin("RP", row.LC_mode), eachrow(Final_table_unique))

RPLC_data = Final_table_unique[indices_RPLC,:]

RPLC_data = remove_missing_identifiers(RPLC_data)

filtered_RPLC_data = remove_outliers_IQR(RPLC_data)

fingerprints = CSV.read("R:\\PHD2024TH-Q6813\\Models and other documents\\RepoRT fingerprints final.csv", DataFrame)

#############################################################################################################
#Train/test split without data leakage

X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(filtered_RPLC_data, 0.1, fingerprints)

#################
#Cross-validation
X_train = Matrix(X_train)
X_test = Matrix(X_test)
y_train
y_test
unique_train = string.(unique(filtered_RPLC_data.InChi[train_indices]))
original_data = string.(filtered_RPLC_data.InChi[train_indices])
Φ = filtered_RPLC_data.Modifier[train_indices]./100
y = assign_labels(Φ)
first_labels = strat_labels(unique_train,original_data,y)

skf = StratifiedKFold(n_splits=3, shuffle=true, random_state=42)

train_fold = []
test_fold = []
unique_indices =  string.(unique(filtered_RPLC_data.InChi[train_indices]))
for (i, (train_index, test_index)) in enumerate(skf.split(unique_indices, first_labels))

    train_index_julia = Int64.(collect(train_index)).+1
    test_index_julia = Int64.(collect(test_index)).+1
    
    push!(train_fold, train_index_julia)
    push!(test_fold, test_index_julia)
end


train_fold_1,test_fold_1 = append_occurrences_folds(unique_train, 1) 

train_fold_2,test_fold_2 = append_occurrences_folds(unique_train, 2) 

train_fold_3,test_fold_3 = append_occurrences_folds(unique_train, 3) 


#### Checking for equal class distribution in each fold
sum(y_train[train_fold_1].=="Very mobile")
sum(y_train[train_fold_1].=="Non-mobile")
sum(y_train[train_fold_1].=="Mobile")

sum(y_train[train_fold_2].=="Very mobile")
sum(y_train[train_fold_2].=="Non-mobile")
sum(y_train[train_fold_2].=="Mobile")

sum(y_train[train_fold_3].=="Very mobile")
sum(y_train[train_fold_2].=="Non-mobile")
sum(y_train[train_fold_2].=="Mobile")



#Setting up cross-validation
function CV_3_RFC(n_estimatorss, max_featuress, max_depths, min_samples_splits, min_samples_leafs)

    rf_cl = RandomForestClassifier(n_estimators = n_estimatorss, max_features = max_featuress,
                                   max_depth = max_depths, min_samples_split = min_samples_splits,
                                   min_samples_leaf = min_samples_leafs,
                                   random_state = 42, class_weight = "balanced", n_jobs = -1)

    #Fold 1
    println("Fold1")
    ScikitLearn.fit!(rf_cl, X_train[train_fold_1,:], y_train[train_fold_1])
    
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

#Initializing the DataFrame
df_results = DataFrame(n_estimators = Int[], max_features = Float64[], max_depth = Int[],
                       min_samples_split = Int[], min_samples_leaf = Int[], 
                       CV_3_train_score = Float64[],
                       CV_3_test_score = Float64[], CV_3_train_F1 = Float64[], CV_3_test_F1 = Float64[])

#Grid search parameters
n_estimatorss = [200]
max_featuress = [0.5, 0.2, 0.1]
max_depths = [10, 25, 50, 100]
min_samples_splits = [2,4,6]
min_samples_leafs = [1,2,4]

#Performing grid search with 3-fold cross-validation
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

#Viewing the results
vscodedisplay(df_results)
#saving the results
CSV.write("R:\\PHD2024TH-Q6813\\Models and other documents\\All RepoRT CV scores mobility data full dataset.csv", df_results)


