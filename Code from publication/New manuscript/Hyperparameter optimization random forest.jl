include("All functions.jl")
import ScikitLearn: RandomForestClassifier
RPLC_pH3_pH26_data = CSV.read("C:\\Users\\uqthulle\\Documents\\RPLC_data_pH3and2.6_updated.csv", DataFrame)
fingerprints = CSV.read("C:\\Users\\uqthulle\\Documents\\pH 3 and 2.6 RepoRT fingerprints final.csv", DataFrame)
fingerprints = Matrix(fingerprints)
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")

#############################################################################################################
#Train/test split without data leakage

X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(filtered_RPLC_data, 0.1, fingerprints)

#################
#Cross-validation
X_train = Matrix(X_train)
X_test = Matrix(X_test)
y_train
y_test
unique_train = string.(unique(RPLC_pH3_pH26_data.InChi[train_indices]))
original_data = string.(RPLC_pH3_pH26_data.InChi[train_indices])
Φ = RPLC_pH3_pH26_data.Modifier[train_indices]./100
y = assign_labels(Φ)
first_labels = strat_labels(unique_train,original_data,y)

skf = StratifiedKFold(n_splits=3, shuffle=true, random_state=42)

train_fold = []
test_fold = []
unique_indices =  string.(unique(RPLC_pH3_pH26_data.InChi[train_indices]))
for (i, (train_index, test_index)) in enumerate(skf.split(unique_indices, first_labels))

    train_index_julia = Int64.(collect(train_index)).+1
    test_index_julia = Int64.(collect(test_index)).+1
    
    push!(train_fold, train_index_julia)
    push!(test_fold, test_index_julia)
end


test_fold_1, train_fold_1 = append_occurrences_folds(unique_train, 1)
###Saving
df_train_fold_1 = DataFrame(indices = train_fold_1)
CSV.write("C:\\Users\\uqthulle\\Documents\\train_fold_1.csv", df_train_fold_1)
df_test_fold_1 = DataFrame(indices = test_fold_1)
CSV.write("C:\\Users\\uqthulle\\Documents\\test_fold_1.csv", df_test_fold_1)

test_fold_2,train_fold_2 = append_occurrences_folds(unique_train, 2) 
#Saving
df_train_fold_2 = DataFrame(indices = train_fold_2)
CSV.write("C:\\Users\\uqthulle\\Documents\\train_fold_2.csv", df_train_fold_2)
df_test_fold_2 = DataFrame(indices = test_fold_2)
CSV.write("C:\\Users\\uqthulle\\Documents\\test_fold_2.csv", df_test_fold_2)

test_fold_3, train_fold_3 = append_occurrences_folds(unique_train, 3) 
#Saving
df_train_fold_3 = DataFrame(indices = train_fold_3)
CSV.write("C:\\Users\\uqthulle\\Documents\\train_fold_3.csv", df_train_fold_3)
df_test_fold_3 = DataFrame(indices = test_fold_3)
CSV.write("C:\\Users\\uqthulle\\Documents\\test_fold_3.csv", df_test_fold_3)

train_fold_1 = vec(CSV.read("C:\\Users\\uqthulle\\Documents\\train_fold_1.csv", DataFrame).indices)
test_fold_1 = vec(CSV.read("C:\\Users\\uqthulle\\Documents\\test_fold_1.csv", DataFrame).indices)
train_fold_2 = vec(CSV.read("C:\\Users\\uqthulle\\Documents\\train_fold_2.csv", DataFrame).indices)
test_fold_2 = vec(CSV.read("C:\\Users\\uqthulle\\Documents\\test_fold_2.csv", DataFrame).indices)
train_fold_3 = vec(CSV.read("C:\\Users\\uqthulle\\Documents\\train_fold_3.csv", DataFrame).indices)
test_fold_3 = vec(CSV.read("C:\\Users\\uqthulle\\Documents\\test_fold_3.csv", DataFrame).indices)

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
n_estimatorss = [100]
max_featuress = [0.1]
max_depths = [25]
min_samples_splits = [2]
min_samples_leafs = [2]

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
CSV.write("C:\\Users\\uqthulle\\Documents\\CV3 results new model.csv", df_results)


