include("All functions.jl")
import ScikitLearn: RandomForestClassifier
merged_data = CSV.read("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_pH3_RPLC_data.csv", DataFrame)
fingerprints = CSV.read("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_entries_FPS.csv", DataFrame)
fingerprints = Matrix(fingerprints)
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")

#############################################################################################################
#Train/test split without data leakage

X_train, X_test, y_train, y_test = train_test_split(fingerprints, merged_data.Class, test_size=0.2, random_state=42, 
stratify = merged_data.Class, shuffle = true)

#################
#Cross-validation

skf = StratifiedKFold(n_splits=3, shuffle=true, random_state=42)

train_fold = []
test_fold = []
for (i, (train_index, test_index)) in enumerate(skf.split(X_train, y_train))

    train_index_julia = Int64.(collect(train_index)).+1
    test_index_julia = Int64.(collect(test_index)).+1
    
    push!(train_fold, train_index_julia)
    push!(test_fold, test_index_julia)
end



#### Checking for equal class distribution in each fold
sum(y_train[train_fold[1]].=="Very mobile")
sum(y_train[train_fold[1]].=="Non-mobile")
sum(y_train[train_fold[1]].=="Mobile")

sum(y_train[train_fold[2]].=="Very mobile")
sum(y_train[train_fold[2]].=="Non-mobile")
sum(y_train[train_fold[2]].=="Mobile")

sum(y_train[train_fold[3]].=="Very mobile")
sum(y_train[train_fold[3]].=="Non-mobile")
sum(y_train[train_fold[3]].=="Mobile")





#Grid search parameters
n_estimators = [50, 100, 200, 500]
max_features = [0.1, 0.5, 0.7, 1]
max_depth = [10, 30, 50 , nothing]
min_samples_split = [2, 5, 10]
min_samples_leaf = [1, 2, 4]

n_tests = 100

#Initializing the DataFrame
df_results = DataFrame(n_estimators = Int[], max_features = Float64[], max_depth = Int[],
                       min_samples_split = Int[], min_samples_leaf = Int[], 
                       CV_3_train_F1 = Float64[], CV_3_test_F1 = Float64[])
for i = 1:n_tests 
    est = rand(n_estimators)
    feat = rand(max_features)
    depth = rand(max_depth)
    split = rand(min_samples_split)
    leaf = rand(min_samples_leaf)
    
    print("Test nr $(i)/$(n_tests)")
    rf_cl = RandomForestClassifier(n_estimators = est, max_features =feat,
                                    min_samples_split = split,
                                    max_depth = depth, 
                                   min_samples_leaf = leaf,
                                   random_state = 42, class_weight = "balanced", n_jobs = -1)
    
    mean_F1_train = zeros(3)
    mean_F1_test = zeros(3)

    for k = 1:3
        labels = ["Very mobile", "Mobile", "Non-mobile"]
        cm = ConfusionMatrix(levels = ["Very mobile", "Mobile", "Non-mobile"])
        ScikitLearn.fit!(rf_cl, X_train[train_fold[k],:], y_train[train_fold[k]])
        y_hat_train_fold = ScikitLearn.predict(rf_cl,X_train[train_fold[k],:])
        y_hat_test_fold = ScikitLearn.predict(rf_cl,X_train[test_fold[k],:])
        
        c_matrix_train = cm(y_hat_train_fold, y_train[train_fold[k]])
        #TPR, FDR and F1-score for all classes for the train set
        results_train = TPR_FDR(c_matrix_train)
        mean_F1_train[k] = mean(results_train[:,4])

        c_matrix_test = cm(y_train[test_fold[k]], y_hat_test_fold)
        #TPR, FDR and F1-score for all classes for the train set
        results_test = TPR_FDR(c_matrix_test)
        mean_F1_test[k] = mean(results_test[:,4])
    end
    mean_F1_train_folds = mean(mean_F1_train)
    mean_F1_test_folds = mean(mean_F1_test)

    println("CV train F1 = $mean_F1_train_folds")
    println("CV test F1 = $mean_F1_test_folds")

    push!(df_results, [est, feat, isnothing(depth) ? -1 : depth, split, leaf, mean_F1_train_folds, mean_F1_test_folds])
end


#Viewing the results
vscodedisplay(df_results)
#saving the results
CSV.write("R:\\PHD2024TH-Q6813\\Research files\\Models and other documents\\CV scores mobility data merged dataset.csv", df_results)
CSV.write("C:\\Users\\uqthulle\\Documents\\CV3 results new model.csv", df_results)


