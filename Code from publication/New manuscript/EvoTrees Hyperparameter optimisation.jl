include("All functions.jl")
using EvoTrees
RPLC_pH3_pH26_data = CSV.read("C:\\Users\\uqthulle\\Documents\\RPLC_data_pH3and2.6_updated.csv", DataFrame)
fingerprints = CSV.read("C:\\Users\\uqthulle\\Documents\\pH 3 and 2.6 RepoRT fingerprints final.csv", DataFrame)
fingerprints = Matrix(fingerprints)
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")

#These are the fingerprint names for each corresponding column
#Splitting data into train and test by avoiding data leakage

shuffle_indices = shuffle(collect(1:length(filtered_RPLC_data[:,1])))
shuffle_indices_test = shuffle(collect(1:length(filtered_RPLC_data_test[:,1])))

X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(RPLC_pH3_pH26_data, 0.1, fingerprints)

X_train, X_eval, y_train, y_eval, train_indices, eval_indices, y = train_test_split_no_leakage_classifier(RPLC_pH3_pH26_data[train_indices,:], 0.111, fingerprints[train_indices,:])

x_train = X_train 

x_eval = X_eval 



function random_hyperparams()
    config = EvoTreeClassifier(
        early_stopping_rounds = 5,
        nrounds = 100,
        eta = rand([0.01, 0.05, 0.1, 0.2, 1]),  
        max_depth = rand([3, 5, 7, 9, 12]),  
        lambda = rand([0.0, 0.1, 1, 10, 100]),  
        gamma = rand([0.0, 0.1, 1, 10, 100]),  
        min_weight = rand([1, 5, 10]),  # Avoiding 10 for now
        rowsample = rand([0.5, 0.6, 0.7, 0.8, 1.0]),  
        colsample = rand([0.5, 0.6, 0.7, 0.8, 1.0]),  
        nbins = rand([16,32, 64, 128]),  # Avoiding small bin counts
        rng = 1  
    )  # Random seed
    return config
end



# Function to evaluate model performance with retries
function evaluate_model(config, x_train, y_train, x_eval, y_eval; max_retries=5)
    retries = 0
    while retries ≤ max_retries
        try
            model = EvoTrees.fit_evotree(config; x_train, y_train, x_eval, y_eval, feature_names=PubChem_keys, print_every_n=1)

            y_hat_train_probas = EvoTrees.predict(model, x_train)
            y_hat_test_probas = EvoTrees.predict(model, x_test)

            labels = ["Mobile", "Non-mobile", "Very mobile"]
            y_hat_train = labels[argmax.(eachrow(y_hat_train_probas))]
            y_hat_test = labels[argmax.(eachrow(y_hat_test_probas))]

            # Scores for train and test
            score_train = sum(y_hat_train .== y_train) / length(y_train)
            score_test = sum(y_hat_test .== y_test) / length(y_test)

            # Confusion matrices
            c_matrix_train = confusion_matrix(y_hat_train, y_train)
            results_train = TPR_FDR(c_matrix_train)
            F1_mean_train = mean(results_train[:, 4])

            c_matrix_test = confusion_matrix(y_hat_test, y_test)
            results_test = TPR_FDR(c_matrix_test)
            F1_mean_test = mean(results_test[:, 4])

            return score_train, score_test, F1_mean_train, F1_mean_test

        catch e
            retries += 1
            println(" ⚠️ Error in evaluation (Attempt $retries/$max_retries): $(e)")
            if retries >= max_retries
                println(" ❌ Skipping this trial after $max_retries failed attempts.\n")
                return nothing, nothing, nothing, nothing  # Return safe values
            end
            
        end
    end
end


####Testing boosted tree with only pH 3

RPLC_pH3_data = CSV.read("C:\\Users\\uqthulle\\Documents\\RPLC_data_pH3.csv", DataFrame)
fingerprints = CSV.read("C:\\Users\\uqthulle\\Documents\\pH 3 RepoRT fingerprints final.csv", DataFrame)
fingerprints = Matrix(fingerprints)

X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(RPLC_pH3_data, 0.1, fingerprints)

X_train, X_eval, y_train, y_eval, train_indices, eval_indices, y = train_test_split_no_leakage_classifier(RPLC_pH3_data[train_indices,:], 0.111, fingerprints[train_indices,:])

x_train = X_train
x_eval = X_eval
x_test = X_test
df_results = DataFrame(
    eta=Float64[], max_depth=Int[], lambda=Float64[], gamma=Float64[],
    min_weight=Int[], rowsample=Float64[], colsample=Float64[], nbins=Int[],
    train_accuracy=Float64[], test_accuracy=Float64[],
    F1_mean_train=Float64[], F1_mean_test=Float64[]
)

# Hyperparameter optimization loop
n_trials = 100
println("Starting hyperparameter optimization with $n_trials trials...\n")

for i in 1:n_trials
    config = random_hyperparams()
    println("Trial $i")

    score_train, score_test, F1_mean_train, F1_mean_test = evaluate_model(config, x_train, y_train, x_eval, y_eval)

    if score_train !== nothing  # Only append successful trials
        println(" -> Train Accuracy: $score_train, Test Accuracy: $score_test, F1_train: $F1_mean_train, F1_test: $F1_mean_test\n")
        
        push!(df_results, [config.eta, config.max_depth, config.lambda, config.gamma, config.min_weight,
                           config.rowsample, config.colsample, config.nbins, score_train, score_test,
                           F1_mean_train, F1_mean_test])
    else
        println(" ❌ Trial $i failed completely. Moving to the next trial.\n")
    end
    if i % 10 == 0
        println("💾 Saving results...")
        CSV.write("C:\\Users\\uqthulle\\Documents\\Evo Hyperparameter pH3 and pH6 only.csv", df_results)
        GC.gc()  # Free up memory periodically
    end
end
 vscodedisplay(df_results)
CSV.write("C:\\Users\\uqthulle\\Documents\\Evo Hyperparameter pH3 and pH6 only.csv", df_results)


##Testing best model
config = EvoTreeClassifier(
        early_stopping_rounds = 5,
        nrounds = 200,
        eta = 0.2,  
        max_depth = 12,
        lambda = 0,
        gamma = 0.1,  
        min_weight = 10,
        rowsample = 0.7,  
        colsample = 0.5,  
        nbins = 64,  # Avoiding small bin counts
        rng = 1  
    )  # Random seed

    model = EvoTrees.fit_evotree(config; x_train, y_train, x_eval, y_eval, feature_names=PubChem_keys, print_every_n=1)

    y_hat_train_probas = EvoTrees.predict(model, x_train)
    y_hat_test_probas = EvoTrees.predict(model, x_test)

    labels = ["Mobile", "Non-mobile", "Very mobile"]
    y_hat_train = labels[argmax.(eachrow(y_hat_train_probas))]
    y_hat_test = labels[argmax.(eachrow(y_hat_test_probas))]

    # Scores for train and test
    score_train = sum(y_hat_train .== y_train) / length(y_train)
    score_test = sum(y_hat_test .== y_test) / length(y_test)

    # Confusion matrices
    c_matrix_train = confusion_matrix(y_hat_train, y_train)
    results_train = TPR_FDR(c_matrix_train)
    F1_mean_train = mean(results_train[:, 4])

    c_matrix_test = confusion_matrix(y_hat_test, y_test)
    results_test = TPR_FDR(c_matrix_test)
    F1_mean_test = mean(results_test[:, 4])