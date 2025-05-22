include("All functions.jl")
import ScikitLearn: RandomForestClassifier
merged_data = CSV.read("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_expanded_RPLC_pH3_data.csv", DataFrame)
fingerprints = CSV.read("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_expanded_entries_FPS.csv", DataFrame)
fingerprints = Matrix(fingerprints)
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")

X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(merged_data, 0.1, fingerprints)

using DecisionTree






n_trees = 50                    # Fewer if tuning
partial_sampling = 1          # Use 30% of samples per tree
max_depth = -1                  # Prevent overly deep trees
n_subfeatures = 30              # sqrt(881) ≈ 30
min_samples_leaf = 20           # Faster & better generalization
min_samples_split = 10          # Reduce overgrowth
min_purity_increase = 0.01 

model = @time build_forest(y_train, X_train,
                     n_subfeatures,
                     n_trees,
                     partial_sampling,
                     max_depth,
                     min_samples_leaf,
                     min_samples_split,
                     min_purity_increase;
                     rng = 42)  # enables parallelization


y_hat_train = apply_forest(model, X_train)
y_hat_test = apply_forest(model, X_test)

c_matrix = confusion_matrix(y_hat_train, y_train)

#TPR, FDR and F1-score for all classes for the train set
results = TPR_FDR(c_matrix)

#Confusion matrix for the test set
c_matrix = confusion_matrix(y_hat_test, y_test)

#TPR, FDR and F1-score for all classes for the test set
results = TPR_FDR(c_matrix)