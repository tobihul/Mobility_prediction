include("All functions.jl")
using EvoTrees
merged_data = CSV.read("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_expanded_RPLC_pH3_data.csv", DataFrame)
fingerprints = CSV.read("R:\\PHD2024TH-Q6813\\Research Files\\Models and other documents\\merged_expanded_entries_FPS.csv", DataFrame)
fingerprints = Matrix(fingerprints)
PubChem_keys = readlines( "C:\\Users\\uqthulle\\OneDrive - The University of Queensland\\Documents\\PubChem keys.txt")

X_train, X_test, y_train, y_test, train_indices, test_indices, y = train_test_split_no_leakage_classifier(merged_data, 0.1, fingerprints)

#Obtaining eval set

X_train, X_eval, y_train, y_eval, train_indices, eval_indices, y = train_test_split_no_leakage_classifier(merged_data[train_indices,:], 0.111, fingerprints[train_indices,:])

y_train = String.(vcat(y_train...))
y_test = String.(vcat(y_test...))
x_train = X_train

x_eval = X_eval

#Build the model with the optimal settings using EvoTrees
using EvoTrees: fit

config = EvoTreeClassifier(
    early_stopping_rounds = 10,
    nrounds = 200,
    eta = 0.1,
    L2 = 0.1,
    lambda = 0,
    gamma = 0.1,
    max_depth = 12,
    min_weight = 10,
    rowsample = 0.7,
    colsample = 1,
    nbins = 32,
    rng = 1,
    device = :cpu

)
y_train_cat = categorical(y_train)
classes = levels(y_train_cat)
@time model = EvoTrees.fit_evotree(config; x_train, y_train, feature_names = PubChem_keys, verbosity=1)

y_hat_train_probas = EvoTrees.predict(model, x_train)
y_hat_test_probas = EvoTrees.predict(model, X_test)

labels = ["Mobile", "Non-mobile", "Very mobile"]

y_hat_train = labels[argmax.(eachrow(y_hat_train_probas))]
y_hat_test = labels[argmax.(eachrow(y_hat_test_probas))]

#Scores of train and test
score_train = sum(y_hat_train .== y_train) / length(y_train)

score_test = sum(y_hat_test .== y_test) / length(y_test)

#Confusion matrix for the train set
c_matrix = confusion_matrix(y_hat_train, y_train)

#TPR, FDR and F1-score for all classes for the train set
results = TPR_FDR(c_matrix)

scatter([1,2,3], results[:,2], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "TPR", shape = :diamond, xlims = (0.5,3.5), grid = :y, markersize = 4,
title = "Training data n = $(length(y_train))", titlefont = font(10),
legend = false)

scatter!([1,2,3], results[:,3], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "FDR", shape = :square, markersize = 4)

p_train = scatter!([1,2,3], results[:,4], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "F1 score", shape = :utriangle, dpi = 300, markersize = 4)

#Confusion matrix for the test set
c_matrix = confusion_matrix(y_hat_test, y_test)

#TPR, FDR and F1-score for all classes for the test set
results = TPR_FDR(c_matrix)

scatter([1,2,3], results[:,2], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "TPR", shape = :diamond, xlims = (0.5,3.5), grid = :y, markersize = 4,
title = "Test data n = $(length(y_test))", titlefont = font(10))

scatter!([1,2,3], results[:,3], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "FDR", shape = :square, markersize = 4)

p_test = scatter!([1,2,3], results[:,4], xticks = ([1,2,3],results[:,1]), ylims = (0,1),
label = "F1 score", shape = :utriangle, dpi = 300, markersize = 4)

plot(p_train, p_test, size = (800,400), dpi = 300)